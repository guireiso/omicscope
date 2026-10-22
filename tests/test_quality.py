"""Sample outliers and batch correction (General/Quality.py) against simulated ground truth."""
import contextlib
import inspect
import io

import matplotlib
import numpy as np
import pandas as pd
import pytest

import omicscope as omics
from omicscope.General.Omicscope import Omicscope
from omicscope.General.Quality import (batch_effect, batch_plot, outlier_plot, remove_batch_effect,
                                       sample_outliers)

matplotlib.use('Agg')


def simulate(seed, n_per_cell=4, batches=('b1', 'b2'), conditions=('Ctrl', 'Treated'), n=2000, de_share=0.15,
             batch_share=0.6, batch_shift=1.5, bio_sd=0.4, outlier=None, missing=True):
    """log2 abundances = baseline + condition effect (DE proteins) + batch effect + biological noise.

    outlier: None, 'noise' (sample 2 gets extra noise, sd 1) or 'shift' (half of its proteins shifted by 1).
    Returns assay (linear, 0 = missing), rdata (truth in simulated_de / simulated_log2fc / simulated_batch)
    and pdata.
    """
    rng = np.random.default_rng(seed)
    baseline = rng.normal(22, 2, n)
    de = rng.random(n) < de_share
    log2fc = np.where(de, rng.choice([-1, 1], n) * rng.uniform(0.8, 2.0, n), 0.0)
    shifted = rng.random(n) < batch_share
    batch_effects = {b: np.where(shifted, rng.normal(0, batch_shift, n), 0.0) for b in batches}
    samples, pdata, bio = {}, [], 0
    for ci, cond in enumerate(conditions):
        for b in batches:
            for _ in range(n_per_cell):
                bio += 1
                log2 = baseline + (log2fc if ci == 1 else 0) + batch_effects[b] + rng.normal(0, bio_sd, n)
                if outlier is not None and bio == 3:
                    log2 = log2 + (rng.normal(0, 1.0, n) if outlier == 'noise' else (rng.random(n) < 0.5) * 1.0)
                value = np.round(2 ** log2, 2)
                if missing:
                    value[rng.random(n) < 0.15 / (1 + np.exp((log2 - 19) * 1.5))] = 0
                name = f'S{bio}'
                samples[name] = value
                pdata.append({'Sample': name, 'Condition': cond, 'Biological': bio, 'Batch': b})
    rdata = pd.DataFrame({'Accession': [f'P{i:05d}' for i in range(n)],
                          'Description': [f'Protein {i} OS=Homo sapiens OX=9606 GN=G{i} PE=1 SV=1' for i in range(n)],
                          'simulated_de': de, 'simulated_log2fc': log2fc.round(3), 'simulated_batch': shifted})
    return pd.DataFrame(samples), rdata, pd.DataFrame(pdata)


def write(path, assay, rdata, pdata):
    with pd.ExcelWriter(path) as writer:
        assay.to_excel(writer, sheet_name='assay', index=False)
        rdata.to_excel(writer, sheet_name='rdata', index=False)
        pdata.to_excel(writer, sheet_name='pdata', index=False)
    return str(path)


def run(path, batch_correction=False, **kw):
    """OmicScope on a General xlsx; batch_correction is passed as a keyword once OmicScope accepts it."""
    with contextlib.redirect_stdout(io.StringIO()):
        if 'batch_correction' in inspect.signature(Omicscope.__init__).parameters:
            return omics.OmicScope(Table=path, Method='General', batch_correction=batch_correction, **kw)
        res = Omicscope.__new__(Omicscope)
        res.batch_correction = batch_correction
        res.__init__(path, 'General', None, **kw)
        return res


def called(res, rdata):
    qd = res.quant_data[res.quant_data['Test'] == 'intensity']
    hits = set(qd.loc[qd['pAdjusted'] <= 0.05, 'Accession'])
    truth = set(rdata.loc[rdata['simulated_de'], 'Accession'])
    tested = truth & set(qd['Accession'])
    return len(hits & tested) / len(tested), len(hits - truth) / max(len(hits), 1)


@pytest.mark.parametrize('kind', ['noise', 'shift'])
def test_injected_outlier_is_flagged_and_clean_samples_are_not(tmp_path, kind):
    for seed in range(3):
        assay, rdata, pdata = simulate(seed, n=1500, batch_share=0, outlier=kind)
        res = run(write(tmp_path / f'out{seed}.xlsx', assay, rdata, pdata.drop(columns='Batch')))
        out = sample_outliers(res)
        out.index = out['Sample'].map(res.pdata.set_index('Sample')['Biological'])
        assert out.loc[3, 'outlier'], (seed, out.loc[3])
        assert not out.drop(index=3)['outlier'].any(), (seed, out[out['outlier']])


def test_no_outlier_in_clean_data_with_unequal_condition_variances():
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table='tests/data/proteins/static_4groups.xlsx', Method='General')
    assert not sample_outliers(res)['outlier'].any()


def test_remove_batch_effect_keeps_condition_effects():
    assay, rdata, pdata = simulate(7, n=2000)
    logm = np.log2(assay.replace(0, np.nan))
    corrected = remove_batch_effect(logm, pdata['Condition'], pdata['Batch'])
    assert corrected.isna().equals(logm.isna())
    b1, b2 = (pdata['Batch'] == 'b1').to_numpy(), (pdata['Batch'] == 'b2').to_numpy()
    ctrl, trt = (pdata['Condition'] == 'Ctrl').to_numpy(), (pdata['Condition'] == 'Treated').to_numpy()

    def batch_gap(m):      # batch difference within conditions, averaged over the two conditions
        m = m.to_numpy()
        return np.nanmean([np.nanmean(m[:, b1 & c], axis=1) - np.nanmean(m[:, b2 & c], axis=1) for c in (ctrl, trt)],
                          axis=0)

    before, after = batch_gap(logm), batch_gap(corrected)
    assert np.nanstd(before) > 1.0 and np.nanstd(after) < 0.05
    # condition effects of the DE proteins are kept
    fc = np.nanmean(corrected.to_numpy()[:, trt], axis=1) - np.nanmean(corrected.to_numpy()[:, ctrl], axis=1)
    de = rdata['simulated_de'].to_numpy()
    ok = de & np.isfinite(fc)
    assert np.corrcoef(fc[ok], rdata['simulated_log2fc'].to_numpy()[ok])[0, 1] > 0.95
    # complete proteins: batch R² drops to ~0, condition R² is kept or rises
    from omicscope.General.Quality import _variance_explained
    complete = logm.notna().all(axis=1).to_numpy()
    rb0, rc0 = _variance_explained(logm.to_numpy()[complete], pdata['Condition'], pdata['Batch'])
    rb1, rc1 = _variance_explained(corrected.to_numpy()[complete], pdata['Condition'], pdata['Batch'])
    shifted = rdata['simulated_batch'].to_numpy()[complete]
    assert np.median(rb0[shifted]) > 0.3 and np.median(rb1) < 0.01
    assert np.median(rc1[de[complete]]) >= np.median(rc0[de[complete]])


def test_remove_batch_effect_is_fast_with_missing_values():
    import time
    rng = np.random.default_rng(1)
    y = rng.normal(20, 1, (3000, 40))
    y[rng.random(y.shape) < 0.1] = np.nan
    cond, batch = np.repeat(['A', 'B'], 20), np.tile(np.repeat(['x', 'y', 'z', 'w'], 5), 2)
    start = time.perf_counter()
    out = remove_batch_effect(y, cond, batch)
    assert time.perf_counter() - start < 2 and np.array_equal(np.isnan(out), np.isnan(y))


@pytest.mark.parametrize('batches', [['b1'] * 4 + ['b2'] * 4, ['b1'] * 2 + ['b2'] * 2 + ['b3'] * 2 + ['b4'] * 2])
def test_confounded_design_raises(batches):
    y = np.random.default_rng(0).normal(size=(50, 8))
    with pytest.raises(ValueError, match='confounded'):
        remove_batch_effect(y, ['A'] * 4 + ['B'] * 4, batches)


def test_batch_correction_in_statistics_recovers_more_de_proteins(tmp_path):
    results = []
    for seed in range(2):
        assay, rdata, pdata = simulate(100 + seed)
        path = write(tmp_path / f'batch{seed}.xlsx', assay, rdata, pdata)
        plain = run(path)
        corrected = run(path, batch_correction=True)
        assert 'Stats_Batch_Correction' in corrected.Params['Params']
        assert all(s.endswith('.' + c) for s, c in zip(corrected.pdata['Sample'], corrected.pdata['Condition']))
        results.append(called(plain, rdata) + called(corrected, rdata))
    power0, fdr0, power1, fdr1 = np.mean(results, axis=0)
    assert power1 > power0 + 0.1 and fdr1 <= 0.07


def test_figures_run(tmp_path):
    assay, rdata, pdata = simulate(3, n=800, outlier='noise')
    res = run(write(tmp_path / 'fig.xlsx', assay, rdata, pdata))
    out = outlier_plot(res, save=str(tmp_path) + '/', vector=False)
    assert out['outlier'].sum() >= 1 and (tmp_path / 'outlier_plot.png').exists()
    result = batch_plot(res, save=str(tmp_path) + '/', vector=False)
    assert (tmp_path / 'batch_plot.png').exists()
    assert result['summary'].loc['after', 'batch'] < result['summary'].loc['before', 'batch']
    assert set(batch_effect(res)) == {'proteins', 'summary', 'before', 'after'}
