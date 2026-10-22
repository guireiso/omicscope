"""Co-expression modules against simulated latent factors (planted modules)."""
import contextlib
import io

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import adjusted_rand_score

import omicscope as omics
from omicscope.General import Modules


def simulate(path, seed=7, n_per=100, n_noise=300, n_bio=8, timecourse=False):
    """Four latent factors drive 100 proteins each (factor 1 follows the condition); 300 proteins are noise.

    rdata column `planted` holds the factor of each protein (0 = noise).
    """
    rng = np.random.default_rng(seed)
    conditions = np.repeat(['Control', 'Treated'], n_bio)
    n_samples = len(conditions)
    factors = rng.normal(0, 1, (4, n_samples))
    factors[0] = 2 * (conditions == 'Treated') - 1 + rng.normal(0, 0.6, n_samples)
    factors = (factors - factors.mean(1, keepdims=True)) / factors.std(1, keepdims=True)
    log2, planted = [], []
    for f in range(4):
        loading = rng.uniform(0.5, 1.0, n_per)
        log2.append(rng.normal(22, 1.5, (n_per, 1)) + loading[:, None] * factors[f]
                    + rng.normal(0, 0.5, (n_per, n_samples)))
        planted += [f + 1] * n_per
    log2.append(rng.normal(22, 1.5, (n_noise, 1)) + rng.normal(0, 0.8, (n_noise, n_samples)))
    planted += [0] * n_noise
    log2 = np.vstack(log2)
    value = np.round(2 ** log2, 2)
    value[rng.random(value.shape) < 0.01] = 0                       # a few missing values
    names = [f'{c}_{i % n_bio + 1}' for i, c in enumerate(conditions)]
    n = len(planted)
    rdata = pd.DataFrame({'Accession': [f'SIM{i:04d}' for i in range(n)], 'gene_name': [f'GENE{i}' for i in range(n)],
                          'planted': planted})
    pdata = pd.DataFrame({'Sample': names, 'Condition': conditions, 'Biological': np.arange(1, n_samples + 1)})
    if timecourse:
        pdata['TimeCourse'] = np.tile(np.arange(n_bio), 2)
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame(value, columns=names).to_excel(writer, sheet_name='assay', index=False)
        rdata.to_excel(writer, sheet_name='rdata', index=False)
        pdata.to_excel(writer, sheet_name='pdata', index=False)
    return rdata.set_index('Accession')['planted']


@pytest.fixture(scope='module')
def planted(tmp_path_factory):
    path = tmp_path_factory.mktemp('modules') / 'planted_modules.xlsx'
    truth = simulate(path)
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=str(path), Method='General')
    return res, truth, Modules.coexpression_modules(res)


def test_planted_modules_are_recovered(planted):
    _, truth, out = planted
    mods = out['modules'].set_index('Accession')
    t = truth.reindex(mods.index)
    real = t > 0
    assert adjusted_rand_score(t[real], mods.loc[real, 'module']) > 0.7
    assert (mods.loc[~real, 'module'] == 0).mean() > 0.75                # noise proteins stay grey
    assert (mods.loc[real, 'module'] > 0).mean() > 0.85
    assert 1 <= out['power'] <= 20
    assert out['eigengenes'].shape == (16, mods['module'].max())
    assigned = mods[mods['module'] > 0]
    assert (assigned['kME'] >= 0.5).all()


def test_condition_module_has_the_strongest_trait_correlation(planted):
    _, truth, out = planted
    mods = out['modules'].set_index('Accession')
    # module holding most proteins of the condition-driven factor
    condition_module = mods.loc[truth[truth == 1].index.intersection(mods.index), 'module'].mode()[0]
    mt = out['module_trait']
    treated = mt[mt['trait'] == 'Treated'].set_index('module')
    strongest = treated['r'].abs().idxmax()
    assert strongest == condition_module
    assert abs(treated.loc[strongest, 'r']) > 0.8 and treated.loc[strongest, 'pvalue'] < 1e-3
    others = treated.drop(strongest)
    assert (others['r'].abs() < 0.6).all()


def test_eigengene_sign_follows_the_module_mean(planted):
    res, _, out = planted
    mods = out['modules']
    logm = np.log2(res.quant_data.set_index('Accession')[out['eigengenes'].index].replace(0, np.nan))
    for col in out['eigengenes'].columns:
        members = mods.loc[mods['module'] == int(col[2:]), 'Accession']
        mean = logm.loc[members].sub(logm.loc[members].mean(axis=1), axis=0).mean()
        assert np.corrcoef(mean.reindex(out['eigengenes'].index), out['eigengenes'][col])[0, 1] > 0.9


def test_timecourse_trait_and_figures(tmp_path):
    path = tmp_path / 'planted_time.xlsx'
    simulate(path, seed=3, timecourse=True)
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=str(path), Method='General')
    out = Modules.coexpression_modules(res, power=8)
    assert out['power'] == 8
    assert 'TimeCourse' in set(out['module_trait']['trait'])
    Modules.module_trait_plot(res, modules=out, save=str(tmp_path) + '/', vector=False)
    Modules.module_eigengene_plot(res, modules=out)
    assert (tmp_path / 'module_trait_plot.png').exists()
    plt.close('all')
