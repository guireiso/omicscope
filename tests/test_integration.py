"""Comparing independent studies: shared universe, directional concordance, effect-size meta-analysis,
per-study diagnostics and the .omics format, against simulated studies with known effects."""
import contextlib
import io
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import omicscope as omics
from omicscope.MultipleData import Integration, Meta
from omicscope.MultipleData.omics_file import omics_text, parse_omics

GENES = np.array([f'G{i}' for i in range(4000)])


def study(rng, effects, coverage, n=4, sigma=0.5, missing_low=False):
    """One study: genes quantified with probability `coverage` (or the most abundant ones), estimated
    differences of log2 means with their SE and two-sided p-values."""
    if missing_low:
        quantified = np.arange(len(effects)) < int(coverage * len(effects))   # low-abundance genes missing
    else:
        quantified = rng.random(len(effects)) < coverage
    se = sigma * np.sqrt(2 / n) * rng.uniform(0.8, 1.25, len(effects))
    est = rng.normal(effects, se)
    p = 2 * stats.norm.sf(np.abs(est / se))
    df = pd.DataFrame({'gene_name': GENES, 'Accession': GENES, 'pvalue': p, 'log2(fc)': est, 'TotalMean': 1.0,
                       'Mean difference': est, 'SE': se})
    return df[quantified].reset_index(drop=True)


def nebula(*studies, names=None):
    names = names or [f'S{i}' for i in range(len(studies))]
    return SimpleNamespace(original=list(studies), groups=names, pvalue='pvalue', pvalue_cutoff=0.05,
                           colors=['#1b9e77', '#d95f02', '#7570b3', '#e7298a'][:len(studies)],
                           enrichment=[None] * len(studies), metadata=[{}] * len(studies))


def test_shared_universe_removes_coverage_bias_from_the_overlap_test():
    """Two unrelated studies; one quantifies only the most abundant genes, where the other also finds most of
    its hits. Against the union of all genes the overlap looks significant; against the shared universe not."""
    old, new = [], []
    for seed in range(20):
        rng = np.random.default_rng(seed)
        e1 = np.where(rng.random(4000) < 0.15, rng.choice([-1, 1], 4000), 0.0)
        e1[2000:] *= 0.1                                              # study 1: hits among the abundant genes
        e2 = np.where(rng.random(4000) < 0.15, rng.choice([-1, 1], 4000), 0.0)   # unrelated to e1
        neb = nebula(study(rng, e1, 1.0), study(rng, e2, 0.5, missing_low=True))
        new.append(Integration.overlap_test(neb, 'S0', 'S1')['pvalue'])
        old.append(Integration.overlap_test(neb, 'S0', 'S1', background=4000)['pvalue'])
    assert np.mean(np.array(old) < 0.05) > 0.8          # the union background calls unrelated studies similar
    assert np.mean(np.array(new) < 0.05) <= 0.15        # the shared universe does not


def test_directional_concordance_survives_low_power_and_rejects_unrelated_studies():
    rng = np.random.default_rng(1)
    effects = np.where(rng.random(4000) < 0.2, rng.choice([-1, 1], 4000) * rng.uniform(0.4, 1.2, 4000), 0.0)
    strong = study(rng, effects, 0.9, n=8)
    weak = study(rng, effects, 0.9, n=2, sigma=0.9)            # same biology, little power
    unrelated = study(rng, np.where(rng.random(4000) < 0.2, rng.choice([-1, 1], 4000), 0.0), 0.9, n=8)
    neb = nebula(strong, weak, unrelated, names=['strong', 'weak', 'unrelated'])
    c = Integration.concordance(neb, 'strong', 'weak')
    assert c['a_in_b']['fraction'] > 0.75 and c['a_in_b']['pvalue'] < 1e-10
    overlap = Integration.overlap_test(neb, 'strong', 'weak')
    assert overlap['overlap'] < 0.5 * overlap['regulated_a']     # few proteins significant in both...
    u = Integration.concordance(neb, 'strong', 'unrelated')
    assert 0.4 < u['a_in_b']['fraction'] < 0.6 and u['a_in_b']['pvalue'] > 0.01     # ...but direction agrees
    assert set(Integration.concordance_table(neb)['source']) == {'strong', 'weak', 'unrelated'}
    Integration.concordance_plot(neb, 'strong', 'weak')
    plt.close('all')


def test_effect_meta_analysis_estimates_effects_and_heterogeneity():
    rng = np.random.default_rng(2)
    base = np.zeros(4000)
    base[:200] = 1.0                                    # consistent effect in every study
    hetero = np.zeros((3, 4000))
    hetero[:, 200:400] = np.array([[1.2], [1.2], [-1.2]])   # opposite in the third study
    studies = [study(rng, base + hetero[i], 0.95, n=5) for i in range(3)]
    neb = nebula(*studies)
    assert Meta.has_effect_sizes(neb)
    meta = Meta.meta_analysis(neb).set_index('gene_name')
    assert meta.attrs['method'] == 'effect'
    consistent, opposite, null = GENES[:200], GENES[200:400], GENES[400:]
    c = meta.loc[meta.index.intersection(consistent)]
    covered = ((c['CI95 low'] <= 1.0) & (1.0 <= c['CI95 high'])).mean()
    assert 0.88 <= covered <= 0.99                               # ~95% coverage of the true effect
    assert (c['FDR'] < 0.05).mean() > 0.9 and c['I2'].median() < 0.2
    o = meta.loc[meta.index.intersection(opposite)]
    assert o['I2'].median() > 0.7                               # disagreement is reported, not hidden
    n = meta.loc[meta.index.intersection(null)]
    assert (n['FDR'] < 0.05).mean() < 0.01
    stouffer = Meta.meta_analysis(neb, method='stouffer')
    assert stouffer.attrs['method'] == 'stouffer'
    Meta.forest_plot(neb, 'G0')
    plt.close('all')


def test_study_diagnostics_estimate_the_share_of_unchanged_proteins():
    rng = np.random.default_rng(3)
    effects = np.where(rng.random(4000) < 0.3, rng.choice([-1, 1], 4000) * 1.0, 0.0)
    neb = nebula(study(rng, effects, 1.0, n=6), study(rng, effects, 0.7, n=3))
    diag = Integration.study_diagnostics(neb).set_index('study')
    assert abs(diag.loc['S0', 'pi0'] - 0.7) < 0.08
    assert diag.loc['S1', 'quantified'] < diag.loc['S0', 'quantified']
    Integration.pvalue_histograms(neb)
    plt.close('all')


def test_omics_file_round_trip_keeps_effect_sizes_and_metadata():
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table='tests/data/proteins/static_2groups_n3.xlsx', Method='General')
    parsed = parse_omics(omics_text(res))
    assert {'Mean difference', 'SE', 'pvalue'} <= set(parsed['original'].columns)
    assert parsed['metadata']['DataType'] == 'proteome'
    assert parsed['metadata']['Replicates'] == {'Control': '3', 'Treated': '3'}
    assert len(parsed['original']) == len(res.quant_data) and parsed['enrichment'] is None
    # files written by earlier versions still read
    old = parse_omics(open('tests/data/MultipleGroups/omics_file/' +
                           __import__('os').listdir('tests/data/MultipleGroups/omics_file')[0], encoding='utf-8').read())
    assert len(old['original']) > 100 and old['metadata'] == {'DataType': 'proteome'}
