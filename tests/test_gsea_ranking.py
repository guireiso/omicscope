"""Pre-ranked GSEA against a simulated ground truth (local gene sets, no network access)."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from omicscope.EnrichmentAnalysis.GSEA import Enrichment, rank_entities


def simulate(seed=0, n=2000, reps=4):
    rng = np.random.default_rng(seed)
    genes = [f'Gene{i}' for i in range(n)]
    sd = rng.uniform(0.15, 0.4, n)
    shift = np.zeros(n)
    shift[:40] = 0.35                     # modest, consistent up-regulation: a real pathway effect
    sd[40:80] = 2.0                        # high-variance proteins with no true change
    ctrl = rng.normal(0, sd[:, None], (n, reps))
    exp = rng.normal(shift[:, None], sd[:, None], (n, reps))
    t, p = stats.ttest_ind(exp, ctrl, axis=1)
    qd = pd.DataFrame({'gene_name': genes, 'log2(fc)': exp.mean(1) - ctrl.mean(1), 'pvalue': p})
    sets = {'TRUE_UP': genes[:40], 'NOISY_NULL': genes[40:80]}
    sets.update({f'NULL_{k}': list(rng.choice(genes[80:], 40, replace=False)) for k in range(20)})
    return SimpleNamespace(quant_data=qd, Conditions=['A', 'B']), sets


def run(obj, sets, rank_by):
    return Enrichment(obj, [sets], padjust_cutoff=1, rank_by=rank_by, permutation_num=1000, threads=1).results.set_index('Term')


def test_ranking_cleans_names_and_duplicates():
    qd = pd.DataFrame({'gene_name': ['A', 'A', None, 'nan', 'B', ' C '], 'log2(fc)': [1, -3, 2, 2, np.nan, 0.5],
                       'pvalue': [0.01, 0.2, 0.01, 0.01, 0.01, 0.5]})
    r = rank_entities(qd)
    assert set(r.index) == {'A', 'C'}
    assert r['A'] == pytest.approx(-np.log10(0.01) + 1e-6, rel=1e-6)   # keeps the strongest of the duplicates


@pytest.mark.parametrize('seed', [0, 1, 2])
def test_signed_p_finds_the_true_set_and_ignores_noise(seed):
    obj, sets = simulate(seed)
    res = run(obj, sets, 'signed_p')
    assert res.loc['TRUE_UP', 'Adjusted P-value'] < 0.05 and res.loc['TRUE_UP', 'NES'] > 0
    assert res.loc['NOISY_NULL', 'Adjusted P-value'] > 0.25
    nulls = res[res.index.str.startswith('NULL_')]['Adjusted P-value']
    assert (nulls < 0.05).sum() <= 1


def test_mixed_case_gene_names_do_not_break_leading_edge():
    obj, sets = simulate(0)
    obj.quant_data['gene_name'] = obj.quant_data['gene_name'].str.capitalize()   # mouse-style symbols
    sets = {k: [g.capitalize() for g in v] for k, v in sets.items()}
    res = run(obj, sets, 'signed_p')
    assert res.loc['TRUE_UP', 'up-regulated'] > 0
    assert np.isfinite(res['-log10(pAdj)']).all()


def test_public_entry_point_forwards_gsea_options():
    import inspect
    import omicscope as omics
    assert 'gsea_options' in inspect.signature(omics.EnrichmentScope).parameters
