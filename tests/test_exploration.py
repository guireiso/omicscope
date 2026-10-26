"""QC, protein map, enrichment extras and Nebula meta-analysis against simulated ground truths."""
import contextlib
import io
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy import stats

import omicscope as omics
from omicscope.EnrichmentAnalysis import EnrichmentExtras
from omicscope.EnrichmentAnalysis.GSEA import Enrichment
from omicscope.MultipleData import Meta


def _load(name):
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=f'tests/data/proteins/{name}.xlsx', Method='General')
    truth = pd.read_excel(f'tests/data/proteins/{name}.xlsx', sheet_name='rdata').set_index('Accession')
    return res, truth


@pytest.fixture(scope='module')
def four_groups():
    return _load('static_4groups')


def test_missingness_follows_the_detection_limit():
    res, _ = _load('static_onoff')
    per_sample, per_protein = res.missing_summary()
    assert len(per_sample) == len(res.pdata) and (per_sample['quantified'] > 0).all()
    rho = stats.spearmanr(per_protein['mean_log2'], per_protein['missing_fraction'])[0]
    assert rho < -0.3                       # simulated censoring below the LOD: low abundance, more missing
    res.missing_values()
    plt.close('all')


def test_cv_ranks_conditions_by_simulated_noise(four_groups):
    res, _ = four_groups
    med = res.cv_table().groupby('Condition')['CV'].median()
    # simulated biological SD: A = B = 0.25 < C = 0.5 < D = 0.8 (log2)
    assert med['A'] < med['C'] < med['D'] and med['B'] < med['C']
    res.cv_plot()
    plt.close('all')


def test_protein_map_keeps_co_regulated_proteins_together(four_groups):
    from sklearn.neighbors import NearestNeighbors
    res, truth = four_groups
    emb = res.protein_embedding(pvalue=0.05)
    shifts = truth[[c for c in truth.columns if c.startswith('simulated_shift_')]]
    pattern = np.sign(shifts.reindex(emb['Accession']).fillna(0)).astype(int).astype(str).agg(''.join, axis=1).to_numpy()
    de = shifts.reindex(emb['Accession']).abs().sum(axis=1).to_numpy() > 0
    xy = emb[['x', 'y']].to_numpy()[de]
    lab = pattern[de]
    idx = NearestNeighbors(n_neighbors=11).fit(xy).kneighbors(xy, return_distance=False)[:, 1:]
    purity = (lab[idx] == lab[:, None]).mean()
    chance = sum((np.sum(lab == u) / len(lab)) ** 2 for u in np.unique(lab))
    assert purity > 0.6 and purity > 3 * chance
    assert emb['cluster'].nunique() >= 2
    res.protein_map(pvalue=0.05)
    plt.close('all')


def _gsea_object():
    import sys
    sys.path.insert(0, 'tests')
    from test_gsea_ranking import simulate
    obj, sets = simulate(0)
    enr = Enrichment(obj, [sets], padjust_cutoff=0.25, permutation_num=200, threads=1)
    return SimpleNamespace(Analysis='GSEA', results=enr.results, curves=enr.curves, ranking=enr.ranking,
                           OmicScope=obj)


def test_gsea_plot_and_term_foldchange():
    e = _gsea_object()
    assert ('custom', 'TRUE_UP') in e.curves
    fig = EnrichmentExtras.gsea_plot(e, 'TRUE_UP')
    plt.close(fig)
    table = EnrichmentExtras.term_foldchange_table(e, top=5)
    assert (table.loc[table['Term'] == 'TRUE_UP', 'log2(fc)'].median()) > 0
    EnrichmentExtras.term_foldchange(e, top=5)
    plt.close('all')


def _studies(seed=0, n=1500, k=4):
    rng = np.random.default_rng(seed)
    genes = [f'P{i}' for i in range(n)]
    effect = np.zeros(n)
    effect[:60] = 2.0                       # same direction in every study, ~50% power in each one
    flip = np.zeros(n, bool)
    flip[60:120] = True                     # strong but opposite between studies
    original = []
    for s in range(k):
        z = rng.normal(effect, 1)
        z[flip] = rng.normal(3.0 * (1 if s % 2 else -1), 1, flip.sum())
        p = 2 * stats.norm.sf(np.abs(z))
        original.append(pd.DataFrame({'gene_name': genes, 'Accession': genes, 'pvalue': p,
                                      'log2(fc)': z * 0.3, 'TotalMean': 1.0}))
    return SimpleNamespace(original=original, groups=[f'S{i}' for i in range(k)], pvalue='pvalue',
                           pvalue_cutoff=0.05)


def test_meta_analysis_rewards_consistency_not_single_hits():
    neb = _studies()
    meta = Meta.meta_analysis(neb).set_index('gene_name')
    hits = set(meta.index[meta['FDR'] <= 0.05])
    consistent = {f'P{i}' for i in range(60)}
    opposite = {f'P{i}' for i in range(60, 120)}
    assert len(hits & consistent) / 60 >= 0.8
    assert len(hits & opposite) == 0
    assert len(hits - consistent) / max(len(hits), 1) <= 0.1
    assert meta.loc['P0', 'consistency'] >= 0.75
    Meta.meta_plot(neb, top=15)
    plt.close('all')


def test_one_extreme_study_cannot_carry_the_meta_analysis():
    neb = _studies()
    neb.original[0].loc[0:9, 'pvalue'] = 1e-80          # P0-P9 extreme in study 0 only ...
    for df in neb.original[1:]:
        df.loc[0:9, 'pvalue'] = 0.9                       # ... and flat elsewhere
        df.loc[0:9, 'log2(fc)'] = -0.01
    meta = Meta.meta_analysis(neb).set_index('gene_name')
    assert (meta.loc[[f'P{i}' for i in range(10)], 'FDR'] > 0.05).all()
    uncapped = Meta.meta_analysis(neb, z_cap=None).set_index('gene_name')
    assert (uncapped.loc[[f'P{i}' for i in range(10)], 'FDR'] < 0.05).all()   # what the cap prevents


def test_cluster_profiles_follow_the_planted_patterns(four_groups):
    res, truth = four_groups
    emb = res.protein_embedding(pvalue=0.05)
    prof = res.protein_profiles(accessions=emb['Accession'])
    n_samples = prof['Sample'].nunique()
    assert len(prof) == len(emb) * n_samples
    assert prof.groupby('Accession')['z'].mean().abs().max() < 1e-9          # z-scored per protein
    first = prof.drop_duplicates('Sample').sort_values('position')['Condition']
    assert list(dict.fromkeys(first)) == list(res.Conditions)                   # samples grouped by condition
    # the mean profile of proteins planted up in D is highest in the D samples
    up_d = truth.index[truth['simulated_shift_D'] > 0.8].intersection(emb['Accession'])
    mean = prof[prof['Accession'].isin(up_d)].groupby('Condition')['z'].mean()
    assert mean.idxmax() == 'D'
    out = res.cluster_profile_plot(embedding=emb)
    assert out['cluster'].notna().all()
    plt.close('all')
