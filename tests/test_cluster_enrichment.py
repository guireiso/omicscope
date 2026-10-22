"""Cluster enrichment, custom gene sets (dict / GMT) and redundant-term collapsing, offline."""
import contextlib
import io
import sys
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import omicscope as omics
from omicscope.EnrichmentAnalysis import ClusterEnrichment as ce
from omicscope.General import Modules

sys.path.insert(0, 'tests')
from test_modules import simulate  # noqa: E402


@pytest.fixture(scope='module')
def experiment(tmp_path_factory):
    """Planted modules (factor 1 follows the condition) and gene sets built from the truth."""
    path = tmp_path_factory.mktemp('cluster_enrichment') / 'planted_modules.xlsx'
    truth = simulate(path)
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=str(path), Method='General')
    gene = res.quant_data.set_index('Accession')['gene_name']
    rng = np.random.default_rng(11)
    all_genes = list(gene)
    sets = {}
    for f in range(1, 5):
        members = list(gene.reindex(truth[truth == f].index).dropna())
        # 80% of the planted proteins plus some unrelated ones, as in a real pathway annotation
        sets[f'SET_{f}'] = list(rng.choice(members, 80, replace=False)) + list(rng.choice(all_genes, 20, replace=False))
    sets.update({f'DECOY_{k}': list(rng.choice(all_genes, 60, replace=False)) for k in range(15)})
    return res, truth, sets


def test_read_gmt_round_trip(tmp_path):
    sets = {'TERM_A': ['TP53', 'MDM2', 'CDKN1A'], 'TERM B (with space)': ['ACTB', 'GAPDH']}
    path = tmp_path / 'my_sets.gmt'
    ce.write_gmt(sets, path)
    assert ce.read_gmt(path) == sets
    path.write_text('X\tdesc\tA\tB\t\tA\nshort\tonly\n', encoding='utf-8')     # empty fields, duplicates, no genes
    assert ce.read_gmt(path) == {'X': ['A', 'B']}
    assert ce.gene_set_library(str(path)) == ({'X': ['A', 'B']}, 'my_sets')
    assert ce.gene_set_library({'X': ['A']})[1] == 'custom'
    assert ce.gene_set_library('KEGG_2021_Human') == ('KEGG_2021_Human', 'KEGG_2021_Human')


def test_each_cluster_finds_its_own_gene_set(experiment):
    res, truth, sets = experiment
    modules = Modules.coexpression_modules(res)
    table = modules['modules'].set_index('Accession')
    factor_of = {m: truth.reindex(g.index).mode()[0] for m, g in table[table['module'] > 0].groupby('module')}
    out = ce.cluster_enrichment(res, modules, dbs=[sets])
    assert set(out['Gene_set']) == {'custom'}
    assert 0 not in set(out['Cluster'])
    for m, f in factor_of.items():
        if f == 0:
            continue
        top = out[out['Cluster'] == m].sort_values('Adjusted P-value').iloc[0]
        assert top['Term'] == f'SET_{f}' and top['Adjusted P-value'] < 1e-10
    assert (out.loc[out['Term'].str.startswith('DECOY'), 'Adjusted P-value'] > 1e-3).all()   # chance overlaps only
    assert (out['GeneRatio'] <= 1).all() and (out['GeneRatio'] > 0).all()
    ce.cluster_dotplot(out, top=3)
    plt.close('all')


def test_cluster_input_formats(experiment):
    res, truth, sets = experiment
    gene = res.quant_data.set_index('Accession')['gene_name']
    planted = {f'factor {f}': list(gene.reindex(truth[truth == f].index).dropna()) for f in (1, 2)}
    from_dict = ce.cluster_enrichment(res, planted, dbs=[sets])
    df = pd.DataFrame([(g, k) for k, v in planted.items() for g in v], columns=['gene_name', 'Cluster'])
    from_frame = ce.cluster_enrichment(res, df, dbs=[sets])
    assert set(from_dict['Cluster']) == set(planted)
    pd.testing.assert_series_equal(from_dict.groupby('Cluster')['Term'].first(),
                                   from_frame.groupby('Cluster')['Term'].first())


def test_ora_and_gsea_accept_dicts_and_gmt_files(experiment, tmp_path):
    res, _, sets = experiment
    gmt = tmp_path / 'planted_sets.gmt'
    ce.write_gmt(sets, gmt)
    ora = omics.EnrichmentScope(res, Analysis='ORA', dbs=[str(gmt), sets], padjust_cutoff=1, background=True)
    assert set(ora.results['Gene_set']) == {'planted_sets', 'custom'}
    for name in ('planted_sets', 'custom'):
        assert ora.results[ora.results['Gene_set'] == name].iloc[0]['Term'] == 'SET_1'   # condition-driven factor
    gsea = omics.EnrichmentScope(res, Analysis='GSEA', dbs=str(gmt), padjust_cutoff=1, permutation_num=200,
                                 threads=1)
    assert set(gsea.results['Gene_set']) == {'planted_sets'}
    assert gsea.results.iloc[0]['Term'] == 'SET_1' and gsea.results.iloc[0]['Adjusted P-value'] < 0.05
    gsea_dict = omics.EnrichmentScope(res, Analysis='GSEA', dbs=[sets], padjust_cutoff=1, permutation_num=200,
                                      threads=1)
    assert set(gsea_dict.results['Gene_set']) == {'custom'} and ('custom', 'SET_1') in gsea_dict.curves


def _results(rows, gsea=False):
    cols = ['Term', 'Adjusted P-value', 'Genes'] + (['NES'] if gsea else [])
    return SimpleNamespace(results=pd.DataFrame(rows, columns=cols), Analysis='GSEA' if gsea else 'ORA')


def test_collapse_terms_keeps_the_most_significant_duplicate():
    g = [f'G{i}' for i in range(300)]
    e = _results([('A copy', 1e-3, g[:36] + g[200:204]),     # near duplicate of A (Jaccard 0.82)
                  ('A', 1e-6, g[:40]),
                  ('B', 1e-4, g[100:140]),                   # distinct
                  ('C', 1e-2, g[:12] + g[250:278])])         # shares 30% with A: stays apart
    out = ce.collapse_terms(e, similarity=0.5)
    rep = out.set_index('Term')['Representative']
    assert rep['A copy'] == 'A' and rep['A'] == 'A'
    assert rep['B'] == 'B' and rep['C'] == 'C'
    assert out.set_index('Term').loc['A', 'Cluster_size'] == 2
    only = ce.collapse_terms(e, similarity=0.5, representative_only=True)
    assert list(only['Term']) == ['A', 'B', 'C']
    # the overlap coefficient also groups a term with its sub-terms
    sub = _results([('parent', 1e-5, g[:60]), ('child', 1e-3, g[:15])])
    assert ce.collapse_terms(sub, metric='jaccard')['Representative'].tolist() == ['parent', 'child']
    assert ce.collapse_terms(sub, metric='overlap')['Representative'].tolist() == ['parent', 'parent']


def test_collapse_terms_keeps_gsea_directions_apart():
    g = [f'G{i}' for i in range(100)]
    e = _results([('up', 1e-4, g[:40], 2.1), ('down', 1e-3, g[:40], -1.9), ('up again', 1e-2, g[:38], 1.5)], gsea=True)
    rep = ce.collapse_terms(e).set_index('Term')['Representative']
    assert rep['up again'] == 'up' and rep['down'] == 'down'


def test_collapse_terms_on_real_ora(experiment):
    res, _, sets = experiment
    dup = dict(sets)
    dup['SET_1_COPY'] = sets['SET_1'][:90] + sets['DECOY_0'][:10]
    ora = omics.EnrichmentScope(res, Analysis='ORA', dbs=[dup], padjust_cutoff=1, background=True)
    out = ce.collapse_terms(ora).set_index('Term')
    assert out.loc['SET_1_COPY', 'Representative'] == 'SET_1'
    assert out.loc['SET_1', 'Representative'] == 'SET_1'
