"""PLS-DA, ROC tables and cross-validated biomarker panels against simulated ground truths.

Two 20 vs 20 datasets are simulated with tests/data/make_static.py: one where 20% of the proteins
change (truth in rdata `simulated_de`) and one with no difference at all. The null dataset is the key
check: a supervised model always separates its own training samples, so only an honest
cross-validation (selection and scaling inside the folds) gives chance-level accuracy / AUC there.
"""
import contextlib
import copy
import importlib.util
import io

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import omicscope as omics
from omicscope.General import Multivariate

_spec = importlib.util.spec_from_file_location('make_static', 'tests/data/make_static.py')
make_static = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(make_static)


def _simulate(path, de_share, n=20, seed=11):
    rng = np.random.default_rng(seed)
    source = (pd.read_excel('tests/data/proteins/general.xlsx', sheet_name='rdata')
              .dropna(subset=['Description']).drop_duplicates('Accession'))
    assay, rdata, pdata = make_static.build(['Control', 'Treated'], n, 1, 2000, de_share, [0.35, 0.35], rng, source)
    with pd.ExcelWriter(path) as writer:
        assay.to_excel(writer, sheet_name='assay', index=False)
        rdata.to_excel(writer, sheet_name='rdata', index=False)
        pdata.to_excel(writer, sheet_name='pdata', index=False)
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=str(path), Method='General')
    return res, rdata.set_index('Accession')


@pytest.fixture(scope='module')
def de(tmp_path_factory):
    return _simulate(tmp_path_factory.mktemp('mv') / 'de.xlsx', de_share=0.2)


@pytest.fixture(scope='module')
def null(tmp_path_factory):
    return _simulate(tmp_path_factory.mktemp('mv') / 'null.xlsx', de_share=0.0, seed=12)


def _shuffled(res, seed=0):
    """Same data, condition labels shuffled across samples."""
    out = copy.copy(res)
    out.pdata = res.pdata.copy()
    out.pdata['Condition'] = np.random.default_rng(seed).permutation(res.pdata['Condition'].to_numpy())
    return out


def test_plsda_separates_real_groups_and_ranks_true_proteins(de):
    res, truth = de
    out = Multivariate.plsda(res)
    assert out['cv_accuracy'] >= 0.95 and out['permutation_p'] <= 0.02
    assert list(out['scores'].columns[:4]) == ['Sample', 'Condition', 'comp1', 'comp2']
    top = out['vip'].head(100)
    assert truth.loc[top['Accession'], 'simulated_de'].mean() >= 0.9      # 20% of all proteins are DE
    assert np.isclose(np.mean(out['vip']['VIP'] ** 2), 1)                  # VIP: mean square is 1 by construction


@pytest.mark.parametrize('case', ['null simulation', 'shuffled labels'])
def test_plsda_does_not_overfit_without_differences(de, null, case):
    res = null[0] if case == 'null simulation' else _shuffled(de[0])
    out = Multivariate.plsda(res)
    # in-sample the groups are always well separated on comp 1...
    comp1 = out['scores'].groupby('Condition')['comp1']
    assert comp1.mean().diff().abs().iloc[-1] > 2 * comp1.std().max()
    # ...but not out of sample
    assert out['cv_accuracy'] <= 0.7 and out['permutation_p'] >= 0.1


def test_auc_of_true_proteins_and_null_confidence_intervals(de, null):
    res, truth = de
    tab = Multivariate.roc_table(res).set_index('Accession')
    is_de = truth['simulated_de'].reindex(tab.index).astype(bool)
    assert tab.loc[is_de, 'AUC'].median() > 0.95 and tab.loc[~is_de, 'AUC'].median() < 0.65
    fc = truth['simulated_log2fc'].reindex(tab.index)
    assert (tab.loc[is_de, 'direction'] == np.where(fc[is_de] > 0, 'up', 'down')).mean() > 0.99
    assert (tab['AUC'] >= 0.5).all() and (tab['CI_low'] <= tab['AUC']).all()
    ntab = Multivariate.roc_table(null[0])
    assert (ntab['CI_low'] <= 0.5).mean() >= 0.9                          # nominal 95%


def test_panel_cross_validation_has_no_leakage(de, null):
    real = Multivariate.biomarker_panel(de[0])
    assert real['auc_mean'] >= 0.95
    noise = Multivariate.biomarker_panel(null[0])
    assert abs(noise['auc_mean'] - 0.5) <= 0.15
    assert len(noise['fpr']) == len(noise['tpr']) and noise['selection']['frequency'].max() <= 1


def test_three_or_more_conditions_take_control_and_case():
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table='tests/data/proteins/static_4groups.xlsx', Method='General')
    truth = pd.read_excel('tests/data/proteins/static_4groups.xlsx', sheet_name='rdata').set_index('Accession')
    tab = Multivariate.roc_table(res, control='A', case='D').set_index('Accession')
    moved = truth['simulated_shift_D'].reindex(tab.index).abs() > 0
    assert tab.loc[moved, 'AUC'].mean() > tab.loc[~moved, 'AUC'].mean() + 0.2
    with pytest.raises(ValueError):
        Multivariate.roc_table(res, control='A', case='A')


def test_figures_run(de, tmp_path):
    res, _ = de
    out = Multivariate.plsda_plot(res, permutations=10, save=str(tmp_path) + '/', vector=False)
    assert (tmp_path / 'plsda.png').exists() and len(out['vip']) > 0
    table, panel = Multivariate.roc_plot(res, save=str(tmp_path) + '/', repeats=2)
    assert (tmp_path / 'roc_plot.svg').exists() and panel['auc'] > 0.9
    Multivariate.roc_plot(res, *table['gene_name'].head(2), panel=False)
    plt.close('all')
