"""Presence/absence layer against simulated on/off entities (tests/data/make_static.py)."""
import io
import contextlib

import numpy as np
import pandas as pd
import pytest

import omicscope as omics


def _run(name, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=f'tests/data/proteins/{name}.xlsx', Method='General', **kw)
    truth = pd.read_excel(f'tests/data/proteins/{name}.xlsx', sheet_name='rdata').set_index('Accession')
    called = set(res.quant_data.loc[res.quant_data['pAdjusted'] <= 0.05, 'Accession'])
    onoff = set(truth.index[truth['simulated_onoff'].fillna('') != ''])
    null = set(truth.index[~truth['simulated_de']])
    return res, len(called & onoff) / len(onoff), len(called & null) / max(len(called), 1)


@pytest.mark.parametrize('name', ['static_onoff', 'static_onoff_3groups'])
def test_on_off_entities_are_recovered_with_fdr_control(name):
    res, found, fdr = _run(name)
    assert found >= 0.9 and fdr <= 0.07
    det = res.quant_data[res.quant_data['Test'] == 'detection']
    assert det['fc_censored'].all() and det['Regulation'].notna().all()


def test_exclusive_label_points_to_the_right_condition():
    res, _, _ = _run('static_onoff')
    truth = pd.read_excel('tests/data/proteins/static_onoff.xlsx', sheet_name='rdata').set_index('Accession')
    exclusive = res.quant_data[res.quant_data['Regulation'].fillna('').str.startswith('exclusive to')]
    absent_in = exclusive['Accession'].map(truth['simulated_onoff'])
    present = exclusive['Regulation'].str.replace('exclusive to ', '')
    known = absent_in.fillna('') != ''
    assert (present[known] != absent_in[known]).mean() > 0.95     # exclusive to the other condition
    # the censored fold change points the right way: positive when exclusive to the experimental group
    assert np.all(np.sign(exclusive['log2(fc)']) == np.where(present == 'Treated', 1, -1))


def test_previous_behaviour_is_available():
    res, found, _ = _run('static_onoff', detection_layer=False)
    assert found < 0.1 and 'detection' not in set(res.quant_data['Test'])


def test_minprob_imputation_and_posthoc_with_detection_rows():
    _, found, fdr = _run('static_onoff', imputation_method='minprob')
    assert found >= 0.9 and fdr <= 0.07
    res, _, _ = _run('static_onoff_3groups', pvalue='pTukey')
    det = res.quant_data[res.quant_data['Test'] == 'detection']
    assert det['pTukey'].notna().all()


@pytest.mark.parametrize('name, welch', [('static_2groups_n3', False), ('static_4groups', True)])
def test_variance_model_is_chosen_automatically(name, welch):
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=f'tests/data/proteins/{name}.xlsx', Method='General')
    assert res._use_welch is welch and res._use_moderated is (not welch)
    assert 'Variance ratio' in res.Params['Params']['Stats_Variance_Model']
