"""The longitudinal example shipped with the docs and the app (tests/data/make_longitudinal_example.py).

Control vs Treated, 4 individuals per condition followed at 0, 6, 12, 24 and 48 h. The 'truth' sheet
(ignored by the General importer) says which proteins were planted with a temporal pattern.
"""
import contextlib
import io

import numpy as np
import pandas as pd
import pytest

import omicscope as omics

np.seterr(divide='ignore', invalid='ignore')
DATA = 'tests/data/proteins/'
EXAMPLE = DATA + 'longitudinal_example.xlsx'


def run(dof=3, **kwargs):
    kwargs.setdefault('Table', EXAMPLE)
    with contextlib.redirect_stdout(io.StringIO()):
        return omics.OmicScope(Method='General', ExperimentalDesign='longitudinal', degrees_of_freedom=dof, **kwargs)


@pytest.fixture(scope='module')
def example():
    truth = pd.read_excel(EXAMPLE, sheet_name='truth')
    return run(), truth


def test_example_layout():
    book = pd.read_excel(EXAMPLE, sheet_name=None)
    assert list(book)[:3] == ['assay', 'rdata', 'pdata']
    pdata = book['pdata']
    assert list(pdata.columns) == ['Sample', 'Condition', 'Biological', 'TimeCourse', 'Individual']
    assert list(book['assay'].columns) == list(pdata['Sample'])
    assert len(book['assay']) == len(book['rdata'])
    # every individual is measured at every time point, and belongs to a single condition
    assert (pdata.groupby('Individual')['TimeCourse'].nunique() == pdata['TimeCourse'].nunique()).all()
    assert (pdata.groupby('Individual')['Condition'].nunique() == 1).all()


def test_example_recovers_planted_patterns(example):
    res, truth = example
    q = res.quant_data.merge(truth, on='Accession')
    assert len(q) > 0.95 * len(truth)
    called, planted = q['pAdjusted'] <= 0.05, q['condition_x_time']
    assert called[planted].mean() >= 0.95
    assert (called & ~planted).sum() / called.sum() <= 0.08
    time_called, time_planted = q['pAdjusted_time'] <= 0.05, q['time_effect']
    assert time_called[q['pattern'] == 'shared_time'].mean() >= 0.95
    assert (time_called & ~time_planted).sum() / time_called.sum() <= 0.05


def test_example_with_separate_pdata(example):
    res, _ = example
    split = run(Table=DATA + 'longitudinal_example_assay.xlsx', pdata=DATA + 'longitudinal_example_pdata.xlsx')
    a = res.quant_data.set_index('Accession')['pvalue']
    b = split.quant_data.set_index('Accession')['pvalue']
    pd.testing.assert_series_equal(a.sort_index(), b.sort_index())
