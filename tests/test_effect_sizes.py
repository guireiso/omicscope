"""Effect sizes and confidence intervals against simulated truth (tests/data/make_static.py)."""
import contextlib
import io

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import omicscope as omics
from omicscope.General import Effects


def _load(name, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=f'tests/data/proteins/{name}.xlsx', Method='General', **kw)
    truth = pd.read_excel(f'tests/data/proteins/{name}.xlsx', sheet_name='rdata').set_index('Accession')
    return res, truth


@pytest.mark.parametrize('kw', [{}, {'moderated': False}, {'welch': True}])
def test_confidence_intervals_cover_the_true_fold_change(kw):
    res, truth = _load('static_2groups_n3', **kw)
    qd = res.quant_data[res.quant_data['Test'] == 'intensity'].set_index('Accession')
    true = truth['simulated_log2fc'].reindex(qd.index)
    covered = (qd['CI95 low'] <= true) & (true <= qd['CI95 high'])
    assert 0.92 <= covered.mean() <= 0.98          # nominal 95%
    # the interval excludes zero exactly when the nominal p-value is below 0.05
    excludes = (qd['CI95 low'] > 0) | (qd['CI95 high'] < 0)
    assert (excludes == (qd['pvalue'] < 0.05)).mean() > 0.999
    assert (np.sign(qd["Hedges' g"]) == np.sign(qd['Mean difference'])).all()


def test_eta_squared_ranks_true_changes_for_several_conditions():
    res, truth = _load('static_4groups')
    qd = res.quant_data.set_index('Accession')
    de = truth['simulated_de'].reindex(qd.index).astype(bool)
    assert qd['Eta squared'].between(0, 1).all()
    assert qd.loc[de, 'Eta squared'].median() > 2 * qd.loc[~de, 'Eta squared'].median()


def test_forest_plot():
    res, _ = _load('static_2groups_n3')
    df = Effects.effect_plot(res, top=10)
    assert len(df) == 10
    plt.close('all')
