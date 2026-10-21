"""Longitudinal statistics against simulated data with a known answer (see tests/data/make_longitudinal.py).

The spline test (Storey, 2005) asks whether the temporal profile differs between conditions, so only the
treatment_* patterns are true positives; flat, time-only and constant-offset entities are true negatives.
"""
import numpy as np
import pytest

import omicscope as omics

np.seterr(divide='ignore', invalid='ignore')
INTERACTION = {'treatment_up', 'treatment_down', 'treatment_late', 'treatment_peak'}


@pytest.mark.parametrize('design', ['independent', 'related'])
@pytest.mark.parametrize('dof', [2, 3])
def test_longitudinal_recovers_interactions(design, dof):
    res = omics.OmicScope(Table=f'tests/data/proteins/longitudinal_{design}.xlsx', Method='General',
                          ExperimentalDesign='longitudinal', pvalue='pAdjusted', degrees_of_freedom=dof)
    q = res.quant_data
    called = q['pAdjusted'] <= 0.05
    truth = q['simulated_pattern'].isin(INTERACTION)
    fdr = (called & ~truth).sum() / max(called.sum(), 1)
    monotonic = q['simulated_pattern'].isin({'treatment_up', 'treatment_down', 'treatment_late'})
    assert fdr <= 0.08, f'false discovery rate {fdr:.2%}'
    assert called[monotonic].mean() >= 0.95
    if dof == 3:  # a transient peak needs the extra spline degree of freedom
        assert called[q['simulated_pattern'] == 'treatment_peak'].mean() >= 0.85


@pytest.mark.parametrize('design', ['independent', 'related'])
def test_moderated_spline_test_and_time_effect(design):
    res = omics.OmicScope(Table=f'tests/data/proteins/longitudinal_{design}.xlsx', Method='General',
                          ExperimentalDesign='longitudinal')          # defaults: moderated, df=2
    q = res.quant_data
    peak = q['simulated_pattern'] == 'treatment_peak'
    assert (q.loc[peak, 'pAdjusted'] <= 0.05).mean() >= 0.9       # moderation recovers the transient peak
    time_only = q['simulated_pattern'].isin({'time_up', 'time_down'})
    no_time = q['simulated_pattern'].isin({'flat', 'condition_offset'})
    assert (q.loc[time_only, 'pAdjusted_time'] <= 0.05).mean() >= 0.95
    called = q['pAdjusted_time'] <= 0.05
    assert (called & no_time).sum() / called.sum() <= 0.06
