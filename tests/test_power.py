"""Power analysis against a Monte-Carlo simulation of t-tests.

The proteins of the simulated dataset have heterogeneous true variances (scaled inverse chi-square,
as assumed by limma's prior); the analytic power averaged over the estimated variances must match
the share of significant t-tests when data are simulated from the true variances.
"""
import contextlib
import io

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy import stats

import omicscope as omics
from omicscope.General import Power

N_PROTEINS, N_PER_GROUP, D0, S0 = 2000, 8, 6, 0.3 ** 2


@pytest.fixture(scope='module')
def sim(tmp_path_factory):
    rng = np.random.default_rng(7)
    sigma = np.sqrt(D0 * S0 / rng.chisquare(D0, N_PROTEINS))
    base = rng.normal(22, 2, N_PROTEINS)
    samples, pdata = {}, []
    for c in ['Control', 'Treated']:
        for b in range(N_PER_GROUP):
            name = f'{c}_{b + 1}'
            samples[name] = np.round(2 ** (base + rng.normal(0, sigma)), 4)
            pdata.append({'Sample': name, 'Condition': c, 'Biological': len(pdata) + 1})
    rdata = pd.DataFrame({'Accession': [f'P{i:05d}' for i in range(N_PROTEINS)],
                          'Description': [f'Protein {i} OS=Homo sapiens GN=G{i} PE=1' for i in range(N_PROTEINS)]})
    path = tmp_path_factory.mktemp('power') / 'power.xlsx'
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame(samples).to_excel(writer, sheet_name='assay', index=False)
        rdata.to_excel(writer, sheet_name='rdata', index=False)
        pd.DataFrame(pdata).to_excel(writer, sheet_name='pdata', index=False)
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table=str(path), Method='General')
    return res, sigma


def _monte_carlo(sigma, fc, n, alpha, reps=20, seed=3):
    rng = np.random.default_rng(seed)
    s = np.repeat(sigma, reps)[:, None]
    a = rng.normal(0, s, (len(s), n))
    b = rng.normal(fc, s, (len(s), n))
    return np.mean(stats.ttest_ind(a, b, axis=1).pvalue <= alpha)


@pytest.mark.parametrize('fc, n, alpha', [(0.5, 3, 0.05), (0.5, 6, 0.05), (1.0, 3, 0.01), (0.3, 10, 0.05)])
def test_power_matches_monte_carlo(sim, fc, n, alpha):
    res, sigma = sim
    pw = Power.power_analysis(res, fold_changes=(fc,), n_range=[n], alpha=alpha)
    assert abs(pw['power'].iloc[0] - _monte_carlo(sigma, fc, n, alpha)) <= 0.04


def test_power_increases_with_n_and_fold_change(sim):
    pw = Power.power_analysis(sim[0])
    assert pw['power'].between(0, 1).all() and pw['power'].notna().all()
    wide = pw.pivot(index='n_per_group', columns='log2FC', values='power')
    assert (wide.diff().dropna() >= -1e-9).all().all()                    # with n
    assert (wide.T.diff().dropna() >= -1e-9).all().all()                  # with fold change
    assert wide.loc[2, 0.5] < 0.3 and wide.loc[20, 1.5] > 0.99


def test_fdr_threshold_is_stricter(sim):
    plain = Power.power_analysis(sim[0], fold_changes=(1.0,), n_range=[4, 8])
    fdr = Power.power_analysis(sim[0], fold_changes=(1.0,), n_range=[4, 8], fdr_share=0.1)
    assert (fdr['alpha'] < 0.05).all() and (fdr['power'].to_numpy() < plain['power'].to_numpy()).all()


def test_power_plot_runs(sim, tmp_path):
    out = Power.power_plot(sim[0], save=str(tmp_path) + '/')
    assert (tmp_path / 'power_plot.svg').exists() and len(out) == 3 * 19
    plt.close('all')
