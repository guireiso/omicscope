"""The vectorised statistics must agree with the reference implementations they replace."""
import numpy as np
import pytest
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.oneway import anova_oneway

from omicscope.General.Stats import fast_stats as fs

rng = np.random.default_rng(7)
LEVELS = ['A', 'B', 'C', 'D']
LABELS = np.repeat(LEVELS, [5, 5, 6, 4])           # unbalanced on purpose
X = rng.normal(0, 1, size=(60, LABELS.size)) + np.repeat(rng.normal(0, 1, (60, 4)), [5, 5, 6, 4], axis=1)
X[:, LABELS == 'D'] *= 2.5                         # unequal variances
XN = X.copy()
XN[rng.random(X.shape) < 0.05] = np.nan            # missing values


@pytest.mark.parametrize('equal_var', [True, False])
def test_ttest_matches_scipy(equal_var):
    a, b = XN[:, LABELS == 'A'], XN[:, LABELS == 'C']
    t, p, _ = fs.ttest(a, b, equal_var=equal_var)
    ref = stats.ttest_ind(a, b, axis=1, equal_var=equal_var, nan_policy='omit')
    np.testing.assert_allclose(t, ref.statistic, rtol=1e-10)
    np.testing.assert_allclose(p, ref.pvalue, rtol=1e-8)


def test_anova_matches_scipy():
    _, p, _, _ = fs.anova(XN, LABELS, LEVELS)
    ref = [stats.f_oneway(*[row[LABELS == g][~np.isnan(row[LABELS == g])] for g in LEVELS]).pvalue for row in XN]
    np.testing.assert_allclose(p, ref, rtol=1e-8)


def test_welch_anova_matches_statsmodels():
    _, p, _, _ = fs.anova(X, LABELS, LEVELS, welch=True)
    ref = [anova_oneway(row, LABELS, use_var='unequal').pvalue for row in X]
    np.testing.assert_allclose(p, ref, rtol=1e-8)


def test_tukey_matches_statsmodels():
    pairs, p = fs.pairwise_posthoc(X, LABELS, LEVELS, method='tukey')
    for r in range(15):
        res = pairwise_tukeyhsd(X[r], LABELS)
        np.testing.assert_allclose(p[r], res.pvalues, atol=2e-6)
    assert pairs == [('A', 'B'), ('A', 'C'), ('A', 'D'), ('B', 'C'), ('B', 'D'), ('C', 'D')]


def test_games_howell_matches_exact_distribution():
    pairs, p = fs.pairwise_posthoc(X[:10], LABELS, LEVELS, method='games-howell')
    for r in range(10):
        for c, (a, b) in enumerate(pairs):
            xa, xb = X[r, LABELS == a], X[r, LABELS == b]
            va, vb = xa.var(ddof=1) / xa.size, xb.var(ddof=1) / xb.size
            q = abs(xa.mean() - xb.mean()) / np.sqrt((va + vb) / 2)
            df = (va + vb) ** 2 / (va ** 2 / (xa.size - 1) + vb ** 2 / (xb.size - 1))
            assert p[r, c] == pytest.approx(stats.studentized_range.sf(q, 4, df), abs=2e-6)


def test_squeeze_var_recovers_prior():
    d0, s0, df = 6.0, 0.4, 4
    sigma2 = s0 * d0 / rng.chisquare(d0, 20000)            # scaled inverse chi-square prior
    s2 = sigma2 * rng.chisquare(df, 20000) / df
    post, d0_hat, s0_hat = fs.squeeze_var(s2, df)
    assert d0_hat == pytest.approx(d0, rel=0.15)
    assert s0_hat == pytest.approx(s0, rel=0.1)
    assert np.all((post - np.minimum(s2, s0_hat) >= -1e-12) & (post - np.maximum(s2, s0_hat) <= 1e-12))


def test_qvalue_estimates_pi0_and_is_not_more_conservative_than_bh():
    p = np.r_[rng.uniform(size=8000), rng.beta(0.2, 20, 2000)]   # 80% nulls
    q, pi0 = fs.qvalue(p)
    assert pi0 == pytest.approx(0.8, abs=0.05)
    assert np.all(q <= multipletests(p, method='fdr_bh')[1] + 1e-12)


def _null_pair_rate(welch):
    import omicscope as omics
    q = omics.OmicScope(Table='tests/data/proteins/static_4groups.xlsx', Method='General', pvalue='pTukey',
                        welch=welch).quant_data
    q = q[q['Comparison'].apply(len) == 2]
    shift = lambda r, c: 0 if c == 'A' else r[f'simulated_shift_{c}']
    same = q.apply(lambda r: shift(r, r.Comparison[0]) == shift(r, r.Comparison[1]), axis=1)
    return q, (q.loc[same, 'pTukey'] <= 0.05).mean()


def test_posthoc_is_a_real_pairwise_test():
    """The old code copied the ANOVA pAdjusted into pTukey whenever an entity had a missing value."""
    q, _ = _null_pair_rate(False)
    assert (q['pTukey'] != q['pAdjusted']).mean() > 0.95


def test_games_howell_protects_unequal_variances_better_than_tukey():
    _, tukey = _null_pair_rate(False)
    _, games_howell = _null_pair_rate(True)
    assert games_howell < tukey
