"""Vectorised statistics: every function works on all entities (rows) at once.

Missing values (NaN) are dropped per entity and per group, as the per-row scipy/statsmodels calls did.

* t-tests (Student / Welch) and one-way ANOVA (classic / Welch)
* post-hoc tests: Tukey-Kramer HSD (equal variances) and Games-Howell (unequal variances)
* empirical-Bayes variance moderation (limma's squeezeVar, Smyth 2004) for t and F statistics
* Storey's q-values (Storey & Tibshirani 2003)

The studentized range distribution is computed by quadrature (the approach of R's ptukey), vectorised
over q and the degrees of freedom: about 500x faster than scipy.stats.studentized_range, absolute error
below 2e-6.
"""
from functools import lru_cache
from itertools import combinations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import special, stats
from scipy.interpolate import CubicSpline

# ---------------------------------------------------------------------------- group summaries


def group_summary(X, labels, levels):
    """Per-row n, mean and variance (ddof=1) of each group. X: (rows, samples); labels: per column."""
    X = np.asarray(X, float)
    labels = np.asarray(labels)
    n = np.empty((X.shape[0], len(levels)))
    mean = np.full_like(n, np.nan)
    var = np.full_like(n, np.nan)
    for j, level in enumerate(levels):
        block = X[:, labels == level]
        n[:, j] = np.sum(~np.isnan(block), axis=1)
        with np.errstate(invalid='ignore', divide='ignore'):
            mean[:, j] = np.nanmean(block, axis=1) if block.shape[1] else np.nan
            var[:, j] = np.nanvar(block, axis=1, ddof=1) if block.shape[1] > 1 else np.nan
    var[n < 2] = np.nan
    mean[n < 1] = np.nan
    return n, mean, var


# ---------------------------------------------------------------------------- empirical Bayes


def _trigamma_inverse(x):
    """Inverse of the trigamma function (Newton iterations, as in limma)."""
    x = np.asarray(x, float)
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(-dif / y < 1e-8):
            break
    return y


def squeeze_var(s2, df):
    """limma::squeezeVar - shrink residual variances towards a common prior.

    Returns (posterior variances, prior degrees of freedom d0, prior variance s0^2)."""
    s2 = np.asarray(s2, float)
    df = np.broadcast_to(np.asarray(df, float), s2.shape)
    ok = np.isfinite(s2) & (s2 > 0) & (df > 0)
    if ok.sum() < 3:
        return s2, 0.0, np.nan
    z = np.log(s2[ok])
    e = z - special.digamma(df[ok] / 2) + np.log(df[ok] / 2)
    emean = e.mean()
    evar = np.sum((e - emean) ** 2) / (ok.sum() - 1) - np.mean(special.polygamma(1, df[ok] / 2))
    if evar > 0:
        d0 = float(2 * _trigamma_inverse(evar))
        s0 = float(np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2)))
        post = (d0 * s0 + df * s2) / (d0 + df)
    else:                                   # no extra variability: every gene shares the prior
        d0, s0 = np.inf, float(np.exp(emean))
        post = np.full_like(s2, s0)
    post = np.where(np.isfinite(s2) | ~np.isfinite(post), post, s0)
    return post, d0, s0


# ---------------------------------------------------------------------------- tests


def ttest(X1, X2, equal_var=True, moderated=False):
    """Two-sample t-test per row (X1 vs X2). Returns (t, p, log2 difference of means)."""
    n, mean, var = group_summary(np.hstack([X1, X2]), np.r_[np.zeros(X1.shape[1]), np.ones(X2.shape[1])], [0, 1])
    diff = mean[:, 0] - mean[:, 1]
    with np.errstate(invalid='ignore', divide='ignore'):
        if equal_var or moderated:
            df = n[:, 0] + n[:, 1] - 2
            s2 = ((n[:, 0] - 1) * var[:, 0] + (n[:, 1] - 1) * var[:, 1]) / df
            if moderated:
                s2, d0, _ = squeeze_var(s2, df)
                df = df + d0
            se = np.sqrt(s2 * (1 / n[:, 0] + 1 / n[:, 1]))
        else:
            a, b = var[:, 0] / n[:, 0], var[:, 1] / n[:, 1]
            se = np.sqrt(a + b)
            df = (a + b) ** 2 / (a ** 2 / (n[:, 0] - 1) + b ** 2 / (n[:, 1] - 1))
        t = diff / se
        p = 2 * stats.t.sf(np.abs(t), np.minimum(df, 1e6))
    return t, p, diff


def anova(X, labels, levels, welch=False, moderated=False):
    """One-way ANOVA per row. Returns (F, p, residual df, residual variance (MSE))."""
    n, mean, var = group_summary(X, labels, levels)
    present = n > 0
    k = present.sum(axis=1)
    N = n.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        grand = np.nansum(n * mean, axis=1) / N
        ss_between = np.nansum(n * (mean - grand[:, None]) ** 2, axis=1)
        ss_within = np.nansum((n - 1) * var, axis=1)
        df1, df2 = k - 1, N - k
        mse = ss_within / df2
        if welch:
            w = n / var
            w[~present] = 0
            wsum = w.sum(axis=1)
            wmean = np.nansum(w * mean, axis=1) / wsum
            a = np.nansum(w * (mean - wmean[:, None]) ** 2, axis=1) / df1
            lam = np.nansum((1 - w / wsum[:, None]) ** 2 / (n - 1), axis=1)
            F = a / (1 + 2 * (k - 2) / (k ** 2 - 1) * lam)
            df2 = (k ** 2 - 1) / (3 * lam)
        else:
            if moderated:
                mse, d0, _ = squeeze_var(mse, df2)
                df2 = df2 + d0
            F = (ss_between / df1) / mse
        p = stats.f.sf(F, df1, np.minimum(df2, 1e6))
    return F, p, df2, mse


# ---------------------------------------------------------------------------- studentized range


_ZN, _ZW = leggauss(200)
_Z, _ZW = 9 * _ZN, 9 * _ZW


@lru_cache(maxsize=None)
def _range_cdf_spline(k):
    """P(range of k iid standard normals <= w) on a fine grid, as a spline."""
    w = np.concatenate([np.linspace(0, 12, 2401), np.linspace(12.05, 40, 200)])
    inner = np.clip(stats.norm.cdf(_Z) - stats.norm.cdf(_Z - w[:, None]), 0, 1)
    g = k * np.sum(_ZW * stats.norm.pdf(_Z) * inner ** (k - 1), axis=-1)
    return CubicSpline(w, np.clip(g, 0, 1)), w[-1]


@lru_cache(maxsize=4)
def _chi_nodes(nodes):
    un, uw = leggauss(nodes)
    return (un + 1) / 2, uw / 2


def srange_sf(q, k, df, nodes=240):
    """Survival function of the studentized range distribution (k groups, df), vectorised over q and df."""
    q = np.atleast_1d(np.asarray(q, float))
    df = np.broadcast_to(np.asarray(df, float), q.shape)
    u, uw = _chi_nodes(nodes)
    spline, wmax = _range_cdf_spline(int(k))
    out = np.full(q.shape, np.nan)
    good = np.isfinite(q) & np.isfinite(df) & (df > 0)
    for d in np.unique(df[good]):
        m = good & (df == d)
        s = np.sqrt(stats.chi2.ppf(u, d) / d) if d < 1e5 else np.ones_like(u)
        w = np.minimum(q[m][:, None] * s[None, :], wmax)
        out[m] = 1 - np.sum(uw * np.clip(spline(w), 0, 1), axis=-1)
    return np.clip(out, 0, 1)


def pairwise_posthoc(X, labels, levels, method='tukey'):
    """All pairwise comparisons per row.

    method='tukey'        Tukey-Kramer HSD (pooled variance, like statsmodels.pairwise_tukeyhsd)
    method='games-howell' Games-Howell (separate variances, Welch df per pair)

    Returns (pairs, p) with pairs as (level_i, level_j) in sorted order and p shaped (rows, pairs)."""
    levels = sorted(levels)
    n, mean, var = group_summary(X, labels, levels)
    k = (n > 0).sum(axis=1)
    pairs = list(combinations(range(len(levels)), 2))
    p = np.full((X.shape[0], len(pairs)), np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        if method == 'tukey':
            df = n.sum(axis=1) - k
            mse = np.nansum((n - 1) * var, axis=1) / df
        for c, (i, j) in enumerate(pairs):
            diff = np.abs(mean[:, i] - mean[:, j])
            if method == 'tukey':
                q = diff / np.sqrt(mse / 2 * (1 / n[:, i] + 1 / n[:, j]))
                dfp = df
            else:
                a, b = var[:, i] / n[:, i], var[:, j] / n[:, j]
                q = diff / np.sqrt((a + b) / 2)
                dfp = (a + b) ** 2 / (a ** 2 / (n[:, i] - 1) + b ** 2 / (n[:, j] - 1))
            for kk in np.unique(k[np.isfinite(q)]):
                m = (k == kk) & np.isfinite(q) & np.isfinite(dfp)
                if kk >= 2 and m.any():
                    p[m, c] = srange_sf(q[m], int(kk), dfp[m])
    return [(levels[i], levels[j]) for i, j in pairs], p


# ---------------------------------------------------------------------------- multiple testing


def qvalue(p, lambdas=np.arange(0.05, 0.96, 0.05)):
    """Storey's q-values with the smoothed pi0 estimate (qvalue's default 'smoother' method)."""
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    pv = p[ok]
    m = pv.size
    if m == 0:
        return p.copy(), np.nan
    pi0_l = np.array([(pv > lam).mean() / (1 - lam) for lam in lambdas])
    coef = np.polyfit(lambdas, pi0_l, 3)
    pi0 = float(np.clip(np.polyval(coef, lambdas[-1]), 0, 1))
    pi0 = pi0 if pi0 > 0 else float(np.clip(pi0_l.min(), 1e-3, 1))
    order = np.argsort(pv)
    ranked = pi0 * m * pv[order] / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.full(p.shape, np.nan)
    tmp = np.empty(m)
    tmp[order] = np.minimum(q, 1)
    out[ok] = tmp
    return out, pi0
