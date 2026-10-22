""" Sample-size planning from the variability observed in an OmicScope object.

power_analysis: expected power of a two-sample t-test for a given log2 fold change and number of
    samples per group, averaged over the protein-wise variances of the current experiment.
power_plot: power against samples per group, one line per fold change.

The variance of each protein is its pooled within-condition variance of the log2 abundances, shrunk
towards the common prior of all proteins (limma's squeezeVar), which removes most of the sampling noise
of variances estimated from few replicates. Only measured values enter the variances: filling missing
values with a protein's minimum would add artificial spread.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .Exploration import _condition_of, _log_matrix, _sample_matrix
from .Stats.fast_stats import group_summary, squeeze_var


def protein_variances(self):
    """Moderated pooled within-condition variance (log2 scale) of each protein.

    Returns:
        pd.Series: variances indexed by Accession (proteins with >= 1 residual degree of freedom).
    """
    mat, _ = _sample_matrix(self, self.quant_data)
    logm = _log_matrix(self, mat)
    cond = _condition_of(self)
    n, _, var = group_summary(logm.to_numpy(float), [cond[s] for s in logm.columns], self.Conditions)
    dfg = np.where(n >= 2, n - 1, 0)
    df = dfg.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        s2 = np.nansum(np.where(dfg > 0, var, 0) * dfg, axis=1) / df
    ok = (df > 0) & np.isfinite(s2) & (s2 > 0)
    post, _, _ = squeeze_var(s2[ok], df[ok])
    return pd.Series(post, index=logm.index[ok], name='variance')


def _power(delta, sigma, n, alpha):
    """Power of a two-sided two-sample t-test with n per group (noncentral t), for each sigma."""
    from scipy.stats import nct, norm, t

    df = 2 * n - 2
    crit = t.ppf(1 - alpha / 2, df)
    ncp = delta / (sigma * np.sqrt(2 / n))
    with np.errstate(invalid='ignore', over='ignore'):
        power = nct.sf(crit, df, ncp) + nct.cdf(-crit, df, ncp)
    # nct fails far in the tail (power ~ 1): fall back to the normal approximation there
    return np.where(np.isfinite(power), power, norm.sf(crit - ncp) + norm.cdf(-crit - ncp))


def power_analysis(self, fold_changes=(0.5, 1.0, 1.5), n_range=range(2, 21), alpha=0.05, fdr_share=None,
                   max_proteins=2000, random_state=0):
    """Expected power of a two-sample t-test per log2 fold change and samples per group.

    For each protein the power at its own variance is computed with the noncentral t distribution and
    the result is averaged over proteins: the share of proteins changing by log2FC that a test with
    n samples per group would detect.

    `alpha` is the per-test significance threshold. With FDR control (Benjamini-Hochberg) the
    threshold actually applied to each p-value is lower than the FDR level, and depends on how many
    proteins change. Give `fdr_share` (expected share of truly changing proteins, e.g. 0.1) to treat
    `alpha` as the FDR level instead: the per-test threshold is then solved from
    alpha_test = q * share * power / ((1 - share) * (1 - q)) (Jung 2005), iterating with the power.

    Args:
        fold_changes (iterable, optional): log2 fold changes. Defaults to (0.5, 1.0, 1.5).
        n_range (iterable, optional): samples per group. Defaults to range(2, 21).
        alpha (float, optional): per-test significance level, or FDR level when fdr_share is set.
            Defaults to 0.05.
        fdr_share (float, optional): expected share of changing proteins (0-1). Defaults to None.
        max_proteins (int, optional): proteins sampled to average over (for speed). Defaults to 2000.
        random_state (int, optional): seed of that sampling. Defaults to 0.

    Returns:
        pd.DataFrame: [log2FC, n_per_group, power, alpha] (alpha = per-test threshold used).
    """
    var = protein_variances(self).to_numpy()
    if len(var) > max_proteins:
        var = np.random.default_rng(random_state).choice(var, max_proteins, replace=False)
    sigma = np.sqrt(var)
    rows = []
    for fc in fold_changes:
        for n in n_range:
            a = alpha
            power = _power(abs(fc), sigma, n, a).mean()
            if fdr_share is not None:
                for _ in range(100):
                    a_new = alpha * fdr_share * power / ((1 - fdr_share) * (1 - alpha))
                    a_new = float(np.clip(a_new, 1e-12, alpha))
                    power = _power(abs(fc), sigma, n, a_new).mean()
                    if abs(a_new - a) <= 1e-6 * a:
                        break
                    a = a_new
                a = a_new
            rows.append({'log2FC': fc, 'n_per_group': n, 'power': float(power), 'alpha': a})
    return pd.DataFrame(rows)


def power_plot(self, fold_changes=(0.5, 1.0, 1.5), n_range=range(2, 21), alpha=0.05, fdr_share=None,
               palette='Dark2', save=None, dpi=300, vector=True):
    """Power against samples per group, one line per log2 fold change.

    The dashed line marks 80% power, the usual design target; the dotted line marks the samples per
    group of the current experiment. Power is averaged over the protein variances of this experiment,
    so it is the share of proteins changing by that fold change that would be detected.

    Args:
        fold_changes (iterable, optional): log2 fold changes. Defaults to (0.5, 1.0, 1.5).
        n_range (iterable, optional): samples per group. Defaults to range(2, 21).
        alpha (float, optional): per-test significance level (FDR level with fdr_share). Defaults to 0.05.
        fdr_share (float, optional): expected share of changing proteins. Defaults to None.
        palette (str, optional): Palette for fold changes. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    power = power_analysis(self, fold_changes=fold_changes, n_range=n_range, alpha=alpha, fdr_share=fdr_share)
    current = int(self.pdata.groupby('Condition')['Sample'].nunique().min())
    fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
    colors = sns.color_palette(palette, len(fold_changes))
    for fc, color in zip(fold_changes, colors):
        sub = power[power['log2FC'] == fc]
        ax.plot(sub['n_per_group'], sub['power'] * 100, marker='o', ms=3.5, color=color, label=f'log2FC = {fc:g}')
        here = sub[sub['n_per_group'] == current]
        ax.scatter(here['n_per_group'], here['power'] * 100, s=70, color=color, edgecolor='black', zorder=3)
    ax.axhline(80, color='#999', ls='--', lw=1)
    ax.axvline(current, color='#555', ls=':', lw=1)
    ax.text(current, 2, f' current n = {current}', fontsize=9, color='#555')
    level = f'FDR {alpha:g}, {fdr_share:.0%} changing' if fdr_share is not None else f'alpha = {alpha:g}'
    ax.set(xlabel='Samples per group', ylabel='Power (%)', ylim=(0, 102), title=f'Two-sample t-test power ({level})')
    ax.legend(frameon=False, loc='lower right')
    sns.despine()
    if save is not None:
        plt.savefig(save + 'power_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return power
