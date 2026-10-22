""" Sample-level quality control of an OmicScope object: outlying samples and batch effects.

sample_outliers: flags samples that deviate from the other replicates of their condition, by a robust
    distance in PCA space and by their correlation with the other replicates.
outlier_plot: both criteria per sample, with the cutoffs and flagged samples labelled.
remove_batch_effect: limma-style removal of batch effects from a log-scale matrix (y ~ condition + batch,
    batch terms subtracted), used by the statistics when OmicScope is run with batch_correction=True.
batch_effect: share of the variance of each protein explained by batch and by condition, before and
    after the correction.
batch_plot: PCA of the samples before and after the correction (colour = condition, marker = batch).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from .Exploration import _condition_of, _log_matrix, _sample_matrix


def _robust_scale(x):
    """Median absolute deviation scaled to the standard deviation of a normal distribution."""
    return stats.median_abs_deviation(x, scale='normal', nan_policy='omit')


def sample_outliers(self, pvalue=1.0, alpha=0.01):
    """Flag samples that deviate from the other replicates of their condition.

    Two criteria, computed on log2 abundances of the proteins quantified in every sample:

    * robust distance: each protein is centred on the median of its condition (so differences between
      conditions are not mistaken for outliers) and scaled to unit variance; the samples are projected on
      the first k principal components (80% of the variance, at most n_samples - 2). Residuals are divided by
      a robust noise level of their condition first, so a noisier condition is not flagged as a whole. With enough samples
      (n >= 5 (k + 1)) the distance is a robust Mahalanobis distance (minimum covariance determinant);
      otherwise it is the distance to the centre rescaled so that its median matches the median of a
      chi-square with k degrees of freedom (per-component MAD estimates are too unstable with few samples).
      A sample beyond the chi-square (k) quantile 1 - alpha is flagged.
    * correlation: median Pearson correlation of the sample with the other samples of its condition.
      Its Fisher z, centred on the median of the condition and divided by the MAD of all samples (at least the
      sampling error of a Fisher z, 1 / sqrt(n_proteins - 3)), is the `correlation_z`; a sample below the
      normal quantile alpha is flagged.

    With three or fewer replicates per condition an outlier also pulls the condition median, so only gross
    outliers can be detected.

    Args:
        pvalue (float, optional): use proteins with self.pvalue <= pvalue. Defaults to 1.0 (all).
        alpha (float, optional): significance level of both cutoffs. Defaults to 0.01.

    Returns:
        pd.DataFrame: [Sample, Condition, robust_distance, distance_cutoff, median_correlation,
        correlation_z, correlation_cutoff, correlation_limit (correlation at the cutoff), outlier], one row
        per sample.
    """
    from sklearn.covariance import MinCovDet
    from sklearn.decomposition import PCA

    qd = self.quant_data
    mat, _ = _sample_matrix(self, qd[qd[self.pvalue] <= pvalue])
    logm = _log_matrix(self, mat).dropna()
    if len(logm) < 10:
        raise ValueError('At least 10 proteins quantified in every sample are needed to look for outliers.')
    cond = _condition_of(self)
    labels = np.array([cond[s] for s in logm.columns])
    n = len(labels)
    if n < 4:
        raise ValueError('At least 4 samples are needed to look for outliers.')
    values = logm.to_numpy(float)

    # robust distance on condition-centred data
    resid = values.copy()
    for c in np.unique(labels):
        block = values[:, labels == c] - np.median(values[:, labels == c], axis=1, keepdims=True)
        resid[:, labels == c] = block / max(_robust_scale(block.ravel()), 1e-12)   # conditions may differ in noise
    sd = resid.std(axis=1, ddof=1)
    resid = resid[sd > 0] / sd[sd > 0, None]
    pca = PCA().fit(resid.T)
    k = int(np.searchsorted(np.cumsum(pca.explained_variance_ratio_), 0.8) + 1)
    k = max(1, min(k, n - 2))
    scores = pca.transform(resid.T)[:, :k]
    if n >= 5 * (k + 1):
        d2 = MinCovDet(random_state=0).fit(scores).mahalanobis(scores)
    else:
        d2 = (scores ** 2).sum(axis=1)
        d2 = d2 * stats.chi2.median(k) / np.median(d2)
    distance_cutoff = np.sqrt(stats.chi2.ppf(1 - alpha, k))

    # correlation with the other replicates of the condition
    corr = np.corrcoef(values.T)
    median_corr = np.array([np.median(corr[i, (labels == labels[i]) & (np.arange(n) != i)])
                            if (labels == labels[i]).sum() > 1 else np.nan for i in range(n)])
    fisher = np.arctanh(np.clip(median_corr, -0.999999, 0.999999))
    centre = pd.Series(fisher).groupby(labels).transform('median').to_numpy()
    # floor: sampling error of a Fisher z correlation, so trivial differences between similar samples are not flagged
    scale = max(_robust_scale(fisher - centre), 1 / np.sqrt(len(values) - 3))
    corr_z = (fisher - centre) / scale
    corr_cutoff = stats.norm.ppf(alpha)

    out = pd.DataFrame({'Sample': logm.columns, 'Condition': labels, 'robust_distance': np.sqrt(d2),
                        'distance_cutoff': distance_cutoff, 'median_correlation': median_corr,
                        'correlation_z': corr_z, 'correlation_cutoff': corr_cutoff,
                        'correlation_limit': np.tanh(centre + corr_cutoff * scale)})
    out['outlier'] = (out['robust_distance'] > distance_cutoff) | (out['correlation_z'] < corr_cutoff)
    self.Params['Params']['Quality_Outliers'] = (
        f'Sample outliers: robust PCA distance ({k} PCs, chi-square cutoff) and within-condition correlation '
        f'(robust z), alpha = {alpha}; flagged: ' + (', '.join(out.loc[out['outlier'], 'Sample']) or 'none'))
    return out


def outlier_plot(self, pvalue=1.0, alpha=0.01, palette='Dark2', save=None, dpi=300, vector=True):
    """Robust distance and within-condition correlation of each sample.

    Bars above the dashed line (left) or correlation points below it (right) are samples that deviate
    from their replicates; they are labelled. Check them (sample preparation, injection, normalisation)
    before excluding anything.

    Args:
        pvalue (float, optional): use proteins with self.pvalue <= pvalue. Defaults to 1.0 (all).
        alpha (float, optional): significance level of the cutoffs. Defaults to 0.01.
        palette (str, optional): Palette for conditions. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    out = sample_outliers(self, pvalue=pvalue, alpha=alpha)
    colors = dict(zip(self.Conditions, sns.color_palette(palette, len(self.Conditions)).as_hex()))
    names = out['Sample'].str.split('.').str[0]
    x = np.arange(len(out))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(8, 0.35 * len(out) + 3), 8), sharex=True, dpi=dpi)
    ax1.bar(x, out['robust_distance'], color=[colors[c] for c in out['Condition']])
    ax1.axhline(out['distance_cutoff'].iloc[0], ls='--', color='#444', lw=1)
    ax1.set(ylabel='Robust distance', title='Sample outliers')
    ax2.scatter(x, out['median_correlation'], c=[colors[c] for c in out['Condition']], s=40, zorder=3)
    ax2.vlines(x, out['median_correlation'].min(), out['median_correlation'], color='#bbb', lw=1)
    for c in self.Conditions:
        sel = (out['Condition'] == c).to_numpy()
        if sel.any():
            ax2.hlines(out.loc[sel, 'correlation_limit'].iloc[0], x[sel].min() - 0.4, x[sel].max() + 0.4,
                       ls='--', color='#444', lw=1)
    ax2.set(ylabel='Median correlation\nwithin condition')
    ax2.set_xticks(x, names, rotation=90, fontsize=8)
    for i in np.where(out['outlier'])[0]:
        ax1.annotate(names.iloc[i], (i, out['robust_distance'].iloc[i]), xytext=(0, 3), textcoords='offset points',
                     ha='center', fontsize=8, color='#d95f02')
        ax2.annotate(names.iloc[i], (i, out['median_correlation'].iloc[i]), xytext=(4, -10),
                     textcoords='offset points', fontsize=8, color='#d95f02')
    ax1.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=v) for v in colors.values()], labels=list(colors),
               frameon=False, bbox_to_anchor=(1, 1), loc='upper left')
    sns.despine(fig)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save + 'outlier_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return out


def _design(conditions, batches):
    """Condition columns (intercept + treatment coding) and batch columns (sum-to-zero coding)."""
    cond = pd.Categorical(np.asarray(conditions, dtype=object))
    bat = pd.Categorical(np.asarray(batches, dtype=object))
    if pd.isna(np.asarray(cond, dtype=object)).any() or pd.isna(np.asarray(bat, dtype=object)).any():
        raise ValueError('Conditions and batches must be defined for every sample.')
    xc = np.column_stack([np.ones(len(cond))] + [(cond == c).astype(float) for c in cond.categories[1:]])
    last = bat.categories[-1]
    xb = np.column_stack([(bat == b).astype(float) - (bat == last).astype(float) for b in bat.categories[:-1]]) \
        if len(bat.categories) > 1 else np.empty((len(bat), 0))
    return xc, xb


def _estimable_columns(x):
    """Columns of x kept in order while they add rank (aliased later columns are dropped, as in a pivoted QR)."""
    keep = []
    for j in range(x.shape[1]):
        if np.linalg.matrix_rank(x[:, keep + [j]]) == len(keep) + 1:
            keep.append(j)
    return keep


def remove_batch_effect(log_matrix, conditions, batches):
    """Remove batch effects from a log-scale matrix, as limma::removeBatchEffect.

    Each row is fitted by least squares as y ~ condition + batch, using only its non-missing values, and the
    batch terms are subtracted. Batches are coded sum-to-zero, so the overall level of each row is kept and
    the condition differences are untouched. Rows whose missing values make a batch term inestimable keep
    that term uncorrected. Missing values stay missing.

    The design must contain samples of at least one condition in more than one batch: if every batch holds a
    single condition (batch nested in, or identical to, condition) batch and condition cannot be told apart
    and a ValueError is raised.

    Args:
        log_matrix (pd.DataFrame or np.ndarray): rows = proteins, columns = samples, log scale (NaN = missing).
        conditions (array-like): condition of each column.
        batches (array-like): batch of each column.

    Returns:
        Same type and shape as log_matrix, with the batch effects removed.
    """
    y = np.asarray(log_matrix, dtype=float)
    if y.ndim != 2 or y.shape[1] != len(conditions) or y.shape[1] != len(batches):
        raise ValueError('conditions and batches need one value per column of log_matrix.')
    xc, xb = _design(conditions, batches)
    if xb.shape[1] == 0:
        corrected = y.copy()
    else:
        x = np.hstack([xc, xb])
        if np.linalg.matrix_rank(x) < x.shape[1]:
            raise ValueError('Batch is confounded with condition (every batch contains a single condition, or the '
                             'batches are the conditions): batch effects cannot be separated from the biological '
                             'differences and cannot be removed. Only designs with conditions measured in more than '
                             'one batch can be corrected.')
        observed = ~np.isnan(y)
        patterns, inverse = np.unique(observed, axis=0, return_inverse=True)
        inverse = inverse.ravel()
        nc = xc.shape[1]
        beta_batch = np.zeros((len(y), xb.shape[1]))
        for p, mask in enumerate(patterns):
            rows = np.where(inverse == p)[0]
            if mask.sum() == 0:
                continue
            xp = x[mask]
            keep = list(range(x.shape[1])) if np.linalg.matrix_rank(xp) == x.shape[1] else _estimable_columns(xp)
            coef = np.linalg.lstsq(xp[:, keep], y[np.ix_(rows, np.where(mask)[0])].T, rcond=None)[0]
            for pos, j in enumerate(keep):
                if j >= nc:
                    beta_batch[rows, j - nc] = coef[pos]
        corrected = y - beta_batch @ xb.T
    if isinstance(log_matrix, pd.DataFrame):
        return pd.DataFrame(corrected, index=log_matrix.index, columns=log_matrix.columns)
    return corrected


def _variance_explained(values, conditions, batches):
    """Per-row share of the variance explained by batch (beyond condition) and by condition (beyond batch)."""
    xc, xb = _design(conditions, batches)
    full = np.hstack([xc, xb])
    only_batch = np.hstack([np.ones((len(xb), 1)), xb])

    def rss(x):
        fitted = x @ np.linalg.lstsq(x, values.T, rcond=None)[0]
        return ((values.T - fitted) ** 2).sum(axis=0)

    tss = ((values - values.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    tss = np.where(tss > 0, tss, np.nan)
    rss_full = rss(full)
    return (rss(xc) - rss_full) / tss, (rss(only_batch) - rss_full) / tss


def _batch_inputs(self, batch_column):
    if batch_column not in self.pdata.columns:
        raise ValueError(f"pdata has no '{batch_column}' column (batch of each sample).")
    expression = self.expression.copy().replace(0, np.nan)
    logm = _log_matrix(self, expression.apply(pd.to_numeric, errors='coerce'))
    info = self.pdata.drop_duplicates('Sample').set_index('Sample').reindex(logm.columns)
    if info['Condition'].isna().any():
        raise ValueError('The samples of the expression matrix are not all described in pdata.')
    return logm, info['Condition'].to_numpy(), info[batch_column].to_numpy()


def batch_effect(self, batch_column='Batch'):
    """Variance explained by batch and by condition for each protein, before and after batch correction.

    Uses the matrix that enters the statistics (self.expression, log2-transformed if logTransform), before
    any correction, and the batch of each sample in pdata[batch_column]. R² values are partial: the share of
    the total variance of the protein explained by batch on top of condition, and by condition on top of
    batch, for proteins quantified in every sample. After correction the batch share drops to ~0 while the
    condition share should be kept (or rise, as batch noise is removed).

    Args:
        batch_column (str, optional): pdata column with the batch of each sample. Defaults to 'Batch'.

    Returns:
        dict: 'proteins' (pd.DataFrame [Accession, R2_batch_before, R2_batch_after, R2_condition_before,
        R2_condition_after]), 'summary' (pd.DataFrame, median R² before/after), 'before' and 'after'
        (log-scale matrices, proteins x samples).
    """
    logm, conditions, batches = _batch_inputs(self, batch_column)
    corrected = remove_batch_effect(logm, conditions, batches)
    complete = logm.notna().all(axis=1)
    b_before, c_before = _variance_explained(logm[complete].to_numpy(float), conditions, batches)
    b_after, c_after = _variance_explained(corrected[complete].to_numpy(float), conditions, batches)
    proteins = pd.DataFrame({'Accession': logm.index[complete], 'R2_batch_before': b_before,
                             'R2_batch_after': b_after, 'R2_condition_before': c_before,
                             'R2_condition_after': c_after})
    summary = pd.DataFrame({'batch': [np.nanmedian(b_before), np.nanmedian(b_after)],
                            'condition': [np.nanmedian(c_before), np.nanmedian(c_after)]},
                           index=['before', 'after'])
    return {'proteins': proteins, 'summary': summary, 'before': logm, 'after': corrected}


def batch_plot(self, batch_column='Batch', palette='Dark2', save=None, dpi=300, vector=True):
    """PCA of the samples before and after batch correction.

    Before correction, samples that group by marker (batch) rather than by colour (condition) show a batch
    effect; after correction they should group by condition. Titles give the median share of variance
    explained by batch and by condition.

    Args:
        batch_column (str, optional): pdata column with the batch of each sample. Defaults to 'Batch'.
        palette (str, optional): Palette for conditions. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    from sklearn.decomposition import PCA

    result = batch_effect(self, batch_column=batch_column)
    _, conditions, batches = _batch_inputs(self, batch_column)
    colors = dict(zip(self.Conditions, sns.color_palette(palette, len(self.Conditions)).as_hex()))
    markers = dict(zip(pd.unique(batches), ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'h', '<', '>']))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=dpi)
    for ax, key in zip(axes, ['before', 'after']):
        data = result[key].dropna()
        data = data[data.std(axis=1) > 0]
        z = (data.sub(data.mean(axis=1), axis=0)).to_numpy(float).T
        pca = PCA(n_components=2).fit(z)
        xy = pca.transform(z)
        for b in markers:
            for c in colors:
                sel = (batches == b) & (conditions == c)
                ax.scatter(xy[sel, 0], xy[sel, 1], color=colors[c], marker=markers[b], s=60, edgecolor='white')
        ratio = pca.explained_variance_ratio_ * 100
        r2 = result['summary'].loc[key]
        ax.set(xlabel=f'PC1 ({ratio[0]:.1f}%)', ylabel=f'PC2 ({ratio[1]:.1f}%)',
               title=f'{key.capitalize()} correction\nmedian R² batch {r2["batch"]:.2f}, '
                     f'condition {r2["condition"]:.2f}')
    handles = ([plt.Line2D([], [], color=v, marker='o', ls='') for v in colors.values()]
               + [plt.Line2D([], [], color='#666', marker=m, ls='') for m in markers.values()])
    axes[1].legend(handles=handles, labels=list(colors) + [str(b) for b in markers], frameon=False,
                   bbox_to_anchor=(1, 1), loc='upper left')
    sns.despine(fig)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save + 'batch_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return result
