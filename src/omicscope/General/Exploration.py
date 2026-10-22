""" Data-quality and exploratory views of an OmicScope object.

missing_values: how many proteins each sample quantified and whether missing values depend on abundance
    (values missing because they fall below the detection limit concentrate on low-abundance proteins).
cv_plot: coefficient of variation of each protein within each condition, a standard measure of
    quantitative precision.
protein_map: t-SNE (or PCA) map of protein abundance profiles; proteins that change together across
    samples end up close together, coloured by k-means trend cluster.

Each figure has a matching function that returns the numbers behind it (missing_summary, cv_table,
protein_embedding).
"""
import copy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def _sample_matrix(self, data=None):
    """Sample columns (Sample.Condition) of quant_data, indexed by Accession, with zeros as missing."""
    df = self.whole_data if data is None else data
    df = df.drop_duplicates('Accession').set_index('Accession')
    samples = [s for s in self.pdata['Sample'] if s in df.columns]
    mat = df[samples].apply(pd.to_numeric, errors='coerce').replace(0, np.nan)
    return mat, df['gene_name'].astype(str)


def _log_matrix(self, mat):
    return np.log2(mat) if getattr(self, 'logTransform', True) else mat


def _condition_of(self):
    return dict(zip(self.pdata['Sample'], self.pdata['Condition']))


def missing_summary(self):
    """Per-sample counts of quantified proteins and per-protein missing fraction against mean abundance.

    Returns:
        (pd.DataFrame, pd.DataFrame): samples [Sample, Condition, quantified, missing] and proteins
        [Accession, gene_name, mean_log2, missing_fraction].
    """
    mat, genes = _sample_matrix(self)
    cond = _condition_of(self)
    per_sample = pd.DataFrame({'Sample': mat.columns, 'Condition': [cond[s] for s in mat.columns],
                               'quantified': mat.notna().sum().to_numpy(), 'missing': mat.isna().sum().to_numpy()})
    logm = _log_matrix(self, mat)
    per_protein = pd.DataFrame({'Accession': mat.index, 'gene_name': genes.reindex(mat.index).to_numpy(),
                                'mean_log2': logm.mean(axis=1).to_numpy(),
                                'missing_fraction': mat.isna().mean(axis=1).to_numpy()})
    return per_sample, per_protein.dropna(subset=['mean_log2'])


def missing_values(self, palette='Dark2', save=None, dpi=300, vector=True):
    """Quantified proteins per sample and missing values against abundance.

    If missingness rises steeply towards low abundance, values are mostly missing not at random
    (below the detection limit): left-censored imputation (MinProb) or the presence/absence layer
    suit these data better than mean/median imputation.

    Args:
        palette (str, optional): Palette for conditions. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    per_sample, per_protein = missing_summary(self)
    colors = dict(zip(self.Conditions, sns.color_palette(palette, len(self.Conditions)).as_hex()))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={'width_ratios': [1.4, 1]}, dpi=dpi)
    ax1.bar(range(len(per_sample)), per_sample['quantified'], color=[colors[c] for c in per_sample['Condition']])
    ax1.set_xticks(range(len(per_sample)), per_sample['Sample'].str.split('.').str[0], rotation=90, fontsize=8)
    ax1.set(ylabel='Quantified proteins', title='Proteins quantified per sample')
    ax1.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=v) for v in colors.values()], labels=list(colors),
               frameon=False, loc='lower right')
    bins = pd.qcut(per_protein['mean_log2'], q=min(20, max(2, len(per_protein) // 20)), duplicates='drop')
    trend = per_protein.groupby(bins, observed=True).agg(x=('mean_log2', 'median'), y=('missing_fraction', 'mean'))
    ax2.scatter(per_protein['mean_log2'], per_protein['missing_fraction'] * 100, s=6, alpha=0.25, color='#555')
    ax2.plot(trend['x'], trend['y'] * 100, color='#d95f02', lw=2.5, label='binned mean')
    ax2.set(xlabel='Mean abundance (log2)', ylabel='Missing values (%)', title='Missing values vs. abundance')
    ax2.legend(frameon=False)
    sns.despine(fig)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save + 'missing_values.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return per_sample, per_protein


def cv_table(self):
    """Coefficient of variation (%) of each protein in each condition, on the linear scale.

    Returns:
        pd.DataFrame: [Accession, gene_name, Condition, CV] for proteins quantified in >= 2 samples.
    """
    mat, genes = _sample_matrix(self, self.quant_data)
    if not getattr(self, 'logTransform', True):
        mat = 2 ** mat
    cond = _condition_of(self)
    out = []
    for c in self.Conditions:
        sub = mat[[s for s in mat.columns if cond[s] == c]]
        ok = sub.notna().sum(axis=1) >= 2
        cv = (sub[ok].std(axis=1, ddof=1) / sub[ok].mean(axis=1) * 100)
        out.append(pd.DataFrame({'Accession': cv.index, 'gene_name': genes.reindex(cv.index).to_numpy(),
                                 'Condition': c, 'CV': cv.to_numpy()}))
    return pd.concat(out, ignore_index=True).dropna()


def cv_plot(self, palette='Dark2', save=None, dpi=300, vector=True):
    """Distribution of protein coefficients of variation per condition.

    Lower CVs mean more precise quantification. For label-free proteomics, median CVs around 10-20%
    between biological replicates are typical; a condition with much higher CVs points to a noisy sample.

    Args:
        palette (str, optional): Palette for conditions. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    cv = cv_table(self)
    fig, ax = plt.subplots(figsize=(1.6 * len(self.Conditions) + 2, 5), dpi=dpi)
    sns.violinplot(data=cv, x='Condition', y='CV', hue='Condition', order=self.Conditions, palette=palette,
                   cut=0, inner='quartile', legend=False, ax=ax)
    upper = cv['CV'].quantile(0.99)
    ax.set_ylim(0, upper * 1.08)
    for i, c in enumerate(self.Conditions):
        med = cv.loc[cv['Condition'] == c, 'CV'].median()
        ax.text(i, upper * 1.02, f'median {med:.1f}%', ha='center', fontsize=10)
    ax.set(xlabel='', ylabel='Coefficient of variation (%)', title='Quantitative precision per condition')
    sns.despine()
    if save is not None:
        plt.savefig(save + 'cv_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return cv


def protein_embedding(self, pvalue=0.05, method='tsne', perplexity=30, k_cluster=None, random_state=0):
    """2-D map of protein abundance profiles (rows z-scored across samples).

    Args:
        pvalue (float, optional): keep proteins with self.pvalue <= pvalue. Defaults to 0.05.
        method (str, optional): 'tsne' or 'pca'. Defaults to 'tsne'.
        perplexity (float, optional): t-SNE perplexity, capped at (n - 1) / 3. Defaults to 30.
        k_cluster (int, optional): k-means clusters; None picks k with the knee of the inertia curve.
        random_state (int, optional): seed for t-SNE and k-means. Defaults to 0.

    Returns:
        pd.DataFrame: [Accession, gene_name, x, y, cluster, log2(fc), <pvalue>].
    """
    from scipy.stats import zscore
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE

    qd = self.quant_data
    sel = qd[qd[self.pvalue] <= pvalue]
    mat, genes = _sample_matrix(self, sel)
    logm = _log_matrix(self, mat)
    logm = logm.apply(lambda r: r.fillna(r.min()), axis=1)   # missing = below detection: lowest value of the protein
    logm = logm[logm.std(axis=1) > 0].dropna()
    if len(logm) < 10:
        raise ValueError('At least 10 proteins are needed for the protein map; relax the p-value cutoff.')
    z = zscore(logm.to_numpy(float), axis=1)
    if method == 'tsne':
        xy = TSNE(n_components=2, perplexity=min(perplexity, (len(z) - 1) / 3), init='pca',
                  learning_rate='auto', random_state=random_state).fit_transform(z)
    elif method == 'pca':
        xy = PCA(n_components=2).fit_transform(z)
    else:
        raise ValueError("method must be 'tsne' or 'pca'")
    if k_cluster is None:
        inertia = [KMeans(n_clusters=k, n_init=4, random_state=random_state).fit(z).inertia_
                   for k in range(1, min(10, len(z)) + 1)]
        from kneed import KneeLocator
        k_cluster = KneeLocator(range(1, len(inertia) + 1), inertia, curve='convex', direction='decreasing').knee or 3
    labels = KMeans(n_clusters=int(k_cluster), n_init=10, random_state=random_state).fit_predict(z)
    info = sel.drop_duplicates('Accession').set_index('Accession')
    return pd.DataFrame({'Accession': logm.index, 'gene_name': genes.reindex(logm.index).to_numpy(),
                         'x': xy[:, 0], 'y': xy[:, 1], 'cluster': labels + 1,
                         'log2(fc)': info['log2(fc)'].reindex(logm.index).to_numpy(),
                         self.pvalue: info[self.pvalue].reindex(logm.index).to_numpy()})


def protein_map(self, *Proteins, pvalue=0.05, method='tsne', perplexity=30, k_cluster=None, palette='Dark2',
                save=None, dpi=300, vector=True):
    """Map of co-regulated proteins.

    Each point is a protein; proteins whose abundance rises and falls together across samples lie close
    together. Colours are k-means clusters of the same profiles, so a cluster forming a compact island
    is a well-separated expression pattern (often a complex or pathway).

    Args:
        *Proteins (str): gene names to label.
        pvalue (float, optional): p-value cutoff. Defaults to 0.05.
        method (str, optional): 'tsne' or 'pca'. Defaults to 'tsne'.
        perplexity (float, optional): t-SNE perplexity. Defaults to 30.
        k_cluster (int, optional): number of k-means clusters (None: knee method).
        palette (str, optional): Palette for clusters. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    emb = protein_embedding(self, pvalue=pvalue, method=method, perplexity=perplexity, k_cluster=k_cluster)
    fig, ax = plt.subplots(figsize=(7.5, 6.5), dpi=dpi)
    sns.scatterplot(data=emb, x='x', y='y', hue='cluster', palette=palette, s=18, linewidth=0, alpha=0.85, ax=ax)
    for _, r in emb[emb['gene_name'].isin(Proteins)].iterrows():
        ax.annotate(r['gene_name'], (r['x'], r['y']), xytext=(4, 4), textcoords='offset points', fontsize=9)
    name = 't-SNE' if method == 'tsne' else 'PC'
    ax.set(xlabel=f'{name} 1', ylabel=f'{name} 2', xticks=[], yticks=[],
           title=f'Protein co-expression map ({len(emb)} proteins, {self.pvalue} <= {pvalue})')
    ax.legend(title='cluster', frameon=False, bbox_to_anchor=(1, 1), loc='upper left')
    sns.despine()
    if save is not None:
        plt.savefig(save + 'protein_map.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return emb
