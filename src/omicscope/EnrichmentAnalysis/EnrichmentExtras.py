""" Extra EnrichmentScope figures.

- gsea_plot: the classic GSEA running enrichment score of one term along the ranked protein list.
- term_foldchange: fold-change distribution of the proteins behind each top term (ORA or GSEA), showing
  whether a term is driven by consistent up- or down-regulation or by a mix of both.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def gsea_plot(self, Term, Gene_set=None, color='#1b9e77', save=None, dpi=300, vector=True):
    """Running enrichment score of a GSEA term.

    Top: running score along the ranked list (peak = enrichment score). Middle: positions of the term's
    proteins in the list. Bottom: ranking metric. Leading-edge proteins are those before the peak
    (positive NES) or after it (negative NES).

    Args:
        Term (str): term name, as in self.results['Term'].
        Gene_set (str, optional): library, needed only if the term appears in more than one.
        color (str, optional): line colour. Defaults to '#1b9e77'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    if self.Analysis != 'GSEA':
        raise ValueError('gsea_plot is only available for GSEA results')
    keys = [k for k in self.curves if k[1] == Term and (Gene_set is None or k[0] == Gene_set)]
    if not keys:
        raise KeyError(f'{Term!r} is not among the significant GSEA terms')
    curve = self.curves[keys[0]]
    row = self.results[(self.results['Term'] == Term) & (self.results['Gene_set'] == keys[0][0])].iloc[0]
    rank = self.ranking.to_numpy()
    x = np.arange(len(rank))
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 5.5), dpi=dpi, sharex=True,
                                        gridspec_kw={'height_ratios': [3, 0.6, 1.4], 'hspace': 0.05})
    ax1.plot(x, curve['RES'], color=color, lw=2)
    ax1.axhline(0, color='#888', lw=0.8)
    peak = int(np.argmax(np.abs(curve['RES'])))
    ax1.axvline(peak, color='#d95f02', ls='--', lw=1)
    ax1.set(ylabel='Running enrichment score')
    # FDR = 0 means no permutation beat the observed score: report the resolution instead
    fdr = f"< {10 ** -row['-log10(pAdj)']:.2g}" if row['Adjusted P-value'] == 0 else f"{row['Adjusted P-value']:.3g}"
    ax1.set_title(f"{Term}\nNES = {row['NES']:.2f}   FDR {'' if fdr[0] == '<' else '= '}{fdr}", loc='left', fontsize=11)
    ax2.vlines(curve['hits'], 0, 1, color='#222', lw=0.6)
    ax2.set(yticks=[])
    ax3.fill_between(x, rank, color='#999', lw=0)
    ax3.axhline(0, color='#555', lw=0.6)
    ax3.set(ylabel='Ranking metric', xlabel=f'Rank in protein list ({len(rank)} proteins)')
    for a in (ax1, ax2, ax3):
        sns.despine(ax=a)
    if save is not None:
        plt.savefig(save + 'gsea_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return fig


def term_foldchange_table(self, top=15, Gene_set=None):
    """log2(fc) of every protein behind the top terms (ordered by adjusted p-value).

    Returns:
        pd.DataFrame: [Term, Gene_set, gene_name, log2(fc), Adjusted P-value].
    """
    res = self.results if Gene_set is None else self.results[self.results['Gene_set'] == Gene_set]
    res = res.sort_values('Adjusted P-value').head(top)
    qd = self.OmicScope.quant_data
    fc = dict(zip(qd['gene_name'].astype(str).str.upper(), pd.to_numeric(qd['log2(fc)'], errors='coerce')))
    rows = [(r['Term'], r['Gene_set'], g, fc[str(g).upper()], r['Adjusted P-value'])
            for _, r in res.iterrows() for g in r['Genes'] if str(g).upper() in fc]
    return pd.DataFrame(rows, columns=['Term', 'Gene_set', 'gene_name', 'log2(fc)', 'Adjusted P-value'])


def term_foldchange(self, top=15, Gene_set=None, palette='RdBu_r', save=None, dpi=300, vector=True):
    """Fold-change distribution of the proteins behind each top term.

    Each point is a protein (coloured by its log2 fold change) and the box summarises the term. Terms
    sitting entirely on one side of zero are consistently up- or down-regulated.

    Args:
        top (int, optional): number of terms. Defaults to 15.
        Gene_set (str, optional): restrict to one library.
        palette (str, optional): diverging palette for fold change. Defaults to 'RdBu_r'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    df = term_foldchange_table(self, top=top, Gene_set=Gene_set)
    if df.empty:
        raise ValueError('No proteins found for the selected terms')
    order = df.groupby('Term', sort=False)['log2(fc)'].median().sort_values().index
    lim = float(np.nanpercentile(np.abs(df['log2(fc)']), 98)) or 1
    fig, ax = plt.subplots(figsize=(8, 0.42 * len(order) + 1.5), dpi=dpi)
    sns.boxplot(data=df, y='Term', x='log2(fc)', order=order, color='white', showfliers=False, width=0.6,
                linewidth=1, ax=ax)
    jitter = np.random.default_rng(0).uniform(-0.18, 0.18, len(df))
    pts = ax.scatter(df['log2(fc)'], df['Term'].map({t: i for i, t in enumerate(order)}) + jitter,
                     c=df['log2(fc)'], cmap=palette, vmin=-lim, vmax=lim, s=16, edgecolor='#333', linewidth=0.3,
                     zorder=3)
    ax.axvline(0, color='#888', lw=0.8, zorder=0)
    ax.set(ylabel='', xlabel='log2(fold change)', title='Fold change of the proteins in each term')
    fig.colorbar(pts, ax=ax, shrink=0.6, label='log2(fc)')
    sns.despine()
    if save is not None:
        plt.savefig(save + 'term_foldchange.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return df
