""" Effect sizes of a two-condition comparison.

Statistical significance says whether a change is detectable; the effect size says how large it is.
For two conditions OmicScope reports, per protein and on the tested (log2) scale, the difference of
means with its 95% confidence interval (same standard error and degrees of freedom as the test, so
moderated, Welch or paired as applicable) and Hedges' g (difference in pooled standard deviations,
small-sample corrected). With more than two conditions, eta squared gives the share of a protein's
variance explained by the conditions.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def effect_table(self, *Proteins, top=20):
    """Effect sizes of the selected proteins (or the `top` most significant ones).

    Returns:
        pd.DataFrame: gene_name, Accession, Mean difference, CI95 low, CI95 high, Hedges' g and p-values.
    """
    qd = self.quant_data
    if 'Mean difference' not in qd.columns:
        raise ValueError('Effect sizes with confidence intervals are available for two-condition static analyses')
    cols = ['gene_name', 'Accession', 'Mean difference', 'CI95 low', 'CI95 high', "Hedges' g", 'pvalue', self.pvalue]
    df = qd[qd['Mean difference'].notna()]
    df = df[df['gene_name'].isin(Proteins)] if Proteins else df.sort_values(self.pvalue).head(top)
    return df[list(dict.fromkeys(cols))].reset_index(drop=True)


def effect_plot(self, *Proteins, top=20, color='#1b9e77', save=None, dpi=300, vector=True):
    """Forest plot: difference of log2 means with its 95% confidence interval.

    Intervals that do not cross zero are significant at the nominal 5% level; their width shows how
    precisely each change is estimated.

    Args:
        *Proteins (str): gene names to show. If none, the `top` most significant proteins are shown.
        top (int, optional): number of proteins when none are given. Defaults to 20.
        color (str, optional): interval colour. Defaults to '#1b9e77'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    df = effect_table(self, *Proteins, top=top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(df) + 1.5), dpi=dpi)
    y = np.arange(len(df))
    ax.hlines(y, df['CI95 low'], df['CI95 high'], color=color, lw=2)
    ax.scatter(df['Mean difference'], y, color=color, s=30, zorder=3, edgecolor='#222', linewidth=0.5)
    ax.axvline(0, color='#888', lw=0.8, ls='--')
    ax.set_yticks(y, df['gene_name'].astype(str))
    ax.set(xlabel=f'log2 difference ({self.experimental[0]} - {self.ctrl}) with 95% CI', ylabel='',
           title='Effect sizes')
    sns.despine()
    if save is not None:
        plt.savefig(save + 'effect_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return df.iloc[::-1].reset_index(drop=True)
