""" Protein-level meta-analysis across Nebula studies.

Each study gives a fold change and a p-value per protein. They are combined with Stouffer's method on
signed z-scores, z = sign(log2FC) * Phi^-1(1 - p/2), Z = sum(z) / sqrt(k): proteins changing in the same
direction in several studies reach strong combined evidence even if no single study is significant, and
opposite changes cancel out. Each study's |z| is capped (default 5, p ~ 6e-7) so that one study with
extreme p-values cannot carry the combination on its own. Combined p-values are corrected with
Benjamini-Hochberg.

Note: .omics files store the p-value column chosen in OmicScope (often pAdjusted); adjusted p-values
make the combination conservative.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


def meta_analysis(self, min_studies=2, z_cap=5.0):
    """Combine each protein's evidence across studies.

    Args:
        min_studies (int, optional): minimum number of studies quantifying the protein. Defaults to 2.
        z_cap (float, optional): maximum |z| contributed by one study (None: no cap). Defaults to 5.

    Returns:
        pd.DataFrame: one row per protein with the number of studies, studies significant at the Nebula
        cutoff, mean log2FC, combined Z, p-value, FDR, direction consistency and per-study log2FC.
    """
    cutoff = getattr(self, 'pvalue_cutoff', 0.05)
    long = []
    for df, group in zip(self.original, self.groups):
        d = df[['gene_name', self.pvalue, 'log2(fc)']].copy()
        d.columns = ['gene_name', 'p', 'fc']
        d['gene_name'] = d['gene_name'].astype(str)
        d = d[(d['gene_name'] != 'nan') & d['p'].notna() & d['fc'].notna()]
        d = d.sort_values('p').drop_duplicates('gene_name')
        d['study'] = group
        long.append(d)
    long = pd.concat(long, ignore_index=True)
    p = long['p'].clip(1e-300, 1)
    z = stats.norm.isf(p / 2)
    long['z'] = np.sign(long['fc']) * (z if z_cap is None else np.minimum(z, z_cap))
    g = long.groupby('gene_name')
    out = pd.DataFrame({'n_studies': g.size(), 'n_significant': g['p'].apply(lambda x: int((x <= cutoff).sum())),
                        'mean log2(fc)': g['fc'].mean(), 'Z': g['z'].sum() / np.sqrt(g.size())})
    out = out[out['n_studies'] >= min_studies].copy()
    if out.empty:
        raise ValueError(f'No protein was quantified in at least {min_studies} studies')
    out['pvalue'] = 2 * stats.norm.sf(np.abs(out['Z']))
    out['FDR'] = multipletests(out['pvalue'], method='fdr_bh')[1]
    agree = np.sign(long['fc']) == np.sign(long['gene_name'].map(out['Z']))
    out['consistency'] = agree.groupby(long['gene_name']).mean().reindex(out.index)
    wide = long.pivot(index='gene_name', columns='study', values='fc').reindex(columns=self.groups)
    wide.columns = [f'log2(fc) {c}' for c in wide.columns]
    out = out.join(wide).sort_values('pvalue')
    return out.reset_index().rename(columns={'index': 'gene_name'})


def meta_plot(self, top=30, min_studies=2, z_cap=5.0, palette='RdBu_r', save=None, dpi=300, vector=True):
    """Top meta-analysis proteins across studies.

    Rows are proteins (best combined p-value first), columns are studies. Colour is the study's log2 fold
    change and size its -log10(p); the left labels give the combined FDR.

    Args:
        top (int, optional): number of proteins. Defaults to 30.
        min_studies (int, optional): minimum number of studies quantifying the protein. Defaults to 2.
        z_cap (float, optional): maximum |z| contributed by one study. Defaults to 5.
        palette (str, optional): diverging palette. Defaults to 'RdBu_r'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    meta = meta_analysis(self, min_studies=min_studies, z_cap=z_cap).head(top)
    genes = list(meta['gene_name'])
    rows = []
    for df, group in zip(self.original, self.groups):
        d = df[df['gene_name'].astype(str).isin(genes)].sort_values(self.pvalue).drop_duplicates('gene_name')
        rows += [(r['gene_name'], group, r['log2(fc)'], r[self.pvalue]) for _, r in d.iterrows()]
    d = pd.DataFrame(rows, columns=['gene_name', 'study', 'fc', 'p'])
    d['size'] = -np.log10(d['p'].clip(1e-300, 1))
    # extreme fold changes (on/off proteins) would wash out the rest: cap the scale (colours clip at the ends)
    lim = min(float(np.nanpercentile(np.abs(d['fc']), 90)) or 1, 4)
    yi = {g: i for i, g in enumerate(genes)}
    xi = {s: i for i, s in enumerate(self.groups)}
    fig, ax = plt.subplots(figsize=(1.1 * len(self.groups) + 3.5, 0.3 * len(genes) + 1.8), dpi=dpi)
    sc = ax.scatter(d['study'].map(xi), d['gene_name'].map(yi), c=d['fc'].clip(-lim, lim), cmap=palette, vmin=-lim, vmax=lim,
                    s=20 + 22 * d['size'].clip(upper=10), edgecolor='#333', linewidth=0.4)
    ax.set_xticks(range(len(self.groups)), self.groups, rotation=45, ha='right')
    ax.set_yticks(range(len(genes)), [f"{g}  (FDR {f:.1e})" for g, f in zip(genes, meta['FDR'])], fontsize=9)
    ax.set_ylim(len(genes) - 0.5, -0.5)
    ax.set_xlim(-0.6, len(self.groups) - 0.4)
    ax.set_title('Proteins consistently changed across studies', loc='left')
    fig.colorbar(sc, ax=ax, shrink=0.5, label='log2(fc)')
    sns.despine(left=True, bottom=True)
    if save is not None:
        plt.savefig(save + 'meta_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return meta
