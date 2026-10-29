""" Protein-level meta-analysis across Nebula studies.

Two ways of combining a protein's results across studies:

- Effect sizes (method='effect', used automatically when every .omics file stores them, as files written by
  OmicScope 2.0 do): random-effects meta-analysis (DerSimonian & Laird, 1986) of each study's difference of
  log2 means, weighted by its precision. It gives the combined effect with its 95% confidence interval and
  I^2, the share of the variation between studies that is beyond chance. A protein that rises in two studies
  and falls in a third gets a small combined effect with a high I^2: real, but context dependent, instead of
  simply disappearing.
- p-values (method='stouffer'): Stouffer's method on signed z-scores, z = sign(log2FC) * Phi^-1(1 - p/2),
  Z = sum(z) / sqrt(k). Each study's abs(z) is capped (default 5, p ~ 6e-7) so that one study with extreme
  p-values cannot carry the combination on its own. Raw p-values are used when the files store them; older
  files only store the calling column (often adjusted), which makes the combination conservative.

Both are restricted to the proteins quantified in at least `min_studies` studies (a protein missing from a
study is unknown there, not unchanged), and combined p-values are corrected with Benjamini-Hochberg.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _long(self):
    """One row per (study, gene): fold change, p-value for Stouffer, effect and SE when stored."""
    cutoff = getattr(self, 'pvalue_cutoff', 0.05)
    rows = []
    for df, group in zip(self.original, self.groups):
        d = df.copy()
        d['gene_name'] = d['gene_name'].astype(str)
        d = d[(d['gene_name'] != 'nan') & d['log2(fc)'].notna()].sort_values(self.pvalue).drop_duplicates('gene_name')
        rows.append(pd.DataFrame({
            'gene_name': d['gene_name'], 'study': group, 'fc': d['log2(fc)'],
            'p': d['pvalue'] if 'pvalue' in d.columns else d[self.pvalue],
            'called': d[self.pvalue] <= cutoff,
            'effect': d['Mean difference'] if 'Mean difference' in d.columns else np.nan,
            'se': d['SE'] if 'SE' in d.columns else np.nan}))
    return pd.concat(rows, ignore_index=True)


def has_effect_sizes(self):
    """True when every study stores effect sizes and standard errors (.omics written by OmicScope 2.0)."""
    return all({'Mean difference', 'SE'} <= set(df.columns) for df in self.original)


def _random_effects(y, se):
    """DerSimonian-Laird random-effects estimate: (effect, SE, Q, I2, tau2)."""
    w = 1 / se ** 2
    fixed = np.sum(w * y) / np.sum(w)
    q = float(np.sum(w * (y - fixed) ** 2))
    k = len(y)
    tau2 = max(0.0, (q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w))) if k > 1 else 0.0
    ws = 1 / (se ** 2 + tau2)
    effect = float(np.sum(ws * y) / np.sum(ws))
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return effect, float(np.sqrt(1 / np.sum(ws))), q, i2, tau2


def meta_analysis(self, min_studies=2, z_cap=5.0, method='auto'):
    """Combine each protein's evidence across studies.

    Args:
        min_studies (int, optional): minimum number of studies quantifying the protein. Defaults to 2.
        z_cap (float, optional): Stouffer only: maximum abs(z) contributed by one study (None: no cap).
            Defaults to 5.
        method (str, optional): 'effect' (random-effects meta-analysis of effect sizes), 'stouffer'
            (combined p-values) or 'auto' (effect sizes when every study stores them). Defaults to 'auto'.

    Returns:
        pd.DataFrame: one row per protein with the number of studies, studies significant at the Nebula cutoff,
        the combined result (effect, SE, CI95, I2 and Q p-value for 'effect'; Z for 'stouffer'), p-value, FDR,
        direction consistency and the log2FC of every study. `result.attrs['method']` says which was used.
    """
    if method == 'auto':
        method = 'effect' if has_effect_sizes(self) else 'stouffer'
    if method not in ('effect', 'stouffer'):
        raise ValueError("method must be 'auto', 'effect' or 'stouffer'")
    long = _long(self)
    if method == 'effect':
        long = long[np.isfinite(long['effect']) & np.isfinite(long['se']) & (long['se'] > 0)]
    g = long.groupby('gene_name')
    out = pd.DataFrame({'n_studies': g.size(), 'n_significant': g['called'].sum().astype(int),
                        'mean log2(fc)': g['fc'].mean()})
    out = out[out['n_studies'] >= min_studies].copy()
    if out.empty:
        raise ValueError(f'No protein was quantified in at least {min_studies} studies')
    if method == 'effect':
        res = {gene: _random_effects(d['effect'].to_numpy(float), d['se'].to_numpy(float))
               for gene, d in long[long['gene_name'].isin(out.index)].groupby('gene_name')}
        res = pd.DataFrame(res, index=['effect', 'SE', 'Q', 'I2', 'tau2']).T.reindex(out.index)
        out = out.join(res)
        out['CI95 low'], out['CI95 high'] = out['effect'] - 1.96 * out['SE'], out['effect'] + 1.96 * out['SE']
        out['pvalue'] = 2 * stats.norm.sf(np.abs(out['effect'] / out['SE']))
        out['Q pvalue'] = stats.chi2.sf(out['Q'], out['n_studies'] - 1)
        direction = np.sign(out['effect'])
    else:
        z = stats.norm.isf(long['p'].clip(1e-300, 1) / 2)
        long = long.assign(z=np.sign(long['fc']) * (z if z_cap is None else np.minimum(z, z_cap)))
        out['Z'] = long.groupby('gene_name')['z'].sum().reindex(out.index) / np.sqrt(out['n_studies'])
        out['pvalue'] = 2 * stats.norm.sf(np.abs(out['Z']))
        direction = np.sign(out['Z'])
    out['FDR'] = multipletests(out['pvalue'], method='fdr_bh')[1]
    agree = np.sign(long['fc']) == np.sign(long['gene_name'].map(direction))
    out['consistency'] = agree.groupby(long['gene_name']).mean().reindex(out.index)
    wide = long.pivot(index='gene_name', columns='study', values='fc').reindex(columns=self.groups)
    wide.columns = [f'log2(fc) {c}' for c in wide.columns]
    out = out.join(wide).sort_values('pvalue').reset_index().rename(columns={'index': 'gene_name'})
    out.attrs['method'] = method
    return out


def forest_plot(self, gene, color='#1b9e77', save=None, dpi=300, vector=True):
    """Effect of one protein in each study (difference of log2 means, 95% CI) and the random-effects estimate.

    Args:
        gene (str): gene name.
        color (str, optional): colour of the study estimates. Defaults to '#1b9e77'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    long = _long(self)
    d = long[(long['gene_name'] == gene) & np.isfinite(long['effect']) & (long['se'] > 0)]
    if len(d) < 2:
        raise ValueError(f'{gene} has effect sizes in fewer than two studies.')
    effect, se, q, i2, _ = _random_effects(d['effect'].to_numpy(float), d['se'].to_numpy(float))
    fig, ax = plt.subplots(figsize=(6.5, 0.55 * len(d) + 1.8), dpi=dpi)
    y = np.arange(len(d))[::-1] + 1
    ax.hlines(y, d['effect'] - 1.96 * d['se'], d['effect'] + 1.96 * d['se'], color=color, lw=2)
    ax.scatter(d['effect'], y, s=30 + 400 * (1 / d['se'] ** 2) / (1 / d['se'] ** 2).max(), color=color, zorder=3,
               edgecolor='#222', linewidth=0.5)
    ax.fill([effect - 1.96 * se, effect, effect + 1.96 * se, effect], [0, 0.25, 0, -0.25], color='#333')
    ax.axvline(0, color='#888', ls='--', lw=0.8)
    ax.set_yticks(list(y) + [0], list(d['study']) + ['Combined'])
    ax.set(xlabel='difference of log2 means (95% CI)', title=f'{gene}: {effect:.2f} [{effect - 1.96 * se:.2f}, '
           f'{effect + 1.96 * se:.2f}], I² = {i2:.0%}')
    sns.despine()
    if save is not None:
        plt.savefig(save + f'forest_{gene}.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return fig


def meta_plot(self, top=30, min_studies=2, z_cap=5.0, method='auto', palette='RdBu_r', save=None, dpi=300,
              vector=True):
    """Top meta-analysis proteins across studies.

    Rows are proteins (best combined p-value first), columns are studies. Colour is the study's log2 fold
    change and size its -log10(p); the left labels give the combined FDR.

    Args:
        top (int, optional): number of proteins. Defaults to 30.
        min_studies (int, optional): minimum number of studies quantifying the protein. Defaults to 2.
        z_cap (float, optional): maximum abs(z) contributed by one study. Defaults to 5.
        palette (str, optional): diverging palette. Defaults to 'RdBu_r'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    meta = meta_analysis(self, min_studies=min_studies, z_cap=z_cap, method=method).head(top)
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
