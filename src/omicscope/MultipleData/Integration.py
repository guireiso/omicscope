""" Comparing independent studies on what they have in common.

Independent experiments differ in platform, sample preparation, coverage and normalisation, so abundances are
never compared across studies: only each study's own effects (log2 fold changes against its control) are. Two
rules follow, and every function here applies them:

1. The shared universe. A protein that was not quantified in study B is unknown there, not unchanged.
   Comparisons between two studies use only the proteins quantified in both.
2. Direction, not just overlap. Being regulated in both studies is weak evidence if the changes go opposite
   ways; concordance is measured on the direction of change.

- shared_universe: proteins quantified in both studies.
- overlap_test: overlap of the regulated proteins of two studies against their shared universe
  (hypergeometric test), split by direction of change.
- overlap_matrix: the overlap p-values of every pair of studies.
- concordance: do the proteins regulated in one study change the same way in the other? (sign test on the
  other study's fold changes, whether significant there or not, in both directions), plus quadrant counts and
  the fold-change correlation of the proteins regulated in either study.
- concordance_table / concordance_plot: all pairs, and the fold-change scatter of one pair with these numbers.
- study_diagnostics / pvalue_histograms: what each study looks like before comparing it with the others
  (proteins quantified and regulated, replicates, estimated share of unchanged proteins pi0 and the p-value
  distribution, which is flat with a peak near zero in a well-calibrated study).
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


def _index(self, study):
    return self.groups.index(study) if not isinstance(study, (int, np.integer)) else int(study)


def _study(self, study):
    """One row per gene: log2(fc), the p-value used for calling, raw p-value, effect and SE when present."""
    df = self.original[_index(self, study)].copy()
    df['gene_name'] = df['gene_name'].astype(str)
    df = df[(df['gene_name'] != 'nan') & df['log2(fc)'].notna()]
    df = df.sort_values(self.pvalue).drop_duplicates('gene_name').set_index('gene_name')
    return df


def _regulated(self, df, cutoff=None):
    return df[df[self.pvalue] <= (getattr(self, 'pvalue_cutoff', 0.05) if cutoff is None else cutoff)]


def shared_universe(self, a, b):
    """Genes quantified in both studies."""
    return set(_study(self, a).index) & set(_study(self, b).index)


def overlap_test(self, a, b, cutoff=None, background=None):
    """Overlap of the regulated proteins of two studies within their shared universe.

    Args:
        a, b (str): study names.
        cutoff (float, optional): p-value cutoff for regulated proteins (default: the Nebula cutoff).
        background (int, optional): fixed number of entities to use instead of the shared universe (then the
            regulated proteins of each study are not restricted to the universe). Not recommended.

    Returns:
        dict: universe (proteins quantified in both), regulated_a/regulated_b (within the universe), overlap,
        expected overlap by chance, enrichment (overlap / expected), pvalue (one-sided hypergeometric) and the
        overlap split into same_direction and opposite_direction.
    """
    da, db = _study(self, a), _study(self, b)
    universe = da.index.intersection(db.index)
    if background is None:
        ra, rb = _regulated(self, da.loc[universe], cutoff), _regulated(self, db.loc[universe], cutoff)
        n = len(universe)
    else:
        ra, rb, n = _regulated(self, da, cutoff), _regulated(self, db, cutoff), int(background)
    both = ra.index.intersection(rb.index)
    na, nb, k = len(ra), len(rb), len(both)
    same = int((np.sign(ra.loc[both, 'log2(fc)']) == np.sign(rb.loc[both, 'log2(fc)'])).sum())
    expected = na * nb / n if n else np.nan
    return {'universe': n, 'regulated_a': na, 'regulated_b': nb, 'overlap': k, 'expected': expected,
            'enrichment': k / expected if expected else np.nan,
            'pvalue': float(stats.hypergeom.sf(k - 1, n, na, nb)) if n else np.nan,
            'same_direction': same, 'opposite_direction': k - same}


def overlap_matrix(self, cutoff=None, background=None):
    """Overlap p-values (overlap_test) of every pair of studies, as a symmetric DataFrame."""
    m = pd.DataFrame(np.zeros((len(self.groups),) * 2), index=self.groups, columns=self.groups)
    for a, b in itertools.combinations(self.groups, 2):
        m.loc[a, b] = m.loc[b, a] = overlap_test(self, a, b, cutoff=cutoff, background=background)['pvalue']
    return m


def _direction(source, target):
    """Sign agreement of `source`'s regulated proteins in `target` (any significance in target)."""
    shared = source.index.intersection(target.index)
    x, y = source.loc[shared, 'log2(fc)'], target.loc[shared, 'log2(fc)']
    keep = (x != 0) & (y != 0) & np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    same = int((np.sign(x) == np.sign(y)).sum())
    n = len(x)
    rho = stats.spearmanr(x, y)[0] if n > 2 else np.nan
    return {'n': n, 'same': same, 'fraction': same / n if n else np.nan,
            'pvalue': float(stats.binomtest(same, n, 0.5, alternative='greater').pvalue) if n else np.nan,
            'spearman': rho}


def _r_ci(r, n, level=0.95):
    if n <= 3 or not np.isfinite(r):
        return np.nan, np.nan
    z, half = np.arctanh(np.clip(r, -0.999999, 0.999999)), stats.norm.ppf(0.5 + level / 2) / np.sqrt(n - 3)
    return float(np.tanh(z - half)), float(np.tanh(z + half))


def concordance(self, a, b):
    """Do two studies change in the same direction?

    Thresholds are applied on one side only: the proteins regulated in study a are looked up in study b
    whatever their significance there (a protein at p = 0.06 in b still counts), and vice versa. This keeps
    a study with less power from hiding a real agreement.

    Returns:
        dict: a_in_b and b_in_a (n proteins, same direction, fraction, one-sided sign-test p-value, Spearman
        rho), quadrants (counts of the proteins regulated in either study by direction in a and b), and the
        Pearson correlation of their fold changes with its 95% CI.
    """
    da, db = _study(self, a), _study(self, b)
    universe = da.index.intersection(db.index)
    da, db = da.loc[universe], db.loc[universe]
    ra, rb = _regulated(self, da), _regulated(self, db)
    either = ra.index.union(rb.index)
    x, y = da.loc[either, 'log2(fc)'], db.loc[either, 'log2(fc)']
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    quadrants = {'up_up': int(((x > 0) & (y > 0)).sum()), 'up_down': int(((x > 0) & (y < 0)).sum()),
                 'down_up': int(((x < 0) & (y > 0)).sum()), 'down_down': int(((x < 0) & (y < 0)).sum())}
    r = stats.pearsonr(x, y)[0] if len(x) > 2 else np.nan
    low, high = _r_ci(r, len(x))
    return {'study_a': a, 'study_b': b, 'universe': len(universe), 'a_in_b': _direction(ra, db),
            'b_in_a': _direction(rb, da), 'quadrants': quadrants, 'pearson': r, 'pearson_ci': (low, high),
            'n_either': len(x)}


def concordance_table(self):
    """Directional concordance of every ordered pair of studies (regulated in `source`, looked up in `target`)."""
    rows = []
    for a, b in itertools.permutations(self.groups, 2):
        c = _direction(_regulated(self, _study(self, a)), _study(self, b))
        rows.append({'source': a, 'target': b, 'proteins': c['n'], 'same direction': c['same'],
                     'fraction': c['fraction'], 'sign-test p': c['pvalue'], 'spearman': c['spearman']})
    return pd.DataFrame(rows)


def concordance_plot(self, a, b, palette='Dark2', save=None, dpi=300, vector=True):
    """Fold changes of two studies on their shared universe, with the concordance numbers.

    Proteins regulated in either study are coloured (in a only, in b only, in both); the corners give the
    quadrant counts and the title the directional agreement in each direction.

    Args:
        a, b (str): study names.
        palette (str, optional): colours of the three groups. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    c = concordance(self, a, b)
    da, db = _study(self, a), _study(self, b)
    universe = da.index.intersection(db.index)
    ra, rb = set(_regulated(self, da.loc[universe]).index), set(_regulated(self, db.loc[universe]).index)
    df = pd.DataFrame({'x': da.loc[universe, 'log2(fc)'], 'y': db.loc[universe, 'log2(fc)']})
    df['group'] = ['both' if g in ra and g in rb else f'{a} only' if g in ra else f'{b} only' if g in rb
                   else 'neither' for g in df.index]
    colors = dict(zip([f'{a} only', f'{b} only', 'both'], sns.color_palette(palette, 3).as_hex()))
    colors['neither'] = '#d5d8dc'
    fig, ax = plt.subplots(figsize=(6.8, 6.4), dpi=dpi)
    for g in ['neither', f'{a} only', f'{b} only', 'both']:
        sub = df[df['group'] == g]
        ax.scatter(sub['x'], sub['y'], s=8 if g == 'neither' else 16, color=colors[g], alpha=0.5 if g == 'neither' else 0.85,
                   lw=0, label=f'{g} ({len(sub)})')
    ax.axhline(0, color='#888', lw=0.8)
    ax.axvline(0, color='#888', lw=0.8)
    q = c['quadrants']
    for (fx, fy, key) in [(0.97, 0.97, 'up_up'), (0.03, 0.97, 'down_up'), (0.03, 0.03, 'down_down'), (0.97, 0.03, 'up_down')]:
        ax.text(fx, fy, str(q[key]), transform=ax.transAxes, ha='right' if fx > 0.5 else 'left',
                va='top' if fy > 0.5 else 'bottom', fontsize=12, fontweight='bold', color='#333')
    ab, ba = c['a_in_b'], c['b_in_a']
    ax.set_title(f"{a} → {b}: {ab['fraction']:.0%} same direction (p = {ab['pvalue']:.1e})\n"
                 f"{b} → {a}: {ba['fraction']:.0%} same direction (p = {ba['pvalue']:.1e});  "
                 f"r = {c['pearson']:.2f} [{c['pearson_ci'][0]:.2f}, {c['pearson_ci'][1]:.2f}]", fontsize=10, loc='left')
    ax.set(xlabel=f'log2(fc) {a}', ylabel=f'log2(fc) {b}')
    ax.legend(frameon=False, fontsize=9, loc='center left', bbox_to_anchor=(1, 0.5))
    sns.despine()
    if save is not None:
        plt.savefig(save + 'concordance_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return c


def _pi0(p):
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    if len(p) < 50:
        return np.nan
    lam = np.arange(0.05, 0.96, 0.05)
    pi = np.array([(p > lmb).mean() / (1 - lmb) for lmb in lam])
    return float(np.clip(np.median(pi[lam >= 0.5]), 0, 1))     # robust, Storey-style estimate


def study_diagnostics(self):
    """One row per study: what was measured and how the statistics behave.

    pi0 is the estimated share of unchanged proteins (from the raw p-values; needs files written by OmicScope
    2.0, which store them, or files whose calling column is the raw p-value). A very low pi0, or a p-value
    histogram with a bump near 1, suggests a miscalibrated test or a normalisation problem in that study.
    """
    rows = []
    meta = getattr(self, 'metadata', [{}] * len(self.groups))
    for study, df, m in zip(self.groups, self.original, meta):
        d = _study(self, study)
        reg = _regulated(self, d)
        raw = d['pvalue'] if 'pvalue' in d.columns else None
        reps = (m or {}).get('Replicates', {})
        rows.append({'study': study, 'data type': (m or {}).get('DataType', 'proteome'),
                     'quantified': len(d), 'regulated': len(reg), 'up': int((reg['log2(fc)'] > 0).sum()),
                     'down': int((reg['log2(fc)'] < 0).sum()), 'regulated (%)': 100 * len(reg) / max(len(d), 1),
                     'pi0': _pi0(raw) if raw is not None else np.nan,
                     'replicates': ', '.join(f'{k}: {v}' for k, v in reps.items()) or 'not recorded',
                     'effect sizes': 'SE' in d.columns})
    return pd.DataFrame(rows)


def pvalue_histograms(self, bins=20, save=None, dpi=300, vector=True):
    """Distribution of the raw p-values of each study (files written by OmicScope 2.0).

    A well-calibrated study shows a flat histogram with a peak near zero (the changed proteins); a U shape or a
    peak near one points to a problem with the test or the normalisation.
    """
    studies = [(s, _study(self, s)) for s in self.groups]
    studies = [(s, d['pvalue']) for s, d in studies if 'pvalue' in d.columns]
    if not studies:
        raise ValueError('The .omics files do not store raw p-values (export them again with OmicScope 2.0).')
    colors = dict(zip(self.groups, self.colors))
    fig, axes = plt.subplots(1, len(studies), figsize=(3.2 * len(studies), 3), dpi=dpi, squeeze=False)
    for ax, (s, p) in zip(axes[0], studies):
        ax.hist(p.dropna(), bins=bins, range=(0, 1), color=colors.get(s, '#555'), edgecolor='white')
        ax.set(title=f'{s} (pi0 = {_pi0(p):.2f})', xlabel='p-value', ylabel='proteins')
        sns.despine(ax=ax)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save + 'pvalue_histograms.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return pd.DataFrame({s: [_pi0(p)] for s, p in studies}, index=['pi0'])
