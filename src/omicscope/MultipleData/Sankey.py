""" Studies -> enriched terms -> proteins (Sankey diagram).

- sankey_table: the flows behind the diagram. A study flows into a term when the term is enriched in that
  study (width: its regulated proteins in the term); a term flows into each of its regulated proteins
  (width: number of studies in which the protein is regulated and the term enriched).
- sankey_plot: the diagram; next to each protein, one cell per study with its log2 fold change in that study
  (grey: not regulated), so proteins shared by several studies show whether they change the same way.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path


def sankey_table(self, *Terms, top=40):
    """Study-term-protein memberships of the selected enriched terms.

    Args:
        *Terms (str): enriched terms (exact names, as in the enrichment results).
        top (int, optional): proteins kept, those regulated in most studies first. Defaults to 40.

    Returns:
        pd.DataFrame: one row per (study, term, protein) with the protein's log2(fc) in that study.
    """
    cutoff = getattr(self, 'pvalue_cutoff', 0.05)
    rows = []
    for study, df, enr in zip(self.groups, self.original, self.enrichment):
        if enr is None:
            continue
        reg = df[df[self.pvalue] <= cutoff].groupby('gene_name')['log2(fc)'].mean()
        for _, r in enr[enr['Term'].isin(Terms)].iterrows():
            rows += [(study, r['Term'], g, float(reg[g])) for g in r['Genes'] if g in reg.index]
    table = pd.DataFrame(rows, columns=['study', 'term', 'gene_name', 'log2(fc)']).drop_duplicates(
        ['study', 'term', 'gene_name'])
    if table.empty:
        raise ValueError('None of the selected terms has regulated proteins in the studies.')
    rank = table.groupby('gene_name').agg(n=('study', 'nunique'), fc=('log2(fc)', 'mean'))
    rank = rank.assign(a=rank['fc'].abs()).sort_values(['n', 'a'], ascending=False)
    return table[table['gene_name'].isin(rank.index[:top])].reset_index(drop=True)


def _ribbon(ax, x0, y0, x1, y1, h0, h1, color, alpha):
    xm = (x0 + x1) / 2
    verts = [(x0, y0), (xm, y0), (xm, y1), (x1, y1), (x1, y1 + h1), (xm, y1 + h1), (xm, y0 + h0), (x0, y0 + h0),
             (x0, y0)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), fc=color, ec='none', alpha=alpha))


def sankey_plot(self, *Terms, top=40, palette='RdBu_r', vmin=-1, vmax=1, save=None, dpi=300, vector=True):
    """Sankey diagram: studies -> selected enriched terms -> their regulated proteins.

    Args:
        *Terms (str): enriched terms (exact names).
        top (int, optional): proteins kept, those regulated in most studies first. Defaults to 40.
        palette (str, optional): diverging palette for the protein fold change. Defaults to 'RdBu_r'.
        vmin, vmax (float, optional): fold-change range of the palette. Defaults to -1 and 1.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.

    Returns:
        pd.DataFrame: sankey_table of the diagram.
    """
    table = sankey_table(self, *Terms, top=top)
    colors = dict(zip(self.groups, self.colors))
    study_term = table.groupby(['study', 'term']).size()
    term_gene = table.groupby(['term', 'gene_name']).size()
    gene_fc = table.groupby('gene_name')['log2(fc)'].mean()
    studies = [s for s in self.groups if s in set(table['study'])]
    terms = [t for t in Terms if t in set(table['term'])]
    # proteins ordered by the mean position of their terms, to limit crossings
    tpos = {t: i for i, t in enumerate(terms)}
    genes = sorted(gene_fc.index, key=lambda g: np.mean([tpos[t] for t in table.loc[table['gene_name'] == g, 'term']]))
    value = {('S', s): study_term[s].sum() for s in studies}
    value.update({('T', t): max(study_term.xs(t, level='term').sum(), term_gene[t].sum()) for t in terms})
    value.update({('P', g): term_gene.xs(g, level='gene_name').sum() for g in genes})
    height = max(6, 0.22 * len(genes) + 2)
    fig, ax = plt.subplots(figsize=(11, height), dpi=dpi)
    cols = {'S': (0.0, studies, 0.04), 'T': (1.0, terms, 0.04), 'P': (2.0, genes, 0.012)}
    # one scale for every column, so widths compare across columns; gaps take their share of the height
    scale = min((0.9 - gap * (len(names) - 1)) / sum(value[(k, n)] for n in names)
                for k, (_, names, gap) in cols.items())
    pos = {}
    for k, (x, names, gap) in cols.items():
        span = sum(value[(k, n)] * scale for n in names) + gap * (len(names) - 1)
        y = (1 - span) / 2
        for n in names[::-1]:
            h = value[(k, n)] * scale
            pos[(k, n)] = [x, y, h, y, y]               # x, y, height, next outgoing y, next incoming y
            y += h + gap
    cmap, norm = plt.get_cmap(palette), plt.Normalize(vmin, vmax)
    for (s, t), w in study_term.items():
        a, b = pos[('S', s)], pos[('T', t)]
        _ribbon(ax, a[0] + 0.03, a[3], b[0], b[4], w * scale, w * scale, colors[s], 0.45)
        a[3] += w * scale
        b[4] += w * scale
    for (t, g), w in term_gene.items():
        a, b = pos[('T', t)], pos[('P', g)]
        _ribbon(ax, a[0] + 0.03, a[3], b[0], b[4], w * scale, w * scale, '#9fb0c2', 0.45)
        a[3] += w * scale
        b[4] += w * scale
    # log2(fc) of each protein in every study (NaN where it is not regulated)
    cutoff = getattr(self, 'pvalue_cutoff', 0.05)
    per_study = {s: df[df[self.pvalue] <= cutoff].groupby('gene_name')['log2(fc)'].mean()
                 for s, df in zip(self.groups, self.original) if s in studies}
    cw = 0.045
    for (k, n), (x, y, h, _, _) in pos.items():
        color = colors.get(n) if k == 'S' else '#34495e' if k == 'T' else '#8a96a3'
        ax.add_patch(Rectangle((x, y), 0.03, h, fc=color, ec='none'))
        shift = 0.045
        if k == 'P':
            cell_h = min(max(h, 0.008), 0.02)
            for i, st in enumerate(studies):
                v = per_study[st].get(n, np.nan)
                ax.add_patch(Rectangle((x + 0.045 + i * cw, y + h / 2 - cell_h / 2), cw * 0.85, cell_h,
                                       fc=cmap(norm(v)) if np.isfinite(v) else '#eceff2', ec='#d6dbe0', lw=0.3))
            shift = 0.06 + cw * len(studies)
        if k == 'S':
            ax.text(x - 0.02, y + h / 2, n, ha='right', va='center', fontsize=11)
        else:
            ax.text(x + shift, y + h / 2, n, ha='left', va='center', fontsize=11 if k == 'T' else 7,
                    fontweight='bold' if k == 'T' else 'normal', zorder=5,
                    bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.5) if k == 'T' else None)
    for x, label in [(0.015, 'Studies'), (1.015, 'Enriched terms')]:
        ax.text(x, 0.945, label, ha='center', va='bottom', fontsize=11, color='#5b6370', fontweight='bold')
    for i, st in enumerate(studies):          # which study each cell column is
        ax.text(2.045 + i * cw + cw * 0.42, 0.945, st, rotation=90, ha='center', va='bottom', fontsize=8,
                color=colors[st])
    ax.text(2.1 + cw * len(studies), 0.945, 'Proteins', ha='left', va='bottom', fontsize=11, color='#5b6370',
            fontweight='bold')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig.colorbar(sm, ax=ax, shrink=0.3, label='protein log2(fc)', pad=0.12)
    ax.set(xlim=(-0.5, 2.9), ylim=(0, 1.08))
    ax.axis('off')
    ax.set_title('Studies, enriched terms and their proteins', loc='left')
    if save is not None:
        plt.savefig(save + 'sankey_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return table
