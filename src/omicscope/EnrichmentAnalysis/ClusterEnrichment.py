""" Enrichment of protein clusters, custom gene sets and redundant terms.

- cluster_enrichment: over-representation analysis (the ORA workflow of EnrichmentScope) of each protein
  group, e.g. co-expression modules or k-means trend clusters, in one long table.
- cluster_dotplot: compareCluster-style dotplot of the top terms of every cluster.
- read_gmt / write_gmt: custom gene sets in the GMT format (term, description, genes, tab separated).
  EnrichmentScope accepts gene-set dicts and .gmt paths in `dbs` next to Enrichr library names.
- collapse_terms: groups enriched terms that share most of their proteins and keeps the most
  significant term of each group as its representative.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def read_gmt(path):
    """Read a GMT file.

    Args:
        path (str): .gmt file; each line is term, description and genes, separated by tabs.

    Returns:
        dict: term -> list of genes (duplicates and empty fields removed, order kept).
    """
    gene_sets = {}
    with open(path, encoding='utf-8') as f:
        for line in f:
            parts = [p.strip() for p in line.rstrip('\r\n').split('\t')]
            if len(parts) < 3 or not parts[0]:
                continue
            gene_sets[parts[0]] = list(dict.fromkeys(g for g in parts[2:] if g))
    return gene_sets


def write_gmt(gene_sets, path, description=''):
    """Write a dict term -> genes as a GMT file."""
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        for term, genes in gene_sets.items():
            f.write('\t'.join([str(term), description] + [str(g) for g in genes]) + '\n')


def gene_set_library(db):
    """Gene sets in the form gseapy expects and the name they get in Gene_set.

    Enrichr library names are kept as they are, dicts are named 'custom' and .gmt files are read and
    named after the file stem.
    """
    if isinstance(db, dict):
        return db, 'custom'
    if isinstance(db, (str, Path)) and str(db).lower().endswith('.gmt'):
        if not Path(db).exists():
            raise FileNotFoundError(f'Gene-set file not found: {db}')
        return read_gmt(db), Path(db).stem
    return db, db


def _cluster_genes(clusters):
    """dict cluster -> unique gene names, from a dict, a DataFrame or a coexpression_modules result."""
    if isinstance(clusters, dict) and isinstance(clusters.get('modules'), pd.DataFrame):
        clusters = clusters['modules']
    if isinstance(clusters, pd.DataFrame):
        col = next((c for c in ('Cluster', 'cluster', 'module', 'Module') if c in clusters.columns), None)
        if col is None or 'gene_name' not in clusters.columns:
            raise ValueError("clusters must have a 'gene_name' column and a 'Cluster', 'cluster' or 'module' column")
        df = clusters
        if col.lower() == 'module':
            df = df[df[col] != 0]      # module 0 holds the unassigned (grey) proteins
        clusters = {k: list(v) for k, v in df.groupby(col, sort=True)['gene_name']}
    out = {}
    for name, genes in clusters.items():
        genes = pd.Series(list(genes), dtype='string').str.strip()
        genes = genes[genes.notna() & (genes != '') & (genes.str.lower() != 'nan')]
        if len(genes):
            out[name] = list(dict.fromkeys(genes))
    return out


def cluster_enrichment(omicscope, clusters, dbs, organism='human', padjust_cutoff=0.05, background=True):
    """Over-representation analysis of each protein cluster.

    Args:
        omicscope (OmicScope): experiment the clusters come from (background and fold changes).
        clusters (dict or pd.DataFrame): cluster name -> gene names, a DataFrame with 'gene_name' and a
            'Cluster'/'cluster'/'module' column, or the result of coexpression_modules (module 0, the
            unassigned proteins, is skipped).
        dbs (list): Enrichr libraries, gene-set dicts (term -> genes) or paths to .gmt files.
        organism (str, optional): organism of the Enrichr libraries. Defaults to 'human'.
        padjust_cutoff (float, optional): adjusted p-value cutoff. Defaults to 0.05.
        background (bool, list, int, optional): True uses every protein quantified in the experiment
            (recommended); None uses all genes of the gene sets; a list or a number as in EnrichmentScope.

    Returns:
        pd.DataFrame: ORA results of all clusters with 'Cluster', 'Cluster_size' (genes in the cluster) and
        'GeneRatio' (genes of the cluster in the term / Cluster_size).
    """
    from .ORA import run_enrichr
    qd = omicscope.quant_data
    if background is True:
        background = list(qd['gene_name'].dropna().astype(str).drop_duplicates())
    foldchange = dict(zip(qd['gene_name'].astype(str).str.upper(), pd.to_numeric(qd['log2(fc)'], errors='coerce')))
    out = []
    for name, genes in _cluster_genes(clusters).items():
        df = run_enrichr(genes, dbs, organism=organism, background=background,
                         padjust_cutoff=padjust_cutoff, foldchange=foldchange)
        if df.empty:
            continue
        df = df.drop(columns='index', errors='ignore')
        df.insert(0, 'Cluster', name)
        df['Cluster_size'] = len(genes)
        df['GeneRatio'] = df['N_Proteins'] / len(genes)
        out.append(df)
    if not out:
        return pd.DataFrame(columns=['Cluster', 'Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value',
                                     'Genes', '-log10(pAdj)', 'N_Proteins', 'Cluster_size', 'GeneRatio'])
    return pd.concat(out, ignore_index=True)


def cluster_dotplot(results, top=5, palette='Reds', save=None, dpi=300, vector=True):
    """compareCluster-style dotplot.

    Columns are clusters and rows the union of the top terms of every cluster; dot size is the gene
    ratio (genes of the cluster in the term / genes in the cluster) and colour -log10(adjusted p).
    Terms enriched in a single cluster point to what makes that cluster distinct; terms in several
    clusters are shared biology.

    Args:
        results (pd.DataFrame): result of cluster_enrichment.
        top (int, optional): terms per cluster. Defaults to 5.
        palette (str, optional): sequential palette for -log10(pAdj). Defaults to 'Reds'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    if results.empty:
        raise ValueError('No enriched term to plot')
    clusters = list(dict.fromkeys(results['Cluster']))
    res = results.sort_values('Adjusted P-value')
    terms = list(dict.fromkeys(t for c in clusters for t in res.loc[res['Cluster'] == c, 'Term'].head(top)))
    df = res[res['Term'].isin(terms)].drop_duplicates(['Cluster', 'Term'])
    x = df['Cluster'].map({c: i for i, c in enumerate(clusters)})
    y = df['Term'].map({t: i for i, t in enumerate(terms)})
    size = 30 + 370 * df['GeneRatio'] / df['GeneRatio'].max()
    fig, ax = plt.subplots(figsize=(0.9 * len(clusters) + 5.5, 0.38 * len(terms) + 1.8), dpi=dpi)
    pts = ax.scatter(x, y, s=size, c=df['-log10(pAdj)'], cmap=palette, edgecolor='#333', linewidth=0.5,
                     vmin=0, zorder=3)
    ax.set_xticks(range(len(clusters)), [str(c) for c in clusters], rotation=45, ha='right')
    ax.set_yticks(range(len(terms)), [t if len(t) <= 60 else t[:57] + '...' for t in terms])
    ax.set(xlim=(-0.6, len(clusters) - 0.4), ylim=(len(terms) - 0.5, -0.5), xlabel='', ylabel='',
           title='Enrichment per cluster')
    ax.grid(color='#e5e5e5', lw=0.6, zorder=0)
    fig.colorbar(pts, ax=ax, shrink=0.5, label='-log10(pAdj)', pad=0.02)
    ratios = np.unique(np.round(np.quantile(df['GeneRatio'], [0.1, 0.5, 1]), 2))
    handles = [ax.scatter([], [], s=30 + 370 * r / df['GeneRatio'].max(), color='#bbb', edgecolor='#333') for r in ratios]
    ax.legend(handles, [f'{r:.2f}' for r in ratios], title='GeneRatio', frameon=False, loc='upper left',
              bbox_to_anchor=(1.25, 1), labelspacing=1.2)
    sns.despine(ax=ax, left=True, bottom=True)
    if save is not None:
        plt.savefig(save + 'cluster_dotplot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return df


def collapse_terms(self, similarity=0.5, metric='jaccard', representative_only=False):
    """Group redundant terms by the proteins they share.

    Terms are visited from the most to the least significant; each term not yet grouped becomes the
    representative of a new group, which takes every other ungrouped term whose proteins overlap it by
    at least `similarity`. For GSEA, terms with opposite NES signs are never grouped.

    Args:
        similarity (float, optional): minimum overlap to join a group. Defaults to 0.5.
        metric (str, optional): 'jaccard' (shared / union) or 'overlap' (shared / smaller set, which
            also groups a term with its sub-terms). Defaults to 'jaccard'.
        representative_only (bool, optional): return only the representative terms. Defaults to False.

    Returns:
        pd.DataFrame: self.results with 'Representative' (term representing the group) and 'Cluster_size'
        (number of terms in the group).
    """
    if metric not in ('jaccard', 'overlap'):
        raise ValueError("metric must be 'jaccard' or 'overlap'")
    res = self.results.copy()
    if res.empty:
        return res.assign(Representative=pd.Series(dtype=object), Cluster_size=pd.Series(dtype=int))
    by, ascending = ['Adjusted P-value'], [True]
    if 'Combined Score' in res.columns:
        by, ascending = by + ['Combined Score'], ascending + [False]
    elif 'NES' in res.columns:
        res['_abs_nes'] = res['NES'].abs()
        by, ascending = by + ['_abs_nes'], ascending + [False]
    res = res.sort_values(by, ascending=ascending, kind='stable').reset_index(drop=True)
    sets = [{str(g).upper() for g in genes} for genes in res['Genes']]
    universe = {g: i for i, g in enumerate(sorted(set().union(*sets)))}
    member = np.zeros((len(sets), len(universe)), np.float32)
    for i, s in enumerate(sets):
        member[i, [universe[g] for g in s]] = 1
    shared = member @ member.T
    size = member.sum(axis=1)
    if metric == 'jaccard':
        sim = shared / np.maximum(size[:, None] + size[None, :] - shared, 1)
    else:
        sim = shared / np.maximum(np.minimum.outer(size, size), 1)
    if 'NES' in res.columns:
        sign = np.sign(pd.to_numeric(res['NES'], errors='coerce').fillna(0).to_numpy())
        sim[sign[:, None] != sign[None, :]] = 0
    group = np.full(len(res), -1)
    for i in range(len(res)):
        if group[i] < 0:
            join = (group < 0) & (sim[i] >= similarity)
            join[i] = True
            group[join] = i
    res['Representative'] = res['Term'].to_numpy()[group]
    res['Cluster_size'] = pd.Series(group).map(pd.Series(group).value_counts()).to_numpy()
    res = res.drop(columns='_abs_nes', errors='ignore')
    if representative_only:
        res = res[group == np.arange(len(res))].reset_index(drop=True)
    return res
