""" Protein co-expression modules (a light WGCNA).

- coexpression_modules: groups proteins whose abundances rise and fall together across samples into
  modules, summarises each module by its eigengene (first principal component) and correlates the
  eigengenes with the experimental conditions (and TimeCourse, when present).
- module_trait_plot: heatmap of module-trait correlations.
- module_eigengene_plot: eigengene of each module in each sample, coloured by condition.
- module_network / module_network_plot: co-expression network of one module (strongest topological
  overlaps between its most connected proteins), with hub proteins highlighted.

Workflow (Langfelder & Horvath, BMC Bioinformatics 2008):

1. log2 abundances of the most variable proteins; missing values are filled with the lowest value of
   the protein (missing = below the detection limit).
2. Signed adjacency a_ij = ((1 + r_ij) / 2) ** power. Signed networks keep positively and negatively
   correlated proteins apart: a protein going up and one going down with the same factor belong to
   different biological programs, and in an unsigned network they would be lumped together.
3. power: the smallest soft threshold whose scale-free fit (signed R^2 of log10 p(k) on log10 k)
   reaches 0.8, or the best-fitting power when none does.
4. Topological overlap (TOM), which also rewards shared neighbours, and average-linkage clustering of
   1 - TOM.
5. Tree cut, top-down from the root: a branch is split in two when both halves have at least
   min_module_size proteins and their eigengenes are not similar (r < 1 - merge_cut); smaller side
   branches are set aside; a branch that cannot be split this way is a module.
6. Modules whose eigengenes correlate above 1 - merge_cut are merged, and every protein is assigned to
   the module it correlates best with (kME), provided kME >= 0.5 and is significant (p < 0.05); the
   others, and modules left with fewer than min_module_size proteins, go to module 0 (grey).
7. Clustering always finds some structure, even in noise. The same steps are run on data whose samples
   are shuffled independently for each protein (null_rounds times); a module is kept only if its
   eigengene explains more variance than every module found in the shuffled data. Likewise a protein
   joins a module only if its kME beats what 95% of shuffled proteins reach against the best of all
   modules. Without this, noise proteins form many small modules.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from .Exploration import _condition_of, _log_matrix, _sample_matrix

POWERS = tuple(range(1, 21))


def _expression_matrix(self, pvalue, max_proteins):
    """log2 matrix (proteins x samples) of the most variable proteins, missing filled with the protein minimum."""
    qd = self.quant_data
    sel = qd[pd.to_numeric(qd[self.pvalue], errors='coerce').fillna(1) <= pvalue]
    mat, genes = _sample_matrix(self, sel)
    logm = _log_matrix(self, mat)
    # filling many missing values with the minimum would create correlation from the missing pattern alone
    logm = logm[logm.notna().mean(axis=1) >= 0.5]
    logm = logm.apply(lambda r: r.fillna(r.min()), axis=1)
    logm = logm[logm.std(axis=1) > 0].dropna()
    logm = logm.loc[logm.var(axis=1).sort_values(ascending=False).index[:max_proteins]]
    if logm.shape[1] < 4:
        raise ValueError('At least 4 samples are needed to build a co-expression network.')
    if len(logm) < 2 * 10:
        raise ValueError('Too few proteins to build a co-expression network; relax the p-value cutoff.')
    return logm, genes


def _scale_free_fit(base, powers=POWERS):
    """Scale-free topology fit of the signed adjacency base ** power for each power."""
    rows = []
    for b in powers:
        adj = base ** b
        k = adj.sum(axis=1) - 1                       # diagonal of the adjacency is 1
        freq, edges = np.histogram(k, bins=10)
        centre = np.array([k[(k >= lo) & (k <= hi)].mean() if f else np.nan
                           for f, lo, hi in zip(freq, edges[:-1], edges[1:])])
        ok = (freq > 0) & (centre > 0)
        x, y = np.log10(centre[ok]), np.log10(freq[ok] / freq.sum())
        slope = np.polyfit(x, y, 1)[0] if ok.sum() > 2 else np.nan
        r2 = np.corrcoef(x, y)[0, 1] ** 2 if ok.sum() > 2 else np.nan
        rows.append({'power': b, 'R2': -np.sign(slope) * r2, 'slope': slope, 'mean_k': k.mean(),
                     'median_k': np.median(k)})
    return pd.DataFrame(rows)


def _tom_dissimilarity(adj):
    """1 - topological overlap of a (dense, symmetric) adjacency matrix."""
    adj = adj.copy()
    np.fill_diagonal(adj, 0)
    k = adj.sum(axis=1)
    tom = adj @ adj
    tom += adj
    tom /= np.minimum.outer(k, k) + 1 - adj
    np.fill_diagonal(tom, 1)
    return 1 - tom


def _eigengene(z):
    """First principal component of z-scored profiles (proteins x samples), z-scored and sign-aligned to their mean."""
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    eg = vt[0]
    if np.corrcoef(eg, z.mean(axis=0))[0, 1] < 0:
        eg = -eg
    return (eg - eg.mean()) / eg.std()


def _cut_tree(z, linkage_matrix, min_module_size, merge_cut):
    """Top-down tree cut: split branches into distinct halves, set aside small side branches."""
    from scipy.cluster.hierarchy import to_tree
    labels = np.zeros(len(z), int)
    stack, n_modules = [to_tree(linkage_matrix)], 0
    while stack:
        node = stack.pop()
        if node.count < min_module_size or node.is_leaf():
            continue
        big = [c for c in (node.left, node.right) if c.count >= min_module_size]
        if len(big) == 2:
            r = np.corrcoef(_eigengene(z[big[0].pre_order()]), _eigengene(z[big[1].pre_order()]))[0, 1]
            if r < 1 - merge_cut:
                stack += big
                continue
        elif len(big) == 1:
            stack.append(big[0])
            continue
        n_modules += 1
        labels[node.pre_order()] = n_modules
    return labels


def _merge_modules(z, labels, merge_cut):
    """Merge modules whose eigengenes correlate above 1 - merge_cut (most similar pair first)."""
    labels = labels.copy()
    while True:
        mods = [m for m in np.unique(labels) if m > 0]
        if len(mods) < 2:
            return labels
        r = np.corrcoef(np.array([_eigengene(z[labels == m]) for m in mods]))
        np.fill_diagonal(r, -1)
        i, j = np.unravel_index(np.argmax(r), r.shape)
        if r[i, j] <= 1 - merge_cut:
            return labels
        labels[labels == mods[j]] = mods[i]


def _assign_by_kme(z, labels, min_module_size, threshold=None):
    """Assign each protein to the module it correlates best with; weak or non-significant kME goes to grey."""
    n = z.shape[1]
    if threshold is None:
        t = stats.t.ppf(0.975, n - 2)
        threshold = max(0.5, t / np.sqrt(n - 2 + t ** 2))      # kME significant at p < 0.05
    for _ in range(3):
        mods = np.array([m for m in np.unique(labels) if m > 0])
        if len(mods) == 0:
            return labels, np.full(len(z), np.nan)
        eg = np.array([_eigengene(z[labels == m]) for m in mods])
        kme = z @ eg.T / n
        best = kme.argmax(axis=1)
        new = np.where(kme.max(axis=1) >= threshold, mods[best], 0)
        size = pd.Series(new).value_counts()
        new[np.isin(new, size.index[size < min_module_size])] = 0
        if np.array_equal(new, labels):
            break
        labels = new
    # modules numbered by size, largest first
    size = pd.Series(labels[labels > 0]).value_counts()
    order = {m: i for i, m in enumerate(size.index, start=1)}
    labels = np.array([order.get(m, 0) for m in labels])
    mods = np.arange(1, len(order) + 1)
    if len(mods) == 0:
        return labels, np.full(len(z), np.nan)
    eg = np.array([_eigengene(z[labels == m]) for m in mods])
    own = np.array([(z[i] @ eg[m - 1]) / n if m > 0 else np.nan for i, m in enumerate(labels)])
    return labels, own


def _variance_explained(z, labels):
    """Share of the variance of each module's (z-scored) proteins captured by its eigengene."""
    out = {}
    for m in np.unique(labels):
        if m > 0:
            sv = np.linalg.svd(z[labels == m], compute_uv=False)
            out[m] = float(sv[0] ** 2 / np.sum(sv ** 2))
    return out


def _modules_from(z, power, min_module_size, merge_cut):
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform
    n = z.shape[1]
    base = ((1 + (z @ z.T) / n) / 2).astype(np.float32)
    np.clip(base, 0, 1, out=base)
    diss = _tom_dissimilarity(base ** power)
    del base
    tree = linkage(squareform(diss.astype(float), checks=False), method='average')
    del diss
    labels = _cut_tree(z, tree, min_module_size, merge_cut)
    labels = _merge_modules(z, labels, merge_cut)
    return _assign_by_kme(z, labels, min_module_size)


def _null_threshold(z, power, min_module_size, merge_cut, rounds, rng):
    """Largest eigengene variance share of any module found in data shuffled within each protein."""
    best = 0.0
    for _ in range(rounds):
        shuffled = np.array([rng.permutation(row) for row in z])
        labels, _ = _modules_from(shuffled, power, min_module_size, merge_cut)
        best = max([best] + list(_variance_explained(shuffled, labels).values()))
    return best


def coexpression_modules(self, pvalue=1.0, power=None, min_module_size=20, merge_cut=0.25, max_proteins=3000,
                         null_rounds=3, random_state=0):
    """Co-expression modules, their eigengenes and module-trait correlations.

    Args:
        pvalue (float, optional): keep proteins with self.pvalue <= pvalue. Defaults to 1.0 (all proteins;
            co-expression is best estimated without pre-filtering on the differential test).
        power (int, optional): soft-thresholding power. None picks it by scale-free fit.
        min_module_size (int, optional): smallest module. Defaults to 20.
        merge_cut (float, optional): modules whose eigengenes correlate above 1 - merge_cut are merged.
            Defaults to 0.25.
        max_proteins (int, optional): network built on the most variable proteins. Defaults to 3000.
        null_rounds (int, optional): shuffled data sets used to discard modules no stronger than noise
            (0 keeps every module). Defaults to 3.
        random_state (int, optional): seed of the shuffles. Defaults to 0.

    Returns:
        dict:
            'modules': [Accession, gene_name, module, kME] (module 0 = unassigned, 'grey').
            'eigengenes': samples x modules (columns ME1, ME2, ...).
            'module_trait': [module, trait, r, pvalue], Pearson correlation of each eigengene with each
                condition (one-hot) and TimeCourse, if present.
            'traits': [Sample, Condition, (TimeCourse)] in the order of the eigengenes.
            'power': soft-thresholding power used.
            'power_table': scale-free fit of each candidate power.
            'null_threshold': eigengene variance share a module must exceed to be kept.
    """
    logm, genes = _expression_matrix(self, pvalue, max_proteins)
    x = logm.to_numpy(float)
    z = (x - x.mean(axis=1, keepdims=True)) / x.std(axis=1, keepdims=True)
    n = z.shape[1]
    base = ((1 + (z @ z.T) / n) / 2).astype(np.float32)
    np.clip(base, 0, 1, out=base)
    power_table = _scale_free_fit(base)
    if power is None:
        fit = power_table[power_table['R2'] >= 0.8]
        power = int(fit['power'].iloc[0] if len(fit) else power_table.loc[power_table['R2'].idxmax(), 'power'])
    del base
    labels, kme = _modules_from(z, power, min_module_size, merge_cut)
    threshold = 0.0
    if null_rounds:
        threshold = _null_threshold(z, power, min_module_size, merge_cut, null_rounds,
                                    np.random.default_rng(random_state))
        weak = [m for m, v in _variance_explained(z, labels).items() if v <= threshold]
        labels = np.where(np.isin(labels, weak), 0, labels)
        mods = [m for m in np.unique(labels) if m > 0]
        if mods:
            # membership cutoff: the best kME that proteins with shuffled samples reach in 95% of cases (it
            # accounts for testing every module and for conditions with very different variances)
            rng = np.random.default_rng(random_state + 1)
            eg = np.array([_eigengene(z[labels == m]) for m in mods])
            shuffled = np.array([rng.permutation(row) for row in z])
            kme_cut = float(np.quantile((shuffled @ eg.T / z.shape[1]).max(axis=1), 0.95))
            labels, kme = _assign_by_kme(z, labels, min_module_size, threshold=max(0.5, kme_cut))
        else:
            kme = np.full(len(z), np.nan)

    modules = pd.DataFrame({'Accession': logm.index, 'gene_name': genes.reindex(logm.index).to_numpy(),
                            'module': labels, 'kME': kme})
    modules = modules.sort_values(['module', 'kME'], ascending=[True, False], ignore_index=True)
    modules = pd.concat([modules[modules['module'] > 0], modules[modules['module'] == 0]], ignore_index=True)
    mods = sorted(m for m in set(labels) if m > 0)
    eigengenes = pd.DataFrame({f'ME{m}': _eigengene(z[labels == m]) for m in mods}, index=logm.columns)
    eigengenes.index.name = 'Sample'

    cond = _condition_of(self)
    traits = pd.DataFrame({'Sample': logm.columns, 'Condition': [cond[s] for s in logm.columns]})
    trait_values = {c: (traits['Condition'] == c).astype(float).to_numpy() for c in self.Conditions}
    if 'TimeCourse' in self.pdata.columns:
        time = dict(zip(self.pdata['Sample'], self.pdata['TimeCourse']))
        traits['TimeCourse'] = [time[s] for s in logm.columns]
        trait_values['TimeCourse'] = pd.to_numeric(traits['TimeCourse'], errors='coerce').to_numpy(float)
    rows = []
    for m in mods:
        for trait, values in trait_values.items():
            ok = np.isfinite(values)
            if np.nanstd(values) == 0:
                continue
            r, p = stats.pearsonr(eigengenes[f'ME{m}'].to_numpy()[ok], values[ok])
            rows.append({'module': m, 'trait': trait, 'r': r, 'pvalue': p})
    module_trait = pd.DataFrame(rows, columns=['module', 'trait', 'r', 'pvalue'])
    return {'modules': modules, 'eigengenes': eigengenes, 'module_trait': module_trait, 'traits': traits,
            'power': power, 'power_table': power_table, 'null_threshold': threshold,
            'settings': {'pvalue': pvalue, 'max_proteins': max_proteins}}


def module_trait_plot(self, modules=None, pvalue=1.0, power=None, min_module_size=20, merge_cut=0.25,
                      max_proteins=3000, palette='RdBu_r', save=None, dpi=300, vector=True):
    """Heatmap of module-trait correlations.

    Each cell is the Pearson correlation between a module eigengene and a trait (membership in a
    condition, or TimeCourse) with its p-value. Modules strongly correlated with a condition are the
    co-regulated protein programs that distinguish it; enrich their proteins (cluster_enrichment) to
    name them.

    Args:
        modules (dict, optional): result of coexpression_modules; computed with the arguments below if None.
        pvalue, power, min_module_size, merge_cut, max_proteins: see coexpression_modules.
        palette (str, optional): diverging palette, centred at r = 0. Defaults to 'RdBu_r'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    if modules is None:
        modules = coexpression_modules(self, pvalue=pvalue, power=power, min_module_size=min_module_size,
                                       merge_cut=merge_cut, max_proteins=max_proteins)
    mt = modules['module_trait']
    if mt.empty:
        raise ValueError('No co-expression module was found.')
    size = modules['modules']['module'].value_counts()
    r = mt.pivot(index='module', columns='trait', values='r')
    p = mt.pivot(index='module', columns='trait', values='pvalue')
    traits = [t for t in list(self.Conditions) + ['TimeCourse'] if t in r.columns]
    r, p = r[traits], p[traits]
    labels = np.array([[f'{rv:.2f}\n({pv:.1e})' for rv, pv in zip(rr, pp)] for rr, pp in zip(r.to_numpy(), p.to_numpy())])
    fig, ax = plt.subplots(figsize=(1.3 * len(traits) + 3, 0.55 * len(r) + 1.5), dpi=dpi)
    sns.heatmap(r, annot=labels, fmt='', cmap=palette, vmin=-1, vmax=1, center=0, linewidths=0.5,
                linecolor='white', annot_kws={'fontsize': 8}, cbar_kws={'label': 'Pearson r', 'shrink': 0.7},
                yticklabels=[f'ME{m} ({size.get(m, 0)})' for m in r.index], ax=ax)
    ax.set(xlabel='', ylabel='', title=f"Module-trait correlation (power = {modules['power']})")
    ax.tick_params(axis='y', rotation=0)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save + 'module_trait_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return mt


def module_eigengene_plot(self, modules=None, pvalue=1.0, power=None, min_module_size=20, merge_cut=0.25,
                          max_proteins=3000, palette='Dark2', save=None, dpi=300, vector=True):
    """Eigengene of each module in each sample, coloured by condition.

    The eigengene is the profile that best represents the proteins of a module: bars consistently
    above zero in one condition mean the module is up-regulated there.

    Args:
        modules (dict, optional): result of coexpression_modules; computed with the arguments below if None.
        pvalue, power, min_module_size, merge_cut, max_proteins: see coexpression_modules.
        palette (str, optional): Palette for conditions. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    if modules is None:
        modules = coexpression_modules(self, pvalue=pvalue, power=power, min_module_size=min_module_size,
                                       merge_cut=merge_cut, max_proteins=max_proteins)
    eg = modules['eigengenes']
    if eg.empty:
        raise ValueError('No co-expression module was found.')
    traits = modules['traits'].copy()
    traits['order'] = traits['Condition'].map({c: i for i, c in enumerate(self.Conditions)})
    traits = traits.sort_values(['order'] + (['TimeCourse'] if 'TimeCourse' in traits.columns else []), kind='stable')
    eg = eg.loc[traits['Sample']]
    colors = dict(zip(self.Conditions, sns.color_palette(palette, len(self.Conditions)).as_hex()))
    size = modules['modules']['module'].value_counts()
    ncol = min(3, eg.shape[1])
    nrow = int(np.ceil(eg.shape[1] / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 2.8 * nrow), dpi=dpi, squeeze=False, sharey=True)
    for ax, col in zip(axes.flat, eg.columns):
        ax.bar(range(len(eg)), eg[col], color=[colors[c] for c in traits['Condition']], width=0.85)
        ax.axhline(0, color='#555', lw=0.6)
        ax.set_title(f'{col} ({size.get(int(col[2:]), 0)} proteins)', fontsize=10)
        ax.set_xticks([])
        sns.despine(ax=ax)
    for ax in axes.flat[eg.shape[1]:]:
        ax.axis('off')
    for ax in axes[:, 0]:
        ax.set_ylabel('Eigengene')
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=v) for v in colors.values()], labels=list(colors),
               frameon=False, loc='upper center', ncol=len(colors), bbox_to_anchor=(0.5, 1.04))
    plt.tight_layout()
    if save is not None:
        plt.savefig(save + 'module_eigengene_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return eg


def module_network(self, module, modules=None, top=150, edges_per_node=3, **module_options):
    """Co-expression network of one module.

    The proteins of the module with the highest kME (at most `top`) are linked by their topological overlap
    (computed within the module with the module's soft-thresholding power); each protein keeps its
    `edges_per_node` strongest links, which keeps the hubs and the structure of the module readable.

    Args:
        module (int): module number (as in coexpression_modules()['modules']['module']).
        modules (dict, optional): result of coexpression_modules (computed with module_options if not given).
        top (int, optional): most connected proteins kept. Defaults to 150.
        edges_per_node (int, optional): strongest links kept per protein. Defaults to 3.

    Returns:
        networkx.Graph: nodes keyed by Accession with gene_name, kME, log2FC, hub (top 10% by kME);
        edges with weight (topological overlap).
    """
    import networkx as nx
    res = coexpression_modules(self, **module_options) if modules is None else modules
    members = res['modules']
    members = members[members['module'] == module].sort_values('kME', ascending=False).head(top)
    if len(members) < 3:
        raise ValueError(f'Module {module} has fewer than three proteins.')
    settings = res.get('settings', {'pvalue': 1.0, 'max_proteins': 3000})
    logm, _ = _expression_matrix(self, settings['pvalue'], settings['max_proteins'])
    x = logm.loc[members['Accession']].to_numpy(float)
    z = (x - x.mean(axis=1, keepdims=True)) / x.std(axis=1, keepdims=True)
    base = np.clip((1 + (z @ z.T) / z.shape[1]) / 2, 0, 1)
    tom = 1 - _tom_dissimilarity(base ** res['power'])
    np.fill_diagonal(tom, 0)
    fc = self.quant_data.drop_duplicates('Accession').set_index('Accession')['log2(fc)']
    hub_cut = members['kME'].quantile(0.9)
    G = nx.Graph()
    for r in members.itertuples(index=False):
        G.add_node(r.Accession, gene_name=str(r.gene_name), kME=float(r.kME),
                   log2FC=float(fc.get(r.Accession, np.nan)), hub=bool(r.kME >= hub_cut))
    acc = list(members['Accession'])
    for i in range(len(acc)):
        for j in np.argsort(tom[i])[::-1][:edges_per_node]:
            G.add_edge(acc[i], acc[j], weight=float(tom[i, j]))
    return G


def module_network_plot(self, module, modules=None, top=150, edges_per_node=3, palette='RdBu_r', save=None,
                        dpi=300, vector=True, **module_options):
    """Co-expression network of one module (see module_network).

    Node size is the module membership (kME), colour the log2 fold change, and hub proteins (top 10% kME) are
    labelled.

    Args:
        module (int): module number.
        modules (dict, optional): result of coexpression_modules.
        top (int, optional): most connected proteins kept. Defaults to 150.
        edges_per_node (int, optional): strongest links kept per protein. Defaults to 3.
        palette (str, optional): diverging palette for fold change. Defaults to 'RdBu_r'.
        save (str, optional): Path to save figure (also writes a GraphML for Cytoscape/Gephi). Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.

    Returns:
        networkx.Graph: the network drawn.
    """
    import networkx as nx
    G = module_network(self, module, modules=modules, top=top, edges_per_node=edges_per_node, **module_options)
    pos = nx.spring_layout(G, weight='weight', seed=0)
    fc = np.array([G.nodes[v]['log2FC'] for v in G])
    lim = float(np.nanpercentile(np.abs(fc), 98)) if np.isfinite(fc).any() else 1
    kme = np.array([G.nodes[v]['kME'] for v in G])
    fig, ax = plt.subplots(figsize=(8, 7), dpi=dpi)
    w = np.array([d['weight'] for _, _, d in G.edges(data=True)])
    nx.draw_networkx_edges(G, pos, ax=ax, width=0.5 + 2.5 * w / max(w.max(), 1e-9), edge_color='#9aa3ad', alpha=0.6)
    nodes = nx.draw_networkx_nodes(G, pos, ax=ax, node_size=30 + 300 * (kme - kme.min()) / max(np.ptp(kme), 1e-9),
                                   node_color=np.nan_to_num(fc), cmap=palette, vmin=-(lim or 1), vmax=lim or 1,
                                   edgecolors='#333', linewidths=0.4)
    nx.draw_networkx_labels(G, pos, {v: G.nodes[v]['gene_name'] for v in G if G.nodes[v]['hub']}, ax=ax, font_size=9)
    fig.colorbar(nodes, ax=ax, shrink=0.6, label='log2(fc)')
    ax.set_title(f'Module {module}: {G.number_of_nodes()} proteins (labels: hubs)')
    ax.axis('off')
    if save is not None:
        plt.savefig(save + f'module_{module}_network.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
        from .GraphExport import write_graphml
        write_graphml(G, save + f'module_{module}_network.graphml', pos)
    plt.show()
    return G
