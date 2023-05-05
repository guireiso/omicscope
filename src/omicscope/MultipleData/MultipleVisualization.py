from copy import copy

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns


def barplot(self, save=None, vector=True, dpi=300):
    """ Barplot
    Bar plot for proteins/genes identified and differentially regulated
    according to each group

    Args:
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg. Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    """
    plt.rcParams['figure.dpi'] = dpi
    data = copy(self)
    conditions = copy(data.groups)
    group_data = data.original
    difreg = data.group_data

    # -Desregulation figures
    whole_proteome = []
    deps = []
    for i, a in zip(group_data, difreg):
        whole_proteome.append(len(i))
        deps.append(len(a))
    conditions.extend(['Total'])

    proteinsIdentified = []
    for i in group_data:
        proteinsIdentified.append(i['gene_name'])
    proteinsIdentified = pd.concat(proteinsIdentified).drop_duplicates()

    proteinsdes = []
    for i in difreg:
        proteinsdes.append(i['gene_name'])
    proteinsdes = pd.concat(proteinsdes).drop_duplicates()

    whole_proteome.extend([len(proteinsIdentified)])
    deps.extend([len(proteinsdes)])
    r = [i for i in range(0, len(conditions))]

    f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
    ax.bar(r, whole_proteome, color='lightgray', edgecolor='black', width=0.5,
           linewidth=1)
    colors = copy(self.colors)
    colors.extend(['gray'])
    ax.bar(r, deps, edgecolor='black', width=0.5, linewidth=1)
    plt.xticks(r,
               fontweight=None, rotation=45, ha='right')

    ax2.bar(r, whole_proteome, color='lightgray', edgecolor='black', width=0.5,
            linewidth=1)
    ax2.bar(r, deps, color=colors,  edgecolor='black', width=0.5, linewidth=1)
    ax.set_ylim(max(whole_proteome[:-1])*1.1, max(whole_proteome)*1.01)  # outliers only
    ax2.set_ylim(0, max(whole_proteome[:-1])*1.005)  # most of the data
    sns.despine()
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(bottom=False)
    plt.xticks(r, conditions,
               fontweight=None, rotation=45, ha='right')
    for i, j, k, c in zip(deps[:-1], r[:-1], whole_proteome[:-1],
                          colors[:-1]):
        ax2.annotate(i, xy=[j, k*1.005], ha='center', color=c,
                     weight='bold')
    ax.annotate(deps[-1], [r[-1], whole_proteome[-1]*1.001], ha='center',
                color=colors[-1], weight='bold')
    if save is not None:
        if vector is True:
            plt.savefig(save + 'barplot.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'barplot.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def diff_reg(self,
             save=None, vector=True, dpi=300):
    """Dotplot

    Dotplot for number of proteins up- and down-regulated in each group.

    Args:
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg.
        Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    """
    plt.rcParams['figure.dpi'] = dpi
    data = copy(self)
    groups = data.groups
    difreg = data.group_data
    up = []
    down = []

    for i in difreg:
        up.append(len(i[i['log2(fc)'] > 0]))
        down.append(-len(i[i['log2(fc)'] < 0]))
    dysregulations = pd.DataFrame(columns=groups,
                                  data=[up, down], index=['Up-regulated', 'Down-regulated']).transpose()

    df = dysregulations

    df = df.reset_index()
    df = df.melt('index')
    df['color'] = np.where(df['value'] > 0, '#e3432d', '#167a9c')
    df['value'] = abs(df['value'])
    M = 2
    N = len(groups)
    fig, ax = plt.subplots()
    fig.set_figwidth(2)
    scatter = ax.scatter(x=df['variable'], y=df['index'], s=df['value'],
                         c=df['color'], ec='black', lw=0.5)
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    sns.despine()
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=45)
    kw = dict(prop="sizes", num=3, color='black', alpha=.6)
    ax.legend(*scatter.legend_elements(**kw), title="# Proteins", handleheight=2,
              bbox_to_anchor=(1, 1), loc="upper left", markerscale=1,
              edgecolor='white')
    plt.margins()
    if save is not None:
        if vector is True:
            plt.savefig(save + 'diff_reg.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'diff_reg.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def whole_network(self, labels=False, save=None, vector=True, dpi=300):
    """Network of entities differentially regulated for each
    group analyzed.

    Args:
        labels (bool, optional): Show graph labels. Defaults to False.
        save (str, optional): Path to save image. Defaults to None.
         Defaults to None.
         dpi (int, optional): Image resolution. Defaults to 300.
    """
    import matplotlib as mpl
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    plt.rcParams['figure.dpi'] = dpi
    data = copy(self)
    network_frame = []
    for group, df, color in zip(data.groups, data.original, data.colors):
        df['Experiment'] = group
        df = df[df[self.pvalue] < 0.05]
        df['Size'] = len(df)
        df['color'] = color
        network_frame.append(df)
    network_frame = pd.concat(network_frame)
    source = pd.DataFrame({'ID': network_frame['Experiment'],
                           'Size': network_frame['Size'],
                           'type': 'Experiment',
                           'color': network_frame['color']})
    source = source.drop_duplicates()
    norm = mpl.colors.TwoSlopeNorm(vmin=min(network_frame['log2(fc)']),
                                   vmax=max(network_frame['log2(fc)']),
                                   vcenter=0)
    cmap = cm.RdYlBu_r
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_hex = [mcolors.to_hex(m.to_rgba(x)) for x in network_frame['log2(fc)']]

    target = pd.DataFrame({'ID': network_frame['gene_name'],
                           'Size':  int(min(network_frame['Size'])*0.5),
                           'type': 'Protein',
                           'color': color_hex})
    target = target.drop_duplicates()
    edgelist = network_frame[['Experiment', 'gene_name']]
    G = nx.from_pandas_edgelist(edgelist,
                                source='Experiment',
                                target='gene_name',
                                create_using=nx.Graph)
    carac = pd.concat([source, target]).reset_index(drop=True)
    carac = carac.drop_duplicates('ID')
    carac = carac.set_index('ID')
    nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
    nx.set_node_attributes(G, dict(zip(carac.index, carac.Size)), name="Size")
    pos = nx.kamada_kawai_layout(G)
    carac = carac.reindex(G.nodes())
    nx.draw(G,
            pos=pos,
            node_color=carac['color'],
            node_size=carac['Size']/20,
            edgecolors='black',
            linewidths=0.4,
            alpha=0.9,
            width=0.2,
            edge_color='gray')
    if labels is True:
        nx.draw_networkx_labels(G, pos, font_size=6)
    if save is not None:
        nx.write_graphml(G, save + 'PPNetwork.graphml', named_key_ids=True)
        if vector is True:
            plt.savefig(save + 'PPNetwork.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'PPNetwork.dpi', dpi=dpi, bbox_inches='tight')
    plt.show()
    return (G)


def dotplot_enrichment(self, *Terms, top=5,  fig_height=None, palette='PuBu', save=None, vector=True,
                       dpi=300):
    """Dotplot Enrichment

    Dotplot to visualize together the enrichment data for each group

    Args:
        top (int, optional): Top N pathway to considered in each group. Defaults to 5.
        palette (str, optional): color palette. Defaults to 'PuBu'.
        fig_height (int, optional): User optionally can define figure height. Defaults to None
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg.
        Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    """
    plt.rcParams['figure.dpi'] = dpi
    if all(enr is None for enr in self.enrichment):
        raise IndexError('There is not Enrichment result in data!')
    genesets = [list(x.Gene_set.drop_duplicates()) for x in self.enrichment]
    genesets = pd.Series(sum(genesets, [])).drop_duplicates()
    for i in genesets:
        data = self.enrichment
        data = [x[x['Gene_set'] == i] for x in data]
        terms = [x.iloc[:top, :] for x in data]
        terms = [list(x['Term']) for x in terms]
        terms = pd.Series(sum(terms, [])).drop_duplicates()
        if len(Terms) > 0:
            terms = Terms
        data = [x[x['Term'].isin(terms)] for x in data]
        data = [x.assign(Group=y) for x, y in zip(data, self.groups)]
        data = pd.concat(data)
        data = data[['Term', 'Overlap', 'Adjusted P-value', 'Group']]
        data['Overlap'] = data.Overlap.str.split('/', regex=False).str[0]
        data['Overlap'] = data.Overlap.astype(int)
        data['-log10(p)'] = -np.log10(data['Adjusted P-value'])
        fig, ax = plt.subplots()
        if fig_height is not None:
            fig.set_figheight(fig_height)
        sns.set_style('white')
        ax = sns.scatterplot(data=data, x='Group', y='Term', size='-log10(p)',
                             hue='-log10(p)', palette=palette, sizes=(40, 280),
                             linewidth=0.5, edgecolor='black'
                             )
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        sns.despine()
        plt.margins()
        plt.ylabel('')
        if save is not None:
            if vector is True:
                plt.savefig(save + '_'+i+'_' + 'dotplot_enrichment.svg', bbox_inches='tight')
            else:
                plt.savefig(save + '_'+i+'_' + 'dotplot_enrichmnet.dpi', dpi=dpi, bbox_inches='tight')
        plt.show()


def protein_overlap(self, min_subset=10, face_color='darkcyan', shad_color="#f0f0f0",
                    edge_color='black', linewidth=1, save=None, vector=True, dpi=300):
    """Upset plot 
    Upset plot to evaluate protein overlap among groups.

    Args:
        min_subset (int, optional): Minimum overlap size to consider for upset plot.
         Defaults to 10.
        face_color (str, optional): Bar and dot colors. Defaults to 'darkcyan'.
        shad_color (str, optional): Shad color in the dot part of the graph.
         Defaults to "#f0f0f0".
        edge_color (str, optional): edge colors. Defaults to 'black'.
        linewidth (int, optional): line widths. Defaults to 1.
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg. Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    """
    plt.rcParams['figure.dpi'] = dpi
    from upsetplot import UpSet
    from upsetplot import from_contents
    plt.style.context('classic')
    plt.rcParams['grid.alpha'] = 0
    plt.rcParams['figure.dpi'] = dpi
    plt.rcParams['patch.linewidth'] = linewidth
    plt.rcParams['patch.edgecolor'] = 'black'
    plt.rcParams['patch.force_edgecolor'] = True

    data = copy(self)
    genes = []
    for i in data.group_data:
        genes.append(i['gene_name'].drop_duplicates())
    dictionary = dict(zip(data.groups, genes))
    upset = from_contents(dictionary)

    figure = UpSet(upset,  facecolor=face_color, shading_color=shad_color,
                   min_subset_size=min_subset, show_counts=True,
                   with_lines=True)
    for i in data.groups:
        figure.style_subsets(present=i, edgecolor=edge_color, linewidth=linewidth)
    figure.plot()
    if save is not None:
        if vector is True:
            plt.savefig(save + 'upset_proteins.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'upset_proteins.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def enrichment_overlap(self,  min_subset=1, face_color='darkcyan', shad_color="#f0f0f0",
                       edge_color='black', linewidth=1, save=None, vector=True, dpi=300):
    """Upset plot 
    Upset plot to evaluate enrichment terms overlap among groups.

    Args:
        min_subset (int, optional): minimum number of overlap size to consider
         for upset plot . Defaults to 10.
        face_color (str, optional): Bar and dot colors. Defaults to 'darkcyan'.
        shad_color (str, optional): Shad color in the dot part of the graph.
         Defaults to "#f0f0f0".
        edge_color (str, optional): edge colors. Defaults to 'black'.
        linewidth (int, optional): line widths. Defaults to 1.
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg. Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    Raises:
        IndexError: If there is no Enrichment data on .omics file.
    """
    from upsetplot import UpSet
    from upsetplot import from_contents
    plt.style.context('classic')
    plt.rcParams['grid.alpha'] = 0
    plt.rcParams['figure.dpi'] = dpi
    plt.rcParams['patch.linewidth'] = linewidth
    plt.rcParams['patch.edgecolor'] = 'black'
    plt.rcParams['patch.force_edgecolor'] = True

    data = copy(self)
    genes = []
    if all(enr is None for enr in self.enrichment):
        raise IndexError('There is not Enrichment result in data!')
    for i in data.enrichment:
        try:
            genes.append(i['Term'].drop_duplicates())
        except KeyError:
            df = pd.DataFrame(columns=['Gene_set', 'Term', 'Overlap', 'Adjusted P-value', 'Genes'])
            genes.append(df['Term'].drop_duplicates())
    dictionary = dict(zip(data.groups, genes))
    upset = from_contents(dictionary)

    figure = UpSet(upset,  facecolor=face_color, shading_color=shad_color, min_subset_size=min_subset, show_counts=True,
                   with_lines=True)
    for i in data.groups:
        figure.style_subsets(present=i, edgecolor=edge_color, linewidth=linewidth)
    figure.plot()
    if save is not None:
        if vector is True:
            plt.savefig(save + 'upset_pathways.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'upset_pathways.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def similarity_network(self, pvalue=1, parameter='TotalMean',
                       metric='jaccard',
                       similarity_cutoff=0.2,
                       save=None, vector=True,
                       dpi=300):
    """Similarity heatmap plot

    Pair-wise similarity comparison heatmap.

    Args:
        pvalue (int, optional): P-value threshold to proteins that
         OmicScope must consider for analysis. Defaults to 1.
        parameter (str, optional): Parameter to take into account in pairwise comparison.
         Defaults to 'TotalMean'. Optionally 'log2(fc)'.
        similarity_cutoff (float, optional): Cuttoff to consider edges between groups.
         Defaults to 0.5.
        metric (str, optional): statistical algorithm to perform pairwise comparison. Defaults to 'correlation'.
         Optionally, user can test other algorithm described in scipy.spatial.distance.
        center(float, optional): number to center the heatmap color gradient.
        palette (str, optional): color palette to plot heatmap.
         Defaults to 'RdYlBu'.
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg. Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    """
    plt.rcParams['figure.dpi'] = dpi
    palette = self.colors
    conditions = self.groups
    pval = self.pvalue
    data = copy(self)
    data1 = data.original
    totalMean = []
    colors = data.colors
    for i in data1:
        df = i.set_index('Accession')
        df = df[df[self.pvalue] < pvalue]
        df = df[[parameter]]
        totalMean.append(df)
    wholedata = pd.concat(totalMean, axis=1, join='outer')
    wholedata.columns = data.groups
    wholedata = wholedata.reset_index()
    wholedata = wholedata.melt(['Accession']).replace(np.nan, 0)
    corr = wholedata.pivot(index='Accession', columns='variable')
    from sklearn.metrics import pairwise_distances
    # Replace -inf to the lowest non-inf value in data
    corr = corr.replace(-np.inf,
                        corr.replace(-np.inf,
                                     np.nan).dropna().min().min())
    corr = 1-pd.DataFrame(pairwise_distances(corr.T.to_numpy(),
                                             metric=metric))
    corr[corr < similarity_cutoff] = 0
    corr.columns, corr.index = data.groups, data.groups
    G = nx.from_pandas_adjacency(corr)
    G.edges(data=True)
    size = [len(set(x[x[pval] < pvalue].gene_name)) for x in self.original]
    carac = pd.DataFrame(zip(conditions, palette, size),
                         columns=['ID', 'color', 'Size'])
    carac = carac.set_index('ID')
    nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
    nx.set_node_attributes(G, dict(zip(carac.index, carac.Size)), name="Size")
    pos = nx.kamada_kawai_layout(G, )
    carac = carac.reindex(G.nodes())
    pos = nx.kamada_kawai_layout(G, weight=None)
    edges = G.edges
    weights = [G[u][v]['weight'] for u, v in edges]
    weights = [round(x, 2) for x in weights]
    norm = [float(i)*5/np.mean(weights) for i in weights]
    G.remove_edges_from(nx.selfloop_edges(G))
    nx.draw(G,
            pos=pos,
            node_color=carac['color'],
            node_size=carac['Size'],
            edgecolors='black',
            linewidths=0.6,
            alpha=0.9,
            width=norm,
            edge_color='gray')
    nx.draw_networkx_labels(G, pos, font_size=6)
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=dict(zip(list(labels.keys()), weights)))
    if save is not None:
        nx.write_graphml(G, save + 'Similarity_network.graphml', named_key_ids=True)
        if vector is True:
            plt.savefig(save + 'Similarity_network.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'Similarity_network.dpi', dpi=dpi, bbox_inches='tight')
    plt.show()


def similarity_heatmap(self, pvalue=1, parameter='TotalMean',
                       metric='correlation',
                       center=0, palette='RdYlBu_r',
                       annotation=True, save=None, vector=True,
                       dpi=300):
    """Similarity heatmap plot

    Pair-wise similarity comparison heatmap.

    Args:
        pvalue (int, optional): P-value threshold to proteins that
         OmicScope must consider for analysis. Defaults to 1.
        parameter (str, optional): Parameter to take into account in pairwise comparison.
        Defaults to 'TotalMean'. Optionally 'log2(fc)'.
        metric (str, optional): statistical algorithm to perform pairwise comparison. Defaults to 'correlation'.
         Optionally, user can test other algorithm described in scipy.spatial.distance.
        center(float, optional): number to center the heatmap color gradient.
        palette (str, optional): color palette to plot heatmap.
         Defaults to 'RdYlBu'.
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg. Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    """
    plt.rcParams['figure.dpi'] = dpi
    data = copy(self)
    data1 = data.original
    totalMean = []
    colors = data.colors
    for i in data1:
        df = i.set_index('Accession')
        df = df[df[self.pvalue] < pvalue]
        df = df[[parameter]]
        totalMean.append(df)
    wholedata = pd.concat(totalMean, axis=1, join='outer')
    wholedata.columns = data.groups
    wholedata = wholedata.reset_index()
    wholedata = wholedata.melt(['Accession']).replace(np.nan, 0)
    corr = wholedata.pivot(index='Accession', columns='variable')
    from sklearn.metrics import pairwise_distances
    # Replace -inf to the lowest non-inf value in data
    corr = corr.replace(-np.inf,
                        corr.replace(-np.inf,
                                     np.nan).dropna().min().min())
    # Calculating distance
    corr = 1-pd.DataFrame(pairwise_distances(corr.T.to_numpy(),
                                             metric=metric))
    corr.columns, corr.index = data.groups, data.groups
    annot = copy(corr)
    if annotation is False:
        annot[annot != -2] = np.nan
    annot[annot == 1] = np.nan
    sns.clustermap(corr, cmap=palette, center=center, annot=annot,
                   mask=annot.isnull(),
                   col_colors=colors, row_colors=colors)
    if save is not None:
        if vector is True:
            plt.savefig(save + 'similarity_heatmap.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'similarity_heatmap.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def overlap_fisher(group1, group2, union):
    """Perform a pair-wise comparison based on hypergeometric
    distribution.

    Args:
        group1 (Series, pandas): Column condition 1
        group2 (Series, pandas): Column condition 2
        union (int): number of whole entities evaluated in the study
        among all conditions

    Returns:
        Pvalue (float): P-value
    """
    import numpy as np
    from scipy.stats import hypergeom

    deps1 = set(group1)
    deps2 = set(group2)
    union = len(union)
    intersection = len(deps1.intersection(deps2))
    distribution = min([len(deps1), len(deps2)])
    [M, n, N] = [union, len(deps1), len(deps2)]
    rv = hypergeom(M, n, N)
    x = np.arange(intersection, distribution)
    pval = sum(rv.pmf(x))
    return pval


def fisher_heatmap(self, palette='Spectral', pvalue=0.05,
                   annotation=True,
                   save=None, vector=True, dpi=300):
    """ Heatmap according to statistical significance

    Perform a pair-wise comparison of all conditions
    based on hypergeometric distribution and plot a heatmap
    with hierarchical clustering. In the plot, it the p-value is 
    labeled in a log10-transformation.

    Args:
        palette (str, optional): color palette. Defaults to 'Spectral'.
        pvalue (float, optional): P-value. Defaults to 0.05.
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg.
        Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.

    Returns:
        DataFrame: Dataframe of pvalues between each condition.
    """
    import numpy as np
    import pandas as pd
    from scipy.spatial.distance import pdist
    from scipy.spatial.distance import squareform
    plt.rcParams['figure.dpi'] = dpi
    conditions = self.groups
    colors = self.colors
    pval = self.pvalue
    union = set(pd.concat(self.original).gene_name)
    sets = [set(x[x[pval] < pvalue].gene_name) for x in self.original]
    sets = pd.DataFrame(sets)
    matrix = pdist(sets, lambda u, v: overlap_fisher(u, v, union=union))
    matrix = squareform(matrix)
    matrix = pd.DataFrame(matrix, columns=conditions, index=conditions)
    annot = -np.log10(matrix)
    if annotation is False:
        annot[annot != np.inf] = np.nan
    annot[annot == np.inf] = np.nan
    sns.clustermap(matrix, cmap=palette, annot=annot,
                   mask=annot.isnull(),
                   col_colors=colors, row_colors=colors)
    if save is not None:
        if vector is True:
            plt.savefig(save + 'overlap_stat.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'overlap_stat.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def fisher_network(self, protein_pvalue=0.05, graph_pvalue=0.1, save=None, vector=True, dpi=300):
    """ Network for pair-wise comparison.

    Network of all groups analysed, linking groups based on
    statistical results from hypergeometric distribution.

    Args:
        protein_pvalue (float, optional): P-value cuttoff for proteins (e.g. differentially regulated).
         Defaults to 0.05.
        graph_pvalue (float, optional): P-value cuttoff for fisher's test and networks. Defaults to 0.1.
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg.
        Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.

    Returns:
        Graph (Networkx.G): Networkx object
    """
    from scipy.spatial.distance import pdist
    from scipy.spatial.distance import squareform
    plt.rcParams['figure.dpi'] = dpi
    palette = self.colors
    conditions = self.groups
    pval = self.pvalue
    union = set(pd.concat(self.original).gene_name)
    sets = [set(x[x[pval] < protein_pvalue].gene_name) for x in self.original]
    sets = pd.DataFrame(sets)
    matrix = pdist(sets, lambda u, v: overlap_fisher(u, v, union=union))
    matrix = squareform(matrix)
    matrix = pd.DataFrame(matrix, columns=conditions, index=conditions)
    matrix[matrix > graph_pvalue] = 0
    G = nx.from_pandas_adjacency(matrix)
    G.edges(data=True)
    size = [len(set(x[x[pval] < protein_pvalue].gene_name)) for x in self.original]
    carac = pd.DataFrame(zip(conditions, palette, size),
                         columns=['ID', 'color', 'Size'])
    carac = carac.set_index('ID')
    nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
    nx.set_node_attributes(G, dict(zip(carac.index, carac.Size)), name="Size")
    pos = nx.kamada_kawai_layout(G, )
    carac = carac.reindex(G.nodes())
    pos = nx.kamada_kawai_layout(G, weight=None)
    edges = G.edges
    weights = [G[u][v]['weight'] for u, v in edges]
    weights = -np.log10(weights)
    weights = [round(x, 2) for x in weights]
    norm = [float(i)*5/np.mean(weights) for i in weights]
    nx.draw(G,
            pos=pos,
            node_color=carac['color'],
            node_size=carac['Size'],
            edgecolors='black',
            linewidths=0.6,
            alpha=0.9,
            width=norm,
            edge_color='gray')
    nx.draw_networkx_labels(G, pos, font_size=6)
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=dict(zip(list(labels.keys()), weights)))
    if save is not None:
        nx.write_graphml(G, save + 'PPNetwork.graphml', named_key_ids=True)
        if vector is True:
            plt.savefig(save + 'groupNetwork.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'groupNetwork.dpi', dpi=dpi, bbox_inches='tight')
    plt.show()


def circular_path(self, Term, protein_cutoff=0.05, save=None, vector=True):
    """Circular Path

    For a determined Term, Circular path links all groups that
    enriched for that term to respective differentially regulated protein.
    To run this code, R and circlize package must be installed in the system.

    Args:
       Term (str): Terms that must be shown. Defaults to 0.05.
       save (str, optional): Path to save image. Defaults to None.
       vector (bool, optional): If image should be export as .svg.
       Defaults to True.

   """
    from .circlize import circlize
    from .circlize import color_matrix
    from .circlize import deps
    from .circlize import deps_matrix
    from .circlize import enrichment_filtering
    if all(enr is None for enr in self.enrichment):
        raise IndexError('There is not Enrichment result in data!')
    data = self
    colors = self.colors
    df = deps(data, pvalue=protein_cutoff)
    terms = enrichment_filtering(data, Term)
    deps = df[df['gene_name'].isin(terms)]
    deps = deps.set_index('gene_name')
    matrix = deps_matrix(deps)
    colmat = color_matrix(deps, colors)
    labels = data.groups
    circlize(matrix, colmat, colors, labels, save=save, vector=vector)
