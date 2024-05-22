from copy import copy

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from pycirclize import Circos


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
    ax.bar(r, deps, color=colors, edgecolor='black', width=0.5, linewidth=1)
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
        df = df[df[self.pvalue] <= 0.05]
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
    enrichments = []
    groups = []
    for g,e in zip(copy(self.groups), copy(self.enrichment)):
        if e is not None:
            enrichments.append(e)
            groups.append(g)
    genesets = [list(x.Gene_set.drop_duplicates()) for x in enrichments]
    genesets = pd.Series(sum(genesets, [])).drop_duplicates()
    for i in genesets:
        data = enrichments
        data = [x[x['Gene_set'] == i] for x in data]
        terms = [x.iloc[:top, :] for x in data]
        terms = [list(x['Term']) for x in terms]
        terms = pd.Series(sum(terms, [])).drop_duplicates()
        if len(Terms) > 0:
            terms = Terms
        data = [x[x['Term'].isin(terms)] for x in data]
        data = [x.assign(Group=y) for x, y in zip(data, groups)]
        data = pd.concat(data)
        data = data[['Term', 'Overlap', 'Adjusted P-value', 'Group']]
        data['Overlap'] = data.Overlap.str.split('/', regex=False).str[0]
        data['Overlap'] = data.Overlap.astype(int)
        data['-log10(p)'] = -np.log10(data['Adjusted P-value'])
        fig, ax = plt.subplots()
        if fig_height is not None:
            fig.set_figheight(fig_height)
        sns.set_style('white')
        sns.scatterplot(data=data, x='Group', y='Term', size='-log10(p)',
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

    data_original = copy(self)
    genes = []
    if all(enr is None for enr in self.enrichment):
        raise IndexError('There is not Enrichment result in data!')
    enrichment = []
    groups = []
    for e, g in zip(data_original.enrichment, data_original.groups):
        if e is not None:
            enrichment.append(e)
            groups.append(g)
    data = copy(self)
    data.enrichment = enrichment
    data.groups = groups
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


def similarity_network(self, pvalue=1, comparison_param='log2(fc)',
                       metric='jaccard',
                       absolute_similarity_cutoff=0.2,
                       save=None, vector=True,
                       dpi=300):
    """Similarity Network plot

        Perform a pairwise correlation analysis and create a graph where groups
        are depicted as nodes, and pairwise similarity indices serve as edges.
        In order to establish a connection between two groups, the function filters
        edges based on an absolute similarity cutoff, excluding edges that fall within
        a specified interval range, for instance, -0.2 to 0.2, when the
        absolute_similarity_cutoff is set to 0.2.

        Furthermore, when utilizing the Jaccard similarity index,
        this function takes into account the shared 'gene_name' between groups.
        In contrast, for the other available options, the function considers either
        'TotalMean' or 'log2(fc)' columns

    Args:
        pvalue (int, optional): P-value threshold to proteins that
         OmicScope must consider for analysis. Defaults to 1.
        comparison_param (str, optional): Parameter/column to take into account in pairwise comparison.
         Defaults to 'log2(fc)'. Optionally 'TotalMean'.
        absolute_similarity_cutoff (float, optional): Cuttoff to consider the links between groups. Since major
         similarity indexes have positive and negative values, the function expect an absolute value to perform cuttof.
         Defaults to 0.2 (which means -0.2 < cutoff< 0.2).
        metric (str, optional): algorithm to perform pairwise comparison. Defaults to 'correlation'.
         Optionally, user can test other algorithm described in scipy.spatial.distance.
        center(float, optional): number to center the heatmap color gradient.
        palette (str, optional): color palette to plot heatmap.
         Defaults to 'RdYlBu'.
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg. Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    """
    from copy import copy
    plt.rcParams['figure.dpi'] = dpi
    palette = self.colors
    conditions = self.groups
    pval = self.pvalue
    data = copy(self)
    data1 = data.original
    totalMean = []
    for i in data1:
        df = i.groupby('gene_name').mean()
        df = df[df[self.pvalue] <= pvalue]
        df = df[[comparison_param]]
        totalMean.append(df)
    wholedata = pd.concat(totalMean, axis=1, join='outer')
    wholedata.columns = data.groups
    corr = wholedata
    from sklearn.metrics import pairwise_distances

    # Replace -inf to the lowest non-inf value in data
    corr = corr.replace(-np.inf,
                        corr.replace(-np.inf,
                                     np.nan).dropna().min().min())

    # Inicializar um DataFrame vazio para armazenar as distâncias
# Inicializar um DataFrame vazio para armazenar as distâncias
    df_dist = pd.DataFrame(index=corr.columns, columns=corr.columns)
    # Calculate the distance between each group
    for col1 in corr.columns:
        for col2 in corr.columns:
            if col1 != col2:
                # Jaccard uses the index for the computation
                if metric == 'jaccard':
                    set1 = set(corr[col1].dropna().index)
                    set2 = set(corr[col2].dropna().index)
                    distance = len(set1.intersection(set2))/len(set1.union(set2))
                    df_dist.at[col1, col2] = 1-distance
                else:
                    # For other methods, we use the non-zero quantitative values from each condition
                    distance = corr[[col1, col2]].dropna(axis=0)
                    distance = pairwise_distances(distance.T.to_numpy(),
                                                  metric=metric, force_all_finite='allow-nan')
                    df_dist.at[col1, col2] = distance[0][1]
            df = 1-df_dist
    corr = df
    # Perform the absoltute_cuttoff for values
    corr = df.applymap(lambda x: 0 if -absolute_similarity_cutoff <= x <= absolute_similarity_cutoff else x)
    corr.columns, corr.index = data.groups, data.groups
    corr = corr.replace(np.nan, 0)
    # Plotting graph
    G = nx.from_pandas_adjacency(corr)
    G.edges(data=True)
    size = [len(set(x[x[pval] <= pvalue].gene_name)) for x in self.original]
    carac = pd.DataFrame(zip(conditions, palette, size),
                         columns=['ID', 'color', 'Size'])
    carac = carac.set_index('ID')
    nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
    nx.set_node_attributes(G, dict(zip(carac.index, carac.Size)), name="Size")
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


def similarity_heatmap(self, pvalue=1, comparison_param='log2(fc)',
                       metric='correlation',
                       center=0, palette='RdYlBu_r',
                       annotation=True, save=None, vector=True,
                       dpi=300):
    """Similarity heatmap plot
        Perform a pair-wise similarity analysis and plot a heatmap.

        When utilizing the Jaccard similarity index,
        this function takes into account the shared 'gene_name' between groups.
        In contrast, for the other available options, the function considers either
        'TotalMean' or 'log2(fc)' columns

    Args:
        pvalue (int, optional): P-value threshold to proteins that
         OmicScope must consider for analysis. Defaults to 1.
        comparison_param (str, optional): Parameter to take into account in pairwise comparison.
        Defaults to 'log2(fc)'. Optionally 'TotalMean'.
        metric (str, optional): algorithm to perform pairwise comparison. Defaults to 'correlation'.
         Optionally, user can test other algorithm described in scipy.spatial.distance.
        center(float, optional): number to center the heatmap color gradient.
        palette (str, optional): color palette to plot heatmap.
         Defaults to 'RdYlBu'.
        save (str, optional): Path to save image. Defaults to None.
        vector (bool, optional): If image should be export as .svg. Defaults to True.
        dpi (int, optional): Image resolution. Defaults to 300.
    """
    from copy import copy
    plt.rcParams['figure.dpi'] = dpi
    data = copy(self)
    data1 = data.original
    totalMean = []
    colors = data.colors
    for i in data1:
        df = i.groupby('gene_name').mean()
        df = df[df[self.pvalue] <= pvalue]
        df = df[[comparison_param]]
        totalMean.append(df)
    wholedata = pd.concat(totalMean, axis=1, join='outer')
    wholedata.columns = data.groups
    corr = wholedata
    from sklearn.metrics import pairwise_distances

    # Replace -inf to the lowest non-inf value in data
    corr = corr.replace(-np.inf,
                        corr.replace(-np.inf,
                                     np.nan).dropna().min().min())

    # Inicializar um DataFrame vazio para armazenar as distâncias
    # Inicializar um DataFrame vazio para armazenar as distâncias
    df_dist = pd.DataFrame(index=corr.columns, columns=corr.columns)
    # Calculate the distance between each group
    for col1 in corr.columns:
        for col2 in corr.columns:
            if col1 != col2:
                # Jaccard uses the index for the computation
                if metric == 'jaccard':
                    set1 = set(corr[col1].dropna().index)
                    set2 = set(corr[col2].dropna().index)
                    distance = len(set1.intersection(set2))/len(set1.union(set2))
                    df_dist.at[col1, col2] = 1-distance
                else:
                    # For other methods, we use the non-zero quantitative values from each condition
                    distance = corr[[col1, col2]].dropna(axis=0)
                    distance = pairwise_distances(distance.T.to_numpy(),
                                                  metric=metric, force_all_finite='allow-nan')
                    df_dist.at[col1, col2] = distance[0][1]
            df = 1-df_dist
    corr = df
    # Perform the absoltute_cuttoff for values
    corr.columns, corr.index = data.groups, data.groups
    corr = corr.replace(np.nan, 0)
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
    from scipy.stats import hypergeom

    deps1 = set(group1)
    deps2 = set(group2)
    intersection = len(deps1.intersection(deps2))
    distribution = min([len(deps1), len(deps2)])
    [M, n, N] = [union, len(deps1), len(deps2)]
    rv = hypergeom(M, n, N)
    x = np.arange(intersection, distribution)
    pval = sum(rv.pmf(x))
    return pval

def distribution_test(self, protein_pvalue, 
                   method):
    """This function performs a statistical analysis on protein data considering overlaps between groups. The function performs different statistical tests depending on the chosen method:
        1. t-test (ttest): This test is used to compare the means of two groups assuming normally distributed data.
        2. Wilcoxon signed-rank test (wilcoxon): This non-parametric test is used to compare two related groups when the data may not be normally distributed.
        3. Kolmogorov-Smirnov test (ks): This test is used to compare the probability distributions of two samples.
    
    Args:
        protein_pvalue (float): The cut-off value for protein p-values.
        method (str):  The statistical method to be used for comparison. Valid options include "ttest" (t-test), "wilcoxon" (Wilcoxon signed-rank test), and "ks" (Kolmogorov-Smirnov test).
        comparison among groups considering t-test (for parametric distributions),
        wilcoxon (for non-parametric distributions), and ks (kolmorov-smirnov test).

    Returns:
        matrix (DataFrame, pandas): P-value
    """
    from scipy.stats import ttest_ind, wilcoxon, kstest
    from scipy.spatial.distance import squareform
    import itertools
    protein_pvalue = protein_pvalue
    stat = method

    conditions = self.groups

    data = [i[i[self.pvalue]<=protein_pvalue] for i in self.original]
    data = [i.set_index('Accession') for i in data]
    data = [i['log2(fc)']for i in data]
    data = [i.replace(-np.inf, np.nan) for i in data]
    data = [i.replace(np.nan, min(i)-1) for i in data]


    pair_data = list(itertools.combinations(data,2))

    overlap = [ set(i[0].index) & set(i[1].index) for i in pair_data]
    pair_data = [ (i[0][i[0].index.isin(j)], i[1][i[1].index.isin(j)]  ) for i,j in zip(pair_data,overlap)]

    if stat == 'ttest':
        stats = [ttest_ind(i[0], i[1])[1] for i in pair_data]
    if stat == 'wilcoxon':
        stats = [wilcoxon(i[0], i[1])[1] for i in pair_data]
    if stat == 'ks':
        stats = []
        for i in pair_data:
            try:
                stats.append(kstest(i[0], i[1])[1])
            except:
                stats.append(0)
    matrix = squareform(stats)
    matrix = pd.DataFrame(matrix, columns=conditions, index=conditions)
    return matrix

def fisher_test(self, protein_pvalue,
                   background_lenght):
    """This function performs a pair-wise statistical analysis using Fisher's exact test.
    Fisher's exact test is a statistical test used to compare two nominal variables from two samples. 
    In this context, it's used to compare the proportions of proteins with significant p-values
    (determined by protein_pvalue) between groups.

    Args:
        protein_pvalue (float): The cut-off value for protein p-values.
        background_lenght (float, optional): The total number of entities in the background set (optional).
        If not provided, all genes from the original data are used as the background (Recommended).
    Returns:
        matrix (DataFrame, pandas): P-value
    """
    from scipy.spatial.distance import squareform, pdist
    conditions = self.groups
    pval = self.pvalue
    if background_lenght is None:
        union = set(pd.concat(self.original).gene_name)
        union = len(union)
    else:
        union = background_lenght
    sets = [set(x[x[pval] <= protein_pvalue].gene_name) for x in self.original]
    sets = pd.DataFrame(sets)
    matrix = pdist(sets, lambda u, v: overlap_fisher(u, v, union=union))
    matrix = squareform(matrix)
    matrix = pd.DataFrame(matrix, columns=conditions, index=conditions)
    return matrix

def stat_matrix(self, method, protein_pvalue, background_lenght):
    """
  Performs a pair-wise statistical comparison between groups based on the 
  chosen method and a protein p-value cut-off.

  Args:
      self: Reference to the class instance where this function is called.
      method (str): The statistical method to be used for the comparison.
          Valid options include:
              * "fisher": Performs Fisher's exact test, suitable for comparing 
                  proportions of significant proteins between groups.
              * "ttest": Performs a t-test, assuming normally distributed data.
              * "wilcoxon": Performs a Wilcoxon signed-rank test, a non-parametric 
                  alternative for comparing related groups when normality cannot 
                  be assumed.
              * "ks": Performs a Kolmogorov-Smirnov test, used to compare the 
                  probability distributions of two samples.
      protein_pvalue (float): The cut-off value for protein p-values. This value 
          determines which proteins are considered significant based on a 
          previous analysis.
      background_lenght (int, optional): The total number of entities in the 
          background set. This argument is only used for the "fisher" method 
          to define the population size. If not provided, all genes from the 
          original data are used as the background.

  Returns:
      pandas.DataFrame: A DataFrame containing the p-values for each pair-wise 
          comparison between groups.
  """
    if method == 'fisher':
        matrix = fisher_test(self, protein_pvalue=protein_pvalue, 
                    background_lenght=background_lenght)
    elif method in ['ttest', 'wilcoxon', 'ks']:
        matrix = distribution_test(self, protein_pvalue=protein_pvalue,
                                   method=method)
    else:
        raise ValueError('''Please, verify if it was selected a valid method ("fisher", "ttest", "wilcoxon", or "ks").''')
    return matrix

def stat_network(self, method='fisher', protein_pvalue=0.05, background_lenght=None, dpi=300,
                 graph_pvalue=0.1,
                 save=None, vector=True):
    """
  Generates a network visualization based on statistical comparisons 
  between groups imported on Nebula

  Args:
      method (str, optional): The statistical method used for comparison 
          between groups. Valid options include:
              * "fisher": Performs Fisher's exact test, suitable for comparing 
                  proportions of significant proteins.
              * "ttest": Performs a t-test, assuming normally distributed data.
              * "wilcoxon": Performs a Wilcoxon signed-rank test, a non-parametric 
                  alternative for related groups when normality cannot be assumed.
              * "ks": Performs a Kolmogorov-Smirnov test, used to compare the 
                  probability distributions of two samples.
          Defaults to "fisher".
      protein_pvalue (float, optional): The cut-off value for protein p-values. 
          This value determines which proteins are considered significant 
          based on a previous analysis. Defaults to 0.05.
      background_lenght (int, optional): The total number of entities in the 
          background set. This argument is only used for the "fisher" method 
          to define the population size. If not provided, all genes from the 
          original data are used as the background. Defaults to None (Recommended).
      dpi (int, optional): The resolution (dots per inch) for the generated 
          plot. Defaults to 300.
      graph_pvalue (float, optional): The threshold for p-values to consider
          the links between network nodes. Edges with 
          p-values greater than (for "fisher") or less than (for other methods) to 
          this value will be excluded from the network.
          Defaults to 0.1.
      save (str, optional): The filename prefix to save the network image 
          (e.g., "groupNetwork"). If provided, the function will save both 
          the GraphML representation of the network and the plot image.
      vector (bool, optional): If True (default), saves the image as an SVG 
          file (scalable vector graphics) suitable for high-quality printing. 
          If False, saves the image as a PNG file.

  Returns:
      None. This function generates a network visualization and potentially 
          saves image files, but it doesn't return any data.
  """
    import matplotlib.pyplot as plt
    plt.rcParams['figure.dpi'] = dpi
    palette = self.colors
    conditions = self.groups
    pval = self.pvalue
    matrix = stat_matrix(self, method=method,
                         protein_pvalue=protein_pvalue,
                         background_lenght=background_lenght)
    if method == 'fisher':
        matrix[matrix >= graph_pvalue] = 0
    elif method in ['ttest', 'wilcoxon','ks']:
        matrix[matrix <= graph_pvalue] = 0
    G = nx.from_pandas_adjacency(matrix)
    G.edges(data=True)
    size = [len(set(x[x[pval] <= protein_pvalue].gene_name)) for x in self.original]
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
    if method == 'fisher':
            weights = -np.log10(weights)
    weights = [round(x, 2) for x in weights]
    G.remove_edges_from(nx.selfloop_edges(G))
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

def stat_heatmap(self, palette='Spectral', method='fisher', pvalue=0.05,
                   background_lenght=None, annotation=True,
                   save=None, vector=True, dpi=300):
    """
    Generates a heatmap visualization to represent the p-values from 
    pair-wise statistical comparisons between groups.

    Args:
        palette (str, optional): The color palette to use for the heatmap. 
            Defaults to "Spectral".
        method (str, optional): The statistical method used for comparison 
            between groups. Valid options include:
                * "fisher": Performs Fisher's exact test, suitable for comparing 
                    proportions of significant proteins.
                * "ttest": Performs a t-test, assuming normally distributed data.
                * "wilcoxon": Performs a Wilcoxon signed-rank test, a non-parametric 
                    alternative for related groups when normality cannot be assumed.
                * "ks": Performs a Kolmogorov-Smirnov test, used to compare the 
                    probability distributions of two samples.
            Defaults to "fisher".
        pvalue (float, optional): The cut-off value for protein p-values. 
            This value determines which proteins are considered significant 
            based on a previous analysis. Defaults to 0.05.
        background_lenght (int, optional): The total number of entities in the 
            background set. This argument is only used for the "fisher" method 
            to define the population size. If not provided, all genes from the 
            original data are used as the background. Defaults to None.
        annotation (bool, optional): If True (default), displays the p-values 
            within each heatmap cell. If False, hides the p-value annotations.
        save (str, optional): The filename prefix to save the heatmap image 
            (e.g., "overlap_stat"). If provided, the function will save the plot 
            image.
        vector (bool, optional): If True (default), saves the image as an SVG 
            file (scalable vector graphics) suitable for high-quality printing. 
            If False, saves the image as a PNG file.
        dpi (int, optional): The resolution (dots per inch) for the generated 
            plot. Defaults to 300.

    Returns:
        None. This function generates a heatmap visualization and potentially 
            saves an image file, but it doesn't return any data.
    """
    plt.rcParams['figure.dpi'] = dpi
    colors = self.colors
    matrix = stat_matrix(self, method=method,
                         protein_pvalue=pvalue,
                         background_lenght=background_lenght)
    annot = matrix.copy()
    if annotation is False:
        annot[annot != np.inf] = np.nan
    annot[annot == 0] = np.nan
    sns.clustermap(matrix, cmap=palette, annot=annot,
                   mask=annot.isnull(),
                   col_colors=colors, row_colors=colors)
    plt.title('Pvalue')
    if save is not None:
        if vector is True:
            plt.savefig(save + 'overlap_stat.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'overlap_stat.png', dpi=dpi, bbox_inches='tight')
    plt.show()



def linkproteins(deps, groups):
    # retrieving overlapped proteins among groups
    overlapped_proteins = pd.concat(deps)
    overlapped_proteins = list(overlapped_proteins[overlapped_proteins.duplicated(
        subset=['gene_name'], keep=False)]['gene_name'].dropna())
    # ordering matrix
    matrixes = [x.assign(duplicated=lambda x: x.gene_name.isin(
        overlapped_proteins)) for x in deps]
    matrixes = [x.sort_values(
        ['duplicated', 'log2(fc)'], ascending=False) for x in matrixes]
    matrixes = [x.reset_index(drop=True) for x in matrixes]
    # Retrieving links
    dataframes = copy(matrixes)
    result = pd.DataFrame(
        columns=['gene_name', 'query_chr', 'query_start', 'ref_chr', 'ref_start'])
    # Perform pairwise comparison between DataFrames
    for i in range(len(dataframes)-1):
        for j in range(i+1, len(dataframes)):
            common_names = set(dataframes[i]['gene_name']).intersection(
                dataframes[j]['gene_name'])
            for name in common_names:
                index1 = dataframes[i][dataframes[i]
                                       ['gene_name'] == name].index[0]
                index2 = dataframes[j][dataframes[j]
                                       ['gene_name'] == name].index[0]
                result = result.append({'gene_name': name, 'query_chr': i, 'query_start': index1, 'ref_chr': j, 'ref_start': index2},
                                       ignore_index=True)
    result['query_end'], result['ref_end'] = result['query_start'] + \
        1, result['ref_start']+1
    group_dict = dict(enumerate(groups))
    result['query_chr'] = result.query_chr.replace(group_dict)
    result['ref_chr'] = result.ref_chr.replace(group_dict)
    return result, matrixes


def linkenrichment(enrichment, groups, number_deps):
    enr_original = copy(enrichment)
    group_original = copy(groups)
    enrichment = []
    group = []
    ndeps = []
    for e, g, n in zip(enr_original, group_original, number_deps):
        if e is not None:
            enrichment.append(e)
            group.append(g)
            ndeps.append(n)
    enrichment = [x.assign(db_term=x['Gene_set']+'-'+x['Term'])
                  for x in enrichment]
    overlapped_enrichment = pd.concat(enrichment)
    overlapped_enrichment['db_term'] = overlapped_enrichment['Gene_set'] + \
        '.'+overlapped_enrichment['Term']
    overlapped_enrichment = list(overlapped_enrichment[overlapped_enrichment.duplicated(
        subset=['db_term'], keep=False)]['db_term'].dropna())
    # ordering matrix
    matrixes = [x.assign(duplicated=lambda x: x.db_term.isin(
        overlapped_enrichment)) for x in enrichment]
    indexes = [list(np.random.randint(0, y, size=len(x)))
               for x, y in zip(matrixes, ndeps)]
    matrixes = [x.set_index(pd.Index(y)) for x, y in zip(matrixes, indexes)]
    # Retrieving links
    dataframes = copy(matrixes)
    result = pd.DataFrame(
        columns=['db_term', 'query_chr', 'query_start', 'ref_chr', 'ref_start'])
    # Perform pairwise comparison between DataFrames
    for i in range(len(dataframes)-1):
        for j in range(i+1, len(dataframes)):
            common_names = set(dataframes[i]['db_term']).intersection(
                dataframes[j]['db_term'])
            for name in common_names:
                index1 = dataframes[i][dataframes[i]
                                       ['db_term'] == name].index[0]
                index2 = dataframes[j][dataframes[j]
                                       ['db_term'] == name].index[0]
                result = result.append({'db_term': name, 'query_chr': i, 'query_start': index1, 'ref_chr': j, 'ref_start': index2},
                                       ignore_index=True)
    result['query_end'], result['ref_end'] = result['query_start'] + \
        1, result['ref_start']+1
    group_dict = dict(enumerate(group))
    result['query_chr'] = result.query_chr.replace(group_dict)
    result['ref_chr'] = result.ref_chr.replace(group_dict)
    return result


def circos_plot(self, vmax=1, vmin=-1, colormap='RdYlBu_r', colorproteins='darkcyan',
                colorenrichment='black',
                linewidth_heatmap=0.1, save=None, vector=True, dpi=300):
    """Circos plot
     This plot offers an overview of proteins differentially regulated
     between groups using circular plots.

    Args:
        vmin (int, optional): minimum value for foldchange. Defaults to -1.
        vmax (int, optional): maximum value for foldchange. Defaults to 1.
        colormap (str, optional): Colormap for heatmap. Defaults to 'RdBu_r'.
        colorproteins (str, optional): Color for protein links. Defaults to 'darkcyan'.
        colorenrichment (str, optional): Color for enrichment links. Defaults to 'black'.
        save (str, optional): Path to save file. Defaults to None.
        vector (bool, optional): Save as svg extension, if False, save as png. Defaults to True.
        dpi (int, optional): Figure resolution. Defaults to 300.
    """
    # Data
    enrichment = copy(self.enrichment)
    groups = copy(self.groups)
    deps = copy(self.group_data)
    deps = [x.dropna().reset_index(drop=True) for x in deps]
    colors = self.colors
    grouplen = [len(x) for x in deps]
    higher_group = max(grouplen)

    # retrieving links and matrixes
    links, matrixes = linkproteins(deps, groups)

    # Mapping heatmap
    matrixes = [y[['log2(fc)']].applymap(
        lambda x: vmax if x > vmax else x) for y in matrixes]
    matrixes = [y[['log2(fc)']].applymap(
        lambda x: vmin if x < vmin else x) for y in matrixes]
    matrixes = [x[['log2(fc)']].T.to_numpy() for x in matrixes]
    # Config circos
    sectors = dict(zip(groups, grouplen))
    sector_colors = dict(zip(groups, colors))
    circos = Circos(sectors, space=5)

    for sector, matrix in zip(circos.sectors, matrixes):
        # Outer name
        # Tamanho da primeira Track, contendo as cores dos grupos
        outer_track = sector.add_track((95, 110))
        # Marcar o nome dos grupos
        outer_track.text(sector.name, color="Black")
        # Outer Track
        # Tamanho da primeira Track, contendo as cores dos grupos
        outer_track = sector.add_track((88, 90))
        outer_track.axis(fc=sector_colors[sector.name])  # cores dos grupos
        outer_track.xticks_by_interval(interval=int(
            higher_group/20), label_orientation="vertical")  # colocar xticks nas tracks
        # foldchange track
        rect_track = sector.add_track((80, 85))
        rect_track.heatmap(matrix, cmap=colormap,
                           rect_kws=dict(ec="black", lw=linewidth_heatmap))

    # drawing enrichment links if applicable
    if len(enrichment) > 1:
        linkenr = linkenrichment(enrichment, groups, grouplen)
        for i in linkenr.to_dict('records'):
            region1 = (i['query_chr'], i['query_start'], i['query_end'])
            region2 = (i['ref_chr'], i['ref_start'], i['ref_end'])
            circos.link(region1, region2, color=colorenrichment)
    # drawing links
    for i in links.to_dict('records'):
        region1 = (i['query_chr'], i['query_start'], i['query_end'])
        region2 = (i['ref_chr'], i['ref_start'], i['ref_end'])
        circos.link(region1, region2, color=colorproteins)
    circos.colorbar(vmin=vmin, vmax=vmax, cmap=colormap,
                    colorbar_kws=dict(label="log2(FoldChange)"))
    fig = circos.plotfig()
    if save is not None:
        if vector is True:
            fig.savefig(save + 'circos.svg')
        else:
            fig.savefig(save + 'circos.png', dpi=dpi)


def circular_term(self, *Terms, pvalue=0.05, vmin=-1, vmax=1, colormap='RdBu_r',
                  label_size=12, save=None, vector=True, dpi=300):
    """Circular term
        Allows the visualization of all proteins related to a pre-specified term.
        This term is extracted from enrichment data.

    Args:
        pvalue (float, optional): Pvalue to consider differentially regulated proteins
         . Defaults to 0.05.
        vmin (int, optional): minimum value for foldchange. Defaults to -1.
        vmax (int, optional): maximum value for foldchange. Defaults to 1.
        colormap (str, optional): Colormap for heatmap. Defaults to 'RdBu_r'.
        save (str, optional): Path to save file. Defaults to None.
        vector (bool, optional): Save as svg extension, if False, save as png. Defaults to True.
        dpi (int, optional): Figure resolution. Defaults to 300.

    Raises:
        TypeError: Term/Terms was/were not found in dataset.
    """
    enrichment = [x for x in self.enrichment if x is not None]
    deps = self.original
    deps = [x[x[self.pvalue] <= pvalue] for x in deps]
    groups = self.groups
    colors = dict(zip(groups, self.colors))
    # Select genes from enriched term.
    enrichTerms = []
    for term in Terms:
        enr = [x[x['Term'].str.contains(term)] for x in enrichment]
        enrichTerms.extend(enr)
    enrichTerms = pd.concat(enrichTerms)
    enrichGenes = list(enrichTerms['Genes'])
    enrichGenes = sum(enrichGenes, [])
    enrichGenes = list(set(enrichGenes))
    # Filtering based on enrichgenes
    data = [x[x['gene_name'].isin(enrichGenes)] for x in deps]
    data = [x[['gene_name', "log2(fc)"]] for x in data]
    data = [x.rename(columns={'log2(fc)': y}) for x, y in zip(data, groups)]
    data = [x.groupby('gene_name').mean() for x in data]
    data = pd.concat(data, axis=1).T
    data = data.sort_index(axis=1)
    matrix = data
    matrix = matrix.notnull().astype(int)
    matrix = matrix.fillna(0)
    row_sums = matrix.sum(axis=1)
    matrix = matrix.drop(index=matrix.index[row_sums == 0])
    heatmaps = data
    heatmaps[heatmaps > vmax] = vmax
    heatmaps[heatmaps < vmin] = vmin
    heatmaps = [np.array([heatmaps[x].dropna()]) for x in heatmaps]
    heatmaps = [np.array([[np.nan]])]*len(matrix) + heatmaps
    heatmaps = [x[:, ::-1] for x in heatmaps]
    if len(matrix.columns) == 0:
        raise TypeError('Matrix length is zero. Term/Terms was/were not found in dataset.')
    circos = Circos.initialize_from_matrix(
        matrix,
        start=-265,
        end=95,
        space=0.3,
        r_lim=(93, 100),
        cmap=colors,
        label_kws=dict(r=101, orientation="vertical", size=label_size),
    )
    for sector, heatmap in zip(circos.sectors, heatmaps):
        # Outer Track
        # foldchange track
        if sector.name not in data.index:
            rect_track = sector.add_track((93, 100))
            rect_track.heatmap(heatmap, cmap=colormap, vmin=vmin, vmax=vmax)

    circos.colorbar(bounds=(1.1, 0.3, 0.02, 0.4), vmin=vmin, vmax=vmax, cmap="RdBu_r",
                    colorbar_kws=dict(label="log2(FoldChange)"))
    fig = circos.plotfig()
    if save is not None:
        if vector is True:
            fig.savefig(save+'Term_circular_plot.svg')
        else:
            fig.savefig(save+'Term_circular_plot.png', dpi=dpi)
