""" Module for EnrichmentScope object visualization

This module allows the user to extract and visualize information from EnrichmentScope object.

Here, it is possible to evaluate enrichment results with 1) dotplots, 2)heatmaps and 3) network analysis.
Dotplots allows the visualization of statistical significance and size of dataset enriched (*number_deps*), and
the number of up- and down-regulated proteins (*number_deps*). Heatmaps can be plotted according statistical
significance and also protein foldchange; for GSEA analysis, the function gsea_heatmap plots the normalization
enrichment score (NES) as color pattern. Finally, the network plots can be used to visualize shared proteins
between pathways (*enrichment_network()*), or perform an enrichment map (*enrichment_map()*).

Additionally, some functions optionally show a specific Term while add the name as Args.

"""

import copy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def dotplot(self, top=10, palette='BuPu', alpha=1, s=10,
            x_size=5, y_size=6, label_wrap=50,
            save=None, dpi=300, vector=True):
    """Dotplot for enriched terms.

    Args:
        top (int, optional): top-N enriched terms to be visualized.
         Defaults to 10.
        palette (str, optional): color map to visualization,
         for more information https://matplotlib.org/stable/tutorials/colors/colormaps.html.
         Defaults to 'BuPu'. For GSEA enrichment, we advise a divergent palette.
        alpha (int, optional): dots transparency. Defaults to 1.
        s (int, optional): dotsize. Defaults to 10.
        x_size (int, optional): Size of horizontal axis. Defaults to 5.
        y_size (int, optional): Size of vertical axis. Defaults to 6.
        label_wrap (int, optional): Label wrap. Defaults to 50.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): Resolution to save figure. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.rcParams["figure.dpi"] = dpi
    dbs = self.dbs
    df = copy.copy(self.results)
    for i in list(dbs):
        df_db = df[df['Gene_set'] == i]
        df_db = df_db.iloc[:top, :]
        xaxis = df_db[['-log10(pAdj)']]
        size = df_db['N_Proteins']
        if self.Analysis == 'GSEA':
            color = df_db['NES']
            xaxis = df_db[['NES']]
        else:
            color = df_db['N_Proteins']
        yaxis = df_db['Term']
        plt.style.use('default')
        plt.figure(figsize=(x_size, y_size))
        from textwrap import wrap
        labels = ['\n'.join(wrap(label, label_wrap)) for label in yaxis]
        plt.scatter(x=xaxis, y=labels, s=size*s, c=color, edgecolors='black',
                    linewidths=0.7, cmap=palette, alpha=alpha)
        plt.gca().invert_yaxis()
        plt.title(i)
        sns.despine()
        if self.Analysis == 'GSEA':
            plt.colorbar().set_label('NES')
            plt.axvline(0, c='black', ls='--')
        else:
            plt.colorbar().set_label('# Proteins')
        plt.xlabel(xaxis.columns[0])
        plt.margins(x=.1, y=0.1)
        if save is not None:
            if vector is True:
                plt.savefig(save + 'dotplot_'+i+'.svg', bbox_inches='tight')
            else:
                plt.savefig(save + 'dotplot_'+i+'.png', dpi=dpi, bbox_inches='tight')
        plt.show(block=True)


def heatmap(self, *Terms, top=5, foldchange=False,
            x_size=5, linewidths=0.01,
            foldchange_range=[-0.5, 0.5], save=None, dpi=300,
            vector=True):
    """Heatmap for enriched terms and proteins

        Args:
            top (int, optional): top N enriched terms to be visualized.
             Defaults to 5.
            foldchange (bool, optional): show the protein fold change.
             Defaults to False.
            foldchange_range (list, optional): Fold change range to plot protein colors,
             such as a heatmap. Defaults to [-0.5, 0.5].
            save (str, optional): Path to save figure. Defaults to None.
            dpi (int, optional): Resolution to save figure. Defaults to 300.
            vector (bool, optional): Save figure in as vector (.svg).
             Defaults to True.
    """
    plt.rcParams["figure.dpi"] = dpi
    omics = copy.copy(self.results)
    foldchange_range = foldchange_range
    if len(Terms) > 0:
        top = len(Terms)
        omics = omics[omics['Term'].isin(Terms)]
    deps = copy.copy(self.OmicScope.deps)
    deps['gene_name'] = deps['gene_name'].str.upper()
    deps.columns = ['Genes', 'Accession', 'pvalue', '-log10(p)', 'log2(fc)']
    for i in omics.Gene_set.drop_duplicates():
        heatmap_path = omics[omics['Gene_set'] == i]
        heatmap_path = heatmap_path.iloc[0:top, :]
        heatmap_path = heatmap_path.explode(['Genes', 'regulation'])
        ordering = list(heatmap_path.Term.drop_duplicates())
        if foldchange is True:
            heatmap = heatmap_path.pivot(values='regulation', index='Genes', columns='Term')
        else:
            heatmap = heatmap_path.pivot(values='-log10(pAdj)', index='Genes', columns='Term')
            heatmap = heatmap[ordering]
            heatmap = heatmap.sort_values(ordering)
        sns.reset_orig()
        px = 1/plt.rcParams['figure.dpi']
        f, ax = plt.subplots(figsize=(x_size,
                                      len(heatmap)*px*20))
        if foldchange is True:
            sns.heatmap(heatmap.astype(float),  cmap='RdYlBu_r',
                        vmin=min(foldchange_range),
                        vmax=max(foldchange_range),
                        yticklabels=True, xticklabels=True,
                        cbar_kws={'label': 'log2(fc)', "shrink": .5},
                        linecolor='black', linewidths=linewidths)
        else:
            sns.heatmap(heatmap,  cmap='Spectral',
                        yticklabels=True, xticklabels=True,
                        cbar_kws={'label': '-log10(p)', "shrink": .5},
                        linecolor='black', linewidths=linewidths)
        plt.xlabel('')
        plt.xticks(rotation=45, ha='right')
        if save is not None:
            if vector is True:
                plt.savefig(save + 'heatmap_'+i+'.svg', bbox_inches='tight')
            else:
                plt.savefig(save + 'heatmap_'+i+'.png', dpi=dpi, bbox_inches='tight')
        plt.show(block=True)


def number_deps(self, *Terms, top=20, palette='RdBu',
                save=None, dpi=300, vector=True):
    """Number of DEPs

    Return number of down- and up-regulated proteins
    in top-N or specified pathways.

    Args:
        top (int, optional): top-N enriched terms to be visualized.
        Defaults to 10.
        palette (str, optional): color map for up- and down-regulated
        representations. Defaults to 'RdBu'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): Resolution to save figure. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.rcParams["figure.dpi"] = dpi
    df = copy.copy(self.results)
    if len(Terms) > 0:
        df = df[df['Term'].isin(Terms)]
    for i in df.Gene_set.drop_duplicates():
        df_db = df[df['Gene_set'] == i]
        df_db = df_db.sort_values(by=['Adjusted P-value', 'N_Proteins'])
        df_db = df_db.iloc[:top, ]
        df_db = df_db.set_index('Term')
        ysize = len(df_db)
        df_db = df_db[['up-regulated', 'down-regulated']]
        df_db = df_db.melt(ignore_index=False)
        df_db = df_db.reset_index()
        df_db = df_db.rename(columns={'variable': 'Regulation',
                                      'value': 'Size'})
        if min(df_db['Size']) == 0:
            sizes = (0, 500)
        else:
            sizes = (10, 500)
        plt.figure(figsize=(1, ysize*0.4))
        sns.scatterplot(data=df_db,
                        x='Regulation',
                        y='Term',
                        size='Size',
                        palette=palette,
                        edgecolor='black',
                        linewidth=0.5,
                        alpha=1,
                        hue='Regulation',
                        sizes=sizes)
        leg = plt.legend(bbox_to_anchor=(1, 1, 0, 0))
        leg.get_frame().set_edgecolor('white')
        plt.xticks(rotation=45, ha='right')
        plt.ylabel('')
        plt.xlabel('')
        plt.margins(x=1, y=1/ysize)
        plt.tight_layout()
        sns.despine()
        if save is not None:
            if vector is True:
                plt.savefig(save + 'number_deps_'+i+'.svg', bbox_inches='tight')
            else:
                plt.savefig(save + 'number_deps_'+i+'.png', dpi=dpi, bbox_inches='tight')
        plt.show(block=True)


def enrichment_network(self, *Terms, top=5, labels=False,
                       term_color='#a1a1a1', foldchange_range=[-0.5, 0.5],
                       save=None, vector=True, dpi=300):
    """EnrichmentTerm-protein network.

    Network visualization to find proteins that are shared by different enriched Terms.

    Args:
        top (int, optional): top N enriched terms to be visualized.
         Defaults to 5.
        labels (bool, optional): Show node labels. Defaults to False.
        term_color (str, optional): Color to plot pathways. Defaults to '#a1a1a1'.
        foldchange_range (list, optional): Fold change range to plot protein colors,
         such as a heatmap. Defaults to [-0.5, 0.5].
        save (str, optional): OmicScope saves networks as figure and Graphml files, which
         can be used in other network software (recommended). Defaults to None.
        vector (bool, optional): Save image as vector (.svg) Defaults to True.
        dpi (int, optional): Figure resolution. Defaults to 300.


    Returns:
        Graph (NetworkX output): Graph
    """
    import matplotlib as mpl
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import networkx as nx

    plt.rcParams["figure.dpi"] = dpi
    omics = copy.copy(self.results)
    if len(Terms) > 0:
        top = len(Terms)
        omics = omics[omics['Term'].isin(Terms)]
    G_objects = []
    for i in omics.Gene_set.drop_duplicates():
        df = omics[omics['Gene_set'] == i]
        df = df.iloc[0:top,]
        df = df.explode(['Genes', 'regulation'])
        if self.Analysis == 'GSEA':
            norm = mpl.colors.TwoSlopeNorm(vmin=omics[omics['Gene_set'] == i]['NES'].min(),
                                           vcenter=0,
                                           vmax=omics[omics['Gene_set'] == i]['NES'].max(),
                                           )
            cmap = cm.RdYlBu_r
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            term_color = [mcolors.to_hex(m.to_rgba(x)) for x in df['NES']]
        source = pd.DataFrame({'ID': df['Term'],
                               'Size': df['-log10(pAdj)']*20,
                               'type': 'pathway',
                               'color': term_color})
        source = source.drop_duplicates()

        norm = mpl.colors.TwoSlopeNorm(vmin=min(foldchange_range),
                                       vmax=max(foldchange_range),
                                       vcenter=0)
        cmap = cm.RdYlBu_r
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        color_hex = [mcolors.to_hex(m.to_rgba(x)) for x in df['regulation']]

        target = pd.DataFrame({'ID': df['Genes'],
                               'Size':  min(df['-log10(pAdj)'])*4,
                               'type': 'Protein',
                               'color': color_hex})
        target = target.drop_duplicates()
        edgelist = df[['Term', 'Genes']]
        G = nx.from_pandas_edgelist(edgelist,
                                    source='Term',
                                    target='Genes',
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
                node_size=carac['Size'],
                edgecolors='black',
                linewidths=0.4,
                alpha=0.9,
                width=0.2,
                edge_color='#a1a1a1')
        if labels is True:
            nx.draw_networkx_labels(G, pos, font_size=6)
        if save is not None:
            nx.write_graphml(G, save + 'PPNetwork_'+i+'.graphml', named_key_ids=True)
            if vector is True:
                plt.savefig(save + 'PPNetwork_'+i+'.svg', bbox_inches='tight')
            else:
                plt.savefig(save + 'PPNetwork_'+i+'.dpi', dpi=300, bbox_inches='tight')
        plt.show()
        G_objects.append(G)
    return G_objects


def enrichment_map(self, *Terms, top=1000, modules=True, labels=False,
                   similarity_cutoff=0.25, metric='jaccard',
                   save=None, vector=True, dpi=300):
    """Enrichment map.

    Since several proteins are presented in more than one pathway, enrichment map
    shows pathway as nodes and the edge thickness is proportional to the amount
    of proteins shared between two terms (similarity score)

    Args:
        top (int, optional): Top terms used to construct network. Defaults to 1000.
        modules (bool, optional): Returns modularity analysis of Terms. Defaults to True.
        labels (bool, optional): Add Term labels to graph. Defaults to False.
        similarity_cutoff (float, optional): similarity score cutoff based on statistical analysis performed.
         Defaults to 0.25.
        metric (str, optional): statistical algorithm to perform pairwise comparison. Defaults to 'jaccard'.
         Optionally, user can test other algorithm described in scipy.spatial.distance.
        save (str, optional): OmicScope saves networks as figure and Graphml files, which
         can be used in other network software (recommended). Defaults to None.
        vector (bool, optional): Save image as vector (.svg) Defaults to True.
        dpi (int, optional): Figure resolution. Defaults to 300.

    Returns:
        Graph (NetworkX output): Graph
    """
    import matplotlib as mpl
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import networkx as nx
    from sklearn.metrics import pairwise_distances

    plt.rcParams["figure.dpi"] = dpi
    omics = copy.copy(self.results)
    # Select terms of interest
    if len(Terms) > 0:
        top = len(Terms)
        omics = omics[omics['Term'].isin(Terms)]
    Analysis = self.Analysis
    G_objects = []
    for i in omics.Gene_set.drop_duplicates():
        original = omics[omics['Gene_set'] == i]
        # Select top enriched terms by sorting values and droping duplicates with < pvalue
        original = original.sort_values('-log10(pAdj)', ascending=False)
        original = original.drop_duplicates('Term')
        original = original.iloc[0:top,]
        # Subset ORA-results in Term, Genes, and regulation columns
        df = original[['Term', 'Genes', 'regulation']]
        df = df.explode(['Genes', 'regulation'])
        # CrossTab and perform Pearson correlation among Terms
        df = pd.crosstab(df.Genes, df.Term)
        colnames = df.columns
        df = pairwise_distances(df.T.to_numpy(), metric=metric)
        df = 1 - pd.DataFrame(df, index=colnames, columns=colnames)
        # Filter correlation (R) to define Adjacency Matrix
        # Filter: R < pearson_cutoff
        df[df < similarity_cutoff] = 0
        # Construct graph
        G = nx.from_pandas_adjacency(df)
        G.edges(data=True)
        # removing self loops
        G.remove_edges_from(nx.selfloop_edges(G))

        if Analysis == 'ORA':
            # Node attributes
            carac = original[['Term', '-log10(pAdj)', 'Combined Score']]
            carac['Label'] = carac.Term
            # Produce color palette according to p-value of each term
            norm = mpl.colors.Normalize(vmin=int(carac['-log10(pAdj)'].drop_duplicates().min()),
                                        vmax=int(carac['-log10(pAdj)'].drop_duplicates().max()))
            cmap = cm.PuBu
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            color_hex = [mcolors.to_hex(m.to_rgba(x)) for x in carac['-log10(pAdj)']]
            carac['color'] = color_hex
        elif Analysis == 'GSEA':
            carac = original[['Term', '-log10(pAdj)', 'NES']]
            carac['Label'] = carac.Term
            norm = mpl.colors.TwoSlopeNorm(vmin=int(carac['NES'].drop_duplicates().min()),
                                           vmax=int(
                carac['NES'].drop_duplicates().max()),
                vcenter=0)
            cmap = cm.RdYlBu_r
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            carac['color'] = [mcolors.to_hex(m.to_rgba(x)) for x in carac['NES']]
        carac = carac.set_index('Term')
        # Set node attributes
        nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
        nx.set_node_attributes(G, dict(zip(carac.index, carac['-log10(pAdj)'])), name="Size")
        # Find Communities
        if modules is True:
            import networkx.algorithms.community as nx_comm
            if Analysis == 'ORA':
                carac = carac[['-log10(pAdj)', 'Combined Score']]
            elif Analysis == 'GSEA':
                carac = carac[['-log10(pAdj)', 'NES']]
            # Find communities based on label propagation
            communities = nx_comm.louvain_communities(G)
            module = []
            color = []
            degree = []
            # Linking Terms, module, colors
            norm = mpl.colors.Normalize(vmin=0,
                                        vmax=len(communities)+1)
            cmap = sns.color_palette('turbo', len(communities)+1, desat=.9)
            for i, g in enumerate(communities):
                subgraph = G.subgraph(g)
                degree_subgraph = subgraph.degree
                degree.extend(degree_subgraph)
                module.extend([i]*len(g))
                if len(g) > 1:
                    color.extend([mcolors.to_hex(cmap[i])]*len(g))
                else:
                    color.extend(['#a6a6a6'])
                # DataFrame
            modularity = pd.DataFrame(zip(module, [x[0] for x in degree], color,
                                          [x[1] for x in degree]), columns=[
                'Module', 'Term', 'color', 'Degree'])
            carac = carac.merge(modularity, on='Term')
            carac = carac.sort_values('-log10(pAdj)', ascending=False)
            Label = carac[['Term', 'Module', 'Degree',
                           '-log10(pAdj)']].sort_values(['Module', 'Degree', '-log10(pAdj)'],
                                                        ascending=False)
            Label['mark'] = Label[['Module']].duplicated()
            carac = carac.merge(Label, how='left', on=['Term',
                                'Module', 'Degree', '-log10(pAdj)'])

            carac['Label'] = np.where(
                carac['mark'] == False,
                carac['Term'], "")
            carac = carac.set_index('Term')
            # Set node attributes to export
            nx.set_node_attributes(G, dict(zip(carac.index, carac.Degree)), name="IntraModuleDegree")
            nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
            nx.set_node_attributes(G, dict(zip(carac.index, carac.Label)), name="Annotation")
            nx.set_node_attributes(G, dict(zip(carac.index, carac.Module)), name="Module")
            carac = carac.reindex(G.nodes())
            # Assign position of each mode
        edges = G.edges
        weights = [G[u][v]['weight'] for u, v in edges]
        weights = [round(x, 2) for x in weights]
        norm = [float(i)/np.mean(weights) for i in weights]
        pos = nx.spring_layout(G, k=1/len(G.nodes)**0.3)
        carac = carac.reindex(G.nodes())
        G.remove_edges_from(nx.selfloop_edges(G))
        # Draw network
        nx.draw(G, pos=pos,
                node_color=carac.color,
                node_size=carac['-log10(pAdj)']*10,
                edgecolors='black',
                linewidths=0.4,
                alpha=0.95,
                width=norm,
                edge_color='gray')
        if labels is True:
            nx.draw_networkx_labels(G, pos, carac.Label, font_size=6)
        if save is not None:
            nx.write_graphml(G, save + 'PathMap.graphml', named_key_ids=True)
            if vector is True:
                plt.savefig(save + 'PathMap.svg', bbox_inches='tight')
            else:
                plt.savefig(save + 'PathMap.dpi', dpi=300, bbox_inches='tight')
        plt.show(block=True)
        G_objects.append(G)
    return G_objects


def gsea_heatmap(self, *Terms, top=5, linewidths=0.01, save=None, dpi=300, vector=True):
    """GSEA Heatmap
    Plot a heatmap with colors based on Normalized Enrichment Score (NES)
    reported by GSEA Analysis.

    Args:
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): Resolution to save figure. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.rcParams['figure.dpi'] = dpi
    omics = copy.copy(self.results)
    omics = omics.sort_values('-log10(pAdj)')
    if len(Terms) > 0:
        top = len(Terms)
        omics = omics[omics['Term'].isin(Terms)]
    for i in omics.Gene_set.drop_duplicates():
        df = omics[omics['Gene_set'] == i]
        df = df.iloc[0:top, :]
    df = df.set_index('Term')
    df = df[['NES']].astype(float)
    df = df.sort_values('NES')
    px = 1/plt.rcParams['figure.dpi']
    f, ax = plt.subplots(figsize=(0.5,
                                  len(df)*px*100))
    sns.heatmap(df, center=0, cmap='RdYlBu_r', vmin=-2, vmax=2,
                linewidths=linewidths, linecolor='black')
    plt.ylabel('')
    plt.xticks(rotation=45, ha='right')
    if save is not None:
        if vector is True:
            plt.savefig(save + 'gsea_heatmap'+i+'.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'gsea_heatmap'+i+'.png', dpi=dpi, bbox_inches='tight')
    plt.show(block=True)
