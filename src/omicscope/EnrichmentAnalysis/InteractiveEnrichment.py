import copy
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import altair as alt
import pandas as pd


def dotplot(self):
    """Dotplot for enriched terms.

    Args:
        OmicScope (_type_): _description_
        top (int, optional): top-N enriched terms to be visualized.
        Defaults to 10.
        palette (str, optional): color map to visualization,
        for more information https://matplotlib.org/stable/tutorials/colors/colormaps.html.
        Defaults to 'Spectral'.
        alpha (int, optional): transparency of dots Defaults to 1.
        s (int, optional): dotsize. Defaults to 10.
        x_size (int, optional): Size of horizontal axis. Defaults to 5.
        y_size (int, optional): Size of vertical axis. Defaults to 6.
        label_wrap (int, optional): Label wrap. Defaults to 50.
        save (str, optional): Path to save figure. Defaults to ''.
        dpi (int, optional): Resolution to save figure. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
        True.
    """
    import altair as alt
    dbs = self.dbs
    df = copy.copy(self.results)
    for i in list(dbs):
       df_db = df[df['Gene_set'] == i]   
       init_value = int(df_db['-log10(pAdj)'].drop_duplicates()[5])
       max_value = int(df_db['-log10(pAdj)'].drop_duplicates()[0])
       slider = alt.binding_range(min= -np.log10(0.05),
         max=max_value+1,
         step=1, name='-log10(pAdj):')
       selector = alt.selection_single(name="SelectorName", fields=['-log10(pAdj)'],
         bind=slider, init={'-log10(pAdj)': init_value})
       base = alt.Chart(df_db).add_selection(
        selector).properties(width=250)
       D = base.transform_filter(
        alt.datum['-log10(pAdj)'] > selector['-log10(pAdj)']).encode(
            y=alt.Y('Term', sort = '-x'),
            
            x=alt.X('-log10(pAdj)'),
            size = '-log10(pAdj)',
            tooltip='Term:N').mark_point(filled = True)
       if self.Analysis == 'GSEA':
            D = D.mark_point(
            color = 'black', strokeWidth = 1).encode(
            fill = alt.Color('NES', scale=alt.Scale(scheme="redyellowblue", reverse = True, domainMid=0),
                title = 'NES'),)
       else:
        D = D.mark_point(
            color = 'black', strokeWidth = 1).encode(
            fill = alt.Color('N_Proteins', scale=alt.Scale(scheme='purplebluegreen')),)
       return D


def heatmap(self, foldchange=False):
    """Heatmap for enriched terms and proteins

        Args:
            top (int, optional): top N enriched terms to be visualized.
                Defaults to 5.
            foldchange (bool, optional): show the protein fold change.
                Defaults to False.
            foldchange_range (list, optional): Foldchange range to plot protein colors,
                such as a heatmap. Defaults to [-0.5, 0.5].
            save (str, optional): Path to save figure. Defaults to ''.
            dpi (int, optional): Resolution to save figure. Defaults to 300.
            vector (bool, optional): Save figure in as vector (.svg). 
                Defaults to True.
         """
    import altair as alt
    dbs = self.dbs
    omics = copy.copy(self.results)
    deps = copy.copy(self.OmicScope.deps)
    deps['gene_name'] = deps['gene_name'].str.upper()
    deps.columns = ['Genes', 'Accession', 'pvalue', '-log10(p)', 'log2(fc)']
    for i in dbs:
      heatmap_path = omics[omics['Gene_set'] == i]
      heatmap_path = heatmap_path.explode(['Genes','regulation'])
      heatmap_path = heatmap_path.sort_values('-log10(pAdj)', ascending = False)
      init_value = int(heatmap_path['-log10(pAdj)'].drop_duplicates()[5])
      max_value = int(heatmap_path['-log10(pAdj)'].drop_duplicates()[0])
      slider = alt.binding_range(min= -np.log10(0.05),
       max=max_value+1,
       step=1, name='-log10(pAdj):')
      selector = alt.selection_single(name="SelectorName", fields=['-log10(pAdj)'],
      bind=slider, init={'-log10(pAdj)': init_value})
      base = alt.Chart(heatmap_path).add_selection(selector).properties(width=250)
      D = base.transform_filter(
            alt.datum['-log10(pAdj)'] > selector['-log10(pAdj)']).encode(
            x=alt.Y('Genes', sort = '-color'),
            y=alt.Y('Term', sort = '-color'),
            tooltip=[
                alt.Tooltip('Genes', title='Proteins'),
                alt.Tooltip('Term', title='Term'),
                alt.Tooltip('-log10(pAdj)', title='-log10(pAdj)'),
                alt.Tooltip('regulation', title='log2(FC)')
            ]
        ).properties(width=550).mark_rect()
      if foldchange is True:  # verify, before was ==
            D = D.encode(
                color=alt.Color('regulation', scale=alt.Scale(scheme="redyellowblue", reverse = True, domainMid=0),
                title = 'log2(FC)'),
            )
      elif self.Analysis == 'GSEA':
            D = D.encode(
                color=alt.Color('NES', scale=alt.Scale(scheme="redyellowblue", reverse = True, domainMid=0),
                title = 'NES'),
                tooltip=[
                alt.Tooltip('Genes', title='Proteins'),
                alt.Tooltip('Term', title='Term'),
                alt.Tooltip('-log10(pAdj)', title='-log10(pAdj)'),
                alt.Tooltip('NES', title='NES')
            ]
            )
      else:
            D = D.encode(
                color=alt.Color('-log10(pAdj)', scale=alt.Scale(scheme="purplebluegreen"),
                title = 'log10(pAdj)')
            )
      return D


def number_deps(self):
    """Return number of down- and up-regulated proteins
    in top-N or specified pathways.

    Args:
        top (int, optional): top-N enriched terms to be visualized.
        Defaults to 10.
        palette (str, optional): color map for up- and down-regulated
        representations. Defaults to 'RdBu'.
        save (str, optional): Path to save figure. Defaults to ''.
        dpi (int, optional): Resolution to save figure. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
        True.
    """
    df = copy.copy(self.results)
    dbs = self.dbs
    for i in dbs:
        df_db = df[df['Gene_set']==i]
        df_db = df_db.set_index(['Term', '-log10(pAdj)'])
        df_db = df_db[['up-regulated', 'down-regulated']]
        df_db = df_db.melt(ignore_index=False)
        df_db = df_db.reset_index()
        df_db = df_db.rename(columns = {'variable':'Regulation',
                                'value':'Size'})
        init_value = int(df_db['-log10(pAdj)'].drop_duplicates()[5])
        max_value = int(df_db['-log10(pAdj)'].drop_duplicates()[0])
        slider = alt.binding_range(min= -np.log10(0.05),
         max=max_value+1,
         step=1, name='-log10(pAdj):')
        selector = alt.selection_single(name="SelectorName", fields=['-log10(pAdj)'],
         bind=slider, init={'-log10(pAdj)': init_value})
        base = alt.Chart(df_db).add_selection(
        selector).properties(width=250)
        D = base.transform_filter(
            alt.datum['-log10(pAdj)'] > selector['-log10(pAdj)']).mark_point(filled = True,
            color = 'black', 
            ).encode(
                y = 'Term',
                x = 'Regulation',
                size = 'Size',
                fill = alt.Color('Regulation', scale=alt.Scale(range = ['#0c90c4', '#d10808']))
            )
    return D


def enrichment_network(self, *Terms, top=5, labels=False,
            path_color = 'gray', foldchange_range = [-0.5, 0.5]):
    """Path-protein network.

    Args:
        top (int, optional): top N enriched terms to be visualized.
                Defaults to 5.
        labels (bool, optional): Show node labels. Defaults to False.
        path_color (str, optional): Color to plot pathways. Defaults to 'gray'.
        foldchange_range (list, optional): Foldchange range to plot protein colors,
                such as a heatmap. Defaults to [-0.5, 0.5].
        save (str, optional): OmicScope saves networks as figure and Graphml files, which
        can be used in other network software (recommended). Defaults to ''.
        vector (bool, optional): Save image as vector (.svg) Defaults to True.
        dpi (int, optional): Figure resolution. Defaults to 300.


    Returns:
        Graph (NetworkX output): Graph
    """
    import matplotlib as mpl
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import networkx as nx
    df = copy.copy(self.results)
    if len(Terms)>0:
       top = len(Terms)
       df = df[df['Term'].isin(Terms)]
    df = df.iloc[0:top,]
    df = df.explode(['Genes', 'regulation'])
            
    source = pd.DataFrame({'ID': df['Term'],
                       'Size': df['-log10(pAdj)']*20,
                       'type': 'pathway',
                       'color': path_color})
    source = source.drop_duplicates()
 
    norm = mpl.colors.Normalize(vmin=min(foldchange_range),
                                vmax=max(foldchange_range))
    cmap = cm.RdYlBu_r
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_hex =  [mcolors.to_hex(m.to_rgba(x)) for x in df['regulation']]
   
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
   
    carac = pd.concat([source, target]).reset_index(drop = True)
    carac = carac.drop_duplicates('ID')
    carac = carac.set_index('ID')
    nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="color")
    nx.set_node_attributes(G, dict(zip(carac.index, carac.Size.astype(int)/10)), name="size")
    carac=carac.reindex(G.nodes())
    from pyvis.network import Network
    nt = Network('600px', '600px')
    nt.from_nx(G)
    nt.set_options(
        """
        const options = {
  "edges": {
    "color": {
      "inherit": "both"
    },
    "selfReferenceSize": null,
    "selfReference": {
      "angle": 0.7853981633974483
    },
    "smooth": {
      "forceDirection": "none"
    }
  },
  "manipulation": {
    "enabled": true
  },
  "physics": {
    "barnesHut": {
      "theta": 0.2,
      "gravitationalConstant": -2050
    },
    "maxVelocity": 40,
    "minVelocity": 0.25
  }
}"""
    )

    nt.show('G.html')
    return G

def enrichment_map(self, *Terms, top=1000, modularity = True, labels=False,
                pearson_cutoff = 0.5, save='', vector = True, dpi=300):
    """Enrichment map.

    Args:
        top (int, optional): Top terms used to construct network. Defaults to 1000.
        modularity (bool, optional): Returns modularity analysis of Terms. Defaults to True.
        labels (bool, optional): Add Term labels to graph. Defaults to False.
        pearson_cutoff (float, optional): Pearson correlation cutoff . Defaults to 0.5.
        save (str, optional): OmicScope saves networks as figure and Graphml files, which
        can be used in other network software (recommended). Defaults to ''.
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
    original = copy.copy(self.results)
    #Select terms of interest
    if len(Terms)>0:
        top = len(Terms)
        original = original[original['Term'].isin(Terms)]
    #Select top enriched terms by sorting values and droping duplicates with < pvalue
    original = original.sort_values('-log10(pAdj)', ascending = False)
    original = original.drop_duplicates('Term')
    original = original.iloc[0:top,]
    #Subset ORA-results in Term, Genes, and regulation columns
    df = original[['Term', 'Genes', 'regulation']]
    df = df.explode(['Genes', 'regulation'])
    # CrossTab and perform Pearson correlation among Terms
    df = pd.crosstab(df.Genes, df.Term)
    df = df.corr()
    # Filter correlation (R) to define Adjacency Matrix
    # Filter: R < pearson_cutoff 
    df[df<pearson_cutoff] = 0
    # Construct graph
    G = nx.from_pandas_adjacency(df)
    #removing self loops
    G.remove_edges_from(nx.selfloop_edges(G))
    # Node attributes 
    carac = original[['Term', '-log10(pAdj)', 'Combined Score']]
    carac['Label'] = carac.Term
    # Produce color palette according to p-value of each term
    norm = mpl.colors.Normalize(vmin=0,
                                vmax=max(carac['-log10(pAdj)']))
    cmap = cm.PuBu
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_hex =  [mcolors.to_hex(m.to_rgba(x)) for x in carac['-log10(pAdj)']]
    carac['color'] = color_hex
    carac = carac.set_index('Term')
    # Set node attributes
    nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
    nx.set_node_attributes(G, dict(zip(carac.index, carac['-log10(pAdj)'].astype(int)*2)), name="size")
    # Find Communities
    if modularity is True:
        import networkx.algorithms.community as nx_comm
        carac = carac[['-log10(pAdj)', 'Combined Score']]
        # Find communities based on label propagation
        communities = nx_comm.louvain_communities(G)
        module = []
        terms = []
        color = []
        # Linking Terms, module, colors
        norm = mpl.colors.Normalize(vmin=0,
                                    vmax=len(communities)+1)
        cmap = sns.color_palette('turbo', len(communities)+1, desat = .9)
        for i, g in enumerate(communities):
            module.extend([i]*len(g))
            terms.extend(list(g))
            if len(g) > 1:
                color.extend([mcolors.to_hex(cmap[i])]*len(g))
            else:
                color.extend(['#a6a6a6'])
            # DataFrame
        modularity = pd.DataFrame(zip(module,terms,color), columns = ['Module', 'Term','color'])
        # Assign best annotation for each module
        # based on following levels: 1) Node Degree (Hub); 2) P-value
        degree = dict(G.degree) #Getting degree
        modularity['Degree'] = modularity.Term.replace(degree)
        ## Merging carac with modularity data
        carac = carac.merge(modularity, on = 'Term')
        label_module = []
        # Assign labels for each module
        for i in range(0, max(carac.Module)+1):
            label = carac[carac.Module == i]
            if len(label[label.Degree == max(label.Degree)]) == 1:
                label['Label'] = np.where(label.Degree == max(label.Degree),
                                        label.Term, ''   )
                label_module.append(label)
            else:
                label['Label'] = np.where(label['-log10(pAdj)'] == max(label['-log10(pAdj)']),
                                            label.Term, '')
                label_module.append(label)
            # Construct feature dataframe
        carac = pd.concat(label_module)
        carac = carac.sort_values('-log10(pAdj)', ascending= False)
        carac = carac.set_index('Term')
        # Set node attributes to export
        nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="color")
        nx.set_node_attributes(G, dict(zip(carac.index, carac.Label)), name="annotation")
        nx.set_node_attributes(G, dict(zip(carac.index, carac.Module)), name="module")
        carac=carac.reindex(G.nodes())
        #Assign position of each mode
    pos = nx.spring_layout(G, seed = 9, k = 1/len(G.nodes)**0.3)
    carac=carac.reindex(G.nodes())
    # Draw network
    nt = Network(height='700px', width='900px')
    nt.from_nx(G)
    nt.set_options(
        """
        const options = {
  "edges": {
    "color": {
      "inherit": "both"
    },
    "selfReferenceSize": null,
    "selfReference": {
      "angle": 0.7853981633974483
    },
    "smooth": {
      "forceDirection": "none"
    }
  },
  "manipulation": {
    "enabled": true
  },
  "physics": {
    "barnesHut": {
      "theta": 0.2,
      "gravitationalConstant": -2050
    },
    "maxVelocity": 40,
    "minVelocity": 0.25
  }
}"""
    )
    nt.show('G.html')
    if save != '':
        nx.write_graphml(G, save+ 'PathMap.graphml', named_key_ids = True)
    return G