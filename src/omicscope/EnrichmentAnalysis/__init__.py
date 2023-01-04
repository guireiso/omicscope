import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class EnrichmentScope():
    def __init__(self, OmicScope, Analysis, dbs=['KEGG_2021_Human'],
                 padjust_cutoff=0.05, organism='Human'):
        """_summary_

        Args:
            OmicScope (_type_): _description_
            Analysis (_type_): _description_
            dbs (list, optional): _description_. Defaults to ['KEGG_2021_Human'].
            padjust_cutoff (float, optional): _description_. Defaults to 0.05.
            organism (str, optional): _description_. Defaults to 'Human'.
            save (str, optional): _description_. Defaults to ''.
        """
        if Analysis == 'ORA':
            from .ORA import Enrichment
        elif Analysis == 'GSEA':
            from .GSEA import Enrichment
        else:
            raise ValueError('You must choose between ORA or GSEA workflow')
        self.OmicScope = copy.copy(OmicScope)
        self.Analysis = Analysis
        self.dbs = copy.copy(dbs)
        self.organism = copy.copy(organism)
        self.padjust_cutoff = copy.copy(padjust_cutoff)
        enrichment = Enrichment(OmicScope=OmicScope, dbs=dbs,
                                padjust_cutoff=padjust_cutoff, organism=organism)
        self.results = enrichment.results

    def libraries(self):
        """Get Libraries from Enrichr
        """
        from gseapy import get_library_name
        libraries = get_library_name()
        return (libraries)

    def savefile(self, Path: str):
        from copy import copy
        save = Path
        data = copy(self)
        string = '-'.join(data.OmicScope.Conditions)
        with open(save + '/' + string + '.omics', 'w') as f:
            expression = data.OmicScope.quant_data[['gene_name', 'Accession',
                                                    data.OmicScope.pvalue, 'log2(fc)', 'TotalMean']].to_csv(sep='\t', index=False)
            if self.Analysis == 'ORA':
                enrichment = data.results[['Gene_set', 'Term', 'Overlap', 'Adjusted P-value', 'Genes']].to_csv(sep='\t', index=False)
            elif self.Analysis == 'GSEA':
                enrichment = data.results[['Gene_set', 'Term', 'NES', 'Overlap', 'Adjusted P-value', 'Genes']].to_csv(sep='\t', index=False)

            f.write("OmicScope v1.0.0" + "\n" +
                    "This file is the output performed by OmicScope pipeline and can be used as input" +
                    " for group comparisons having the controling group used as used according to OmicScope." +
                    "Please, cite: Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: an Comprehensive Python " +
                    "library for Systems Biology Visualization" +
                    '\nControlGroup:' + '\t' + data.OmicScope.ctrl + '\n' +
                    'Experimental:' + '\t' + '\t'.join(data.OmicScope.experimental) + '\n' +
                    'Statistics:' + '\t' + data.OmicScope.pvalue + '\n' +
                    'Expression:\n' + '-------\n' +
                    expression + '\n' +
                    'Enrichment Analysis:\n' + '-------\n' +
                    enrichment)

    def dotplot(self, top=10, palette='Spectral', alpha=1, s=10,
                x_size=5, y_size=6, label_wrap=50,
                save='', dpi=300, vector=True):
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
        plt.rcParams["figure.dpi"] = dpi
        dbs = self.dbs
        df = copy.copy(self.results)
        for i in list(dbs):
            df_db = df[df['Gene_set'] == i]
            df_db = df_db.iloc[:top, :]
            xaxis = df_db['-log10(pAdj)']
            size = df_db['N_Proteins']
            if self.Analysis == 'GSEA':
                color = df_db['NES']
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
            else:
                plt.colorbar().set_label('# Proteins')
            plt.xlabel('-log(p-Adjusted)')
            plt.margins(x=.1, y=0.1)
            if save != '':
                if vector is True:
                    plt.savefig(save + 'dotplot.svg')
                else:
                    plt.savefig(save + 'dotplot.png', dpi=dpi)
            plt.show(block=True)

    def heatmap(self, *Terms, top=5, foldchange=False,
                x_size=5, y_size=6, linewidths=0.01,
                foldchange_range=[-0.5, 0.5], save='', dpi=300,
                vector=True):
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
        plt.rcParams["figure.dpi"] = dpi
        dbs = self.dbs
        omics = copy.copy(self.results)
        foldchange_range = foldchange_range
        if len(Terms) > 0:
            top = len(Terms)
            omics = omics[omics['Term'].isin(Terms)]
        deps = copy.copy(self.OmicScope.deps)
        deps['gene_name'] = deps['gene_name'].str.upper()
        deps.columns = ['Genes', 'Accession', 'pvalue', '-log10(p)', 'log2(fc)']
        for i in dbs:
            heatmap_path = omics[omics['Gene_set'] == i]
            heatmap_path = heatmap_path.iloc[0:top, :]
            heatmap_path = heatmap_path.explode(['Genes', 'regulation'])
            ordering = list(heatmap_path.Term.drop_duplicates())
            if foldchange is True:  # verify, before was ==
                heatmap = heatmap_path.pivot(values='regulation', index='Genes', columns='Term')
            else:
                heatmap = heatmap_path.pivot(values='-log10(pAdj)', index='Genes', columns='Term')
                heatmap = heatmap[ordering]
                heatmap = heatmap.sort_values(ordering)
                sns.reset_orig()
                f, ax = plt.subplots(figsize=(x_size, y_size))
                if foldchange is True:  # verify, before was ==
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
                if save != '':
                    if vector is True:
                        plt.savefig(save + 'heatmap.svg')
                    else:
                        plt.savefig(save + 'heatmap.png', dpi=dpi)
                plt.show(block=True)

    def number_deps(self, *Terms, top=20, palette='RdBu',
                    save='', dpi=300, vector=True):
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
        if len(Terms) > 0:
            df = df[df['Term'].isin(Terms)]
        for i in dbs:
            df_db = df[df['Gene_set'] == i]
            df_db = df_db.sort_values(by=['Adjusted P-value', 'N_Proteins'])
            df_db = df_db.iloc[:top, ]
            df_db = df_db.set_index('Term')
            ysize_graph = len(df_db)
            df_db = df_db[['up-regulated', 'down-regulated']]
            df_db = df_db.melt(ignore_index=False)
            df_db = df_db.reset_index()
            df_db = df_db.rename(columns={'variable': 'Regulation',
                                          'value': 'Size'})
            plt.figure(figsize=(1, ysize_graph*0.4))
            sns.scatterplot(data=df_db,
                            x='Regulation',
                            y='Term',
                            size='Size',
                            palette=palette,
                            edgecolor='black',
                            linewidth=0.5,
                            alpha=1,
                            hue='Regulation',
                            sizes=(0, 500))
            leg = plt.legend(bbox_to_anchor=(1, 1, 0.1, 0))
            leg.get_frame().set_edgecolor('white')
            plt.xticks(rotation=45, ha='right')
            plt.ylabel('')
            plt.xlabel('')
            plt.margins(x=1,)
            sns.despine()
            if save != '':
                if vector is True:
                    plt.savefig(save + 'number_deps.svg')
                else:
                    plt.savefig(save + 'number_deps.png', dpi=dpi)
            plt.show(block=True)

    def enrichment_network(self, *Terms, top=5, labels=False,
                           path_color='gray', foldchange_range=[-0.5, 0.5],
                           save='', vector=True, dpi=300):
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

        plt.rcParams["figure.dpi"] = dpi
        df = copy.copy(self.results)
        if len(Terms) > 0:
            top = len(Terms)
            df = df[df['Term'].isin(Terms)]
        df = df.iloc[0:top,]
        df = df.explode(['Genes', 'regulation'])
        if self.Analysis == 'GSEA':
            norm = mpl.colors.TwoSlopeNorm(vmin=min(df['NES']),
                                           vmax=max(df['NES']),
                                           vcenter=0)
            cmap = cm.RdYlBu_r
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            path_color = [mcolors.to_hex(m.to_rgba(x)) for x in df['NES']]
        source = pd.DataFrame({'ID': df['Term'],
                               'Size': df['-log10(pAdj)']*20,
                               'type': 'pathway',
                               'color': path_color})
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
                edge_color='gray')
        if labels is True:
            nx.draw_networkx_labels(G, pos, font_size=6)
        if save != '':
            nx.write_graphml(G, save + 'PPNetwork.graphml', named_key_ids=True)
            if vector is True:
                plt.savefig(save + 'PPNetwork.svg')
            else:
                plt.savefig(save + 'PPNetwork.dpi', dpi=300)
        plt.show()
        return G

    def enrichment_map(self, *Terms, top=1000, modularity=True, labels=False,
                       pearson_cutoff=0.5, save='', vector=True, dpi=300):
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
        # Select terms of interest
        if len(Terms) > 0:
            top = len(Terms)
            original = original[original['Term'].isin(Terms)]
        # Select top enriched terms by sorting values and droping duplicates with < pvalue
        original = original.sort_values('-log10(pAdj)', ascending=False)
        original = original.drop_duplicates('Term')
        original = original.iloc[0:top,]
        # Subset ORA-results in Term, Genes, and regulation columns
        df = original[['Term', 'Genes', 'regulation']]
        df = df.explode(['Genes', 'regulation'])
        # CrossTab and perform Pearson correlation among Terms
        df = pd.crosstab(df.Genes, df.Term)
        df = df.corr()
        # Filter correlation (R) to define Adjacency Matrix
        # Filter: R < pearson_cutoff
        df[df < pearson_cutoff] = 0
        # Construct graph
        G = nx.from_pandas_adjacency(df)
        # removing self loops
        G.remove_edges_from(nx.selfloop_edges(G))

        if self.Analysis == 'ORA':
            # Node attributes
            carac = original[['Term', '-log10(pAdj)', 'Combined Score']]
            carac['Label'] = carac.Term
            # Produce color palette according to p-value of each term
            norm = mpl.colors.Normalize(vmin=0,
                                        vmax=max(carac['-log10(pAdj)']))
            cmap = cm.PuBu
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            color_hex = [mcolors.to_hex(m.to_rgba(x)) for x in carac['-log10(pAdj)']]
            carac['color'] = color_hex
        elif self.Analysis == 'GSEA':
            carac = original[['Term', '-log10(pAdj)', 'NES']]
            carac['Label'] = carac.Term
            norm = mpl.colors.TwoSlopeNorm(vmin=min(carac['NES']),
                                           vmax=max(carac['NES']),
                                           vcenter=0)
            cmap = cm.RdYlBu_r
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            carac['color'] = [mcolors.to_hex(m.to_rgba(x)) for x in carac['NES']]
        carac = carac.set_index('Term')
        # Set node attributes
        nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
        nx.set_node_attributes(G, dict(zip(carac.index, carac['-log10(pAdj)'])), name="Size")
        # Find Communities
        if modularity is True:
            import networkx.algorithms.community as nx_comm
            if self.Analysis == 'ORA':
                carac = carac[['-log10(pAdj)', 'Combined Score']]
            elif self.Analysis == 'GSEA':
                carac = carac[['-log10(pAdj)', 'NES']]
            # Find communities based on label propagation
            communities = nx_comm.louvain_communities(G)
            module = []
            terms = []
            color = []
            # Linking Terms, module, colors
            norm = mpl.colors.Normalize(vmin=0,
                                        vmax=len(communities)+1)
            cmap = sns.color_palette('turbo', len(communities)+1, desat=.9)
            for i, g in enumerate(communities):
                module.extend([i]*len(g))
                terms.extend(list(g))
                if len(g) > 1:
                    color.extend([mcolors.to_hex(cmap[i])]*len(g))
                else:
                    color.extend(['#a6a6a6'])
                # DataFrame
            modularity = pd.DataFrame(zip(module, terms, color), columns=['Module', 'Term', 'color'])
            # Assign best annotation for each module
            # based on following levels: 1) Node Degree (Hub); 2) P-value
            degree = dict(G.degree)  # Getting degree
            modularity['Degree'] = modularity.Term.replace(degree)
            # Merging carac with modularity data
            carac = carac.merge(modularity, on='Term')
            label_module = []
            # Assign labels for each module
            for i in range(0, max(carac.Module)+1):
                label = carac[carac.Module == i]
                if len(label[label.Degree == max(label.Degree)]) == 1:
                    label['Label'] = np.where(label.Degree == max(label.Degree),
                                              label.Term, '')
                    label_module.append(label)
                else:
                    label['Label'] = np.where(label['-log10(pAdj)'] == max(label['-log10(pAdj)']),
                                              label.Term, '')
                    label_module.append(label)
                # Construct feature dataframe
            carac = pd.concat(label_module)
            carac = carac.sort_values('-log10(pAdj)', ascending=False)
            carac = carac.set_index('Term')
            # Set node attributes to export
            nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
            nx.set_node_attributes(G, dict(zip(carac.index, carac.Label)), name="Annotation")
            nx.set_node_attributes(G, dict(zip(carac.index, carac.Module)), name="Module")
            carac = carac.reindex(G.nodes())
            # Assign position of each mode
        pos = nx.spring_layout(G, seed=9, k=1/len(G.nodes)**0.3)
        carac = carac.reindex(G.nodes())
        # Draw network
        nx.draw(G, pos=pos,
                node_color=carac.color,
                node_size=carac['-log10(pAdj)']*10,
                edgecolors='black',
                linewidths=0.4,
                alpha=0.95,
                width=0.2,
                edge_color='gray')
        if labels is True:
            nx.draw_networkx_labels(G, pos, carac.Label, font_size=6)
        if save != '':
            nx.write_graphml(G, save + 'PathMap.graphml', named_key_ids=True)
            if vector is True:
                plt.savefig(save + 'PathMap.svg')
            else:
                plt.savefig(save + 'PathMap.dpi', dpi=300)
        plt.show(block=True)
        return G

    def gsea_heatmap(self, save=''):
        if self.Analysis == 'ORA':
            raise Exception('The function gsea_heatmap cannot be run with ORA analysis.')
        df = copy.copy(self.results)
        df = df[df['Adjusted P-value'] <= 1]
        df = df.set_index('Term')
        df = df[['NES']].astype(float)
        sns.heatmap(df, center=0, cmap='RdYlBu_r', vmin=-2, vmax=2,
                    linewidths=0.8, linecolor='black')
        plt.ylabel('')
        plt.xticks(rotation=45, ha='right')
        if save != '':
            plt.savefig(save + '_gsea_heatmap.png', dpi=300)
        plt.show()
