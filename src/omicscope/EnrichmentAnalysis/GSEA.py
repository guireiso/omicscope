""" Module for Functional Enrichment Analysis.

This OmicScope module is totally dedicated to Functional Enrichment
Analysis of Omics Experiment.

In this workflow we used GSEApy algorithm to
get two main analysis: Over-Representation Analysis and Gene-Set Enrichment
Analysis. For both cases, GSEApy runs the analysis against Enrichr databases.

Additionally to prepare data to perform enrichment analysis,
OmicScope perform several analysis associated to enrichment analysis and
also can produces ready-to-publish figures.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""

import copy
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

class GSEA:
    def __init__(self, OmicScope, dbs=['KEGG_2021_Human', 'Reactome_2016',
                                           'GO_Molecular_Function_2018',
                                           'GO_Cellular_Component_2018',
                                           'GO_Biological_Process_2018'],
                 organism='Human', pvalue=0.05):
        import copy
        import pandas as pd
        self.Analysis = 'GSEA'
        self.dbs = copy.copy(dbs)
        self.OmicScope = copy.copy(OmicScope)
        self.organism = copy.copy(organism)
        self.pvalue = copy.copy(pvalue)
        dfs = []
        for i in dbs:
            df = self.gsea(db=i)
            dfs.append(df)
        results = pd.concat(dfs)
        self.results = results.sort_values('Adjusted P-value', ignore_index=True)

    def gsea(self, db):
            """GSEA workflow according to gseapy
            """
            from gseapy import gsea
            organism = self.organism
            pvalue_cutoff = self.OmicScope.PValue_cutoff
            foldchange_cutoff = self.OmicScope.FoldChange_cutoff
            omics = self.OmicScope.quant_data
            # Filtering data based on Fold Change and P-value
            omics = omics.loc[(omics['log2(fc)']<=-foldchange_cutoff) | (omics['log2(fc)']>=foldchange_cutoff)]
            omics = omics[omics['pvalue'] < pvalue_cutoff]
            foldchange = dict(zip(omics.gene_name.str.upper(), omics['log2(fc)']))
            
            gsea_input = omics.set_index('gene_name')
            gsea_input = gsea_input.iloc[:, gsea_input.columns.str.contains('\.')]
            gsea_input = np.log2(gsea_input)
            gsea_input = gsea_input.replace(-np.inf, 0)
            groups = gsea_input.columns.str.split('\.').str[1]
            gsea_input = gsea_input.reset_index()
            gsea_input['gene_name'] = gsea_input['gene_name'].str.upper()

            gsea_result = gsea(gsea_input, gene_sets=db,
                        cls=groups, no_plot=True, outdir=None)
            gsea_result = gsea_result.res2d
            gsea_result = gsea_result.reset_index()
            gsea_result['DB'] = db
            gsea_result = gsea_result.rename( {'DB' : 'Gene_set',
                                    'genes':'Genes',
                                    'es' : 'ES',
                                    'nes': 'NES',
                                    'fdr' : 'Adjusted P-value',
                                    'pval':'P-value',
                                    'matched_size':'N_Proteins'},
                                    axis = 'columns')
            gsea_result['Overlap'] = gsea_result.N_Proteins.astype(str) + '/' + gsea_result.geneset_size.astype(str)
            gsea_result.Genes = gsea_result.Genes.str.split(';')
            gsea_result['regulation'] = gsea_result['Genes'].apply(lambda x: [foldchange[i] for i in x])
            gsea_result['down-regulated'] = gsea_result['regulation'].apply(lambda x: len([i for i in x if i < 0]))
            gsea_result['up-regulated'] = gsea_result['regulation'].apply(lambda x: len([i for i in x if i > 0]))
            gsea_result['-log10(pAdj)'] = -np.log10(gsea_result['Adjusted P-value'])
            gsea_result = gsea_result[['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value',
                                    'NES', 'ES', 'Genes', 'N_Proteins', '-log10(pAdj)', 'regulation',
                                    'down-regulated', 'up-regulated']]
            return gsea_result


    def dotplot(self, top=10, palette='Spectral', alpha=1, s=10,
                x_size = 5, y_size = 6, label_wrap = 50,
                save = '', dpi = 300, vector = True):
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
           yaxis = df_db['Term']
           plt.style.use('default')
           plt.figure(figsize=(x_size, y_size))
           from textwrap import wrap
           labels = [ '\n'.join(wrap(l, label_wrap)) for l in yaxis ]
           plt.scatter(x=xaxis, y=labels, s=size*s, c=df_db['NES'], edgecolors='black',
                       linewidths=0.7, cmap=palette, alpha=alpha)
           plt.gca().invert_yaxis()
           plt.title(i)
           sns.despine()
           plt.colorbar().set_label('NES')
           plt.xlabel('-log(p-Adjusted)')
           plt.margins(x=.1, y=0.1)
           if save != '':
               if vector == True:
                   plt.savefig(save + 'dotplot.svg')
               else:
                   plt.savefig(save + 'dotplot.png', dpi=dpi)
           plt.show(block=True)

    def number_deps(self, *Terms, top=20, palette = 'RdBu',
                    save = '', dpi=300, vector = True):
        """Return number of down- and up-regulated proteins
        in top-N or specified pathways.

        Args:
            top (int, optional): top-N enriched terms to be visualized.
            Defaults to 10.
            palette (str, optional): color map for up- and down-regulated
            representations. Defaults to 'RdBu_r'.
            save (str, optional): Path to save figure. Defaults to ''.
            dpi (int, optional): Resolution to save figure. Defaults to 300.
            vector (bool, optional): Save figure in as vector (.svg). Defaults to
            True.
        """
        df = copy.copy(self.results)
        dbs = self.dbs
        if len(Terms)>0:
            df=df[df['Term'].isin(Terms)]
        for i in dbs:
            df_db = df[df['Gene_set']==i]
            df_db = df_db.sort_values(by=['Adjusted P-value', 'N_Proteins'])
            df_db = df_db.iloc[:top, ]
            df_db = df_db.set_index('Term')
            ysize_graph = len(df_db)
            df_db = df_db[['up-regulated', 'down-regulated']]
            df_db = df_db.melt(ignore_index=False)
            df_db = df_db.reset_index()
            df_db = df_db.rename(columns = {'variable':'Regulation',
                                    'value':'Size'})
            plt.figure(figsize = (1,ysize_graph*0.4))
            sns.scatterplot(data=df_db, 
                            x = 'Regulation',
                            y = 'Term',
                            size = 'Size', 
                            palette = palette,
                            edgecolor = 'black',
                            linewidth = 0.5,
                            alpha = 1,
                            hue= 'Regulation',
                            sizes = (0,500))    
            leg = plt.legend(bbox_to_anchor=(1, 1, 0.1, 0))
            leg.get_frame().set_edgecolor('white')
            plt.xticks(rotation = 45, ha = 'right')
            plt.ylabel('')
            plt.xlabel('')
            plt.margins(x =1,)
            sns.despine()
            if save != '':
                if vector == True:
                    plt.savefig(save + 'heatmap.svg')
                else:
                    plt.savefig(save + 'heatmap.png', dpi=dpi)
            plt.show(block=True)

    def pathprotein_network(self, *Terms, top=5, labels=False,
                path_color = 'gray', foldchange_range = [-0.5, 0.5],
                save='', dpi=300):
        """Path-protein network.

        Args:
            top (int, optional): top N enriched terms to be visualized.
                    Defaults to 5.
            labels (bool, optional): Show node labels. Defaults to False.
            path_color (str, optional): Color to plot pathways. Defaults to 'gray'.
            foldchange_range (list, optional): Foldchange range to plot protein colors,
                    such as a heatmap. Defaults to [-0.5, 0.5].
            save (str, optional): _description_. Defaults to ''.
            dpi (int, optional): _description_. Defaults to 300.

        Returns:
            _type_: _description_
        """
        import matplotlib as mpl
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        import networkx as nx
        
        plt.rcParams["figure.dpi"] = dpi
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
        nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
        nx.set_node_attributes(G, dict(zip(carac.index, carac.Size)), name="Size")
        pos = nx.kamada_kawai_layout(G)
        carac=carac.reindex(G.nodes())
        nx.draw(G,
               pos = pos,
               node_color=carac['color'],
               node_size= carac['Size'],
               edgecolors='black',
               linewidths=0.4,
               alpha=0.9,
               width=0.2,
               edge_color='gray')
        if labels is True:
            nx.draw_networkx_labels(G,pos,font_size = 6)
        if save != '':
            nx.write_graphml(G, save+ 'PPNetwork.graphml', named_key_ids = True)
            plt.savefig(save + 'PPNetwork.svg', dpi=dpi)
        plt.show()
        return G


    def gseaheatmap(self, fdr=0.05, save=''):
        import copy
        import matplotlib.pyplot as plt
        import seaborn as sns
        df = copy.copy(self.results)
        df = df[df['fdr'] <= fdr]
        df = df[['nes']]
        sns.heatmap(df, center=0, cmap='RdYlBu_r', vmin=-2, vmax=2)
        plt.ylabel('')
        plt.xticks(rotation=45, ha='right')
        if save != '':
            plt.savefig(save + '_gsea_heatmap.png', dpi=300)
        plt.show()
        plt.show(block=True)
