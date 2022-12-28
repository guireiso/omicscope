# -*- coding: utf-8 -*-
"""
@author: gui_d
"""
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
from matplotlib.collections import PatchCollection
import warnings

def barplot(self, palette = 'Spectral', save = '', vector = True):
    data = copy(self)
    conditions = copy(data.groups)
    group_data = data.original
    difreg = data.group_data

    #-Figura das desregulações
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
            linewidth = 1)
    colors = sns.color_palette(palette, as_cmap=False, n_colors= len(conditions))
    ax.bar(r, deps, edgecolor = 'black', width = 0.5, linewidth = 1)
    plt.xticks(r,
                fontweight= None, rotation = 45, ha = 'right')
    
    ax2.bar(r, whole_proteome, color='lightgray', edgecolor='black', width=0.5,
            linewidth = 1)
    ax2.bar(r, deps,color = colors,  edgecolor = 'black', width = 0.5, linewidth = 1)
    ax.set_ylim(max(whole_proteome[:-1])*1.1, max(whole_proteome)*1.01)  # outliers only
    ax2.set_ylim(0, max(whole_proteome[:-1])*1.005)  # most of the data
    sns.despine()
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(bottom = False)
    plt.xticks(r, conditions,
            fontweight= None, rotation = 45, ha = 'right')
    for i,j,k, c in zip (deps[:-1], r[:-1], whole_proteome[:-1],
                      colors[:-1]):
        ax2.annotate(i, xy = [j,k*1.005], ha = 'center',color = c,
                     weight='bold')
    ax.annotate(deps[-1], [r[-1], whole_proteome[-1]*1.001], ha='center',
                color = colors[-1], weight='bold')
    if save != '':
        if vector == True:
            plt.savefig(save + 'bar_deps.svg')
        else:
            plt.savefig(save + 'bar_deps.png', dpi = 600)

def protein_overlap(self, dpi = 600, min_subset = 10, face_color = 'darkcyan', shad_color="#f0f0f0",
                      edge_color = 'black', linewidth = 1, save = '', vector = True,):
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_contents
    plt.style.context('classic')
    plt.rcParams['grid.alpha']= 0
    plt.rcParams['figure.dpi']=dpi
    plt.rcParams['patch.linewidth'] = linewidth
    plt.rcParams['patch.edgecolor'] = 'black'
    plt.rcParams['patch.force_edgecolor'] = True
    
    data = copy(self)
    genes = []
    for i in data.group_data:
        genes.append(i['gene_name'].drop_duplicates())
    dictionary = dict(zip(data.groups, genes))
    upset = from_contents(dictionary)

    figure = UpSet(upset,  facecolor= face_color, shading_color= shad_color, min_subset_size= min_subset, show_counts=True,
        with_lines=True)
    #for i in data.groups:
    #    figure.style_subsets(present = i, edgecolor=edge_color, linewidth=linewidth)
    figure.plot()
    if save != '':
        if vector == True:
            plt.savefig(save + 'upset_proteins.svg')
        else:
            plt.savefig(save + 'upset_proteins.png', dpi = 600)

def enrichment_overlap(self, dpi = 600, min_subset = 1, face_color = 'darkcyan', shad_color="#f0f0f0",
                      edge_color = 'black', linewidth = 1,save = '', vector = True):
    from upsetplot import UpSet, from_contents
    plt.style.context('classic')
    plt.rcParams['grid.alpha']= 0
    plt.rcParams['figure.dpi']=dpi
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
        except:
            df = pd.DataFrame(columns = ['Gene_set', 'Term', 'Overlap', 'Adjusted P-value', 'Genes'])
            genes.append(df['Term'].drop_duplicates())
    dictionary = dict(zip(data.groups, genes))
    upset = from_contents(dictionary)

    figure = UpSet(upset,  facecolor= face_color, shading_color= shad_color, min_subset_size= min_subset, show_counts=True,
        with_lines=True)
    for i in data.groups:
        figure.style_subsets(present = i, edgecolor=edge_color, linewidth=linewidth)
    figure.plot()
    if save != '':
        if vector == True:
            plt.savefig(save + 'upset_pathways.svg')
        else:
            plt.savefig(save + 'upset_pathways.png', dpi = dpi)
    
def group_pearson(multiples_output, pvalue = 1, palette = 'Spectral', 
                  save = '', vector = True, dpi = 600):
    data = copy(multiples_output)
    conditions = data.groups
    data1 = data.original
    totalMean = []
    colors = sns.color_palette(palette, as_cmap=False, n_colors= len(conditions))
    for i in data1:
        df = i.set_index('Accession')
        df = df[df['pvalue'] < pvalue]
        df = df[['TotalMean']]
        totalMean.append(df)
    wholedata = pd.concat(totalMean, axis = 1, join = 'outer')
    corr = wholedata.corr()
    corr.columns = data.groups
    sns.clustermap(corr, cmap = 'RdYlBu', center = 0.9,
                   col_colors=colors, row_colors = colors)
    
    if save != '':
        if vector == True:
            plt.savefig(save + 'clustermap.svg')
        else:
            plt.savefig(save + 'clustermap.png', dpi = dpi)
    plt.plot()

def Differentially_Regulated(multiples_output,
                             save = '', vector = True, dpi = 600):
    data = copy(multiples_output)
    groups =  data.groups
    difreg = data.group_data
    up = []
    down = []
    
    for i in difreg:
        up.append(len(i[i['log2(fc)']>0]))
        down.append(-len(i[i['log2(fc)']<0]))
    dysregulations = pd.DataFrame(columns = groups,
                                  data=[up, down], index = ['Up-regulated', 'Down-regulated']).transpose()
    
    df = dysregulations

    df = df.reset_index()
    df = df.melt('index')
    df['color'] = np.where(df['value']>0, '#e3432d', '#167a9c')
    df['value'] = abs(df['value'])
    M = 2
    N = len(groups)
    fig, ax = plt.subplots(figsize=(1.15,5/9*len(df)/2))
    scatter = ax.scatter(x = df['variable'], y = df['index'], s = df['value'],
               c = df['color'], ec = 'black', lw = 0.5)
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    sns.despine()
    plt.xticks(rotation = 45, ha = 'right')
    plt.yticks(rotation = 45)
    kw = dict(prop="sizes", num=3, color = 'black', alpha = .6)
    ax.legend(*scatter.legend_elements(**kw), title="# Proteins",handleheight = 2,
                        bbox_to_anchor=(1, 1),loc="upper left", markerscale = 1,
                        edgecolor = 'white')
    plt.margins(x =1, y = 0.1)
    
    if save != '':
        if vector == True:
            plt.savefig(save + 'clustermap.svg')
        else:
            plt.savefig(save + 'clustermap.png', dpi = dpi)
    plt.show()

def network(self, labels = False, save = '',):
    import matplotlib as mpl
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import networkx as nx

    palette = sns.color_palette("Set3", 12).as_hex()
    data = copy(self)
    palette = palette[:len(data.groups)]
    network_frame = []
    for group, df, color in zip(data.groups, data.original, palette):
        df['Experiment'] = group
        df = df[df['pvalue']<0.05]
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
                            vcenter = 0)
    cmap = cm.RdYlBu_r
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_hex =  [mcolors.to_hex(m.to_rgba(x)) for x in network_frame['log2(fc)']]
   
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
    carac = pd.concat([source, target]).reset_index(drop = True)
    carac = carac.drop_duplicates('ID')
    carac = carac.set_index('ID')
    nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="Color")
    nx.set_node_attributes(G, dict(zip(carac.index, carac.Size)), name="Size")
    pos= nx.kamada_kawai_layout(G)
    carac=carac.reindex(G.nodes())
    nx.draw(G,
           pos = pos,
           node_color=carac['color'],
           node_size= carac['Size']/20,
           edgecolors='black',
           linewidths=0.4,
           alpha=0.9,
           width=0.2,
           edge_color='gray')  
    if labels is True:
        nx.draw_networkx_labels(G,pos,font_size = 6)
    if save != '':
        nx.write_graphml(G, save+ 'PPNetwork.graphml', named_key_ids = True)
        if vector == True:
            plt.savefig(save + 'PPNetwork.svg')
        else:
            plt.savefig(save + 'PPNetwork.dpi', dpi=300)
    plt.show()
    return(G)

