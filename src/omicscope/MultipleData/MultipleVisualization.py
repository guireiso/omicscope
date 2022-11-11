# -*- coding: utf-8 -*-
"""
@author: gui_d
"""

import omicscope as omics
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
    cell_data = data.original
    difreg = data.cell_data

    #-Figura das desregulações
    whole_proteome = []
    deps = []
    for i, a in zip(cell_data, difreg):
        whole_proteome.append(len(i))
        deps.append(len(a))
    conditions.extend(['Total'])
    
    proteinsIdentified = []
    for i in cell_data:
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
    for i in data.cell_data:
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
    difreg = data.cell_data
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