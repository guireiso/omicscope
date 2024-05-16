""" Module for OmicScope Object Visualization

This module allows the user to extract and visualize information from OmicScope object.
Here, it is possible to evaluate data normalization (MA-Plot, Volcano Plot, Dynamic range plot),
individual protein abundance (barplot, boxplot), search for protein-protein interactions,
and perform Principal Component Analysis (PCA), Hierarchical clustering analysis (heatmap, pearson
correlation plot) and K-means clustering (k-trend).

Some functions below allow user to choose protein to be highlighted and/or plotted.
For that, user must write protein 'gene_name' (See examples in OmicScope Object tab).
Additionally, colors and color palettes follows the matplotlib and seaborn libraries options.

"""

import copy
import itertools
import random

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import zscore
from sklearn import preprocessing
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


def bar_ident(self, logscale=False, col='darkcyan', save=None, dpi=300,
              vector=True):
    """Show the amount of entities identified and differentially regulated
    in the study.

    Args:
        logscale (bool, optional): Y-axis log-scaled. Defaults to True.
        col (str, optional): Color. Defaults to 'darkcyan'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): Resolution to save figure. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to True.

    Returns:
        ax [matplotlib object]: Barplot
    """
    # Define plt parameters
    plt.style.use('default')
    sns.set(rc={'figure.figsize': (11.7, 8.27)})
    sns.set_style("ticks")
    plt.rcParams["figure.dpi"] = dpi
    OmicScope = copy.copy(self)
    df = OmicScope.quant_data
    # Get number of identified proteins
    identified = df.Accession.count()
    # Get number of quantified proteins
    quantified = df.dropna(axis=0, subset=[OmicScope.pvalue]).Accession.count()
    # Get number of differentially regulated proteins
    deps = OmicScope.deps['Accession'].count()
    if identified != quantified:
        df = ['Identified', 'Quantified', 'Differentially\nregulated']
        protein_number = [identified, quantified, deps]
    else:
        df = ['Quantified', 'Differentially\nregulated']
        protein_number = [quantified, deps]
    df = pd.DataFrame(protein_number, df)
    # Log scaling y-axis
    if logscale is True:
        ax = df.plot(kind='bar', color=col, rot=0, edgecolor="black", log=True,
                     figsize=(2.4, 3.5))
    else:
        ax = df.plot(kind='bar', color=col, rot=0, edgecolor="black",
                     figsize=(2.4, 3.5))
    # Add label upside bar

    def autolabel(rects, ax):
        for rect in rects:
            x = rect.get_x() + rect.get_width() / 2.
            y = rect.get_height()
            ax.annotate("{}".format(y), (x, y), xytext=(0, 5), textcoords="offset points",
                        ha='center', va='bottom')
    autolabel(ax.patches, ax)
    ax.margins(y=0.1)
    ax.legend().remove()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.ylabel('#Proteins')
    plt.title(label=OmicScope.ctrl + ' vs ' + '-'.join(OmicScope.experimental), loc='left')
    plt.grid(b=False)
    if save is not None:
        if vector is True:
            plt.savefig(save + 'barplot.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'barplot.png', dpi=dpi, bbox_inches='tight')
    plt.show()
    return ax

def normalization_boxplot(self, palette='Dark2', dpi=300,
                          save=None, vector=True):
        """
        Generates a boxplot of log2-transformed expression data for each sample to 
        evaluate normalization.

        Parameters:
        - palette (str): The color palette to use for the boxplot. Default is 'Dark2'.
        - dpi (int): resolution of the saved figure. Default is 300.
        - save (str): Filepath to save the generated plot. If None, the plot is not saved. Default is None.
        - vector (bool): If True, save the plot in vector format (SVG). If False, save in raster format (PNG). Default is True.

        """
        plt.rcParams['figure.dpi']=dpi
        df = self.expression.copy()
        df = np.log2(df)
        df = df.replace([-np.inf,np.inf],0)
        conditions = self.pdata.Condition
        colors = sns.color_palette(palette=palette,
                                as_cmap=False, n_colors=len(conditions.drop_duplicates()))
        color_dict = {i:j for i,j in zip(conditions.drop_duplicates(),colors)}
        bplot = plt.boxplot(df, vert=True,  # vertical box alignment
                            patch_artist=True, labels=df.columns,
                            notch=True, medianprops=dict(color='black'),
                            boxprops=dict(linewidth=0.8))
        plt.xticks(rotation=90)
        plt.title('Normalization Boxplot')
        for patch, condition in zip(bplot['boxes'], conditions):
            patch.set_facecolor(color_dict[condition])
        # Create legend for colors
        legend_labels = [mpatches.Patch(color=color_dict[i], label=i) for i in color_dict]
        plt.legend(handles=legend_labels, loc='center', bbox_to_anchor=(1.1, .5))
        sns.despine()
        plt.tight_layout()
        if save is not None:
            if vector is True:
                plt.savefig(save + 'boxplot_normalization.svg', bbox_inches='tight')
            else:
                plt.savefig(save + 'boxplot_normalization.png', dpi=dpi, bbox_inches='tight')
        plt.show()

def volcano_Multicond(self, *Proteins, pvalue=0.05, palette='viridis',
                      bcol='#962558', non_regulated='#606060',
                      save=None, dpi=300, vector=True):
    """Creates a volcano plot for multiple conditions.

    In general, volcano plots are designed to compare 2 conditions.
    Here, we aim to see the distribution of quantified proteins' p-value and
    fold changes among multiple conditions.

    Args:
        pvalue (float, optional): p-value threshold. Defaults to 0.05.
        palette (str, optional): Color palette to differentiate dots. Defaults to 'viridis'.
        bcol (str, optional): color for density plot. Defaults to '#962558'.
        non_regulated (str, optional): Proteins not differentially regulated. Defaults to '#606060'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to True.
    """
    plt.style.use('default')
    plt.rcParams["figure.dpi"] = dpi
    OmicScope = copy.copy(self)
    FoldChange_cutoff = OmicScope.FoldChange_cutoff
    # Definitions for the axes
    df_initial = OmicScope.quant_data
    dropped_inf = copy.copy(df_initial)
    dropped_inf = dropped_inf.replace([np.inf, -np.inf], np.nan)
    dropped_inf = dropped_inf.dropna()
    fc = df_initial['log2(fc)']
    fc = fc.replace(np.inf, dropped_inf['log2(fc)'].max() + 1)
    fc = fc.replace(-np.inf, dropped_inf['log2(fc)'].min() - 1)
    pval = df_initial[f'-log10({OmicScope.pvalue})']
    pval = pval.replace(
        np.inf, dropped_inf[f'-log10({OmicScope.pvalue})'].max() + 1)
    df_initial[f'-log10({OmicScope.pvalue})'] = pval
    df_initial['log2(fc)'] = fc

    # colors per condition
    comparisons = df_initial.Comparison
    comparisons = comparisons.str[-1] + '-' + comparisons.str[0]
    number_of_comparison = len(comparisons.drop_duplicates())
    color_per_comparison = sns.color_palette(palette=palette,
                                             n_colors=number_of_comparison).as_hex()
    color_comp_dict = dict(zip(comparisons, color_per_comparison))
    print(color_comp_dict)
    col = []
    comparison = []
    for pv, comp in zip(pval, comparisons):
        if pv >= -np.log10(pvalue):
            comparison.append(comp)
            col.append(comp)
        else:
            col.append(non_regulated)
            comparison.append('Non-regulated')

    col = pd.Series(col)
    comparison = pd.Series(comparison)
    df = pd.DataFrame(data=zip(pval, fc, col, comparison))
    df.columns = [OmicScope.pvalue, 'fc', 'col', 'comparison']
    df[OmicScope.pvalue] = df[OmicScope.pvalue].astype(float)
    df['fc'] = df['fc'].astype(float)
    df['col'] = df.col.replace(color_comp_dict)
    # annotation if it is on args

    if len(Proteins) > 0:
        ls = pd.DataFrame(Proteins)
        ls.columns = ['gene_name']
        ls_final = ls.merge(df_initial)
        ls_final = ls_final[['gene_name', f'-log10({OmicScope.pvalue})', 'log2(fc)']]
        ls_final = ls_final.sort_values('log2(fc)', ascending=False)

    # dimensions
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    # plot
    plt.style.use('bmh')
    plt.figure(figsize=(8, 8))
    # Scatter plot
    ax_scatter = plt.axes(rect_scatter)
    sns.scatterplot(x=df.fc, y=df[OmicScope.pvalue], alpha=0.7,
                    hue=df.comparison)
    sns.despine(left=False)
    plt.xlabel("log2( FC)")
    plt.ylabel(f'-log10({OmicScope.pvalue})')
    ax_scatter.tick_params(direction='in', top=True, right=True)
    plt.grid(b=None)
    ax_scatter.plot(alpha=0.5)
    # Density plot on x-axis (upper)
    ax_histx = plt.axes(rect_histx)
    plt.ylabel("Density")
    plt.title(label=OmicScope.ctrl + ' vs ' + '-'.join(OmicScope.experimental), loc='left')
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histx.set_facecolor('white')
    ax_histx.grid(b=None)
    sns.kdeplot(fc, fill=True, color=bcol, edgecolor='black')
    sns.despine(bottom=True, top=True)
    # Density plot on y-axis (right)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)
    ax_histy.set_facecolor('white')
    ax_histy.grid(b=None)
    sns.kdeplot(y=pval, fill=True, color=bcol, edgecolor='black')
    sns.despine(top=True, left=True)
    plt.xlabel("Density")
    ax_scatter.set_facecolor('white')
    # Annotation positions
    if len(Proteins) > 0:
        from adjustText import adjust_text
        texts = []
        for a, b, c in zip(ls_final['log2(fc)'], ls_final[f'-log10({OmicScope.pvalue})'],
                           ls_final['gene_name']):
            texts.append(ax_scatter.text(a, b, c, size=8))
        adjust_text(texts, ax=ax_scatter, force_points=0.25, force_text=0.25,
                    expand_points=(1.5, 1.5), expand_text=(1.5, 1.5),
                    arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    ax_scatter.axhline(y=-np.log10(0.05), color='gray', linestyle=':')
    ax_scatter.legend()
    if FoldChange_cutoff == 0:
        ax_scatter.axvline(x=0, color='gray', linestyle=':')
    else:
        ax_scatter.axvline(x=FoldChange_cutoff, color='gray', linestyle=':')
    limh = round(fc.max(), 2)
    liml = round(fc.min(), 2)
    limp = round(pval.max() + 0.5, 0)
    ax_scatter.set_xlim((liml - (limh - liml) * 0.10, limh + ((limh - liml) * 0.10)))
    ax_scatter.set_ylim((pval.min(), limp + limp * .1))
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    if save is not None:
        if vector is True:
            plt.savefig(save + 'volcano.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'volcano.png', dpi=dpi, bbox_inches='tight')
    plt.show()
    plt.show(block=True)


def volcano_2cond(self, *Proteins, pvalue=0.05, bcol='darkcyan',
                  non_regulated='#606060', up_regulated='#E4001B', down_regulated='#6194BC',
                  save=None, dpi=300, vector=True):
    """Plot a conventional volcano plot

    Args:
        pvalue (float, optional): p-value threshold. Defaults to 0.05.
        bcol (str, optional): color for density plot. Defaults to 'darkcyan'.
        non_regulated (str, optional): Proteins not differentially regulated.
         Defaults to '#606060'.
        up_regulated (str, optional): Proteins up-regulated in relation
         to Control group. Defaults to '#E4001B'.
        down_regulated (str, optional): Proteins down-regulated in relation
         to Control group. Defaults to '#6194BC'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.style.use('default')
    plt.rcParams["figure.dpi"] = dpi
    OmicScope = copy.copy(self)
    FoldChange_cutoff = OmicScope.FoldChange_cutoff
    # Definitions for the axes
    df_initial = OmicScope.quant_data
    dropped_inf = copy.copy(df_initial)
    dropped_inf = dropped_inf.replace([np.inf, -np.inf], np.nan)
    dropped_inf = dropped_inf.dropna()
    fc = df_initial['log2(fc)']
    fc = fc.replace(np.inf, dropped_inf['log2(fc)'].max() + 1)
    fc = fc.replace(-np.inf, dropped_inf['log2(fc)'].min() - 1)
    pval = df_initial[f'-log10({OmicScope.pvalue})']
    pval = pval.replace(
        np.inf, dropped_inf[f'-log10({OmicScope.pvalue})'].max() + 1)
    df_initial[f'-log10({OmicScope.pvalue})'] = pval
    df_initial['log2(fc)'] = fc
    # Defining colors for dots
    col = []
    comparison = []
    for i, j in zip(fc, pval):
        if j <= -np.log10(pvalue):
            col.append(non_regulated)
            comparison.append('non-regulated')
        else:
            if i >= FoldChange_cutoff:
                col.append(up_regulated)
                comparison.append('up-regulated')
            elif i <= -FoldChange_cutoff:
                col.append(down_regulated)
                comparison.append('down-regulated')
            else:
                col.append(non_regulated)
                comparison.append('non-regulated')
    col, comparison = pd.Series(col), pd.Series(comparison)
    df = pd.DataFrame(zip(pval, fc, col, comparison))
    df.columns = ['pval', 'fc', 'col', 'Regulation']
    # annotation
    if len(Proteins) > 0:
        ls = pd.DataFrame(Proteins)
        ls.columns = ['gene_name']
        ls_final = ls.merge(df_initial)
        ls_final = ls_final[['gene_name', f'-log10({OmicScope.pvalue})', 'log2(fc)']]
        ls_final = ls_final.sort_values('log2(fc)', ascending=False)
    # plot dimensions
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    # Plot
    plt.style.use('bmh')
    plt.figure(figsize=(8, 8))
    ax_scatter = plt.axes(rect_scatter)
    # Scatter plot
    sns.scatterplot(x=df.fc, y=df.pval, alpha=0.5, linewidth=0.5,
                    hue=df.Regulation,
                    palette=list(df.col.drop_duplicates().dropna()))
    plt.xlabel("log2( FC)")
    plt.ylabel(f'-log10({OmicScope.pvalue})')
    ax_scatter.tick_params(direction='in', top=True, right=True)
    plt.grid(b=None)
    ax_scatter.plot(alpha=0.5)
    # Density Plot in x-axis (upper)
    ax_histx = plt.axes(rect_histx)
    plt.ylabel("Density")
    plt.title(label=OmicScope.experimental[0] + ' - ' + OmicScope.ctrl, loc='left')
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histx.set_facecolor('white')
    ax_histx.grid(b=None)
    sns.kdeplot(fc, fill=True, color=bcol, alpha=0.8, edgecolor='black')
    sns.despine(bottom=True, top=True)
    # Density plot in y axis (right)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)
    ax_histy.set_facecolor('white')
    ax_histy.grid(b=None)
    sns.kdeplot(y=pval, fill=True, color=bcol, alpha=0.8, edgecolor='black')
    sns.despine(top=True, left=True)
    plt.xlabel("Density")
    ax_scatter.set_facecolor('white')
    # Annotation positions
    if len(Proteins) > 0:
        from adjustText import adjust_text
        texts = []
        for a, b, c in zip(ls_final['log2(fc)'], ls_final[f'-log10({OmicScope.pvalue})'],
                           ls_final['gene_name']):
            texts.append(ax_scatter.text(a, b, c, size=8))
        adjust_text(texts, ax=ax_scatter,
                    arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    ax_scatter.axhline(y=-np.log10(0.05), color='gray', linestyle=':')

    # Fold change cutoff
    if FoldChange_cutoff != 0:
        ax_scatter.axvline(x=-FoldChange_cutoff,
                           color='gray', linestyle=':')
        ax_scatter.axvline(x=FoldChange_cutoff,
                           color='gray', linestyle=':')
    else:
        ax_scatter.axvline(x=FoldChange_cutoff,
                           color='gray', linestyle=':')
    # Figure limits
    limh = round(fc.max(), 2)
    liml = round(fc.min(), 2)
    limp = round(pval.max() + 0.5, 0)
    ax_scatter.set_xlim((liml - ((limh - liml) * 0.10),
                         limh + ((limh - liml) * 0.10)))
    ax_scatter.set_ylim((pval.min(), limp))
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    # save figure and how to save
    if save is not None:
        if vector is True:
            plt.savefig(save + 'volcano.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'volcano.png', dpi=dpi, bbox_inches='tight')
    plt.show()
    plt.show(block=True)


def volcano(self, *Proteins,
            pvalue=0.05, bcol='#962558', palette='viridis',
            non_regulated='#606060', up_regulated='#E4001B', down_regulated='#6194BC',
            save=None, dpi=300, vector=True):
    """Plot volcano plot.

    Args:
        pvalue (float, optional): p-value threshold. Defaults to 0.05.
        bcol (str, optional): Density plot color. Defaults to '#962558'.
        palette (str, optional): Palette for Multiconditions volcano plot.
         Defaults to 'viridis'.
        non_regulated (str, optional): Color of non-differentially regulated
         proteins. Defaults to '#606060'.
        up_regulated (str, optional): Color of up-regulated proteins for volcano
         with 2 conditions. Defaults to '#E4001B'.
        down_regulated (str, optional): Color of down-regulated proteins for volcano
         with 2 conditions. Defaults to '#6194BC'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): . Resolution to save figure. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    OmicScope = copy.copy(self)
    if len(OmicScope.Conditions) == 2:
        volcano_2cond(OmicScope, *Proteins, pvalue=pvalue, bcol=bcol,
                      non_regulated=non_regulated, up_regulated=up_regulated,
                      down_regulated=down_regulated,
                      save=save, dpi=dpi, vector=vector)
    if len(OmicScope.Conditions) > 2:
        volcano_Multicond(OmicScope, *Proteins,
                          pvalue=pvalue, palette=palette, bcol=bcol,
                          save=None, dpi=dpi, vector=vector)


def heatmap(self, *Proteins, pvalue=0.05, c_cluster=True,
            clust_metric="euclidean", clust_method='average',
            palette='RdYlBu_r', linewidth=0.01, color_groups='tab20',
            sample_label=False,
            save=None, dpi=300, vector=True):
    """ Heatmap with hierarchical clustering

    Args:
        pvalue (float, optional): p-value threshold. Defaults to 0.05.
        c_cluster (bool, optional): Applies Hierarchical clustering for
         columns. Defaults to True.
        clust_metric (str, optional): The distance metric to use. Optionally: braycurtis, canberra,
         chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard,
         jensenshannon, kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto, russellrao,
         seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule.
        clust_method (str, optional): Linkage method to use for calculating clusters.
         Optionally: single, complete, average, weighted, centroid, median, or ward.
        color_groups (str, optional): Palette for group colors.
         Defaults to 'tab20'.
        palette (str, optional): Palette for protein abundance. Defaults to 'RdYlBu_r'.
        linewidth (float, optional): Line width. Defaults to 0.01.
        sample_label (bool, optional): insert biological sample label and condition.
        Defaults to False.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.style.use('default')
    plt.rcParams["figure.dpi"] = dpi
    OmicScope = copy.copy(self)
    FoldChange_cutoff = OmicScope.FoldChange_cutoff
    df = OmicScope.quant_data
    pdata = OmicScope.pdata
    sample_condition = dict(zip(pdata.Sample, pdata.Condition))
    sort_columns = list(pdata.sort_values(['Condition'])['Sample'])
    if 'TimeCourse' in pdata.columns:
        sort_columns = list(pdata.sort_values(
            ['Condition', 'TimeCourse'])['Sample'])
        times = pdata.TimeCourse.drop_duplicates()
        ntimes = len(times)
        pal = sns.color_palette(palette='BuPu', as_cmap=False, n_colors=ntimes)
        dic = {}
        for C, c in zip(times, pal):
            dic.update({C: c})
        replacer = dic.get
        time_colors = [replacer(n, n) for n in list(
            pdata.sort_values(['Condition', 'TimeCourse']).TimeCourse)]
    else:
        time_colors = []
    Conditions = OmicScope.Conditions
    # Creating Heatmap Matrix
    heatmap = df[df.loc[:, OmicScope.pvalue] <= pvalue]
    heatmap = heatmap.loc[(heatmap['log2(fc)'] <= -FoldChange_cutoff)
                          | (heatmap['log2(fc)'] >= FoldChange_cutoff)]
    heatmap = heatmap.dropna()
    # Filtering for specific proteins
    if len(Proteins) > 0:
        heatmap = heatmap[heatmap['gene_name'].isin(Proteins)]
    heatmap = heatmap.set_index('gene_name')
    heatmap = heatmap.loc[:, heatmap.columns.str.contains('.', regex=False)]
    # Log2 transform
    heatmap = np.log2(heatmap).replace([-np.inf], int(0))
    heatmap = heatmap[sort_columns]
    heatmap_original_label = heatmap.copy()
    heatmap = heatmap.rename(columns=sample_condition)
    Col = list(heatmap.columns)
    # Creating matrix for group colors
    ngroup = len(Conditions)
    pal = sns.color_palette(palette=color_groups,
                            as_cmap=False, n_colors=ngroup)
    dic = {}
    for C, c in zip(Conditions, pal):
        dic.update({C: c})
    replacer = dic.get
    colcolors = [replacer(n, n) for n in Col]
    colors = [colcolors] + [time_colors]
    if colors[-1] == []:
        colors = colcolors
    # Title
    title = 'Heatmap - ' + OmicScope.ctrl + \
        ' vs ' + '-'.join(OmicScope.experimental)
    # Plot
    heatmap.columns.name = 'Samples'
    heatmap.index.name = 'Gene name'
    if sample_label is True:
        heatmap = heatmap_original_label
    sns.clustermap(heatmap,
                   cmap=palette, z_score=0, linewidths=linewidth, linecolor='black',
                   col_colors=colors, metric=clust_metric, method=clust_method,
                   col_cluster=c_cluster, center=0).fig.suptitle(title, y=1.02)
    
    # Create legend for colors
    if type(colors[0]) is list:
        col_dict_Time = {i:j for i,j in zip(pdata.TimeCourse, colors[1])}
        legend_labels_time = [mpatches.Patch(color=col_dict_Time[i], label=i) for i in col_dict_Time]
        legend1 = plt.legend(handles=legend_labels_time, bbox_to_anchor=(1, 1),
                              bbox_transform=plt.gcf().transFigure, loc='upper right',
                              title='Time Point')

        col_dict = {i:j for i,j in zip(pdata.Condition, colors[0])}
        legend_labels = [mpatches.Patch(color=col_dict[i], label=i) for i in col_dict]
        legend2 = plt.legend(handles=legend_labels, bbox_to_anchor=(1, 0.88),
                              bbox_transform=plt.gcf().transFigure, loc='upper right',
                              title='Condition')
        plt.gca().add_artist(legend1)

    else:
        col_dict = {i:j for i,j in zip(heatmap.columns, colors)}
        legend_labels = [mpatches.Patch(color=col_dict[i], label=i) for i in col_dict]
        plt.legend(handles=legend_labels, bbox_to_anchor=(1, 1),
                    bbox_transform=plt.gcf().transFigure, loc='upper right',
                    title='Condition')

    # Save
    if save is not None:
        if vector is True:
            plt.savefig(save + 'heatmap.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'heatmap.png', dpi=dpi, bbox_inches='tight')
    plt.show()

def correlation(self, *Proteins, pvalue=1.0,
                sample_method='pearson',
                clust_metric='euclidean',
                clust_method='average', palette='RdYlBu_r', linewidth=0.005,
                color_groups='tab20', sample_label=False,
                save=None, dpi=300, vector=True):
    """Pairwise correlation plot for samples.

    Args:
        pvalue (float, optional): p-value threshold. Defaults to 1.0.
        sample_method (str, optional): Method to compute pair-wise correlation.
         Options: pearson, kendal, spearman. Defaults to 'pearson'.
        clust_metric (str, optional): The distance metric to use. Optionally: braycurtis, canberra,
         chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard,
         jensenshannon, kulczynski1, mahalanobis, matching, minkowski, rogerstanimoto, russellrao,
         seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule.
        clust_method (str, optional): Linkage method to use for calculating clusters.
         Optionally: single, complete, average, weighted, centroid, median, or ward.
        palette (str, optional): Palette for R-distribution. Defaults to 'RdYlBu_r'.
        linewidth (float, optional): Line width. Defaults to 0.005.
        sample_label (bool, optional): insert biological sample label and condition.
        Defaults to False.
        color_groups (str, optional): Color of each group. Defaults to 'tab20'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.style.use('default')
    plt.rcParams["figure.dpi"] = dpi
    OmicScope = copy.copy(self)
    FoldChange_cutoff = OmicScope.FoldChange_cutoff
    df = OmicScope.quant_data
    pdata = OmicScope.pdata
    sample_condition = dict(zip(pdata.Sample, pdata.Condition))
    sort_columns = list(pdata.sort_values(['Condition'])['Sample'])
    if 'TimeCourse' in pdata.columns:
        sort_columns = list(pdata.sort_values(
            ['Condition', 'TimeCourse'])['Sample'])
        times = pdata.TimeCourse.drop_duplicates()
        ntimes = len(times)
        pal = sns.color_palette(palette='BuPu', as_cmap=False, n_colors=ntimes)
        dic = {}
        for C, c in zip(times, pal):
            dic.update({C: c})
        replacer = dic.get
        time_colors = [replacer(n, n) for n in list(
            pdata.sort_values(['Condition', 'TimeCourse']).TimeCourse)]
    else:
        time_colors = []
    Conditions = OmicScope.Conditions
    # Selecting Matrix for correlation
    pearson = df[df.loc[:, OmicScope.pvalue] <= pvalue]
    pearson = pearson.loc[(pearson['log2(fc)'] <= -FoldChange_cutoff)
                          | (pearson['log2(fc)'] >= FoldChange_cutoff)]
    pearson = pearson.dropna()
    # Filtering for specific proteins
    if len(Proteins) > 0:
        pearson = pearson[pearson['gene_name'].isin(Proteins)]
    pearson = pearson.set_index('gene_name')
    pearson = pearson.loc[:, pearson.columns.str.contains('.', regex=False)]
    # log2 transform
    pearson = np.log2(pearson).replace([-np.inf], int(0))
    # Performing Pearson's Correlation
    corr_matrix_raw = pearson.corr(method=sample_method)
    # Creating matrix for group colors
    corr_matrix_raw = corr_matrix_raw[sort_columns]
    corr_matrix = corr_matrix_raw.rename(columns=sample_condition)
    Col = list(corr_matrix.columns)
    # Creating matrix for group colors
    ngroup = len(Conditions)
    pal = sns.color_palette(palette=color_groups, as_cmap=False, n_colors=ngroup)
    dic = {}
    for C, c in zip(Conditions, pal):
        dic.update({C: c})
    replacer = dic.get
    colcolors = [replacer(n, n) for n in Col]
    colors = [colcolors] + [time_colors]
    if colors[-1] == []:
        colors = colcolors
    # Title
    title = 'Heatmap - ' + OmicScope.ctrl + ' vs ' + '-'.join(OmicScope.experimental)
    # Plot
    corr_matrix.columns.name = 'Samples'
    corr_matrix.index.name = 'Samples'
    if sample_label is False:
        sns.clustermap(corr_matrix, metric=clust_metric, method=clust_method,
                       xticklabels=corr_matrix.columns, row_colors=colors,
                       yticklabels=corr_matrix.columns, col_colors=colors,
                       cmap=palette, linewidths=linewidth, linecolor='black').fig.suptitle(title, y=1.02)
    else:
        sns.clustermap(corr_matrix_raw, metric=clust_metric, method=clust_method,
                       xticklabels=corr_matrix_raw.columns, row_colors=colors,
                       yticklabels=corr_matrix_raw.columns, col_colors=colors,
                       cmap=palette, linewidths=linewidth, linecolor='black').fig.suptitle(title, y=1.02)
        # Create legend for colors
    if type(colors[0]) is list:
        col_dict_Time = {i:j for i,j in zip(pdata.TimeCourse, colors[1])}
        legend_labels_time = [mpatches.Patch(color=col_dict_Time[i], label=i) for i in col_dict_Time]
        legend1 = plt.legend(handles=legend_labels_time, bbox_to_anchor=(1, 1),
                              bbox_transform=plt.gcf().transFigure, loc='upper right',
                              title='Time Point')

        col_dict = {i:j for i,j in zip(pdata.Condition, colors[0])}
        legend_labels = [mpatches.Patch(color=col_dict[i], label=i) for i in col_dict]
        legend2 = plt.legend(handles=legend_labels, bbox_to_anchor=(1, 0.88),
                              bbox_transform=plt.gcf().transFigure, loc='upper right',
                              title='Condition')
        plt.gca().add_artist(legend1)
    else:
        col_dict = {i:j for i,j in zip(pdata.Condition, colors)}
        legend_labels = [mpatches.Patch(color=col_dict[i], label=i) for i in col_dict]
        plt.legend(handles=legend_labels, bbox_to_anchor=(1, 1),
                    bbox_transform=plt.gcf().transFigure, loc='upper right',
                    title='Condition')

    if save is not None:
        if vector is True:
            plt.savefig(save + 'pearson.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'pearson.png', dpi=dpi, bbox_inches='tight')
    plt.xlabel('')
    plt.ylabel('')
    plt.show()
    plt.show(block=True)


def DynamicRange(self, *Proteins, color='#565059',
                 protein_color='orange', max_min=False,
                 min_color='#18ab75', max_color='#ab4e18', dpi=300, save=None,
                 vector=True):
    """Dynamic range plot

    Args:
        mean_color (str, optional): default color for dots (mean abundance).
         Defaults to '#565059'.
        protein_color (str, optional): Color for specific proteins (Args).
         Defaults to 'orange'.
        max_min (bool, optional): Plot the maximum and minimum abundance
         value for each protein. Defaults to False.
        min_color (str, optional): Color of minimum values. Defaults to '#1daec7'.
        max_color (str, optional): Color of maximum values. Defaults to '#f7463d'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """

    plt.style.use('default')
    plt.rcParams["figure.dpi"] = dpi
    OmicScope = copy.copy(self)
    df = OmicScope.quant_data
    # Dictionary for Accessions
    accession = dict(zip(df.Accession, df.gene_name))
    df = df.set_index('Accession')
    df = df.loc[:, df.columns.str.contains('.', regex=False)]
    df = np.log10(df).transpose()
    # Getting Mean, Max and Min values for each protein
    df = df.describe().transpose()
    df = df.dropna()
    df = df.reset_index()
    # Applying Dictionary to get gene name
    df['gene_name'] = df.Accession.replace(accession)
    # Ranking entities
    df['rank'] = df['mean'].rank(method='first')
    # Return variation for each protein
    ax_scatter = plt.axes()
    if max_min is True:
        plt.scatter(y=df['rank'], x=df['min'],
                    alpha=0.7, s=2, c=min_color)
        plt.scatter(y=df['rank'], x=df['max'],
                    alpha=0.7, s=2, c=max_color)
    # Plot mean
    plt.scatter(y=df['rank'], x=df['mean'],
                c=color, s=10, alpha=0.5,
                linewidths=0.5)
    # Highlight specific proteins
    if len(Proteins) > 0:
        df = df[df.gene_name.isin(Proteins)]
        plt.scatter(y=df['rank'], x=df['mean'],
                    c=protein_color, s=10, alpha=1, edgecolors='black',
                    linewidths=0.5)
        ls_final = df
        ls_final = ls_final[['gene_name', 'rank', 'mean']]
        from adjustText import adjust_text
        texts = []
        for a, b, c in zip(ls_final['mean'], ls_final['rank'],
                           ls_final['gene_name']):
            texts.append(ax_scatter.text(a, b, c, size=8))
        adjust_text(texts, ax=ax_scatter, force_points=0.25, force_text=0.25,
                    expand_points=(1.5, 1.5), expand_text=(1.5, 1.5),
                    arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    plt.grid(b=False)
    plt.xlabel('log10(Abundance)')
    plt.ylabel('Rank')
    plt.title('Dynamic Range')
    if save is not None:
        if vector is True:
            plt.savefig(save + 'DynamicRange.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'DynamicRange.png', dpi=dpi, bbox_inches='tight')
    sns.despine()
    plt.show()


def pca(self, pvalue=1.00, scree_color='#900C3F',
        sample_label=False, palette='tab20', FoldChange_cutoff=0,
        save=None, dpi=300, vector=True):
    """ Perform Principal Component Analysis.

    Args:
        pvalue (float, optional): p-value threshold. Defaults to 1.00.
        scree_color (str, optional): Color of Scree plot. Defaults to '#900C3F'.
        sample_label (bool, optional): Insert biological sample label. Defaults to False.
        palette (str, optional): Palette for groups. Defaults to 'tab20'.
         FoldChange_cutoff (int, optional): Fold change threshold. Defaults to 0.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.style.use('default')
    plt.rcParams["figure.dpi"] = dpi
    OmicScope = copy.copy(self)
    df = OmicScope.quant_data
    FoldChange_cutoff = OmicScope.FoldChange_cutoff
    # Filtering P-value and Fold Change
    df = df[df[OmicScope.pvalue] <= pvalue]
    if len(OmicScope.Conditions) == 2:
        df = df.loc[(df['log2(fc)'] <= -FoldChange_cutoff) | (df['log2(fc)'] >= FoldChange_cutoff)]
    df = df.loc[:, df.columns.str.contains('.', regex=False)]
    samples = df.columns
    # Getting Conditions
    Conditions = OmicScope.Conditions
    group = []
    for i in Conditions:
        ncond = df.columns.str.contains(i).sum()
        nlist = list(itertools.repeat(i, ncond))
        group.extend(nlist)
    # Center and Scale data
    scaled_data = preprocessing.scale(df.T)
    pca = PCA()  # PCA object
    pca.fit(scaled_data)
    # get PCA coordinates for scaled_data
    pca_data = pca.transform(scaled_data)
    # Scree plot
    per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
    labels = ['PC' + str(x) for x in range(1, len(per_var) + 1)]
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_figwidth(16)
    fig.set_figheight(7)
    ax1.bar(x=range(1, len(per_var) + 1), height=per_var, tick_label=labels,
            color=scree_color, edgecolor="black", linewidth=1)
    ax1.set(ylabel='Percentage of Explained Variance',
            xlabel='Principal Component', title='Scree plot')
    ax1.tick_params(labelrotation=45)
    sns.despine()
    ax1.grid(b=False)
    # PC1 vs PC2 plot
    pca_df = pd.DataFrame(pca_data, index=samples, columns=labels)
    pdata = OmicScope.pdata
    conditions = pdata.set_index('Sample')[['Condition']]
    if 'TimeCourse' in pdata.columns:
        conditions = pdata.set_index('Sample')[['Condition', 'TimeCourse']]
        pca_df = pca_df.merge(conditions, right_index=True, left_index=True)
        pca_df = pca_df.set_index(['Condition', 'TimeCourse'])
    else:
        pca_df = pca_df.merge(conditions, right_index=True, left_index=True)
        pca_df = pca_df.set_index(['Condition'])

    ax_scatter = sns.scatterplot(x="PC1", y="PC2", data=pca_df,
                                 hue=list(pca_df.index), palette=palette,
                                 edgecolor='black', linewidth=1, s=100, alpha=0.7)
    plt.title('Principal Component Analysis ' + '(pvalue=' + str(pvalue) + ')')
    plt.xlabel('PC1 - {0}%'.format(per_var[0]))
    plt.ylabel('PC2 - {0}%'.format(per_var[1]))
    plt.grid(b=False)
    sns.despine()
    # Annotate samples in PCA plot
    if sample_label is True:
        from adjustText import adjust_text
        texts = []
        for a, b, c in zip(pca_df['PC1'], pca_df['PC2'],
                           conditions.index):
            texts.append(ax_scatter.text(a, b, c.split('.')[0], size=8))
        adjust_text(texts, ax=ax_scatter, force_points=0.25, force_text=0.25,
                    expand_points=(1.5, 1.5), expand_text=(1.5, 1.5),
                    arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    if save is not None:
        if vector is True:
            plt.savefig(save + 'pca.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'pca.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def color_scheme(df, palette):
    """Generate colors to barplot and

    Args:
        df (DataFrame): dataframe
        palette (str): palette

    Returns:
        color (dict): dictionary assign each condition with respective color.
    """
    color = df.index.to_frame().reset_index(drop=True)
    color['variable'] = color.astype(str).apply(lambda x: '_'.join(x), axis=1)
    color = color.drop_duplicates().sort_values('variable')
    colors_scheme = ['red', 'blue', 'orange', 'green', 'darkcyan', 'magenta']
    if 'TimeCourse' in color.columns:
        timecourse = len(color.TimeCourse.drop_duplicates())
    else:
        timecourse = 1
    random.seed(10)
    colors_random = random.sample(colors_scheme, color.Condition.nunique())
    palettes = []
    if timecourse > 1:
        for i in colors_random:
            palette = sns.light_palette(i, n_colors=timecourse).as_hex()
            palettes.append(palette)
    else:
        ncolors = color.Condition.nunique()
        palettes = sns.color_palette(palette=palette, n_colors=ncolors).as_hex()
    color['colors'] = list(pd.Series(palettes).explode())
    color = dict(zip(color.variable, color.colors))
    return color


def bar_protein(self, *Proteins, logscale=True,
                palette='Spectral', save=None, dpi=300,
                vector=True):
    """Bar plot to show protein abundance in each condition

    Args:
        logscale (bool, optional): Apply log-transformed abundance.
        Defaults to True.
        palette (str, optional): Palette for groups. Defaults to 'Spectral'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.rcParams['figure.dpi'] = dpi
    OmicScope = copy.copy(self)
    df = copy.copy(OmicScope.quant_data)
    # Proteins to plot
    df = df[df['gene_name'].isin(Proteins)]
    # Get conditions
    ctrl = [copy.copy(OmicScope.ctrl)]
    conditions = copy.copy(OmicScope.Conditions)
    conditions.remove(ctrl[0])
    # Get protein abundance for each condition
    df = df.set_index('gene_name')
    df = df.iloc[:, df.columns.str.contains('.', regex=False)]
    df = df.melt(var_name='variable', ignore_index=False)
    df = df.reset_index()
    pdata = copy.copy(OmicScope.pdata)
    pdata = pdata.set_index('Sample')
    df = df.set_index('variable')
    df = pdata.merge(df, left_index=True, right_index=True)
    if 'TimeCourse' in df.columns:
        df = df[['Condition',  'TimeCourse', 'gene_name', 'value']]
        df = df.set_index(['Condition', 'TimeCourse'])
    else:
        df = df[['Condition', 'gene_name', 'value']]
        df = df.set_index('Condition')
    # Apply log transformation
    if logscale is True:
        df['value'] = np.log2(df['value'])
        df['value'] = df['value'].replace(-np.inf, np.nan)

    color = color_scheme(df, palette=palette)
    df_index = df.index.to_frame().reset_index(drop=True)
    df_index = df_index.astype(str).apply(lambda x: '_'.join(x), axis=1)
    df.index = df_index
    # Size of figure
    sns.set(rc={'figure.figsize': (1.5*len(Proteins), 5)})
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    # Plot
    if len(Proteins) == 1:
        sns.barplot(x=df.index, y='value',
                    data=df, errwidth=1, capsize=0.07,
                    palette=color,
                    errorbar='se',
                    edgecolor='black', linewidth=1, dodge=False)
    else:
        sns.barplot(x='gene_name', y='value', hue=list(df.index),
                    data=df, errwidth=1,
                    capsize=0.07, edgecolor='black',
                    errorbar='se',
                    palette=color,
                    linewidth=1)
    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
    sns.despine()
    plt.title('Abundance - ' + ' and '.join(Proteins))
    plt.xlabel('')
    plt.xticks(rotation=70)
    plt.ylabel('log2(Abundance)')
    if save is not None:
        if vector is True:
            plt.savefig(save + 'barplot_' + '_'.join(Proteins) + '.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'barplot_' + '_'.join(Proteins) + '.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def boxplot_protein(self, *Proteins, logscale=True,
                    palette='Spectral', save=None, dpi=300,
                    vector=True):
    """Boxplot to show protein abundance in each condition

    Args:
        logscale (bool, optional): Apply abundance log-transformed.
        Defaults to True.
        palette (str, optional): Palette for groups. Defaults to 'Spectral'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
        True.
    """
    plt.rcParams['figure.dpi'] = dpi
    OmicScope = copy.copy(self)
    df = copy.copy(OmicScope.quant_data)
    # Proteins to plot
    df = df[df['gene_name'].isin(Proteins)]
    # Get conditions
    ctrl = [copy.copy(OmicScope.ctrl)]
    conditions = copy.copy(OmicScope.Conditions)
    conditions.remove(ctrl[0])
    # Get protein abundance for each condition
    df = df.set_index('gene_name')
    df = df.iloc[:, df.columns.str.contains('.', regex=False)]
    df = df.melt(var_name='variable', ignore_index=False)
    df = df.reset_index()
    pdata = copy.copy(OmicScope.pdata)
    pdata = pdata.set_index('Sample')
    df = df.set_index('variable')
    df = pdata.merge(df, left_index=True, right_index=True)
    if 'TimeCourse' in df.columns:
        df = df[['Condition',  'TimeCourse', 'gene_name', 'value']]
        df = df.set_index(['Condition', 'TimeCourse'])
    else:
        df = df[['Condition', 'gene_name', 'value']]
        df = df.set_index('Condition')
    # Apply log transformation
    if logscale is True:
        df['value'] = np.log2(df['value'])
        df['value'] = df['value'].replace(-np.inf, np.nan)

    color = color_scheme(df, palette=palette)
    df_index = df.index.to_frame().reset_index(drop=True)
    df_index = df_index.astype(str).apply(lambda x: '_'.join(x), axis=1)
    df.index = df_index
    # Size of figure
    sns.set(rc={'figure.figsize': (1.5*len(Proteins), 5)})
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    # Plot
    if len(Proteins) == 1:
        sns.boxplot(x=df.index, y='value',
                    data=df,
                    palette=color,
                    linewidth=1, dodge=False)
    else:
        sns.boxplot(x='gene_name', y='value', hue=list(df.index),
                    data=df,
                    palette=color,
                    linewidth=1)
    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
    sns.despine()
    plt.title('Abundance - ' + ' and '.join(Proteins))
    plt.xlabel('')
    plt.xticks(rotation=70)
    plt.ylabel('log2(Abundance)')
    if save is not None:
        if vector is True:
            plt.savefig(save + 'barplot_' + '_'.join(Proteins) + '.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'barplot_' + '_'.join(Proteins) + '.png', dpi=dpi, bbox_inches='tight')
    plt.show()


def MAplot(self, *Proteins,
           pvalue=0.05, non_regulated='#606060', up_regulated='#E4001B',
           down_regulated='#6194BC', FoldChange_cutoff=0,
           save=None, dpi=300, vector=True):
    """MA plot

    Args:
        pvalue (float, optional): p-value threshold. Defaults to 0.05.
        non_regulated (str, optional): color for non-regulated proteins.
         Defaults to '#606060'.
        up_regulated (str, optional): color for up-regulated proteins.
         Defaults to '#E4001B'.
        down_regulated (str, optional): color for down-regulated proteins.
         Defaults to '#6194BC'.
        FoldChange_cutoff (int, optional): Foldchange threshold. Defaults to 0.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.
    """
    plt.rcParams['figure.dpi'] = dpi
    OmicScope = copy.copy(self)
    df = copy.copy(OmicScope.quant_data)
    df['TotalMean'] = np.log2(df['TotalMean'])
    df['TotalMean'] = df['TotalMean'].replace(-np.inf, 0.01)
    # Defining axis
    y = df[f'-log10({OmicScope.pvalue})']
    x = df['log2(fc)']
    # Defining colors
    col = np.where(y >= -np.log10(pvalue), np.where(x >= -FoldChange_cutoff,
                                                   np.where(x <= FoldChange_cutoff, non_regulated, up_regulated),
                                                   down_regulated), non_regulated)
    regulation = np.where(y >= -np.log10(pvalue), np.where(x >= -FoldChange_cutoff,
                                                          np.where(x <= FoldChange_cutoff, 'non-regulated',
                                                                   'up-regulated'), 'down-regulated'), 'non-regulated')
    df['col'] = col
    df['Regulation'] = regulation
    # Plot
    sns.set_style("whitegrid")
    ax_scatter = sns.scatterplot(x='TotalMean', y='log2(fc)', data=df, hue='Regulation',
                                 palette=list(df.col.drop_duplicates()))
    if len(Proteins) > 0:
        df = df[df.gene_name.isin(Proteins)]
        ls_final = df
        ls_final = ls_final[['gene_name', 'log2(fc)', 'TotalMean']]
        from adjustText import adjust_text
        texts = []
        for a, b, c in zip(ls_final['TotalMean'], ls_final['log2(fc)'],
                           ls_final['gene_name']):
            texts.append(ax_scatter.text(a, b, c, size=8))
        adjust_text(texts, ax=ax_scatter, force_points=0.25, force_text=0.25,
                    expand_points=(1.5, 1.5), expand_text=(1.5, 1.5),
                    arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    sns.despine()
    plt.axhline(y=0, color='gray', linestyle=':')
    plt.grid(b=False)
    plt.xlabel('log2(Mean)')
    plt.title('MA plot')
    if save is not None:
        if vector is True:
            plt.savefig(save + 'MAPlot.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'MAPlot.png', dpi=dpi, bbox_inches='tight')
    plt.show(block=True)


def find_k(df):
    """Find the optimum k clusters"""

    df_norm = df
    sse = {}

    for k in range(1, 21):
        kmeans = KMeans(n_clusters=k, random_state=1)
        kmeans.fit(df_norm)
        sse[k] = kmeans.inertia_

    kn = KneeLocator(x=list(sse.keys()),
                     y=list(sse.values()),
                     curve='convex',
                     direction='decreasing')
    k = kn.knee
    print('KneeLocator identifies: ' + str(k) + ' clusters')
    return k


def k_trend(self, pvalue=0.05, k_cluster=None,
            save=None, dpi=300, vector=True):
    """Perform a K-mean algorithm

    BigTrend apply k-mean algorithm to identify co-expressed
    proteins/genes. For longitudinal analysis, k-means can help users to
    visualize the trends of proteins in the evaluated timecourse.

    Optionally, the user can define the number of clusters that k-means will cluster proteins and samples.
    By default, OmicScope run KneeLocator algorithm to suggest an optimal number of clusters.

    Args:
        OmicScope (_type_): OmicScope object.
        pvalue (float, optional): _description_. Defaults to 0.05.
        k_cluster (int, optional): Number of cluster to perform k-means.
         If None, OmicScope defines k-clusters based on KneeLocator algorithm.
         Defaults to None.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.

    Returns:
        _type_: _description_
    """
    plt.rcParams['figure.dpi'] = dpi
    omics = copy.copy(self)
    pdata = omics.pdata
    try:
        pdata['sample'] = pdata[['Condition', 'TimeCourse']].astype(str).apply(lambda x: '-'.join(x), axis=1)
    except KeyError:
        pdata['sample'] = pdata[['Condition', 'Biological']].astype(str).apply(lambda x: '-'.join(x), axis=1)
    data = omics.quant_data
    protein_dictionary = dict(zip(data['Accession'], data['gene_name']))
    data = data[data[omics.pvalue] <= pvalue]
    data = data.set_index('Accession')
    data = data.iloc[:, data.columns.str.contains('.', regex=False)]
    data = data.replace(0, 0.01)
    data = np.log2(data)
    # Scale protein abundance according to their mean (z-score)
    zscored_data = zscore(data, axis=1)
    if k_cluster is None:
        n_clusters = find_k(zscored_data)
    else:
        n_clusters = k_cluster
    kmeans = KMeans(n_clusters=n_clusters, init="k-means++",
                    max_iter=500)
    protein_cluster_assig = kmeans.fit(zscored_data)
    k_data = zscored_data
    k_data.columns = pdata['sample']
    k_data['cluster'] = protein_cluster_assig.labels_
    k_data_protein = k_data.reset_index().melt(id_vars=['Accession', 'cluster'],
                                               var_name='sample')
    dictionary = dict(zip(pdata['sample'], pdata['Condition']))
    k_data_protein['Condition'] = k_data_protein[['sample']].replace(dictionary)
    k_data_protein['gene_name'] = k_data_protein[['Accession']].replace(protein_dictionary)
    g = sns.catplot(
        data=k_data_protein, x="sample", y="value", hue="Condition", col="cluster",
        capsize=.2, palette="viridis",
        kind="point", col_wrap=2
    )
    g.despine()
    (g.map(plt.axhline, y=0, color="black", dashes=(2, 1), zorder=0)
     .set_ylabels('z-score', clear_inner=True)
     .set_xlabels('Conditions', clear_inner=True)
     .set(xticklabels=[]))
    if save is not None:
        if vector is True:
            plt.savefig(save + 'MAPlot.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'MAPlot.png', dpi=dpi, bbox_inches='tight')
    plt.show()
    k_data_protein = k_data_protein.groupby('cluster')['gene_name'].apply(set).reset_index()
    return k_data_protein


def PPInteractions(self, *Proteins,
                   network_k=None,
                   network_iterations=50,
                   score_threshold=0.6, labels=False, modules=False,
                   module_palette='Paired', species=9606,
                   pvalue=0.05, network_type='functional', save=None, dpi=300, vector=True):
    """Protein-Protein interaction

    Using String API, OmicScope search or protein-protein interactions among
    the proteins from users's dataset.

    Args:
        network_k (int, optional): Optimal distance between nodes. If None the distance is set to
          1/sqrt(n) where n is the number of nodes. Defauts to None.
        network_iterations (int, optional): Maximum number of iterations taken. Defaults to 50.
        score_threshold (float, optional): String-score threshold. Defaults to 0.6.
        labels (bool, optional): Show protein names. Defaults to False.
        modules (bool, optional): Perform Louvain Method to find communities/modules.
         Defaults to False.
        module_palette (str, optional): Color palette to assign modules. Defaults to 'Paired'.
        species (int, optional): String requires the definition of organism to search
         for proteins name, according to NCBI identifier. Defaults to 9606 (Human).
         Mus musculus = 10090; Rattus norvegicus = 10116.
        pvalue (float, optional): p-value threshold to proteins to be accepted.
         Defaults to 0.05 (differentially regulated).
        network_type (str, optional): Interactions can be defined as 'functional' or
         'physical'. Defaults to 'functional'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure in as vector (.svg). Defaults to
         True.

    Returns:
        Networkx object
    """
    # Set parameters to import data
    plt.rcParams['figure.dpi'] = dpi
    string_api_url = "https://version12.string-db.org/api"
    output_format = "json"
    method = "network"
    data = copy.copy(self)
    data = data.quant_data
    data = data[data[self.pvalue] <= pvalue]
    data = data.sort_values(self.pvalue)
    data = data.reset_index(drop=True)
    if len(data) > 2000:
        data = data.sort_values(self.pvalue)
        data = data.iloc[0:1999, :]
    data['TopTerms'] = data.index+1
    foldchange_range = data['log2(fc)']
    genes = data['gene_name'].dropna().drop_duplicates()
    genes = list(genes.astype(str))

    request_url = "/".join([string_api_url, output_format, method])

    params = {

        "identifiers": "%0d".join(genes),
        "species": species,  # species NCBI identifier
        "caller_identity": "Omicscope",
        'show_query_node_labels': 1,
        'network_type': network_type

    }
    response = requests.post(request_url, data=params)
    string = pd.read_json(response.text)
    # Filter string results based on score
    string = string[string['score'] >= score_threshold]

    norm = mpl.colors.TwoSlopeNorm(vmin=min(foldchange_range),
                                   vmax=max(foldchange_range),
                                   vcenter=0)
    cmap = cm.RdBu_r
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_hex = [mcolors.to_hex(m.to_rgba(x)) for x in data['log2(fc)']]
    source = pd.DataFrame({'ID': data['gene_name'],
                           'Size': 10,
                           'type': 'protein',
                           'color': color_hex,
                           'regulation': data['log2(fc)']})
    source = source.drop_duplicates()
    source = source.dropna()
    string = string.rename(columns={'score': 'weight'})
    G = nx.from_pandas_edgelist(string,
                                source='preferredName_A',
                                target='preferredName_B',
                                edge_attr='weight',
                                create_using=nx.Graph)

    carac = source.set_index('ID')
    carac['color_edge'] = 'black'
    linewidths = 0.8
    nx.set_node_attributes(G, dict(zip(carac.index, carac.color)), name="color")
    nx.set_node_attributes(G, dict(zip(carac.index, carac.Size)), name="Size")
    # widths = nx.get_edge_attributes(G, 'score')
    G.edges(data=True)
    # Find Communities
    if modules is True:
        import networkx.algorithms.community as nx_comm

        # Find communities based on label propagation
        communities = nx_comm.louvain_communities(G)
        module = []
        terms = []
        color = []
        degree = []
        # Linking Terms, module, colors
        cmap = sns.color_palette(module_palette, len(communities), desat=.9)
        for i, g in enumerate(communities):
            subgraph = G.subgraph(g)
            degree_subgraph = subgraph.degree
            degree.extend(degree_subgraph)
            module.extend([i]*len(g))
            terms.extend(list(g))
            if len(g) > len(source)*0.05:
                color.extend([mcolors.to_hex(cmap[i])]*len(g))
            else:
                color.extend(['#666666']*len(g))
            # DataFrame
        modularity = pd.DataFrame(zip(module, terms, color, degree), columns=[
                                  'Module', 'ID', 'color_edge', 'Degree'])
        # Merging carac with modularity data
        carac = carac.merge(modularity, on='ID', suffixes=('_x', ''))
        carac = carac.set_index('ID')
        nx.set_node_attributes(G, dict(zip(carac.index, carac.Module)), name="Module")
        carac = carac.reindex(G.nodes())
        linewidths = 1.5
    pos = nx.spring_layout(G,
                           k=network_k,
                           iterations=network_iterations,
                           threshold=0.0001, weight='weight',
                           center=None,)
    carac = carac.reindex(G.nodes())
    weights = nx.get_edge_attributes(G, 'weight')
    nx.draw(G,
            pos=pos,
            node_color=carac['color'],
            node_size=carac['Size']*10,
            edgecolors=carac['color_edge'],
            linewidths=linewidths,
            alpha=0.9,
            width=[i**4 for i in list(weights.values())],
            edge_color='#666666')
    if labels is True:
        nx.draw_networkx_labels(G, pos, font_size=6)
    if save is not None:
        nx.write_graphml(G, save + 'PPinteraction.graphml', named_key_ids=True)
        if vector is True:
            plt.savefig(save + 'PPinteractions.svg', bbox_inches='tight')
        else:
            plt.savefig(save + 'PPinteractions.dpi', dpi=dpi, bbox_inches='tight')
        plt.show(block=True)
    return G
