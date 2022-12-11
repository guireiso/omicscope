
from copy import copy

import omicscope as omics

teste = omics.Omicscope(Table = 'C:/Users/Guilherme/Desktop/ALL.csv',
                          Method = 'Progenesis',
                          ControlGroup= None,
                          pdata = 'C:/Users/Guilherme/Desktop/pdata.xls',
                          ExperimentalDesign= 'longitudinal', 
                          pvalue = 'pvalue')

# omics.heatmap(teste, c_cluster = False)
# import matplotlib.pyplot as plt
# import pandas as pd
# import numpy as np
# import seaborn as sns
# import copy
# import itertools
# from sklearn import preprocessing
# from sklearn.decomposition import PCA

# def bar_protein(OmicScope, *Proteins, logscale=True,
#                 palette='Spectral', save='', dpi=300,
#                 vector=True):
#     """Bar plot to show protein abundance in each condition

#     Args:
#         OmicScope (OmicScope object): OmicScope Experiment
#         logscale (bool, optional): Apply abundance log-transformed.
#         Defaults to True.
#         palette (str, optional): Palette for groups. Defaults to 'Spectral'.
#         save (str, optional): Path to save figure. Defaults to ''.
#         dpi (int, optional): figure resolution. Defaults to 300.
#         vector (bool, optional): Save figure in as vector (.svg). Defaults to
#         True.
#     """
#     plt.rcParams['figure.dpi']=dpi
#     df = copy.copy(OmicScope.quant_data)
#     # Proteins to plot
#     df = df[df['gene_name'].isin(Proteins)]
#     # Get conditions
#     ctrl = [copy.copy(OmicScope.ctrl)]
#     conditions = copy.copy(OmicScope.Conditions)
#     conditions.remove(ctrl[0])
#     # Get protein abundance for each condition 
#     df = df.set_index('gene_name')
#     df = df.iloc[:,df.columns.str.contains('.', regex = False)]
#     df = df.melt(ignore_index=False)
#     df = df.reset_index()
#     pdata = copy.copy(OmicScope.pdata)
#     pdata = pdata.set_index('Sample')
#     df = df.set_index('variable')
#     df = pdata.merge(df, left_index = True, right_index = True)
#     if 'TimeCourse' in df.columns:
#         df = df[[ 'Condition',  'TimeCourse','gene_name', 'value']]
#         df = df.set_index(['Condition', 'TimeCourse'])
#     else:
#         df = df[['Condition', 'gene_name', 'value']]
#         df = df.set_index('Condition')
#     # Apply log transformation
#     if logscale == True:
#         df['value'] = np.log2(df['value'])
#         df['value'] = df['value'].replace(-np.inf, np.nan)
#     # Size of figure
#     sns.set(rc={'figure.figsize': (1.5*len(Proteins),5)})
#     custom_params = {"axes.spines.right": False, "axes.spines.top": False}
#     sns.set_theme(style="ticks", rc=custom_params)
#     # Plot
#     if len(Proteins) == 1:
#         sns.barplot(x='Condition', y='value',
#                    data=df, errwidth=1, capsize=0.07,
#                    palette=palette, order=ctrl + conditions,
#                    edgecolor='black', linewidth=1, dodge=False)
#     else:
#         sns.barplot(x='gene_name', y='value', hue=list(df.index),
#                     data=df, errwidth=1, capsize=0.07,
#                     palette=palette,
#                     edgecolor='black', linewidth=1)
#         plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.)
#     sns.despine()
#     plt.title('Abundance - ' + ' and '.join(Proteins))
#     plt.xlabel('')
#     plt.ylabel('log2(Abundance)')
#     if save != '':
#         if vector == True:
#             plt.savefig(save + 'barplot_' + '_'.join(Proteins) + '.svg')
#         else:
#             plt.savefig(save + 'barplot_' + '_'.join(Proteins) + '.png', dpi=dpi)
#     plt.show()

# bar_protein(teste, 'GAPDH', 'ACTN4')