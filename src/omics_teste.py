
from copy import copy

import omicscope as omics
# path = '..\\tests\\data\\MultipleGroups\\original\\'

# for i in ['astrocytes.csv', 'HB_VCC_2021.xls', 'neurons.csv', 'sh.csv']:
#      df = omics.Omicscope(path + i, 'CTRL', 'Progenesis', 'static', 'pAdjusted')
#      enrichment = omics.EnrichmentScope(df, 'ORA')
#      enrichment.savefile(Path = '..\\tests\\data\\MultipleGroups\\omics_file\\')
#      print('done')



path = '..\\tests\\data\\MultipleGroups\\omics_file\\'

df = omics.multiples(path)

teste = df.circos_plot()

# class circos():
#     def __init__(self, multiples_output):
#         self.group_data = multiples_output.group_data
#         self.enrichment = multiples_output.enrichment
#         self.groups = multiples_output.groups
#         self.labels = multiples_output.labels
#         self.original = multiples_output.original
#         self.original_path = multiples_output.original_path
#         self.foldering()
#         self.karyotype()
#         self.heatmap()
#         self.linkProteins()
#         if any(enr is not None for enr in self.enrichment):
#             self.test = self.linkEnrichment()
#         # self.changing_folders()
#         # self.perl()
#         # self.copycircos()
#         # self.plot()

#     def foldering(self):
#         string = os.path.dirname(os.path.abspath(__file__))
#         string = os.path.normpath(string)
#         string = string.split('\\')
#         string = '\\'.join(string[:-1])
#         try:
#             os.makedirs(string + '\\circos\\newData\\data')
#         except FileExistsError:
#             shutil.rmtree(string + '\\circos\\newData')
#             os.makedirs(string + '\\circos\\newData\\data')
#         self.newFolder = string + '\\circos\\newData\\data\\'

#     def karyotype(self):
#         labels = copy(self.labels)
#         group_data = copy(self.group_data)
#         colors = ['set3-12-qual-'+str(x) for x in range(1,13)]
#         colors = colors[:len(labels)]
#         karyotype = pd.DataFrame({'type': repeat('chr', len(labels)),
#                                 'chr': repeat('-', len(labels)),
#                                 'name': labels,
#                                 'label': labels,
#                                 'START': repeat(0, len(labels)),
#                                 'END': [len(group_data[i]) for i in range(0, len(labels))],
#                                 'COLOR': colors})
#         karyotype.to_csv(self.newFolder + 'karyotype.txt', sep='\t', header=False,
#                           index=False)

#     def heatmap(self):
#         from copy import copy

#         import pandas as pd
#         group_data = copy(self.group_data)
#         labels = copy(self.labels)
#         whole = pd.concat(group_data)
#         whole['dup'] = whole['gene_name'].duplicated(keep=False)
#         whole = whole.sort_values(by=['group', 'dup'], ascending=False)
#         self.whole = whole

#         def transforming(whole, group):
#             df = whole[whole['group'] == group]
#             df['type'] = 'band'
#             df['chr'] = group
#             df['name'] = df['gene_name']
#             df['label'] = df['gene_name']
#             df['START'] = range(0, len(df))
#             df['END'] = range(1, len(df) + 1)
#             df['COLOR'] = df['log2(fc)']
#             df = df[['type', 'chr', 'name', 'label', 'START', 'END', 'COLOR']]
#             return(df)

#         transformed = []
#         for i in labels:
#             transformed.append(transforming(whole, i))
#         self.transformed = transformed
#         heatmap = pd.concat(transformed)[['chr', 'START', 'END', 'COLOR']]
#         heatmap.to_csv(self.newFolder + 'heatmap.txt', sep="\t",
#                             index=False)

#     def linkProteins(self):
#         from copy import copy

#         import pandas as pd
#         transformed = copy(self.transformed)
#         labels = copy(self.labels)

#         links = pd.concat(transformed)[['name', 'chr', 'START', 'END']]
#         links['dup'] = links['name'].duplicated(keep=False)
#         links = links[links['dup'] == True].iloc[:, :-1]
#         unique_names = pd.DataFrame(links['name'].unique())
#         unique_names.columns = ['name']

#         modified = []

#         for i, group in enumerate(labels):
#             for x in range(i + 1, len(labels)):
#                 modified.append(unique_names.merge(links[links['chr'] == labels[i]], how='left', on='name').merge(
#                     links[links['chr'] == labels[x]], how='left', on='name').dropna().iloc[:, 1:])

#         joined = pd.concat(modified)
#         joined.columns = ['chr.x', 'START.x', 'END.x', 'chr.y', 'START.y', 'END.y']
#         joined[['START.x', 'END.x', 'START.y', 'END.y']] = joined[[
#             'START.x', 'END.x', 'START.y', 'END.y']].astype(int)
#         joined.to_csv(self.newFolder + 'links.txt', sep="\t",
#                           index=False)

#     def linkEnrichment(self):
#         import pandas as pd
#         enrichment = self.enrichment
#         groups = self.groups
#         def enrichment_links(original, group):
#             df = original
#             df = df.reset_index()
#             databases = list(df.Gene_set.drop_duplicates())
#             pal = 'color=paired'
#             n = len(databases)
#             t = 'qual'
#             strings = []
#             for m in range(1, n + 1):
#                 string = '-'.join([pal, str(n), t, str(m)])
#                 strings.append(string)
#             di = dict(zip(databases, strings))
#             df = df[df['Adjusted P-value'] < 0.05][['Term', 'Gene_set']]
#             df['group'] = group
#             df = df.replace({"Gene_set": di}).reset_index(drop=True)
#             df.columns = ['pathway', 'color', 'group']
#             return(df)
#         whole_enr = []
#         for i, g in zip(enrichment, groups):
#             whole_enr.append(enrichment_links(i, g))
#         whole_enr = pd.concat(whole_enr)
#         whole_enr['dup'] = whole_enr['pathway'].duplicated(keep=False)
#         whole_enr = whole_enr.sort_values(by=['group', 'dup'], ascending=False)

#         def transforming_enr(whole_enr, whole_proteome, group):
#             import random

#             import numpy as np
#             random.seed(123)
#             chrlen = len(whole_proteome[whole_proteome['group'] == group])
#             df = whole_enr[whole_enr['group'] == group]
#             df['type'] = 'band'
#             df['chr'] = group
#             df['name'] = df['pathway']
#             df['label'] = df['pathway']
#             df['START'] = np.random.randint(0, chrlen, len(df))
#             df['END'] = df['START'] + 1
#             df['COLOR'] = df['color']
#             df = df[['type', 'chr', 'name', 'label', 'START', 'END', 'COLOR']]
#             return(df)
#         links = []
#         for g in groups:
#             links.append(transforming_enr(whole_enr, self.whole, g))

#         links_enr = pd.concat(links)[['name', 'chr', 'START', 'END', 'COLOR']]
#         links_enr['dup'] = links_enr['name'].duplicated(keep=False)
#         links_enr = links_enr[links_enr['dup'] == True].iloc[:, :-1]
#         unique_names_enr = pd.DataFrame(links_enr['name'].unique())
#         unique_names_enr.columns = ['name']

#         modified_enr = []

#         for i, group in enumerate(groups):
#             for x in range(i + 1, len(groups)):
#                 modified_enr.append(unique_names_enr.merge(links_enr[links_enr['chr'] == groups[i]], how='left', on='name').merge(
#                     links_enr[links_enr['chr'] == groups[x]], how='left', on='name').dropna().iloc[:, 1:])

#         joined_enr = pd.concat(modified_enr)
#         return joined_enr
#         joined_enr = joined_enr[['chr_x', 'START_x',
#             'END_x', 'chr_y', 'START_y', 'END_y', 'COLOR_y']]
#         joined_enr.columns = ['chr.x', 'START.x', 'END.x', 'chr.y', 'START.y', 'END.y', 'COLOR']
#         joined_enr[['START.x', 'END.x', 'START.y', 'END.y']] = joined_enr[[
#             'START.x', 'END.x', 'START.y', 'END.y']].astype(int)
#         joined_enr.to_csv(self.newFolder + "links_path.txt", sep="\t",
#                           index=False)

# import glob
# import os
# import random
# import shutil
# from copy import copy
# from itertools import repeat

# import pandas as pd

# teste = circos(df)
# teste.test