
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