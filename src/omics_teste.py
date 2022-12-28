
from copy import copy

import omicscope as omics
# path = '..\\tests\\data\\MultipleGroups\\original\\'

# for i in ['astrocytes.csv', 'HB_VCC_2021.xls', 'neurons.csv', 'sh.csv']:
#     df = omics.Omicscope(path + i, 'CTRL', 'Progenesis', 'static', 'pAdjusted')
#     enrichment = omics.EnrichmentScope(df, 'ORA')
#     enrichment.savefile(Path = '..\\tests\\data\\MultipleGroups\\omics_file\\')
#     print('done')

# print('terminou caraio')

path = '..\\tests\\data\\MultipleGroups\\omics_file\\'

df = omics.multiples(path)

#sum(dhyper(t:b, a, n - a, b))
# sum(dhyper(10:30, 70, 130, 30))
#0.6561562

def overlap_fisher(group1, group2, union):
    from scipy.stats import hypergeom
    import numpy as np
    deps1 = set(group1)
    deps2 = set(group2)
    union = len(union)
    intersection = len(deps1.intersection(deps2))
    distribution = min([len(deps1), len(deps2)])
    [M, n, N] = [union, len(deps1), len(deps2)]
    rv = hypergeom(M, n, N)
    x = np.arange(intersection, distribution)
    pval = sum(rv.pmf(x))
    return pval

def overlap_stat(self):
    import pandas as pd
    from scipy.spatial.distance import pdist, squareform
    pval = self.pvalue
    union = set(pd.concat(self.original).gene_name)
    sets = [set(x[x[pval] < 0.05].gene_name) for x in self.original]
    sets = pd.DataFrame(sets)
    pvalue = pdist(sets, lambda u, v: overlap_fisher(u, v, union = union))
    pvalue = squareform(pvalue)
    pvalue = pd.DataFrame(pvalue, columns = self.groups, index = self.groups)
    return pvalue

# pdist(X, lambda u, v: np.sqrt(((u-v)**2).sum()))
teste = overlap_groups(df)
import seaborn as sns
sns.clustermap(teste, )