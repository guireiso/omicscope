
from copy import copy

import omicscope as omics

longitudinal = omics.Omicscope(Table = 'C:/Users/Guilherme/Desktop/ALL.csv',
                                 Method = 'Progenesis',
                                 ControlGroup= None,
                                 pdata = 'C:/Users/Guilherme/Desktop/pdata.xls',
                                 ExperimentalDesign= 'longitudinal', 
                                 pvalue = 'pvalue')

static = omics.Omicscope(Table = 'C:/Users/Guilherme/Desktop/progenesis.csv',
                                Method = 'Progenesis',
                                ControlGroup= None,
                                ExperimentalDesign= 'static', 
                                pvalue = 'pvalue')
