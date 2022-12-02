
from copy import copy

import omicscope as omics

teste = omics.Omicscope(Table = 'C:/Users/Guilherme/Desktop/ALL.csv',
                        Method = 'Progenesis',
                        ControlGroup= 'VEH',
                        pdata = 'C:/Users/Guilherme/Desktop/pdata.xls',
                        ExperimentalDesign= 'longitudinal', 
                        pvalue = 'pvalue')

