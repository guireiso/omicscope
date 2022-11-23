
from copy import copy

import omicscope as omics

teste = omics.Omicscope(Table = 'C:/Users/Guilherme/Desktop/ALL.csv',
                        Method = 'Progenesis',
                        ControlGroup= None,
                        pdata = 'C:/Users/Guilherme/Desktop/pdata.xls',
                        ExperimentalDesign= 'longitudinal', 
                        pvalue = 'pvalue')

teste.quant_data