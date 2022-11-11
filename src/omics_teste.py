import omicscope as omics


teste = omics.multiples('C:\\Users\\Guilherme\\omicscope\\omicscope\\tests\\data\\MultipleGroups\\omics_file')

from omicscope.MultipleData.MultipleVisualisation import *


general_data(teste)

upsetplot_proteins(teste)

group_pearson(teste) # TODO

Differentially_Regulated(teste) # TODO