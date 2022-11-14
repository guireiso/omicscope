import omicscope as omics

teste = omics.Omicscope('C:\\Users\\Guilherme\\omicscope\\omicscope\\tests\\data\\proteins\\progenesis.csv', 
        'WT', 'Progenesis')  

teste2 = omics.EnrichmentScope(teste, 'GSEA')


teste2.output('C:/Users/Guilherme/Desktop/')