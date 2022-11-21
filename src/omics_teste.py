#import omicscope as omics

#teste = omics.Omicscope('C:\\Users\\Guilherme\\Desktop\\ALL.csv', ControlGroup= 'VEH1', Method = 'Progenesis')



from omicscope.Input import Input
patternLabNone = Input('..\\tests\\data\\proteins\\patternlab.xlsx', ControlGroup = None,
            Method = 'PatternLab', filtering_method = 70)