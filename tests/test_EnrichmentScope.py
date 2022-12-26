from src import omicscope as omics

def test_omicscope_ORA(self):
    df = omics.Omicscope('tests//data//proteins//progenesis.csv',
    Method = 'Progenesis', ControlGroup= 'WT')
    enr = omics.EnrichmentScope(OmicScope = df, Analysis = 'ORA')
        