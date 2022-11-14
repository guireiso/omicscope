from src import omicscope as omics


class TestEnrichment(object):
    def test_omicscope_progenesis(self):
        df = omics.Omicscope('tests//data//proteins//progenesis.csv',
            Method = 'Progenesis', ControlGroup= 'WT')
        self.df = df
        assert len(df.quant_data) == 1594

    
    def test_omicscope_ORA(self):
        df = omics.Omicscope('tests//data//proteins//progenesis.csv',
            Method = 'Progenesis', ControlGroup= 'WT')
        enr = omics.EnrichmentScope(OmicScope = df, Analysis = 'ORA')
        