class Test_Enrichment(object):
    def test_omicscope_ORA(self):
        from src import omicscope as omics
        progenesis = omics.Omicscope(Table = 'tests//data//proteins//progenesis.csv',
                        Method = 'Progenesis',
                        ControlGroup= 'WT',
                        pvalue = 'pvalue')
        enr = omics.EnrichmentScope(OmicScope = progenesis, Analysis = 'ORA')
        assert enr.Analysis == 'ORA'