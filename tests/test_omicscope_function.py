from src.omicscope import Omicscope



class TestOmicScope(object):
    def test_user_input_stat(self):
        progenesis = Omicscope(Table = 'tests//data//proteins//progenesis.csv',
                        Method = 'Progenesis',
                        ControlGroup= 'WT',
                        pvalue = 'pvalue')
        general = Omicscope(Table = 'tests//data//proteins//general.xls',
                        Method = 'General',
                        ControlGroup= None,
                        pvalue = 'pvalue')
        assert progenesis.ctrl == general.experimental[0]
        assert general.ctrl == progenesis.experimental[0]
        assert len(progenesis.Conditions) == 2
    
    # def test_static_workflow():
    #     df = Omicscope('tests//data//proteins//progenesis.csv',
    #     Method = 'Progenesis', ControlGroup= 'KO')
    #     assert df.ctrl == 'KO'
