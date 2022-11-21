from src.omicscope.cli import main


def test_main():
    assert main([]) == 0

class TestInput(object):
    def test_Progenesis(self):
        from src.omicscope.Input import Input
        ctrlWT = Input('tests//data//proteins//progenesis.csv', Method = 'Progenesis', ControlGroup='WT')
        ctrlNone = Input('tests//data//proteins//progenesis.csv', Method = 'Progenesis', ControlGroup=None)
        assert len(ctrlWT.Conditions) == len(ctrlNone.Conditions)
        assert ctrlWT.ControlGroup == 'WT'
        assert ctrlNone.ControlGroup == 'KO'
    
    def test_General(self):
        from src.omicscope.Input import Input
        ctrlWT = Input('tests//data//proteins//general.xls', Method = 'General', ControlGroup='WT')
        ctrlNone = Input('tests//data//proteins//general.xls', Method = 'General', ControlGroup=None)
        assert len(ctrlWT.Conditions) == len(ctrlNone.Conditions)
        assert ctrlWT.ControlGroup == 'WT'
        assert ctrlNone.ControlGroup == 'KO'
    
    def test_MaxQuant(self):
        from src.omicscope.Input import Input
        pdata = 'tests//data//proteins//MQ_pdata.xlsx'
        maxquantNone = Input('tests//data//proteins//MQ.txt', ControlGroup = None,
         Method = 'MaxQuant', pdata = pdata, filtering_method = 70)
        maxquantWT = Input('tests//data//proteins//MQ.txt', ControlGroup = 'CoV',
         Method = 'MaxQuant', pdata = pdata, filtering_method = 70)
        assert len(maxquantWT.Conditions) == len(maxquantNone.Conditions)
        assert maxquantNone.ControlGroup == maxquantWT.ControlGroup

    def test_PatternLab(self):
        from src.omicscope.Input import Input
        patternLabNone = Input('tests//data//proteins//patternlab.xlsx', ControlGroup = None,
         Method = 'PatternLab', filtering_method = 70)
        patternLabCN = Input('tests//data//proteins//patternlab.xlsx', ControlGroup = 'CN',
         Method = 'PatternLab', filtering_method = 70)
        assert patternLabNone.ControlGroup == patternLabCN.ControlGroup