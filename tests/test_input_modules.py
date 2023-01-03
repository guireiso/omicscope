from omicscope.cli import main


def test_main():
    assert main([]) == 0

class TestInput(object):
    def test_Progenesis(self):
        from omicscope.Input import Input
        ctrlWT = Input('tests//data//proteins//progenesis.csv', Method = 'Progenesis')
        ctrlNone = Input('tests//data//proteins//progenesis.csv', Method = 'Progenesis')
        assert len(ctrlWT.Conditions) == len(ctrlNone.Conditions)
        assert 'Sample' in ctrlWT.pdata.columns
        assert 'Condition' in ctrlWT.pdata.columns
    
    def test_General(self):
        from omicscope.Input import Input
        ctrlWT = Input('tests//data//proteins//general.xls', Method = 'General')
        assert 'Sample' in ctrlWT.pdata.columns
        assert 'Condition' in ctrlWT.pdata.columns
    
    def test_MaxQuant(self):
        from omicscope.Input import Input
        pdata = 'tests//data//proteins//MQ_pdata.xlsx'
        maxquantNone = Input('tests//data//proteins//MQ.txt',
         Method = 'MaxQuant', pdata = pdata, filtering_method = 70)
        maxquantWT = Input('tests//data//proteins//MQ.txt',
         Method = 'MaxQuant', pdata = pdata, filtering_method = 70)
        assert len(maxquantWT.Conditions) == len(maxquantNone.Conditions)
        assert 'Sample' in maxquantNone.pdata.columns
        assert 'Condition' in maxquantNone.pdata.columns

        

    def test_PatternLab(self):
        from omicscope.Input import Input
        patternLabNone = Input('tests//data//proteins//patternlab.xlsx',
         Method = 'PatternLab', filtering_method = 70)
        assert 'Sample' in patternLabNone.pdata.columns
        assert 'Condition' in patternLabNone.pdata.columns