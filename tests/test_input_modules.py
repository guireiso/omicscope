from omicscope.General.Input import Input


class TestInput(object):
    def test_Progenesis(self):
        ctrlWT = Input('tests/data/proteins/progenesis.xls', Method='Progenesis')
        ctrlNone = Input('tests/data/proteins/progenesis.xls', Method='Progenesis')
        assert len(ctrlWT.Conditions) == len(ctrlNone.Conditions)
        assert 'Sample' in ctrlWT.pdata.columns
        assert 'Condition' in ctrlWT.pdata.columns

    def test_General(self):
        ctrlWT = Input('tests//data//proteins//general.xlsx', Method='General')
        assert 'Sample' in ctrlWT.pdata.columns
        assert 'Condition' in ctrlWT.pdata.columns

    def test_MaxQuant(self):
        pdata = 'tests//data//proteins//MQ_pdata.xlsx'
        maxquantNone = Input('tests//data//proteins//MQ.txt',
                             Method='MaxQuant', pdata=pdata)
        assert 'Sample' in maxquantNone.pdata.columns
        assert 'Condition' in maxquantNone.pdata.columns

    def test_PatternLab(self):
        patternLabNone = Input('tests//data//proteins//patternlab.plp',
                               Method='PatternLab')
        assert 'Sample' in patternLabNone.pdata.columns
        assert 'Condition' in patternLabNone.pdata.columns
