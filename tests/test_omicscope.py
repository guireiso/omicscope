from src.omicscope.cli import main


def test_main():
    assert main([]) == 0

class TestGeneral(object):
    def test_size_data(self):
        from src.omicscope.Input.General import Input
        df = Input('tests//data//proteins//general.xls', 'WT')
        assert len(df.assay.columns) == len(df.pdata)
        assert len(df.assay.index) == len(df.rdata)

    def test_assay(self):
        from src.omicscope.Input.General import Input
        df = Input('tests//data//proteins//general.xls', 'WT')
        assay = df.assay
        assert assay.size == 52192 - len(assay.columns)

    def test_pdata(self):
        from src.omicscope.Input.General import Input
        df = Input('tests//data//proteins//general.xls', 'WT')
        pdata = df.pdata
        Conditions = list(pdata.Condition.drop_duplicates())
        assert Conditions == df.Conditions
    
class TestProgenesis(object):
    def test_size_data(self):
        from src.omicscope.Input.Progenesis import Input
        df = Input('tests//data//proteins//progenesis.csv', 'WT')
        assert len(df.assay.columns) == len(df.pdata)
        assert len(df.assay.index) == len(df.rdata)

    def test_assay(self):
        from src.omicscope.Input.Progenesis import Input
        df = Input('tests//data//proteins//progenesis.csv', 'WT')
        assay = df.assay
        assert assay.size == 52192 - len(assay.columns)

    def test_pdata(self):
        from src.omicscope.Input.Progenesis import Input
        df = Input('tests//data//proteins//progenesis.csv', 'WT')
        pdata = df.pdata
        Conditions = list(pdata.Condition.drop_duplicates())
        assert Conditions == df.Conditions