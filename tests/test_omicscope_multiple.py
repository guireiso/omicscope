from src.omicscope import *

class Test(object):
    def test_input_multiple(self):
        df = multiples('tests//data//MultipleGroups//omics_file')
        assert len(df.groups) == 4
        