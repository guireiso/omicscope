from omicscope import Nebula


class Test(object):
    def test_input_nebula(self):
        df = Nebula('tests//data//MultipleGroups//omics_file')
        assert len(df.groups) == 4
