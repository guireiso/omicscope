import omicscope as omics


class Test(object):
    def test_input_nebula(self):
        df = omics.Nebula('tests//data//MultipleGroups//omics_file')
        assert len(df.groups) == 4
