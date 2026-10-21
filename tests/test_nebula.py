import omicscope as omics


class Test(object):
    def test_input_nebula(self):
        df = omics.Nebula('tests//data//MultipleGroups//omics_file')
        assert len(df.groups) == 4

    def test_nebula_figures(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        df = omics.Nebula('tests//data//MultipleGroups//omics_file')
        calls = [
            lambda: df.barplot(),
            lambda: df.diff_reg(),
            lambda: df.similarity_heatmap(),
            lambda: df.similarity_network(),
            lambda: df.stat_heatmap(),
            lambda: df.circos_plot(),
        ]
        for call in calls:
            call()
            plt.close('all')
