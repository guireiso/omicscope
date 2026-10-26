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


def test_sankey_follows_the_enrichment_tables(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    df = omics.Nebula('tests//data//MultipleGroups//omics_file')
    enr = next(e for e in df.enrichment if e is not None)
    terms = list(enr['Term'].head(2))
    table = df.sankey_table(*terms, top=25)
    assert set(table['term']) <= set(terms) and table['gene_name'].nunique() <= 25
    for study, term, gene in table[['study', 'term', 'gene_name']].itertuples(index=False):
        e = df.enrichment[df.groups.index(study)]
        assert gene in e.loc[e['Term'] == term, 'Genes'].iloc[0]       # each flow is backed by the study's results
    df.sankey_plot(*terms, top=25, save=str(tmp_path) + '/')
    plt.close('all')
    assert (tmp_path / 'sankey_plot.svg').exists()
