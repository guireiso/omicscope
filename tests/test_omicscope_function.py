import numpy as np

import omicscope as omics

np.seterr(divide='ignore', invalid='ignore')


class TestOmicScope(object):
    def test_user_input_stat(self):
        progenesis = omics.OmicScope(Table='tests/data/proteins/progenesis.xls',
                                     Method='Progenesis',
                                     ControlGroup='CTRL',
                                     pvalue='pvalue')
        general = omics.OmicScope(Table='tests//data//proteins//general.xlsx',
                                  Method='General',
                                  ControlGroup=None,
                                  pvalue='pvalue')
        assert progenesis.ctrl == 'CTRL'
        assert general.ctrl == 'COVID'
        assert len(progenesis.Conditions) == 2

    def test_snapshot(self):
        snapshot = omics.OmicScope(Table='tests//data//proteins//snapshot.xlsx',
                                   Method='Snapshot')
        assert snapshot.ControlGroup == ' CTRL'

    def test_longitudinal(self, tmp_path):
        import pandas as pd
        rng = np.random.default_rng(0)
        pdata = pd.DataFrame({'Sample': [f'S{i}' for i in range(24)],
                              'Condition': ['Control'] * 12 + ['Treated'] * 12,
                              'TimeCourse': ([1] * 4 + [2] * 4 + [3] * 4) * 2,
                              'Biological': list(range(1, 5)) * 6})
        assay = pd.DataFrame(rng.lognormal(10, 1, (60, 24)), columns=pdata.Sample)
        rdata = pd.DataFrame({'Accession': [f'P{i}' for i in range(60)],
                              'gene_name': [f'G{i}' for i in range(60)],
                              'Description': 'protein'})
        path = tmp_path / 'longitudinal.xlsx'
        with pd.ExcelWriter(path) as writer:
            assay.to_excel(writer, sheet_name='assay', index=False)
            rdata.to_excel(writer, sheet_name='rdata', index=False)
            pdata.to_excel(writer, sheet_name='pdata', index=False)
        result = omics.OmicScope(Table=str(path), Method='General', ExperimentalDesign='longitudinal',
                                 pvalue='pvalue')
        assert len(result.quant_data) > 0
