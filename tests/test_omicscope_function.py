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
