import numpy as np

import omicscope as omics

np.seterr(divide='ignore', invalid='ignore')


class TestEnrichmentScope(object):
    def test_user(self):
        data = omics.OmicScope(Table='tests/data/proteins/progenesis.xls',
                                     Method='Progenesis',
                                     ControlGroup='CTRL',
                                     pvalue='pvalue')
        enrichment = omics.EnrichmentScope(data)
        assert len(enrichment.results) > 0
