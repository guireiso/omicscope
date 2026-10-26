
from copy import copy
import numpy as np
import pandas as pd

from .utils import read_table, require_pdata


class Input:

    def __init__(self, Table, pdata=None):
        """ RPKM for transcriptomics analysis output
        RPKM method is composed by a tab-separated file and pdata table. 
        Following Transcriptomics guidelines, we suggest the use of RPKM/FPKM/TPM 
        values to perform differential analysis. To avoid NA-values the guideline 
        suggest to perform RPKM log-transformation using log2(value+1).

        The tab-separated file must contain:
         1. 'gene_name' column: Gene Symbols
         2. 'Accession' column: specific gene cone (eg. 'NM_130786.4')
         3. Other columns: RPKM/TPM/FPKM values.

        Args:
            Table (str): Path to the Matrix.txt file (tab-separated file).
            pdata (str): Path to the phenotype data, which is a .xlsx or .xls
            file that contains all the features associated to the samples

        """
        self.Table = Table
        data = read_table(copy(self.Table))
        assay = data.drop(['gene_name', 'Accession'], axis='columns')
        assay = np.log2(assay+1)
        rdata = data[['gene_name', 'Accession']]
        self.assay = assay
        self.rdata = rdata
        self.pdata = require_pdata(pdata, 'The RPKM table', 'the column names of the RPKM table')
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())
