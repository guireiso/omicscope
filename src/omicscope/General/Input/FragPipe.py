
from copy import copy
import sys

import pandas as pd


class Input:

    def __init__(self, Table, pdata):
        """ FragPipe output
        FragPipe output is composed by a combined_protein.tsv file and it contains
        several parameters that can be used to filter the data.

        Args:
            Table (str): Path to the combined_protein.tsv file
            pdata (str): Path to the phenotype data, which is a .xlsx or .xls
            file that contains all the features associated to the samples
        """
        self.Table = Table
        self.rdata_assay()
        try:
            pdata = pd.read_excel(pdata)
        except ValueError:
            pdata = pd.read_csv(pdata)
            if pdata.shape[1] == 1:
                raise KeyError('OmicScope does not import pdata correctly. Please, confirm if data is xlsx, xls or csv.')
        self.pdata = pdata
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())

    def rdata_assay(self):
        import pandas as pd

        data = pd.read_table(self.Table)
        if data.shape[1] == 1:
            data = pd.read_csv(self.Table, sep=';')
            if data.shape[1] == 1:
                raise KeyError('OmicScope does not import data correctly. Please, confirm if data is tab-separated values (tsv).')
        data = data.rename(columns={'Protein ID':'Accession',
                      'Gene': 'gene_name'},)
        data = data[~data['Protein'].str.contains('contam|rever')]
        data = data.set_index(['Accession'])
        assay = data.iloc[:,data.columns.str.contains('MaxLFQ')]
        assay.columns = assay.columns.str.split(' MaxLFQ').str[0]
        if assay.size == 0:
            assay = data.iloc[:,data.columns.str.contains('Intensity')]
            assay.columns = assay.columns.str.split(' Intensity').str[0]

        rdata = data.iloc[:,~data.columns.str.contains('Intensity')]
        self.rdata = copy(rdata.reset_index())
        self.assay = copy(assay)
