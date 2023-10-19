
import numpy as np
import pandas as pd


class Input:

    def __init__(self, Table, pdata):
        """ MaxQuant output
        MaxQuant output is composed by a ProteinGroup.txt file and it contains
        several parameters that can be used to filter the data.

        Args:
            Table (str): Path to the ProteinGroup.txt file
            pdata (str): Path to the phenotype data, which is a .xlsx or .xls
            file that contains all the features associated to the samples
            quant_strategy (str, optional): Quantification strategy that will be performed
            to the proteomic data. Defaults to 'LFQ intensity'. Also allowed
            'Intensity', and 'iBAQ' strategies.
            Defaults to 70.
        """
        self.Table = Table
        self.assay, self.rdata = self.assay_and_rdata()
        try:
            self.pdata = pd.read_excel(pdata)
        except ValueError:
            self.pdata = pd.read_csv(pdata)
        self.Conditions = []

    def assay_and_rdata(self):
        table = pd.read_csv(self.Table, sep='\t')
        table = table.rename(columns={
            'Protein.Group': 'Accession',
            'Protein.Ids': 'ProteinIds',
            'Protein.Names': 'ProteinNames',
            'Genes': 'gene_name',
            'First.Protein.Description': 'FirstProteinDescription',
            'Run': 'Sample'
        })
        table.columns = table.columns.str.replace('.', '', regex=False)
        assay = pd.pivot_table(table, values='PGMaxLFQ', index='Accession', columns='Sample')
        rdata = table[['Accession', 'ProteinIds', 'ProteinNames',
                       'gene_name', 'FirstProteinDescription']]
        rdata = rdata[rdata['Accession'].isin(assay.index)]
        rdata = rdata.drop_duplicates('Accession')

        rdata['gene_name'] = np.where(rdata['gene_name'] == np.nan, rdata['Accession'], rdata['gene_name'])
        rdata['gene_name'] = rdata.gene_name.str.split(';').str[0]
        rdata['Description'] = rdata['FirstProteinDescription'] + ' GN=' + rdata['gene_name'] + ' PE'
        rdata = rdata.sort_values('Accession')
        return assay, rdata
