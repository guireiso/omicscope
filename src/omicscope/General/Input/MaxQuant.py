
from copy import copy

import pandas as pd


class Input:

    def __init__(self, Table, pdata,
                 quant_strategy='LFQ intensity'):
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
        self.quant_strategy = quant_strategy
        self.rdata = self.rdata()
        try:
            self.pdata = pd.read_excel(pdata)
        except ValueError:
            self.pdata = pd.read_csv(pdata)
        self.assay = self.assay()
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())

    def rdata(self):
        df = pd.read_csv(copy(self.Table), sep='\t')
        # Filter Contaminants and Reverses
        df = df[~df['Protein IDs'].str.contains('REV__|CON__')]
        # Filter Andromeda Score >0
        df = df[df['Score'] > 0]
        # Columns of interest
        columns = ['Majority protein IDs', 'Peptide counts (all)',
                   'Peptide counts (razor+unique)', 'Peptide counts (unique)',
                   'Gene names', 'Fasta headers', 'Number of proteins', 'Peptides',
                   'Razor + unique peptides', 'Unique peptides', 'Mol. weight [kDa]',
                   'Sequence length', 'Score']
        df = df.iloc[:,df.columns.isin(columns)]
        df = df.rename(columns={'Majority protein IDs': "Accession",
                                'Fasta headers': 'Description',
                                'Gene names': 'gene_name',
                                'Mol. weight [kDa]': "Mol weight [kDa]"})
        return df

    def assay(self):
        QUANTIFICATION_STRATEGIES = ['Intensity', 'LFQ intensity', 'iBAQ', ]
        quant_strategy = self.quant_strategy
        assay = []
        if quant_strategy in QUANTIFICATION_STRATEGIES:
            assay_raw = pd.read_csv(copy(self.Table), sep='\t')
            assay_raw = assay_raw.rename(columns={'Majority protein IDs': "Accession"})
            assay_raw = assay_raw[assay_raw['Accession'].isin(self.rdata['Accession'])]
            assay_raw = assay_raw.set_index('Accession')
            assay = assay_raw.iloc[:, assay_raw.columns.str.startswith(quant_strategy)]
            assay.columns = assay.columns.str.removeprefix(quant_strategy + ' ')
            if len(assay.columns)==0:
                assay = assay_raw.iloc[:, assay_raw.columns.str.startswith('Intensity')]
                assay.columns = assay.columns.str.removeprefix('Intensity' + ' ')
        else:
            raise ValueError(f'{quant_strategy} is not a MaxQuant quantification strategy.\n' +
                  f'User should try one of: {", ".join(QUANTIFICATION_STRATEGIES)}.')
        return assay

    def filtering_data(self):
        assay = copy(self.assay)
        rdata = copy(self.rdata)
        not_null_matrix = assay[assay > 0]
        # If user use percentage of null values, filter according to percentage
        if isinstance(self.filtering_method, int):
            ncol = len(not_null_matrix.columns)
            nan_count = not_null_matrix.apply(lambda x: (1 -
                                                         x.isna().sum()/ncol)*100,
                                              axis=1)
            nan_count = nan_count[nan_count > self.filtering_method]
            self.assay = assay[assay.index.isin(nan_count.index)]
            self.rdata = rdata[rdata['Accession'].isin(nan_count.index)]

        elif self.filtering_method == 'minimum':
            pdata = copy(self.pdata)
            conditions = dict(zip(pdata.Sample, pdata.Condition))
            valid_values = assay.rename(columns=conditions)
            valid_values = not_null_matrix.rename(columns=conditions)
            valid_values[valid_values.notna()] = 1
            valid_values = valid_values.groupby(by=valid_values.columns,
                                                axis=1).sum()
            valid_values = valid_values[valid_values > 1].dropna()
            self.assay = assay[assay.index.isin(valid_values.index)]
            self.rdata = rdata[rdata['Accession'].isin(valid_values.index)]
