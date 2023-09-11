"""   Import PatternLab Output

PatternLab V allows user to export Excel file with filtered identified
proteins, sample information, biological and technical replicates.

Here, we advise users to export data with at least 1 unique peptide/protein
and protein abundance normalized by XIC.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""
from copy import copy

import numpy as np
import pandas as pd


class Input:
    def __init__(self, Table):
        """  PatternLab V output for OmicScope input

        Args:
            Table (str): Path to PatternLab V output (.xlsx)
        """
        self.Table = Table
        patternlab = self.PatternLab()
        assay = patternlab[0]  # Expression data
        assay.columns = patternlab[1].Sample
        norm_factor = patternlab[3]
        assay = assay.divide(norm_factor, axis='columns')
        assay = assay.replace(0, np.nan)
        self.assay = assay
        self.pdata = patternlab[1]  # Phenotype data
        self.rdata = patternlab[2]  # row/gene data
        # Extract analyzed conditions
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())

    def PatternLab(self):
        """ Extract all PatternLab V information
        """

        df = self.Table
        PLV_output = pd.read_excel(df)
        # filtering contaminants
        PLV_output = PLV_output[~PLV_output.Locus.str.startswith('contaminant')]
        PLV_output = PLV_output[~PLV_output.Locus.str.startswith('Reverse')]
        PLV_output = PLV_output.reset_index(drop=True)

        # Getting the Conditions used for PLV experiment
        labels = pd.read_excel(df, sheet_name='Class Description')
        labels = dict(zip(labels['Label'].astype(str), labels['Description']))

        # Defining assay
        # As suggested by PLV developers, for quantitation use XIC values
        assay = PLV_output.iloc[:, PLV_output.columns.str.contains('XIC')]

        # Defining pdata
        # PLV output columns have several information regarding
        # phenotype of each sample
        pdata = pd.DataFrame(assay.columns, columns=['Data'])
        pdata['Data'] = pdata['Data'].str.split('\n').str[1:]
        pdata['Sample'] = pdata.Data.str[0].str.split(': ').str[1]
        pdata['Condition'] = pdata.Data.str[2].str.split(': ').str[1]
        pdata['Condition'] = pdata.Condition.replace(labels, regex=True)
        pdata['Biological'] = pdata.Data.str[1].str.split(': ').str[1]
        pdata = pdata.iloc[:, 1:]

        # Defining rdata
        rdata = PLV_output[['Locus', 'Description']]
        Accession = rdata.Locus
        gene_name = rdata['Description'].str.split('GN=').str[-1]
        gene_name = gene_name.str.split(' ').str[0]
        rdata = rdata.assign(Accession=Accession,
                             gene_name=gene_name)
        # Normalization factors
        norm_factor = pd.read_excel(df, sheet_name='Normalization Factors')
        norm_factor = dict(zip(norm_factor['File Name'], norm_factor['Normalization Factor']))
        return (assay, pdata, rdata, norm_factor)

    def filtering_data(self):
        assay = copy(self.assay)
        rdata = copy(self.rdata)
        assay.index = rdata.Locus
        rdata.index = rdata.Locus
        not_null_matrix = assay[assay > 0]
        # If user use percentage of null values, filter according to percentage
        if isinstance(self.filtering_method, int):
            ncol = len(not_null_matrix.columns)
            nan_count = not_null_matrix.apply(lambda x: (1 -
                                                         x.isna().sum()/ncol)*100,
                                              axis=1)
            nan_count = nan_count[nan_count > self.filtering_method]
            self.assay = assay[assay.index.isin(nan_count.index)]
            self.rdata = rdata[rdata['Locus'].isin(nan_count.index)]
        # If user use "minimum", proteins will be select if at least 2 samples
        # present identified protein
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
            self.rdata = rdata[rdata['Locus'].isin(valid_values.index)]
