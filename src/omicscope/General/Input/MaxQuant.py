
from copy import copy

import pandas as pd

from .utils import drop_decoys_contaminants, read_table, require_pdata


class Input:

    def __init__(self, Table, pdata=None,
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
            'Intensity', 'iBAQ', 'Reporter intensity corrected' (TMT/iTRAQ) and
            'Reporter intensity'. If the chosen columns are absent, OmicScope uses
            'Reporter intensity corrected' (isobaric labelling) or 'Intensity'.
        """
        self.Table = Table
        self.quant_strategy = quant_strategy
        self.table = read_table(Table)
        self.rdata = self.rdata()
        self.pdata = require_pdata(pdata, 'MaxQuant output',
                                   'the experiment names used in the quantitative columns of '
                                   'proteinGroups.txt (e.g. "CoV1" for "LFQ intensity CoV1")')
        self.assay = self.assay()
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())

    def rdata(self):
        df = copy(self.table)
        # Filter Contaminants and Reverses
        df = drop_decoys_contaminants(df, 'Protein IDs')
        for column in ['Reverse', 'Potential contaminant', 'Only identified by site']:
            if column in df.columns:
                df = df[df[column] != '+']
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
        # First gene of each group; missing genes are taken from the FASTA header
        if 'gene_name' not in df.columns:
            df['gene_name'] = None
        if 'Description' in df.columns:
            from_header = df['Description'].astype(str).str.extract(r'GN=(\S+)')[0]
            df['gene_name'] = df['gene_name'].fillna(from_header)
        df['gene_name'] = df['gene_name'].fillna(df['Accession']).astype(str).str.split(';').str[0]
        return df

    def assay(self):
        QUANTIFICATION_STRATEGIES = ['Intensity', 'LFQ intensity', 'iBAQ',
                                     'Reporter intensity corrected', 'Reporter intensity']
        quant_strategy = self.quant_strategy
        if quant_strategy not in QUANTIFICATION_STRATEGIES:
            raise ValueError(f'{quant_strategy} is not a MaxQuant quantification strategy.\n' +
                  f'User should try one of: {", ".join(QUANTIFICATION_STRATEGIES)}.')
        assay_raw = copy(self.table)
        assay_raw = assay_raw.rename(columns={'Majority protein IDs': "Accession"})
        assay_raw = assay_raw[assay_raw['Accession'].isin(self.rdata['Accession'])]
        assay_raw = assay_raw.set_index('Accession')
        for strategy in [quant_strategy, 'Reporter intensity corrected', 'Intensity']:
            assay = self.quant_columns(assay_raw, strategy)
            if len(assay.columns) > 0:
                break
        self.quant_strategy = strategy
        return assay

    @staticmethod
    def quant_columns(assay_raw, strategy):
        """Per-sample columns '<strategy> <sample>' (totals and counts are skipped)."""
        columns = pd.Series(assay_raw.columns, dtype=str)
        samples = columns.str.extract('^' + strategy + r' (?!count|corrected|peptides)(.+)$')[0]
        if strategy.startswith('Reporter'):
            # Several experiments: '<strategy> 1 Exp1'; the summed '<strategy> 1' is dropped
            per_experiment = samples.str.contains(r'^\d+ .+', na=False)
            if per_experiment.any():
                samples[~per_experiment] = None
        keep = samples.notna().values
        assay = assay_raw.loc[:, keep].apply(pd.to_numeric, errors='coerce')
        assay.columns = list(samples[keep])
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
