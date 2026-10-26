"""   Import Spectronaut Output

Spectronaut (Biognosys) exports protein-group quantities in two layouts:

1. Report (long format, e.g. a protein-level report built from the
   "BGS Factory Report" schema). One row per protein group and run, with
   R.Condition, R.FileName, R.Replicate, PG.ProteinGroups, PG.Genes and
   PG.Quantity columns. Conditions are read from R.Condition, so pdata
   is optional.
2. Pivot report (wide format, e.g. "Protein Quant (Pivot)"). One row per
   protein group and one "[n] <run>.PG.Quantity" column per run. It does
   not contain the conditions, so a pdata file is required.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""

import re

import numpy as np
import pandas as pd

from .utils import (clean_run_columns, drop_decoys_contaminants, read_pdata,
                    read_table, require_pdata)

RDATA = {'PG.ProteinGroups': 'Accession',
         'PG.Genes': 'gene_name',
         'PG.ProteinDescriptions': 'ProteinDescription',
         'PG.ProteinNames': 'ProteinNames',
         'PG.ProteinAccessions': 'ProteinAccessions',
         'PG.Organisms': 'Organisms',
         'PG.MolecularWeight': 'MolecularWeight',
         'PG.NrOfStrippedSequencesIdentified (Experiment-wide)': 'NrOfStrippedSequences',
         'PG.NrOfPrecursorsIdentified (Experiment-wide)': 'NrOfPrecursors'}

PIVOT_COLUMN = re.compile(r'^\[\d+\]\s*(.+)\.PG\.Quantity$')


def _numeric(values):
    """Spectronaut writes 'Filtered'/'NaN' for missing values and may use
    decimal commas depending on the system locale."""
    if values.dtype == object:
        values = values.astype(str).str.replace(',', '.', regex=False)
    return pd.to_numeric(values, errors='coerce')


class Input:

    def __init__(self, Table, pdata=None):
        """ Spectronaut output

        Args:
            Table (str): Path to the Spectronaut report (.tsv/.xls/.csv).
            pdata (str, optional): Path to the phenotype data. Required for
            pivot reports; sample names are the run names (R.FileName).
        """
        self.Table = Table
        data = read_table(Table)
        if {'R.FileName', 'PG.Quantity'}.issubset(data.columns):
            self.assay, self.rdata, self.pdata = self.report(data)
            if pdata is not None:
                self.pdata = read_pdata(pdata)
        else:
            self.assay, self.rdata = self.pivot_report(data)
            self.pdata = require_pdata(pdata, 'The Spectronaut pivot report',
                                       'the run names (R.FileName), i.e. the text between "[n] " and '
                                       '".PG.Quantity" in the column headers')
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())

    @staticmethod
    def _rdata(data):
        if 'PG.ProteinGroups' not in data.columns:
            raise KeyError('OmicScope did not find the "PG.ProteinGroups" column in the Spectronaut report. '
                           'Please, include PG.ProteinGroups, PG.Genes and PG.Quantity in the report schema.')
        rdata = data.iloc[:, data.columns.isin(RDATA.keys())].rename(columns=RDATA)
        rdata = rdata.drop_duplicates('Accession')
        if 'gene_name' not in rdata.columns:
            rdata['gene_name'] = np.nan
        rdata['gene_name'] = rdata['gene_name'].fillna(rdata['Accession']).astype(str).str.split(';').str[0]
        if 'ProteinDescription' in rdata.columns:
            rdata['Description'] = (rdata['ProteinDescription'].fillna('').astype(str).str.split(';').str[0] +
                                    ' GN=' + rdata['gene_name'] + ' PE')
        return rdata

    def report(self, data):
        """Long report: one row per protein group (or precursor) and run."""
        for decoy in ['PG.IsDecoy', 'EG.IsDecoy']:
            if decoy in data.columns:
                data = data[~data[decoy].astype(str).str.lower().isin(['true', '1'])]
        data = drop_decoys_contaminants(data, 'PG.ProteinGroups')
        data = data.assign(**{'PG.Quantity': _numeric(data['PG.Quantity']),
                              'R.FileName': clean_run_columns(data['R.FileName'].astype(str))})
        assay = pd.pivot_table(data, values='PG.Quantity', index='PG.ProteinGroups',
                               columns='R.FileName', aggfunc='max')
        assay.columns.name, assay.index.name = None, 'Accession'
        if 'R.Condition' not in data.columns:
            raise KeyError('The Spectronaut report does not contain "R.Condition". '
                           'Add it to the report schema or provide a pdata file.')
        pdata = data[['R.FileName', 'R.Condition'] +
                     (['R.Replicate'] if 'R.Replicate' in data.columns else [])].drop_duplicates('R.FileName')
        pdata = pdata.rename(columns={'R.FileName': 'Sample', 'R.Condition': 'Condition',
                                      'R.Replicate': 'Biological'})
        pdata['Condition'] = pdata['Condition'].astype(str)
        if 'Biological' not in pdata.columns:
            pdata['Biological'] = pdata.groupby('Condition').cumcount() + 1
        pdata = pdata.sort_values(['Condition', 'Biological']).reset_index(drop=True)
        rdata = self._rdata(data)
        rdata = rdata[rdata['Accession'].isin(assay.index)].sort_values('Accession').reset_index(drop=True)
        assay = assay.loc[rdata['Accession']]
        return assay, rdata, pdata

    def pivot_report(self, data):
        """Pivot report: one '[n] <run>.PG.Quantity' column per run."""
        quant = [c for c in data.columns if PIVOT_COLUMN.match(str(c))]
        if len(quant) == 0:
            raise KeyError('OmicScope did not find "[n] <run>.PG.Quantity" columns (pivot report) nor the '
                           'R.FileName/PG.Quantity columns (long report) in the Spectronaut file.')
        data = drop_decoys_contaminants(data, 'PG.ProteinGroups')
        rdata = self._rdata(data).reset_index(drop=True)
        assay = data.set_index('PG.ProteinGroups')[quant].apply(_numeric)
        assay.index.name = 'Accession'
        assay.columns = clean_run_columns([PIVOT_COLUMN.match(c).group(1) for c in quant])
        assay = assay[~assay.index.duplicated()].loc[rdata['Accession']]
        return assay, rdata

