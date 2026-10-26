"""   Import DIA-NN Output

DIA-NN reports protein-group quantities in two ways:

1. Protein-group matrix (report.pg_matrix.tsv; recommended). One row per
   protein group and one column per run (the column names are the raw file
   paths). Written by DIA-NN 1.8, 1.9 and 2.x.
2. Main report (report.tsv in DIA-NN <= 1.9; report.parquet in DIA-NN 2.x).
   One row per precursor and run; OmicScope uses the PG.MaxLFQ column and
   keeps protein groups with PG.Q.Value <= 0.01.

DIA-NN does not store the sample conditions, so a pdata file is required.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""

import numpy as np
import pandas as pd

from .utils import (clean_run_columns, drop_decoys_contaminants, read_table,
                    require_pdata)

# Annotation columns of DIA-NN matrices (with and without dots, as some
# users/tools strip them).
ANNOTATION = {'Protein.Group': 'Accession', 'ProteinGroup': 'Accession',
              'Protein.Ids': 'ProteinIds', 'ProteinIds': 'ProteinIds',
              'Protein.Names': 'ProteinNames', 'ProteinNames': 'ProteinNames',
              'Genes': 'gene_name', 'gene_name': 'gene_name',
              'First.Protein.Description': 'FirstProteinDescription',
              'FirstProteinDescription': 'FirstProteinDescription',
              'N.Sequences': 'NSequences', 'N.Proteotypic.Sequences': 'NProteotypicSequences'}


class Input:

    def __init__(self, Table, pdata=None, pg_qvalue=0.01):
        """ DIA-NN output

        Args:
            Table (str): Path to report.pg_matrix.tsv (recommended) or to the
            main report (report.tsv, or report.parquet if pyarrow is installed).
            pdata (str): Path to the phenotype data (.xlsx, .xls or .csv).
            Sample names are the run names (raw file names; folders and
            extensions such as .raw, .d or .mzML are ignored).
            pg_qvalue (float, optional): Protein-group q-value cutoff applied
            to the main report. Defaults to 0.01.
        """
        self.Table = Table
        self.pg_qvalue = pg_qvalue
        self.pdata = require_pdata(pdata, 'DIA-NN output',
                                   'the run names (raw file names, as in the "Run" column '
                                   'of report.tsv or the column names of report.pg_matrix.tsv)')
        self.assay, self.rdata = self.assay_and_rdata()
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())

    def assay_and_rdata(self):
        table = read_table(self.Table)
        if 'Run' in table.columns:
            assay, rdata = self.main_report(table)
        else:
            assay, rdata = self.pg_matrix(table)
        rdata = rdata.drop_duplicates('Accession')
        rdata = rdata[rdata['Accession'].isin(assay.index)]
        assay = assay[assay.index.isin(rdata['Accession'])]
        if 'gene_name' not in rdata.columns:
            rdata['gene_name'] = np.nan
        rdata['gene_name'] = rdata['gene_name'].fillna(rdata['Accession']).astype(str)
        rdata['gene_name'] = rdata['gene_name'].str.split(';').str[0]
        if 'FirstProteinDescription' in rdata.columns and rdata['FirstProteinDescription'].notna().any():
            rdata['Description'] = (rdata['FirstProteinDescription'].fillna('').astype(str) +
                                    ' GN=' + rdata['gene_name'] + ' PE')
        rdata = rdata.sort_values('Accession').reset_index(drop=True)
        assay = assay.loc[rdata['Accession']]
        return assay, rdata

    def main_report(self, table):
        """Precursor-level report (report.tsv / report.parquet)."""
        # Older OmicScope examples stripped the dots from column names
        table.columns = [c if c in ANNOTATION or c == 'Run' else c.replace('.', '')
                         for c in table.columns]
        table = table.rename(columns=ANNOTATION)
        quant = [q for q in ['PGMaxLFQ', 'PGNormalised', 'PGQuantity'] if q in table.columns]
        if len(quant) == 0:
            raise KeyError('OmicScope did not find the "PG.MaxLFQ" column in the DIA-NN report. '
                           'Please, use the report.pg_matrix.tsv file instead.')
        if 'Decoy' in table.columns:
            table = table[table['Decoy'] == 0]
        if 'QValue' in table.columns:
            table = table[table['QValue'] <= 0.01]
        if 'PGQValue' in table.columns:
            table = table[table['PGQValue'] <= self.pg_qvalue]
        table = drop_decoys_contaminants(table, 'Accession')
        table = table.assign(Run=clean_run_columns(table['Run'].astype(str)))
        assay = pd.pivot_table(table, values=quant[0], index='Accession',
                               columns='Run', aggfunc='max')
        assay = assay.replace(0, np.nan)
        assay.columns.name = None
        rdata = table.iloc[:, table.columns.isin(set(ANNOTATION.values()))]
        return assay, rdata

    def pg_matrix(self, table):
        """Protein-group matrix (report.pg_matrix.tsv)."""
        table = table.rename(columns=ANNOTATION)
        if 'Accession' not in table.columns:
            raise KeyError('OmicScope did not find the "Protein.Group" column. Please, confirm the file '
                           'is a DIA-NN report.pg_matrix.tsv or report.tsv.')
        table = drop_decoys_contaminants(table, 'Accession')
        annotation = list(table.columns[table.columns.isin(set(ANNOTATION.values()))])
        rdata = table[annotation]
        assay = table.set_index('Accession').drop(columns=[c for c in annotation if c != 'Accession'])
        assay = assay.apply(pd.to_numeric, errors='coerce')
        assay.columns = clean_run_columns(assay.columns)
        return assay, rdata
