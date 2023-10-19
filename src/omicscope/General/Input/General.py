"""  Import a General data.

Import data from any source, even transcriptomics/ metabolomics/etc.

To construct its own data, user should prepare an Excel file with 3 sheets.

1. The 1st sheet must be the assay data, in which the first row is the
sample name and furthers are the expression levels of each
protein/metabolite/transcripts.

2. The 2nd sheet must be row-data, which is the information regarding each row
(aka protein/metabolite/transcripts) contained in assay data. It is crucial to
have at least two columns: 'Accession' and 'Gene_name'/'Description'.

*Accession. It must be a unique identifier to each row (e.g. Uniprot Accession
for proteomics, Gene ID (Entrez) for transcriptomics).

*Gene_name or Description. The Gene name (or symbol) that makes Accession more
readable; Alternatively, for proteomic results, users can also input
Description column, which must contain 'GN=' as exported from fasta UniprotKB.

3.  The 3rd sheet must be phenotype-data, which is the information regarding
each sample. Here the user must specify for each sample the biological
Condition ('Condition'), biological replicates ('Biological'), and Technical
replicates ('TechRep') that was designed for experiment.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""


class Input:
    def __init__(self, Table, **kwargs):
        """  Format general data for OmicScope input

        Args:
            Table (str): excel file.
        """

        self.Table = Table
        self.assay = self.readxl(sheetno=0)
        if ['gene_name'] is self.readxl(sheetno=1).columns:
            self.rdata = self.readxl(sheetno=1)
        else:
            rdata = self.readxl(sheetno=1)
            rdata.columns = rdata.columns.str.replace('.', '', regex=False)
            self.rdata = rdata
            if 'gene_name' not in self.rdata.columns:
                self.rdata['gene_name'] = self.rdata['Description'].str.split(
                    'GN=').str[1].str.split(' ').str[0]
        self.pdata = self.readxl(sheetno=2)
        self.Conditions = list(self.pdata['Condition'].drop_duplicates())
        self.assay.index = self.rdata.Accession

    def readxl(self, sheetno):
        """  Read excel file

        Args:
            sheetno (int): position of sheet to be imported

        """
        import pandas as pd
        excel = self.Table
        df = pd.read_excel(excel,
                           sheet_name=sheetno)
        return (df)
