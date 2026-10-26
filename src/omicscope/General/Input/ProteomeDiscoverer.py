"""   Import Proteome Discoverer Output

Proteome Discoverer (PD) allows user to export proteomic
quantitative data as .xlsx.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""

import pandas as pd

# Sample types that Proteome Discoverer writes before the study factors
SAMPLE_TYPES = {'Sample', 'Control', 'Blank', 'Standard', 'Identification Only', 'Reference'}


def pd_pdata(columns):
    """Build pdata from abundance columns such as
    'Abundances (Normalized): F17: Sample, SCZ_5, CBD (SCZ)' (label-free) or
    'Abundances (Normalized): F1: 126, Sample, CTRL' (TMT).
    The last study factor is the Condition; a second-to-last factor is used
    as biological replicate, otherwise each file/channel is one replicate.
    """
    rows = []
    for column in columns:
        label = column.split(': ', 1)[-1]
        tokens = [t.strip() for t in label.split(',')]
        file_channel = tokens[0]
        first = file_channel.split(':')
        if first[-1].strip() in SAMPLE_TYPES:
            file_channel = ':'.join(first[:-1]).strip()
        factors = [t for t in tokens[1:] if t not in SAMPLE_TYPES]
        condition = factors[-1] if factors else file_channel
        biological = factors[-2] if len(factors) >= 2 else file_channel
        rows.append([column, condition, biological])
    return pd.DataFrame(rows, columns=['Sample', 'Condition', 'Biological'])


class Input:
    def __init__(self, Table, **kwargs):
        """   Proteome Discoverer output for OmicScope input

        Args:
            Table (str): Path to Proteome Discoverer output (.xslx)
            pdata (DataFrame): Phenotype data. Default: None.
            If user desires it is possible to defines externally
            the sample conditions.
        """
        self.Table = Table
        try: 
            data = self.ProteomeDiscovererV2()
        except AttributeError:
            data = self.ProteomeDiscovererV3()
        self.assay = data[0]
        self.rdata = data[1]
        self.pdata = data[2]
        self.assay.index = self.rdata['Accession']

    def ProteomeDiscovererV2(self):
        """ Extract rdata and assay from Proteome Discoverer output
        """
        # Read .csv or .xlsx/.xls files
        Table = self.Table
        try:
            df = pd.read_csv(Table, header=0)
        except UnicodeDecodeError:
            df = pd.read_excel(Table, header=0)
        
        df = df.dropna(subset=['Checked']) #checking if it is protein
        df = self.filter_proteins(df)
        # rdata
        rdata = df.copy()
        rdata_columns = ['# AAs','# PSMs', '# Peptides', '# Protein Groups', '# Razor Peptides',
                            '# Unique Peptides', 'Accession', 'Coverage [%]', 'Description',
                            'Ensembl Gene ID', 'Entrez Gene ID', 'Gene ID', 'Gene Symbol', 'MW [kDa]',
                            'pvalue', 'pAdjusted']
        rdata = df.iloc[:,df.columns.isin(rdata_columns)]
        if not {'Accession', 'Description'}.issubset(rdata.columns):
                raise ValueError(
                    """OmicScope did not find "Accession" or "Description" columns in the dataset.
                    Possible causes:
                    - Incorrect export from Proteome Discoverer.
                    - Missing variable assignments during export.
                    Please check your Proteome Discoverer export settings and ensure the data contains the necessary abundance information.
                    """)
        rdata['gene_name'] = rdata['Description'].str.split(
                'GN=').str[1].str.split(' ').str[0]
        # assay
        df = df.set_index('Accession')
        assay = df.iloc[:,df.columns.str.contains('Normalized')]
        if len(assay.columns) == 0:
            assay = df.iloc[:,df.columns.str.contains('Abundance:', regex=False)]
            if len(assay.columns) == 0:
                raise ValueError(
                    """OmicScope did not find "Abundances (Normalized)" or "Abundance" columns in the dataset.
                    Possible causes:
                    - Incorrect export from Proteome Discoverer.
                    - Missing variable assignments during export.
                    - Incompatible experimental design.
                    Please check your Proteome Discoverer export settings and ensure the data contains the necessary abundance information.
                    """)
        # pdata
        pdata = pd_pdata(assay.columns)
        self.Conditions = list(pdata.Condition.drop_duplicates())
        return (assay, rdata, pdata)

    def ProteomeDiscovererV3(self):
        """ Extract rdata and assay from Proteome Discoverer output
        """
        # Read .csv or .xlsx/.xls files
        Table = self.Table
        try:
            df = pd.read_csv(Table, header=0)
        except UnicodeDecodeError:
            df = pd.read_excel(Table, header=0)
        df = df.dropna(subset=['Checked']) #checking if it is protein
        df = self.filter_proteins(df)
        # rdata
        rdata = df.copy()
        rdata_columns = list(rdata.columns)
        AbundanceRatios = [x.startswith('AbundanceRatio') for x in rdata_columns]
        PvalueRatios = [x.startswith('AbundanceRatioP-Value') for x in rdata_columns]
        rdata_columns = ['Accession','Description','Sum PEP Score',
                        'Coverage [%]','# PSMs','# Unique Peptides',
                        'MW [kDa]','calc. pI','Gene Symbol','Gene ID'	
                        'pvalue', 'pAdjusted']
        rdata_columns.extend(AbundanceRatios)
        rdata_columns.extend(PvalueRatios)
        rdata = df.iloc[:,df.columns.isin(rdata_columns)]
        rdata = rdata.rename(columns={'Gene Symbol':'gene_name'})

        if not {'Accession', 'Description', 'gene_name'}.issubset(rdata.columns):
                raise ValueError(
                    """OmicScope did not find "Accession", "Description", 'Gene Symbol' columns in the dataset.
                    Possible causes:
                    - Incorrect export from Proteome Discoverer.
                    - Missing variable assignments during export.
                    Please check your Proteome Discoverer export settings and ensure the data contains the necessary abundance information.
                    """)
        # assay
        df = df.set_index('Accession')
        assay = df.iloc[:,df.columns.str.contains('Normalized')]
        if len(assay.columns) == 0:
            assay = df.iloc[:,df.columns.str.contains('Abundance:', regex=False)]
            if len(assay.columns) == 0:
                raise ValueError(
                    """OmicScope did not find "Abundances (Normalized)" or "Abundance" columns in the dataset.
                    Possible causes:
                    - Incorrect export from Proteome Discoverer.
                    - Missing variable assignments during export.
                    - Incompatible experimental design.
                    Please check your Proteome Discoverer export settings and ensure the data contains the necessary abundance information.
                    """)
        # pdata
        pdata = pd_pdata(assay.columns)
        self.Conditions = list(pdata.Condition.drop_duplicates())
        return (assay, rdata, pdata)

    @staticmethod
    def filter_proteins(df):
        """Keep master proteins and drop contaminants when PD exports these columns."""
        if 'Master' in df.columns:
            df = df[df['Master'].isna() | (df['Master'] == 'Master Protein')]
        if 'Contaminant' in df.columns:
            df = df[~df['Contaminant'].astype(str).str.lower().isin(['true', '1'])]
        return df