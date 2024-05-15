"""   Import Proteome Discoverer Output

Proteome Discoverer (PD) allows user to export proteomic
quantitative data as .xlsx.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""

import numpy as np
import pandas as pd


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
        data = self.ProteomeDiscoverer(**kwargs)
        self.assay = data[0]
        self.rdata = data[1]
        self.pdata = data[2]
        self.rdata['gene_name'] = self.rdata['Description'].str.split(
            'GN=').str[1].str.split(' ').str[0]
        self.assay.index = self.rdata['Accession']

    def ProteomeDiscoverer(self):
        """ Extract rdata and assay from Proteome Discoverer output
        """
        # Read .csv or .xlsx/.xls files
        Table = self.Table
        try:
            df = pd.read_csv(Table, header=0)
        except UnicodeDecodeError:
            df = pd.read_excel(Table, header=0)
        
        df = df.dropna(subset=['Checked']) #checking if it is protein
        # rdata
        rdata = df.copy()
        rdata_columns = ['# AAs','# PSMs', '# Peptides', '# Protein Groups', '# Razor Peptides',
                            '# Unique Peptides', 'Accession', 'Coverage [%]', 'Description',
                            'Ensembl Gene ID', 'Entrez Gene ID', 'Gene ID', 'Gene Symbol', 'MW [kDa]',
                            'pvalue', 'pAdjusted']
        rdata.columns = rdata.columns.str.split(':', regex=False).str[0]
        rdata = rdata.rename(columns={'Abundance Ratio P-Value': 'pvalue',
                                'Abundance Ratio Adj. P-Value':'pAdjusted'})
        rdata = df.iloc[:,df.columns.isin(rdata_columns)]
        if not {'Accession', 'Description'}.issubset(rdata.columns):
                raise ValueError(
                    """OmicScope did not find "Accession" or "Description" columns in the dataset.
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
        pdata = pd.DataFrame({'Sample':assay.columns})
        pdata['Condition'] = pdata.Sample.str.split(', ').str[-1]
        self.Conditions = pdata.Condition
        pdata['Biological'] = pdata.Sample.str.split(', ').str[-2]
        return (assay, rdata, pdata)
