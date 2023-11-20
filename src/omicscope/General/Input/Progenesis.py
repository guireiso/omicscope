"""   Import Progenesis Output

Progenesis QI for proteomics allows user to export proteomic
quantitative data as .csv.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""

import numpy as np
import pandas as pd


class Input:
    def __init__(self, Table, **kwargs):
        """   Progenesis output for OmicScope input

        Args:
            Table (str): Path to Progenesis output (.xslx or .csv)
            pdata (DataFrame): Phenotype data. Default: None.
            If user desires it is possible to defines externally
            the sample conditions.
        """
        self.Table = Table
        data = self.progenesis(**kwargs)
        self.assay = data[0]
        self.rdata = data[1]
        self.rdata['gene_name'] = self.rdata['Description'].str.split(
            'GN=').str[1].str.split(' ').str[0]
        self.pdata = self.pdata()
        self.assay.columns = self.assay.columns.str.split('.', regex=False).str[0]
        self.assay.index = self.rdata['Accession']

    def progenesis(self, **kwargs):
        """ Extract rdata and assay from Progenesis output
        """

        # Read .csv or .xlsx/.xls files
        Table = self.Table
        try:
            df = pd.read_csv(Table, header=0)
        except UnicodeDecodeError:
            df = pd.read_excel(Table, header=0)
        df = df.replace("Best identification", np.nan)
        df.iloc[1] = df.iloc[1].astype(str)

        # Extracting conditions
        self.Conditions = list(df.iloc[0, :].dropna().drop_duplicates())

        # Adjusting dataframe
        try:
            df = df.drop(df.loc[:, 'Raw abundance':].columns, axis=1)
        except KeyError:
            pass
        finally:
            df.iloc[0] = df.iloc[1].str.cat(df.iloc[0]
                                            .where(df.iloc[0].notnull() & df.iloc[0].ne(''))
                                            .ffill(), '.').fillna(df.iloc[1])
        df = df.drop(1).reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.drop(df.index[0])
        df = df.apply(pd.to_numeric, errors='ignore')
        df['Description'] = df.Description.astype(str)

        # Filtering unique peptides
        if 'UniquePeptides' in kwargs:
            df = df[df['Unique peptides'] >= kwargs['UniquePeptides']]
        df = df.reset_index(drop=True)
        df = df.rename(columns={'Anova (p)': 'pvalue',
                                'q Value': 'pAdjusted'})
        # Defining assay
        assay = df.iloc[:, df.columns.str.contains(".", regex=False)]
        # Defining rdata
        rdata = df.iloc[:, ~df.columns.str.contains(".", regex=False)]
        return (assay, rdata)

    def pdata(self):
        """Defining pdata
        """
        pdata = pd.Series(self.assay.columns)
        samples = pdata.str.split('.').str[0]
        pdata = pdata.str.split('.').str[-1]
        pdata = pd.DataFrame([samples, pdata]).T
        pdata.columns = ['Sample', 'Condition']
        # Defining Biological replicates
        Biological = []
        for i in self.Conditions:
            df = pdata
            df = df[df.Condition == i]
            for j in range(1, len(df) + 1):
                Biological.append(j)
        pdata['Biological'] = Biological
        return (pdata)
