"""  Statistic Module comprises all statistical tests performed
by OmicScope for differential expression analysis.

    Returns:
        DataFrame: DataFrame with recommended statistical
        analysis performed
"""
import os

import pandas as pd


def perform_static_stat(self):
    """log2 transformation for expression
    """
    from copy import copy

    import numpy as np
    expression = copy(self.expression)
    log = copy(self.logTransformed)
    rdata = copy(self.rdata)
    pvalue = copy(self.pvalue)
    # Log-normalize data if it was not
    if log is False:
        expression = expression.replace(0, np.nan)
        expression = np.log2(expression)
        # Apply t-test if len(conditions) == 2
    if len(self.Conditions) == 2:
        from .Static_Statistics import ttest
        params = [self.ind_variables, self.ctrl, self.experimental[0], expression, rdata,
                  pvalue]
        data = ttest(params=params)
        data = params[4].merge(data, on='Accession')

    # Apply ANOVA if len(conditions) > 2
    elif len(self.Conditions) > 2:
        from .Static_Statistics import anova
        params = [expression, rdata, self.Conditions, pvalue]
        data = anova(params=params)
        data = params[1].merge(data, on='Accession')
    data = data.sort_values('pvalue')
    data = data.reset_index(drop=True)
    # Filtering Keratin
    if self.ExcludeContaminants is True:
        path = os.path.dirname(os.path.abspath(__file__))
        contaminants = pd.read_csv(path+'/contaminants.csv')[['Accession', 'gene_name']]
        data = data[~data['Accession'].isin(contaminants['Accession'])]
        data = data[~data['gene_name'].isin(contaminants['gene_name'].str.split(' ').explode())]
    data = data.reset_index(drop=True)
    data = copy(data)
    data.iloc[:, data.columns.str.contains(
        '.', regex=False)] = data.iloc[:, data.columns.str.contains('.', regex=False)].replace(np.nan, 0)
    return (data)


def perform_longitudinal_stat(self):
    """log2 transformation for expression
    """
    from copy import copy

    import numpy as np
    degrees_of_freedom = copy(self.degrees_of_freedom)
    expression = copy(self.expression)
    log = copy(self.logTransformed)
    rdata = copy(self.rdata)
    pdata = copy(self.pdata)
    pvalue = copy(self.pvalue)
    ctrl = copy(self.ctrl)
    # Log-normalize data if it was not
    if log is False:
        expression = expression.replace(0, 0.01)
        expression = np.log2(expression)
    from .Longitudinal_Stat import Longitudinal_Stats
    data = Longitudinal_Stats(assay=expression, pdata=pdata,
                              degrees_of_freedom=degrees_of_freedom, pvalue=pvalue, ctrl=ctrl)
    data = rdata.merge(data, on='Accession')
    # longitudinal modules
    data = data.sort_values('pvalue')
    data = data.reset_index(drop=True)
    # # Filtering Keratin
    if self.ExcludeContaminants is True:
        path = os.path.dirname(os.path.abspath(__file__))
        contaminants = pd.read_csv(path+'/contaminants.csv')[['Accession', 'gene_name']]
        data = data[~data['Accession'].isin(contaminants['Accession'])]
        data = data[~data['gene_name'].isin(contaminants['gene_name'].str.split(' ').explode())]
    data = data.reset_index(drop=True)
    data = copy(data)
    data.iloc[:, data.columns.str.contains(
        '.', regex=False)] = data.iloc[:, data.columns.str.contains('.', regex=False)].replace(np.nan, 0)
    return (data)
