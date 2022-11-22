
from copy import copy

import omicscope as omics
general = omics.Omicscope(Table = 'C:/Users/Guilherme/Desktop/general.xls', ControlGroup = None, Method = 'General', pvalue = 'pAdjusted')

general.quant_data


def tukey_correction(df):
    """Apply tukey hsd algorithm developed by
    Statmodels

    Args:
        df (DataFrame): Dataframe 
        i (Index): [description]

    Returns:
       padj [list]: Adjusted p-value for each differentially regulated gene and
       comparison
       comparison [list]: Comparison for which p-value was evaluated
    """
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    df = copy(df)
    df = df.reset_index()
    df.columns = ['Comparison', 'abundance']
    df['Comparison'] = df.Comparison.str.split('.').str[-1]
    df = pairwise_tukeyhsd(endog=df['abundance'],
                            groups=df['Comparison'],
                            alpha=0.05)
    df = pd.DataFrame(data=df._results_table.data[1:], columns=df._results_table.data[0])
    df['Comparison'] = df.group1 + '-' + df.group2
    padj = list(df['p-adj'])
    comparison = list(df['Comparison'])
    return padj, comparison


def Tukey_hsd(quant_data):
    """Perform Tukey HSD test for each differentially regulated entity
    found according ANOVA's test.

    Args:
        quant_data (DataFrame): Normalized dataframe.

    Returns:
       padj [list]: Adjusted p-value for each differentially regulated entity
       in each comparison. Entities that were not differentially regulated
       returned 2
       comparison [list]: Comparison for which p-value was evaluated. Entities 
       that were not differentially regulated returned 2
    """
    df = copy(quant_data)
    df = df[df['pvalue'] < 0.05] #TODO: Trocar para pAdjust
    df = df.iloc[:, df.columns.str.contains('.', regex = False)]
    # Perform Tukey's post-hoc correction for each of differentially
    # regulated entity according ANOVA's test
    df = df.apply(lambda x: tukey_correction(x), axis = 1)
    df = pd.DataFrame(df.to_list(), columns=['pTukey',
                                              'Comparison'], index = df.index)
    return df
import numpy as np
teste = copy(general.quant_data)
teste