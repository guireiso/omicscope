
from copy import copy

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def ttest(self, params):
    """Perform a T-test on the data

    Args:
        params (list): list containing control and experimental groups,
        normalized dataframe and if the variables are dependent or related.
    """
    ind_variables = copy(params[0])
    ControlGroup = copy(params[1])
    Experimental = copy(params[2])
    quant_data = copy(params[3])
    pvalue = copy(params[5])

    if ind_variables is True:
        from scipy.stats import ttest_ind

        #  Perform an independent t-test
        result = ttest_ind(copy(quant_data.iloc[:, quant_data.columns.str.endswith(Experimental)]),
                           copy(
                               quant_data.iloc[:, quant_data.columns.str.endswith(ControlGroup)]),
                           axis=1, nan_policy='omit')
        self.Params['Params']['Stats_Workflow_3'] = 'OmicScope performed Independent T-test'

        #  Perform a sample-related t-test
    elif ind_variables is False:
        from scipy.stats import ttest_rel
        result = ttest_rel(copy(quant_data.iloc[:, quant_data.columns.str.endswith(Experimental)]),
                           copy(
                               quant_data.iloc[:, quant_data.columns.str.endswith(ControlGroup)]),
                           axis=1, nan_policy='omit')
        self.Params['Params']['Stats_Workflow_3'] = 'OmicScope performed Independent T-test'

    quant_data['pvalue'] = result[1]
    quant_data['pvalue'] = quant_data['pvalue'].replace(np.nan, 1)
    self.Params['Params']['Stats_Workflow_4'] = 'OmicScope replaceD p-value = NaN to 1.0'
    quant_data = quant_data.sort_values('pvalue')

    #  Correcting multilple hypothesis test according to fdr_bh
    quant_data['pAdjusted'] = multipletests(quant_data['pvalue'],alpha=params[6],
                                            method='fdr_tsbh')[1]
    self.Params['Params']['Stats_Workflow_5'] = 'OmicScope performed two-stage Benjamini-Hochberg correction'

    #  Mean abundance for each protein among conditions
    quant_data.loc[:, quant_data.columns.str.endswith(ControlGroup)] = np.exp2(
        copy(quant_data.loc[:, quant_data.columns.str.endswith(ControlGroup)]))
    quant_data['mean ' + ControlGroup] = quant_data.loc[:,
                                                        quant_data.columns.str.endswith(ControlGroup)].mean(axis=1)

    quant_data.loc[:, quant_data.columns.str.endswith(Experimental)] = np.exp2(
        copy(quant_data.loc[:, quant_data.columns.str.endswith(Experimental)]))
    quant_data['mean ' + Experimental] = quant_data.loc[:,
                                                        quant_data.columns.str.endswith(Experimental)].mean(axis=1)
    #  Mean abundance for each protein
    quant_data['TotalMean'] = quant_data.loc[:,
                                             quant_data.columns.str.contains('.')].mean(axis=1)
    #  Protein Fold change (Experimental/Control)
    quant_data['fc'] = quant_data['mean ' + Experimental] / \
        quant_data['mean ' + ControlGroup]
    #  Log2(FC)
    quant_data['log2(fc)'] = np.log2(quant_data['fc'])
    #  -log10(pvalue)
    quant_data[f'-log10({pvalue})'] = -np.log10(quant_data[pvalue])
    quant_data = quant_data.reset_index()
    return quant_data


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
    df = pd.DataFrame(
        data=df._results_table.data[1:], columns=df._results_table.data[0])
    df['Comparison'] = df.group1 + '-vs-' + df.group2
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
    df = df[df['pAdjusted'] < 0.05]
    df = df.iloc[:, df.columns.str.contains('.', regex=False)]
    # Perform Tukey's post-hoc correction for each of differentially
    # regulated entity according ANOVA's test
    df = df.apply(lambda x: tukey_correction(x), axis=1)
    df = pd.DataFrame(df.to_list(), columns=['pTukey',
                                             'Comparison'], index=df.index)
    return df


def anova(self, params):
    """Test for ANOVA for parametric and multiple comparisons

    Args:
        params ([type]): [description]
    """
    from scipy.stats import f_oneway
    self.Params['Params']['Stats_Workflow_3'] = 'OmicScope performed One-way ANOVA'
    quant_data = params[0]
    Conditions = params[2]
    pvalue = params[3]

    # Perform ANOVA test according to SCI-PY algorithm
    expression = []  # List of abundances for each group
    for i in Conditions:
        df = quant_data.iloc[:, quant_data.columns.str.endswith(i)]
        df = df.to_numpy()
        df = [x[~np.isnan(x)] for x in df]
        expression.append(df)
    stat = pd.DataFrame(expression).T
    stat = stat.apply(lambda x: f_oneway(*x)[1], axis=1)
    # Multiple hypothesis correction
    stat = pd.Series(stat)
    stat.index = quant_data.index
    quant_data['pvalue'] = stat
    quant_data['pvalue'] = quant_data['pvalue'].replace(np.nan, 1)
    self.Params['Params']['Stats_Workflow_4'] = 'OmicScope replaceD p-value = NaN to 1.0'
    quant_data = quant_data.sort_values('pvalue')
    # BH-correction
    quant_data['pAdjusted'] = multipletests(quant_data['pvalue'], alpha=params[4],
                                            method='fdr_tsbh', is_sorted=False, returnsorted=False)[1]
    self.Params['Params']['Stats_Workflow_5'] = 'OmicScope performed two-stage Benjamini-Hochberg correction'
    # Perform Tukey's Post-hoc correction
    Tukey = Tukey_hsd(quant_data)
    self.Params['Params']['Stats_Workflow_6'] = 'OmicScope performed Tukey HSD correction with pAdjusted < 0.05'
    quant_data = quant_data.join(Tukey)
    quant_data = quant_data.explode(['pTukey', 'Comparison'])
    quant_data['pTukey'] = quant_data.pTukey.astype(float)
    quant_data = quant_data.sort_values(['pTukey', 'pAdjusted'])
    # # Replace '2' with repeated values from ANOVA for
    # # non diffenrentially regulated entities.
    quant_data['pTukey'] = quant_data['pTukey'].replace(np.nan, 2)
    quant_data['pTukey'] = np.where(quant_data['pTukey'] == 2,
                                    quant_data['pAdjusted'],
                                    quant_data['pTukey'])
    # Replace '2' with right comparisons performed
    quant_data['Comparison'] = quant_data['Comparison'].replace(np.nan,
                                                                '-vs-'.join(Conditions))
    comp = copy(quant_data)
    comp['Comparison'] = comp.Comparison.str.split('-vs-')

    quant_data['Condition1'] = comp.apply(lambda x: x[x.index.str.endswith('.' + x.Comparison[0])].mean(),
                                          axis=1)
    quant_data['Condition2'] = comp.apply(lambda x: x[x.index.str.endswith('.' + x.Comparison[-1])].mean(),
                                          axis=1)
    quant_data['Condition_All'] = comp.apply(lambda x: x[x.index.str.contains('.', regex=False)].mean(),
                                             axis=1)
    quant_data['Condition2'] = np.where(
        quant_data['pvalue'] > 0.05, quant_data['Condition_All'], quant_data['Condition2'])
    quant_data = quant_data.sort_values(['pTukey', 'pAdjusted', 'pvalue'])

    # #  Mean abundance for each protein among conditions
    quant_data.iloc[:, quant_data.columns.str.contains('.', regex=False)] = np.exp2(
        quant_data.iloc[:, quant_data.columns.str.contains('.', regex=False)])
    #  Mean abundance for each protein
    quant_data['TotalMean'] = quant_data.loc[:,
                                             quant_data.columns.str.contains('.', regex=False)].mean(axis=1)
    #  Log2(FC)
    quant_data['log2(fc)'] = quant_data['Condition2'] - \
        quant_data['Condition1']
    #  -log10(pvalue)
    quant_data[f'-log10({pvalue})'] = -np.log10(quant_data[pvalue])
    quant_data = quant_data.iloc[:, ~
                                 quant_data.columns.isin(['Condition_All'])]
    quant_data = quant_data.reset_index()
    return (quant_data)
