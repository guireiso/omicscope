
from copy import copy

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from . import fast_stats

FDR_NAMES = {'fdr_tsbh': 'two-stage Benjamini-Hochberg', 'fdr_bh': 'Benjamini-Hochberg', 'qvalue': 'Storey q-value'}


def adjust_pvalues(self, pvalues, alpha):
    """Multiple-testing correction chosen by the user (self.fdr)."""
    method = getattr(self, 'fdr', 'fdr_tsbh')
    self.Params['Params']['Stats_Workflow_5'] = f'OmicScope performed {FDR_NAMES[method]} correction'
    if method == 'qvalue':
        q, pi0 = fast_stats.qvalue(pvalues)
        self.Params['Params']['Stats_Workflow_5'] += f' (estimated proportion of null entities = {pi0:.2f})'
        return q
    return multipletests(pvalues, alpha=alpha, method=method)[1]


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
        welch, moderated = getattr(self, '_use_welch', False), getattr(self, '_use_moderated', False)
        _, p, _ = fast_stats.ttest(quant_data.iloc[:, quant_data.columns.str.endswith(Experimental)].to_numpy(float),
                                   quant_data.iloc[:, quant_data.columns.str.endswith(ControlGroup)].to_numpy(float),
                                   equal_var=not welch, moderated=moderated)
        result = (None, p)
        self.Params['Params']['Stats_Workflow_3'] = ('OmicScope performed moderated (empirical Bayes) T-test' if moderated
                                                     else 'OmicScope performed Welch T-test' if welch
                                                     else 'OmicScope performed Independent T-test')

        #  Perform a sample-related t-test
    elif ind_variables is False:
        from scipy.stats import ttest_rel
        result = ttest_rel(copy(quant_data.iloc[:, quant_data.columns.str.endswith(Experimental)]),
                           copy(
                               quant_data.iloc[:, quant_data.columns.str.endswith(ControlGroup)]),
                           axis=1, nan_policy='omit')
        self.Params['Params']['Stats_Workflow_3'] = 'OmicScope performed Paired T-test'

    quant_data['pvalue'] = result[1]
    quant_data['pvalue'] = quant_data['pvalue'].replace(np.nan, 1)
    self.Params['Params']['Stats_Workflow_4'] = 'OmicScope replaceD p-value = NaN to 1.0'
    # Effect size on the tested (log) scale: mean difference with its 95% CI and Hedges' g
    X1 = quant_data.iloc[:, quant_data.columns.str.endswith(Experimental)].to_numpy(float)
    X2 = quant_data.iloc[:, quant_data.columns.str.endswith(ControlGroup)].to_numpy(float)
    diff, low, high, g, se = fast_stats.effect_sizes(
        X1, X2, equal_var=not getattr(self, '_use_welch', False), moderated=getattr(self, '_use_moderated', False),
        paired=ind_variables is False)
    # added after the per-condition means: sample columns are found by suffix
    effects = pd.DataFrame({'Mean difference': diff, 'CI95 low': low, 'CI95 high': high, "Hedges' g": g, 'SE': se},
                           index=quant_data.index)
    quant_data = quant_data.sort_values('pvalue')

    #  Correcting multiple hypothesis test
    quant_data['pAdjusted'] = adjust_pvalues(self, quant_data['pvalue'].to_numpy(), params[6])

    #  Mean abundance for each protein among conditions
    if self.logTransform is True:
        quant_data.loc[:, quant_data.columns.str.endswith(ControlGroup)] = np.exp2(
            copy(quant_data.loc[:, quant_data.columns.str.endswith(ControlGroup)]))
    quant_data['mean ' + ControlGroup] = quant_data.loc[:,
                                                        quant_data.columns.str.endswith(ControlGroup)].mean(axis=1)
    if self.logTransform is True:
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
    quant_data = quant_data.join(effects)
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


def Tukey_hsd(quant_data, cutoff=0.05, method='tukey'):
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
    df = quant_data[quant_data['pAdjusted'] <= cutoff]
    df = df.iloc[:, df.columns.str.contains('.', regex=False)]
    labels = df.columns.str.split('.').str[-1].to_numpy()
    pairs, p = fast_stats.pairwise_posthoc(df.to_numpy(float), labels, sorted(set(labels)), method=method)
    comparison = [f'{a}-vs-{b}' for a, b in pairs]
    return pd.DataFrame({'pTukey': [list(row) for row in p], 'Comparison': [comparison] * len(df)}, index=df.index)


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

    welch, moderated = getattr(self, '_use_welch', False), getattr(self, '_use_moderated', False)
    samples = quant_data.columns[quant_data.columns.str.contains('.', regex=False)]
    labels = np.array([next((c for c in Conditions if s.endswith(c)), None) for s in samples], dtype=object)
    _, stat, _, _ = fast_stats.anova(quant_data[samples].to_numpy(float), labels, list(Conditions),
                                     welch=welch, moderated=moderated)
    if welch:
        self.Params['Params']['Stats_Workflow_3'] = 'OmicScope performed Welch ANOVA'
    elif moderated:
        self.Params['Params']['Stats_Workflow_3'] = 'OmicScope performed moderated (empirical Bayes) ANOVA'
    # Multiple hypothesis correction
    stat = pd.Series(stat)
    stat.index = quant_data.index
    quant_data['pvalue'] = stat
    eta2 = pd.Series(fast_stats.eta_squared(quant_data[samples].to_numpy(float), labels, list(Conditions)),
                     index=quant_data.index, name='Eta squared')   # share of variance explained by the conditions
    quant_data['pvalue'] = quant_data['pvalue'].replace(np.nan, 1)
    self.Params['Params']['Stats_Workflow_4'] = 'OmicScope replaceD p-value = NaN to 1.0'
    quant_data = quant_data.sort_values('pvalue')
    # Multiple-testing correction
    quant_data['pAdjusted'] = adjust_pvalues(self, quant_data['pvalue'].to_numpy(), params[4])
    # Post-hoc test on the entities significant in the ANOVA (pAdjusted <= cutoff)
    if pvalue == 'pTukey':
        method = 'games-howell' if welch else 'tukey'
        Tukey = Tukey_hsd(quant_data, cutoff=params[4], method=method)
        self.Params['Params']['Stats_Workflow_6'] = (f"OmicScope performed {'Games-Howell' if welch else 'Tukey HSD'} "
                                                     f"post-hoc test on entities with pAdjusted <= {params[4]}")
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
        if self.logTransform is True:
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
        quant_data['Comparison'] = quant_data.Comparison.str.split('-vs-') 
    else:
        quant_data = quant_data.sort_values('pAdjusted') 
        MeanAbundance = {}
        for x in Conditions:
             MeanAbundance[x] = quant_data.iloc[:,quant_data.columns.str.contains('.'+x, regex=False)].mean(axis=1)
        MeanAbundance = pd.DataFrame(MeanAbundance)

        HighestGroup = MeanAbundance.apply(lambda x: x.idxmax(), axis=1)
        LowestGroup = MeanAbundance.apply(lambda x: x.idxmin(), axis=1)
        MeanAbundance['HighestMean'] = MeanAbundance.max(axis=1, numeric_only=True)
        MeanAbundance['LowestMean'] = MeanAbundance.min(axis=1, numeric_only=True)
        MeanAbundance['HighestGroup'] = HighestGroup
        MeanAbundance['LowestGroup'] = LowestGroup
        quant_data = quant_data.merge(MeanAbundance, how='left', left_index=True, right_index=True)
        quant_data['log2(fc)'] = MeanAbundance['HighestMean']-MeanAbundance['LowestMean']
        quant_data['Comparison'] = MeanAbundance['HighestGroup']+'-vs-'+MeanAbundance['LowestGroup']
        quant_data['Comparison'] = quant_data.Comparison.str.split('-vs-')
        # #  Mean abundance for each protein among conditions
        if self.logTransform is True:
            quant_data.iloc[:, quant_data.columns.str.contains('.', regex=False)] = np.exp2(
                quant_data.iloc[:, quant_data.columns.str.contains('.', regex=False)])
        #  Mean abundance for each protein
        quant_data['TotalMean'] = quant_data.loc[:,
                                                quant_data.columns.str.contains('.', regex=False)].mean(axis=1)
        #  -log10(pvalue)
        quant_data[f'-log10({pvalue})'] = -np.log10(quant_data[pvalue])
        quant_data = quant_data.iloc[:, ~
                                    quant_data.columns.isin(['Condition_All'])]
        
    quant_data = quant_data.join(eta2)
    quant_data = quant_data.reset_index()
    return (quant_data)
