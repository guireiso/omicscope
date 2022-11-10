"""  Statistic Module comprises all statistical tests performed
by OmicScope for differential expression analysis.

    Returns:
        DataFrame: DataFrame with recommended statistical
        analysis performed
"""

from copy import copy

import numpy as np
import pandas as pd


def imported_stat(self):
    """  If user have already imported statistics from other
    softwares, this function 1) joins pdata, rdata and assay;
    2) calculates mean of each protein abundance;
    3) calculates mean of each protein per condition;
    4) calculates the fold-change in relation to ControlGroup;
    5) calculates log2(fc);
    6) calculates -log10(p-value);
    7) If ExcludeKeratins = True, they are excluded.

        Returns:
        DataFrame: Quantitative data
    """
    for i in ['Anova (p)','Anova', 'pvalue', 'p-value']:
        if i in self.rdata.columns:
            quant_data = copy(self.expression)
            quant_data['TotalMean'] = quant_data.mean(axis=1)
            rdata = copy(self.rdata)
            quant_data = rdata.merge(quant_data, on='Accession')
            quant_data = quant_data.sort_values(i)
            for g in self.Conditions:
                quant_data['mean ' + g] = quant_data.loc[:,
                                                         quant_data.columns.str.endswith(g)].mean(axis=1)
            if len(self.experimental) == 1:
                exp = self.experimental
                quant_data['fc'] = quant_data['mean ' + exp[0]] / quant_data['mean ' + self.ctrl]
                quant_data['log2(fc)'] = np.log2(quant_data['fc'])
            quant_data = quant_data.rename(columns={i: 'pvalue'})
            quant_data['-log10(p)'] = -np.log10(quant_data['pvalue'])
            if self.ExcludeKeratins is True:
                quant_data = quant_data[~quant_data['Description'].str.contains('Krt|KRT|krt')]
            break
    return(quant_data)


def ttest(params):
    """Perform a T-test on the data

    Args:
        params (list): list containing control and experimental groups,
        normalized dataframe and if the variables are dependent or related.
    """
    from statsmodels.stats.multitest import multipletests
    ind_variables = params[0]
    ControlGroup = params[1]
    Experimental = params[2]
    quant_data = params[3]

    if ind_variables == 'True':
        from scipy.stats import ttest_ind

        #  Perform an independent t-test
        result = ttest_ind(quant_data.iloc[:, quant_data.columns.str.endswith(Experimental)],
                           quant_data.iloc[:, quant_data.columns.str.endswith(ControlGroup)], axis=1)
        quant_data['pvalue'] = result[1]
        quant_data = quant_data.sort_values('pvalue')
        print('Independent T-test was carried out!')
        #  Perform a sample-related t-test
    else:
        from scipy.stats import ttest_rel
        result = ttest_rel(quant_data.iloc[:, quant_data.columns.str.endswith(Experimental)],
                           quant_data.iloc[:, quant_data.columns.str.endswith(ControlGroup)], axis=1)
        quant_data['pvalue'] = result[1]
        quant_data = quant_data.sort_values('pvalue')
        print('T-test of related variables was carried out!')

    #  Correcting multilple hypothesis test according to fdr_bh
    quant_data['q Value'] = multipletests(quant_data['pvalue'], alpha=0.1,
                                          method='fdr_tsbh', is_sorted=False, returnsorted=False)[1]

    #  Mean abundance for each protein among conditions
    quant_data.loc[:, quant_data.columns.str.endswith(ControlGroup)] = np.exp2(
        quant_data.loc[:, quant_data.columns.str.endswith(ControlGroup)])
    quant_data['mean ' + ControlGroup] = quant_data.loc[:,
                                                        quant_data.columns.str.endswith(ControlGroup)].mean(axis=1)

    quant_data.loc[:, quant_data.columns.str.endswith(Experimental)] = np.exp2(
        quant_data.loc[:, quant_data.columns.str.endswith(Experimental)])
    quant_data['mean ' + Experimental] = quant_data.loc[:,
                                                        quant_data.columns.str.endswith(Experimental)].mean(axis=1)
    #  Mean abundance for each protein
    quant_data['TotalMean'] = quant_data.loc[:, quant_data.columns.str.contains('\.')].mean(axis=1)
    #  Protein Fold change (Experimental/Control)
    quant_data['fc'] = quant_data['mean ' + Experimental] / quant_data['mean ' + ControlGroup]
    #  Log2(FC)
    quant_data['log2(fc)'] = np.log2(quant_data['fc'])
    #  -log10(pvalue)
    quant_data['-log10(p)'] = -np.log10(quant_data['pvalue'])
    quant_data = quant_data.reset_index()
    return(quant_data)


def tukey_correction(df, i):
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
    df.columns = df.columns.str.split('\.').str[-1]
    df = pd.DataFrame(df.iloc[i, :])
    df = df.replace(0, 0.01)
    df = df.reset_index()
    df.columns = ['groups', 'abundance']
    df = pairwise_tukeyhsd(endog=df['abundance'],
                           groups=df['groups'],
                           alpha=0.05)
    df = pd.DataFrame(data=df._results_table.data[1:], columns=df._results_table.data[0])
    df['Comparison'] = df.group1 + '-' + df.group2
    padj = df['p-adj']
    comparison = df.Comparison
    return list(padj), list(comparison)


def Tukey_hsd(df):
    """Perform Tukey HSD test for each differentially regulated entity
    found according ANOVA's test.

    Args:
        df (DataFrame): Normalized dataframe.

    Returns:
       padj [list]: Adjusted p-value for each differentially regulated entity
       in each comparison. Entities that were not differentially regulated
       returned 2
       comparison [list]: Comparison for which p-value was evaluated. Entities 
       that were not differentially regulated returned 2
    """
    df_original = copy(df)
    df = copy(df)
    df = df[df['pvalue'] < 0.05]
    df = df.iloc[:, df.columns.str.contains('\.')]
    padj = []
    comparison = []
    # Perform Tukey's post-hoc correction for each of differentially
    # regulated entity according ANOVA's test
    for i in range(0, len(df)):
        result = tukey_correction(df, i)
        padj.append(result[0])
        comparison.append(result[1])
    padj.extend([2] * (len(df_original) - len(padj)))
    comparison.extend([2] * (len(df_original) - len(comparison)))
    return padj, comparison


def anova(params):
    """Test for ANOVA for parametric and multiple comparisons

    Args:
        params ([type]): [description]
    """
    from scipy.stats import f_oneway
    quant_data = params[0]
    Conditions = params[2]

    # Perform ANOVA test according to SCI-PY algorithm
    expression = []
    for i in Conditions:
        df = quant_data.iloc[:, quant_data.columns.str.endswith(i)]
        expression.append(df)
    result = f_oneway(*expression, axis=1)

    # Perform Tukey's Post-hoc correction
    # TODO: Needs major comparisons with edgeR
    quant_data['pvalue'] = result[1]
    quant_data = quant_data.sort_values('pvalue')
    quant_data['padj'], quant_data['Comparison'] = Tukey_hsd(quant_data)
    quant_data = quant_data.explode(['padj', 'Comparison'])
    quant_data['padj'] = quant_data.padj.astype(float)
    quant_data = quant_data.sort_values('padj')
    # Replace '2' with repeated values from ANOVA for
    # non diffenrentially regulated entities.
    quant_data['padj'] = np.where(quant_data['padj'] == 2, quant_data['pvalue'],
                                  quant_data['padj'])
    # Replace '2' with right comparisons performed
    quant_data['Comparison'] = np.where(
        quant_data['Comparison'] == 2, '-'.join(Conditions),
        quant_data['Comparison'])
    comp = quant_data
    comp['Comparison'] = comp.Comparison.str.split('-')
    comparisons = comp.Comparison
    comp = comp.iloc[:, comp.columns.str.contains('\.')]
    quant_data['Condition1'] = [comp.iloc[z, comp.columns.str.endswith(i[0])].mean()
                                for z, i in enumerate(comparisons)]
    quant_data['Condition2'] = [comp.iloc[z, comp.columns.str.endswith(i[-1])].mean()
                                for z, i in enumerate(comparisons)]
    quant_data['Condition_All'] = [comp.iloc[z, ~comp.columns.str.endswith(i[0])].mean()
                                   for z, i in enumerate(comparisons)]
    quant_data['Condition2'] = np.where(
        quant_data['pvalue'] > 0.05, quant_data['Condition_All'], quant_data['Condition2'])
    quant_data = quant_data.sort_values('pvalue')
    print('Anova test was performed!')

    #  Mean abundance for each protein among conditions
    quant_data.iloc[:, quant_data.columns.str.contains('\.')] = np.exp2(
        quant_data.iloc[:, quant_data.columns.str.contains('\.')])
    #  Mean abundance for each protein
    quant_data['TotalMean'] = quant_data.loc[:, quant_data.columns.str.contains('\.')].mean(axis=1)
    #  Log2(FC)
    quant_data['log2(fc)'] = quant_data['Condition2'] - quant_data['Condition1']
    quant_data['pvalue'] = quant_data['padj']
    #  -log10(pvalue)
    quant_data['-log10(p)'] = -np.log10(quant_data['pvalue'])
    quant_data = quant_data.iloc[:, ~quant_data.columns.isin(['Condition_All', 'padj'])]
    quant_data = quant_data.reset_index()
    return(quant_data)


def perform_stat(self):
    """log2 transformation for expression
    """
    from copy import copy

    import numpy as np
    expression = copy(self.expression)
    log = copy(self.logTransformed)
    rdata = copy(self.rdata)
    # Log-normalize data if it was not
    if log is False:
        expression = expression.replace(0, 0.01)
        expression = np.log2(expression)
    # Apply t-test if len(conditions) == 2
    if len(self.Conditions) == 2:
        # Are variables independent?
        # If true (default) independent t-test is applied.
        # If false related t-test is applied
        self.ind_variables = (input('Are variables independent? [Default: True]: ') or 'True')
        params = [self.ind_variables, self.ctrl, self.experimental[0], expression, rdata]
        data = ttest(params=params)
        data = params[4].merge(data, on='Accession')
    # Apply ANOVA if len(conditions) > 2
    elif len(self.Conditions) > 2:
        params = [expression, rdata, self.Conditions]
        data = anova(params=params)
        data = params[1].merge(data, on='Accession')
    data = data.sort_values('pvalue')
    data = data.reset_index(drop=True)
    # Filter Unique peptides if it was set.
    try:
        data = data[data['Unique peptides'] >= self.uniquepeptides]
    except Exception:
        pass
    data['Description'] = data['Description'].astype(str)
    # Filter Keratin proteins
    data = data[~data['Description'].str.contains('Krt|KRT|krt')]
    data = data.reset_index(drop=True)
    return(data)
