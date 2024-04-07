from copy import copy

import numpy as np
import pandas as pd
from patsy import PatsyError
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests


def Spline_Model_Full(self, pdata, df):
    """
    Natural Cubic Spline
    This algorithm perform the Natural Cubic Spline regression

    Args:
        pdata (list): Time Series.
        expression (list): Expression/Abundance values.
        df (int): Degrees of freedom.
        plot (bool): False

    Returns:
        F_stat (Values): DESCRIPTION.

    """
    TimeCourse = pdata.TimeCourse  # Get timecourse
    dictTC = {'TimeCourse': TimeCourse}
    adjust_variables = list(pdata.iloc[:, pdata.columns.str.startswith('feature')].columns)
    variables = list(pdata.iloc[:, ~pdata.columns.str.startswith('feature')].columns)
    variables.remove('TimeCourse')
    if len(pdata.Condition.drop_duplicates()) == 1:
        variables = []
    if len(variables) > 0:
        variables = '+'.join(variables)
        dictVAR = dict(zip(variables.split('+'), [pdata[x] for x in variables.split('+')]))
        dictionary = dictTC | dictVAR
        sentence = (f"~ {variables} + cr((TimeCourse), df ={df}, constraints='center') + ({variables}):cr((TimeCourse), df ={df}, "
                    "constraints='center')")
    elif len(adjust_variables) > 0:
        adjust_variables = '+'.join(adjust_variables)
        dictionary = dict(zip(adjust_variables.split('+'), [pdata[x] for x in adjust_variables.split('+')]))
        sentence = f"~ {adjust_variables} + cr((TimeCourse), df ={df}, constraints='center')"
    elif len(adjust_variables) == 0:
        variables = 1
        dictionary = dictTC
        sentence = f"~ cr((TimeCourse), df ={df}, constraints='center')"
    transformed_x = dmatrix(sentence,
                            dictionary, return_type='dataframe')
    transformed_x = transformed_x.reset_index(drop=True)
    self.Params['Params']['Stats_Workflow_4'] = "FullModel: "+sentence
    return transformed_x


def Spline_Model_Null(self, pdata, df):
    """
    Natural Cubic Spline
    This algorithm perform the Natural Cubic Spline regression

    Args:
        pdata (list): Time Series.
        expression (list): Expression/Abundance values.
        df (int): Degrees of freedom.
        plot (bool): False

    Returns:
        F_stat (Values): DESCRIPTION.

    """
    TimeCourse = pdata.TimeCourse  # Get timecourse
    dictTC = {'TimeCourse': TimeCourse}
    adjust_variables = list(pdata.iloc[:, pdata.columns.str.startswith('feature')].columns)
    variables = list(pdata.iloc[:, ~pdata.columns.str.startswith('feature')].columns)
    variables.remove('TimeCourse')
    if len(pdata.Condition.drop_duplicates()) == 1:
        variables = []
    if len(variables) > 0:
        variables = '+'.join(variables)
        dictVAR = dict(zip(variables.split('+'), [pdata[x] for x in variables.split('+')]))
        dictionary = dictTC | dictVAR
        sentence = f"~ {variables} + cr((TimeCourse), df ={df}, constraints='center')"
    elif len(adjust_variables) > 0:
        adjust_variables = '+'.join(adjust_variables)
        dictionary = dict(zip(adjust_variables.split('+'), [pdata[x] for x in adjust_variables.split('+')]))
        sentence = f"~ {adjust_variables}"
    elif len(adjust_variables) == 0:
        variables = 1
        dictionary = dictTC
        sentence = "~ 1"
    try:
        transformed_x = dmatrix(sentence,
                                dictionary, return_type='dataframe')
    except PatsyError:
        transformed_x = pd.DataFrame([1] * len(pdata), columns=['Intercept'])
    self.Params['Params']['Stats_Workflow_5'] = "NullModel: "+sentence
    transformed_x = transformed_x.reset_index(drop=True)
    return transformed_x


def rm_zero_cols(matrix):
    matrix = np.asarray(matrix)
    module = abs(matrix)
    matrix = matrix[:, module.sum(axis=0) > 10e-12]
    return matrix


def projMatrix(matrix):
    ginv = np.linalg.pinv(matrix.transpose() @ matrix)
    H = matrix @ ginv @ matrix.transpose()
    return H


def Individuals(expression, pdata):
    expression = np.asmatrix(expression)
    ind = pd.Series(pdata.Individual, dtype="category")
    ind_matrix = dmatrix('~ -1 + individuals',
                         {'individuals': ind})
    Hi = projMatrix(ind_matrix)
    fitInd = Hi @ expression.transpose()
    fitInd = fitInd.transpose()
    expression = expression - fitInd
    return expression, Hi


def Model_Adjustments(full_model, null_model, Individuals):
    expression = Individuals[0]
    Hi = Individuals[1]
    # Full_model
    full_model = np.asmatrix(full_model)
    full_model = full_model - (Hi@full_model)
    full_model = rm_zero_cols(full_model)
    # Null model
    null_model_original = copy(null_model)
    null_model_adj = np.asmatrix(copy(null_model))
    null_model_adj = null_model_adj - (Hi@null_model_adj)
    null_model_adj = rm_zero_cols(null_model_adj)
    if null_model_adj.size == 0:
        null_model = null_model_original
    else:
        null_model = null_model_adj
    return null_model, full_model, expression


def longitudinal(expression, pdata, full_model, null_model):
    individuals = Individuals(expression, pdata)
    model_adjustment = Model_Adjustments(full_model=full_model,
                                         null_model=null_model,
                                         Individuals=individuals)
    return (model_adjustment)


def Spline_Normalization(expression, full_model, null_model):
    import statsmodels.api as sm
    expression = expression.reset_index(drop=True)
    # Full Model
    full = sm.GLM(expression, full_model).fit()
    # Null Model
    null = sm.GLM(expression, null_model).fit()
    F_stat = (null.deviance - full.deviance)/full.deviance
    return F_stat


def Longitudinal_pval(self, assay, pdata, df, ctrl):
    """
    Perform Longitudinal Statistics based on Natural Cubic Spline Regression

    Args:
        assay (pandas dataframe): data containing protein/gene abundance for each
        biological replicate.
        pdata (pandas datafreme/Series): variables associated with time course data.
        df (int): degrees of freedom

    Returns:
        pvalue (pandas series): P-value for each gene.
    """
    ctrl = ctrl
    phenotypedata = copy(pdata)
    expression = copy(assay)
    # Calculating full and null models
    # While data from the same individuals was collected overtime
    if 'Individual' in pdata.columns:
        self.Params['Params']['Stats_Workflow_3'] = 'Related sampling detected'
        pdata2 = pdata.loc[:, pdata.columns != 'Individual']
        pdata2 = pdata2.loc[:, pdata2.columns != 'Biological']
        full_model = Spline_Model_Full(self, pdata=pdata2, df=df)
        full_model_len = len(full_model.columns)
        null_model = Spline_Model_Null(self, pdata=pdata2, df=df)
        longitudinal_workflow = longitudinal(assay, pdata, full_model, null_model)
        full_model = longitudinal_workflow[1]
        null_model = longitudinal_workflow[0]
        assay = longitudinal_workflow[2]
        assay = pd.DataFrame(assay)
        df1 = df
        df2 = len(assay.columns) - full_model_len

    else:
        # While data from the different individuals was collected overtime
        self.Params['Params']['Stats_Workflow_3'] = 'Independent sampling detected'
        pdata = pdata.loc[:, pdata.columns != 'Biological']
        full_model = Spline_Model_Full(self, pdata=pdata, df=df)
        null_model = Spline_Model_Null(self, pdata=pdata, df=df)
        df1 = df
        df2 = len(assay.columns) - len(full_model.columns)
    # Fi for each gene
    Ftest = assay.apply(lambda x:
                        Spline_Normalization(expression=x,
                                             full_model=full_model,
                                             null_model=null_model),
                        axis=1)
    stat = Ftest*df2/df1
    from scipy.stats import f
    pvalue = stat.apply(lambda x: 1-f.cdf(x, df1, df2))
    pvalue = pvalue
    # getting mean expression levels and log2fc
    expression.columns = pd.MultiIndex.from_frame(phenotypedata)
    Mean_Conditions = expression.groupby('Condition', axis=1).mean(numeric_only=True)
    control = Mean_Conditions.loc[:, Mean_Conditions.columns == ctrl].mean(axis=1)
    othergroups = Mean_Conditions.loc[:, Mean_Conditions.columns != ctrl].mean(axis=1)
    log2fc = othergroups.subtract(control)
    if log2fc.isnull().sum() > 0.5*len(log2fc):
        first_TC = pdata.TimeCourse.min()
        last_TC = pdata.TimeCourse.max()
        Mean_Conditions = expression.groupby('TimeCourse', axis=1).mean(numeric_only=True)
        first_TC = Mean_Conditions.loc[:, Mean_Conditions.columns == first_TC].mean(axis=1)
        last_TC = Mean_Conditions.loc[:, Mean_Conditions.columns == last_TC].mean(axis=1)
        log2fc = last_TC.subtract(first_TC)
    return pvalue, Mean_Conditions, log2fc


def Longitudinal_Stats(self, assay, pdata, degrees_of_freedom, pvalue, ctrl, PValue_cutoff):
    quant_data = copy(assay)
    pdata = copy(pdata)
    pdata = pdata.set_index(['Biological', 'Sample'])
    stat = Longitudinal_pval(self, assay=assay, pdata=pdata, df=degrees_of_freedom,
                             ctrl=ctrl)
    pval = stat[0]
    quant_data['pvalue'] = list(pval)
    #  Correcting multilple hypothesis test according to fdr_bh
    pAdjusted = multipletests(pval, alpha=PValue_cutoff,
                              method='fdr_tsbh', is_sorted=False, returnsorted=False)[1]

    quant_data['pAdjusted'] = pAdjusted
    #  Mean abundance for each protein among conditions
    quant_data.loc[:, quant_data.columns.str.contains('.', regex=False)] = np.exp2(
        quant_data.loc[:, quant_data.columns.str.contains('.', regex=False)])
    #  Mean abundance for each protein
    quant_data['TotalMean'] = quant_data.loc[:, quant_data.columns.str.contains('.', regex=False)].mean(axis=1)
    # #  Protein Fold change (Experimental/Control)
    #  Log2(FC)
    quant_data['log2(fc)'] = list(stat[2])
    #  -log10(pvalue)
    quant_data[f'-log10({pvalue})'] = -np.log10(quant_data[pvalue])
    quant_data = quant_data.reset_index()
    return (quant_data)
