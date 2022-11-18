import numpy as np
import pandas as pd
import pandas as pd
import statsmodels.api as sm
from scipy.stats import f
from statsmodels.stats.multitest import multipletests
from copy import copy

def Spline_Model_Full(pdata, df):
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
    from patsy import dmatrix
    TimeCourse = pdata.TimeCourse #Get timecourse
    dictTC = {'TimeCourse':TimeCourse}
    adjust_variables = list(pdata.iloc[:,pdata.columns.str.startswith('feature')].columns)
    variables = list(pdata.iloc[:,~pdata.columns.str.startswith('feature')].columns)
    variables.remove('TimeCourse')
    if len(variables)>0:
        variables = '+'.join(variables)
        dictVAR = dict(zip(variables.split('+'), [pdata[x] for x in variables.split('+')]))
        dictionary = dictTC | dictVAR
        sentence = f"~ {variables} + cr((TimeCourse), df ={df}, constraints='center') + ({variables}):cr((TimeCourse), df ={df}, constraints='center')"
    elif len(adjust_variables) > 0:
        adjust_variables = '+'.join(adjust_variables)
        dictionary = dict(zip(adjust_variables.split('+'), [pdata[x] for x in adjust_variables.split('+')]))
        sentence = f"~ {adjust_variables} + cr((TimeCourse), df ={df}, constraints='center')"
    elif len(adjust_variables) == 0:
        variables = 1
        dictionary = dictTC
        sentence = f"~ cr((TimeCourse), df ={df}, constraints='center')"
    transformed_x = dmatrix(sentence,
                            dictionary,return_type='dataframe')
    transformed_x =  transformed_x.reset_index(drop = True)
    return transformed_x

def Spline_Model_Null(pdata, df):
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
    from patsy import dmatrix
    TimeCourse = pdata.TimeCourse #Get timecourse
    dictTC = {'TimeCourse':TimeCourse}
    adjust_variables = list(pdata.iloc[:,pdata.columns.str.startswith('feature')].columns)
    variables = list(pdata.iloc[:,~pdata.columns.str.startswith('feature')].columns)
    variables.remove('TimeCourse')
    if len(variables)>0:
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
        sentence = f"~ 1"
    try:
        transformed_x = dmatrix(sentence,
                                dictionary,return_type='dataframe')
    except:
        transformed_x = pd.DataFrame([1] * len(pdata), columns = ['Intercept'])
    transformed_x =  transformed_x.reset_index(drop = True)
    return transformed_x

def rm_zero_cols(matrix):
    matrix = np.asarray(matrix)
    module = abs(matrix)
    matrix = matrix[:, module.sum(axis = 0)> 10e-12]
    return matrix

def projMatrix(matrix):
    ginv = np.linalg.pinv(matrix.transpose() @ matrix)
    H = matrix @ ginv @ matrix.transpose()
    return H

def Individuals(expression, pdata):
    expression = np.asmatrix(expression)
    ind = pd.Series(pdata.ind, dtype="category")
    ind_matrix = dmatrix(f'~ -1 + individuals',
                         {'individuals':ind})
    Hi = projMatrix(ind_matrix)
    fitInd = Hi @ expression.transpose()
    fitInd = fitInd.transpose()
    expression = expression - fitInd
    return expression, Hi

def Model_Adjustments(full_model, null_model, Individuals):
    expression = Individuals[0]
    Hi = Individuals[1]
    #Full_model
    full_model = np.asmatrix(full_model)
    full_model = full_model - (Hi@full_model)
    full_model = rm_zero_cols(full_model)
    #Null model    
    null_model = np.asmatrix(null_model)
    null_model = null_model - (Hi@null_model)
    null_model = rm_zero_cols(null_model)
    return null_model, full_model, expression

def longitudinal(expression, pdata, full_model, null_model):
    individuals = Individuals(expression, pdata)
    model_adjustment = Model_Adjustments(full_model = full_model,
                                         null_model = null_model,
                                         Individuals = individuals)
    return(model_adjustment)


def Spline_Normalization(expression, full_model, null_model, plot = False):
    import numpy as np
    import statsmodels.api as sm
    import matplotlib.pyplot as plt
    expression = expression.reset_index(drop = True)
    # Full Model
    full = sm.GLM(expression, full_model).fit()
    # Null Model
    null = sm.GLM(expression, null_model).fit()       
    F_stat = (null.deviance - full.deviance)/full.deviance  
    return F_stat


def Longitudinal_pval(assay, pdata, df):
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
    #Calculating full and null models
    if 'ind' in pdata.columns:
        pdata2 = pdata.loc[:, pdata.columns!='ind']
        full_model = Spline_Model_Full(pdata2, df)
        null_model = Spline_Model_Null(pdata2, df)
    
        full_model = longitudinal(assay, pdata, full_model, null_model)[1]
        null_model = longitudinal(assay, pdata, full_model, null_model)[0]
        assay = longitudinal(assay, pdata, full_model, null_model)[2]
        assay = pd.DataFrame(assay)
        df1 = df
        df2 = len(assay.columns) - len(full_model[0]) - 2
    else:
        full_model = Spline_Model_Full(pdata, df)
        null_model = Spline_Model_Null(pdata, df)
        df1 = df
        df2 = len(assay.columns) - len(full_model.columns)
    #Fi for each gene
    Ftest = assay.apply(lambda x:
                        Spline_Normalization(expression = x,
                                             full_model = full_model,
                                             null_model = null_model),
                        axis = 1)

    stat = Ftest*df2/df1
    from scipy.stats import f
    pvalue = stat.apply(lambda x: 1-f.cdf(x, df1, df2))
    pvalue = pvalue
    return pvalue

def Longitudinal_Stats(assay, pdata, df):
    quant_data = copy(assay)
    pval = Longitudinal_pval(assay = assay, pdata = pdata, df = df)
    #  Correcting multilple hypothesis test according to fdr_bh
    pAdjusted = multipletests(pval, alpha=0.1,
                                          method='fdr_tsbh', is_sorted=False, returnsorted=False)[1]
    quant_data['pvalue'] = pval
    quant_data['pAdjusted'] = pAdjusted
    return quant_data
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
    quant_data['TotalMean'] = quant_data.loc[:, quant_data.columns.str.contains('.')].mean(axis=1)
    #  Protein Fold change (Experimental/Control)
    quant_data['fc'] = quant_data['mean ' + Experimental] / quant_data['mean ' + ControlGroup]
    #  Log2(FC)
    quant_data['log2(fc)'] = np.log2(quant_data['fc'])
    #  -log10(pvalue)
    quant_data['-log10(p)'] = -np.log10(quant_data['pvalue'])
    quant_data = quant_data.reset_index()
    return(quant_data)