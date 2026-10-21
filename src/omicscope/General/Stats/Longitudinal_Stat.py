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


def _rss(Y, X):
    """Residual sum of squares of every row of Y regressed on X (least squares via the pseudo-inverse)."""
    H = X @ np.linalg.pinv(X)
    resid = Y - Y @ H.T
    return np.sum(resid ** 2, axis=1)


def spline_f_test(Y, full, null, individuals=None, moderated=False):
    """Nested-model F test for every row of Y (entities x samples), vectorised per missing-value pattern.

    Missing values (NaN) are left out of each entity's fit. With `moderated`, the residual variances are
    shrunk by empirical Bayes (fast_stats.squeeze_var) and the denominator df grow by the prior df.
    Returns (p-values, F statistics)."""
    from scipy.stats import f as fdist

    from .fast_stats import squeeze_var
    Y = np.asarray(Y, float)
    if individuals is not None:
        full, null = np.hstack([individuals, full]), np.hstack([individuals, null])
    rss_full = np.full(len(Y), np.nan)
    rss_null = np.full(len(Y), np.nan)
    df1 = np.zeros(len(Y))
    df2 = np.zeros(len(Y))
    observed = ~np.isnan(Y)
    patterns = {}
    for i, m in enumerate(observed):
        patterns.setdefault(m.tobytes(), []).append(i)
    for rows in patterns.values():
        m = observed[rows[0]]
        if m.sum() < 3:
            continue
        rows = np.array(rows)
        Xf, Xn = full[m], null[m]
        rf, rn = np.linalg.matrix_rank(Xf), np.linalg.matrix_rank(Xn)
        rss_full[rows] = _rss(Y[np.ix_(rows, m)], Xf)
        rss_null[rows] = _rss(Y[np.ix_(rows, m)], Xn)
        df1[rows], df2[rows] = rf - rn, m.sum() - rf
    ok = (df1 > 0) & (df2 > 0) & np.isfinite(rss_full)
    s2 = np.full(len(Y), np.nan)
    s2[ok] = rss_full[ok] / df2[ok]
    dfd = df2.copy()
    if moderated and ok.sum() > 2:
        post, d0, _ = squeeze_var(s2[ok], df2[ok])
        s2[ok] = post
        dfd[ok] = df2[ok] + d0
    with np.errstate(invalid='ignore', divide='ignore'):
        F = ((rss_null - rss_full) / df1) / s2
        p = np.where(ok, fdist.sf(F, np.maximum(df1, 1), np.minimum(np.maximum(dfd, 1), 1e6)), np.nan)
    return p, F


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
    related = 'Individual' in pdata.columns
    self.Params['Params']['Stats_Workflow_3'] = 'Related sampling detected' if related else 'Independent sampling detected'
    design = pdata.loc[:, ~pdata.columns.isin(['Individual', 'Biological'])]
    full_model = Spline_Model_Full(self, pdata=design, df=df)
    null_model = Spline_Model_Null(self, pdata=design, df=df)
    # Individuals enter both models as fixed effects (equivalent to projecting them out, but it stays
    # exact when an entity has missing samples and the residual df come out right)
    individuals = (dmatrix('~ -1 + C(individuals)', {'individuals': pd.Series(pdata.Individual).astype(str).to_numpy()},
                           return_type='dataframe').to_numpy() if related else None)
    moderated = getattr(self, 'moderated', False)
    pvalue, _ = spline_f_test(assay.to_numpy(float), full_model.to_numpy(float), null_model.to_numpy(float),
                              individuals, moderated)
    pvalue = pd.Series(pvalue, index=assay.index)
    if moderated:
        self.Params['Params']['Stats_Workflow_6'] = 'Moderated (empirical Bayes) F statistics'
    # time effect common to the conditions: condition-only model vs condition + time spline
    conditions = [c for c in design.columns if c != 'TimeCourse' and not c.startswith('feature')]
    if len(design.Condition.drop_duplicates()) > 1 and conditions:
        condition_only = dmatrix('~ ' + '+'.join(conditions), {c: design[c].to_numpy() for c in conditions},
                                 return_type='dataframe').to_numpy(float)
        ptime, _ = spline_f_test(assay.to_numpy(float), null_model.to_numpy(float), condition_only, individuals, moderated)
        self.time_pvalues = pd.Series(ptime, index=assay.index)
    else:
        self.time_pvalues = None
    # getting mean expression levels and log2fc
    expression.columns = pd.MultiIndex.from_frame(phenotypedata)
    Mean_Conditions = expression.T.groupby('Condition').mean(numeric_only=True).T
    control = Mean_Conditions.loc[:, Mean_Conditions.columns == ctrl].mean(axis=1)
    othergroups = Mean_Conditions.loc[:, Mean_Conditions.columns != ctrl].mean(axis=1)
    log2fc = othergroups.subtract(control)
    if log2fc.isnull().sum() > 0.5*len(log2fc):
        first_TC = pdata.TimeCourse.min()
        last_TC = pdata.TimeCourse.max()
        Mean_Conditions = expression.T.groupby('TimeCourse').mean(numeric_only=True).T
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
    from .Static_Statistics import adjust_pvalues
    pval = stat[0].fillna(1)
    quant_data['pvalue'] = list(pval)
    quant_data['pAdjusted'] = adjust_pvalues(self, pval.to_numpy(), PValue_cutoff)
    if getattr(self, 'time_pvalues', None) is not None:
        ptime = self.time_pvalues.fillna(1)
        quant_data['pTime'] = list(ptime)
        message = self.Params['Params'].get('Stats_Workflow_5')
        quant_data['pAdjusted_time'] = adjust_pvalues(self, ptime.to_numpy(), PValue_cutoff)
        self.Params['Params']['Stats_Workflow_5'] = message
        self.Params['Params']['Stats_Workflow_8'] = ('pTime: time effect shared by the conditions '
                                                     '(condition-only model vs condition + time spline)')
    #  Mean abundance for each protein among conditions
    if self.logTransform is True:
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
