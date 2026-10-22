"""  Statistic Module comprises all statistical tests performed
by OmicScope for differential expression analysis.

    Returns:
        DataFrame: DataFrame with recommended statistical
        analysis performed
"""
import os

import pandas as pd


def choose_variance_model(self, expression, ratio_cutoff=2.0, min_replicates=4):
    """Decide once, for the whole dataset, between pooled-variance (optionally moderated) and Welch tests.

    welch='auto' uses Welch/Games-Howell when the median per-entity variance differs more than
    `ratio_cutoff`-fold between conditions and every condition has >= `min_replicates` biological
    replicates (with fewer, Welch has almost no power and the moderated pooled test is kept)."""
    import numpy as np
    labels = expression.columns.str.split('.').str[-1]
    sizes = {c: int((labels == c).sum()) for c in self.Conditions}
    medians = {c: float(np.nanmedian(expression.loc[:, labels == c].var(axis=1, ddof=1))) for c in self.Conditions}
    ratio = max(medians.values()) / max(min(medians.values()), 1e-12)
    self.variance_ratio = ratio
    welch = self.welch
    if welch == 'auto':
        heteroscedastic = ratio > ratio_cutoff
        welch = heteroscedastic and min(sizes.values()) >= min_replicates
        note = (f'Variance ratio between conditions = {ratio:.2f} (median per-entity variance, log scale); '
                + ('Welch tests selected' if welch else
                   'heteroscedastic, but fewer than 4 replicates in a condition: pooled moderated test kept'
                   if heteroscedastic else 'pooled-variance tests selected'))
        self.Params['Params']['Stats_Variance_Model'] = note
    self._use_welch = bool(welch)
    self._use_moderated = bool(self.moderated) and not self._use_welch


def condition_last_sample_names(self):
    """Make every biological sample name end with '.<Condition>'.

    The static tests find the samples of each condition by that suffix. expression_function joins every
    pdata column into the name in pdata order (e.g. 'BioRep_1.Ctrl-b1' when pdata has a Batch column),
    so extra columns are moved before the dot: 'BioRep_1-b1.Ctrl'. pdata, expression, partial_expression
    and excluded_values are renamed together."""
    pdata = self.pdata
    if all(str(s).endswith('.' + str(c)) for s, c in zip(pdata['Sample'], pdata['Condition'])):
        return
    extra = sorted(c for c in pdata.columns if c not in ('Sample', 'Condition', 'Biological', 'technical'))
    new = pdata.apply(lambda r: 'BioRep_' + str(r['Biological'])
                      + ''.join('-' + str(r[c]).replace('.', '_') for c in extra) + '.' + str(r['Condition']), axis=1)
    names = dict(zip(pdata['Sample'], new))
    self.pdata = pdata.assign(Sample=new)
    self.expression = self.expression.rename(columns=names)
    for attr in ('partial_expression', 'excluded_values'):
        if getattr(self, attr, None) is not None:
            setattr(self, attr, getattr(self, attr).rename(columns=names))


def correct_batch(self, expression):
    """Remove batch effects from the (log-scale) expression matrix before testing (batch_correction=True).

    limma-style removeBatchEffect (Quality.remove_batch_effect): each entity is fitted as
    y ~ condition + batch and the batch terms are subtracted. The tests then run on the corrected values
    without counting the degrees of freedom spent on the batch terms, which slightly overstates them (p-values
    a little optimistic with few samples per batch). Designs in which batch is confounded with condition
    cannot be corrected and raise a ValueError. Missing values (0 or NaN) stay missing."""
    import numpy as np

    from ..Quality import remove_batch_effect
    column = getattr(self, 'batch_column', 'Batch')
    if column not in self.pdata.columns:
        raise ValueError(f"batch_correction=True needs a '{column}' column in pdata with the batch of each sample.")
    info = self.pdata.drop_duplicates('Sample').set_index('Sample').reindex(expression.columns)
    missing = expression.isna() | (expression == 0) if not self.logTransform else expression.isna()
    corrected = remove_batch_effect(expression.mask(missing), info['Condition'].to_numpy(), info[column].to_numpy())
    if not self.logTransform:
        corrected = corrected.mask(missing, expression)
    n_batches = info[column].nunique()
    self.Params['Params']['Stats_Batch_Correction'] = (
        f"Batch effects ('{column}', {n_batches} batches) removed before testing by fitting each entity as "
        "condition + batch and subtracting the batch terms (limma removeBatchEffect)")
    return corrected


def perform_static_stat(self):
    """log2 transformation for expression, optional batch correction (self.batch_correction, see
    correct_batch: it slightly overstates the degrees of freedom of the tests, and confounded designs
    cannot be corrected) and differential expression tests.
    """
    from copy import copy

    import numpy as np
    self.Params['Params']['Stats_Workflow_1'] = 'OmicScope performed Static Workflow'
    condition_last_sample_names(self)
    expression = copy(self.expression)
    log = copy(self.logTransform)
    rdata = copy(self.rdata)
    pvalue = copy(self.pvalue)
    PValue_cutoff = copy(self.PValue_cutoff)
    # Log-normalize data
    if log is True:
        expression = expression.replace(0, np.nan)
        expression = np.log2(expression)
        self.Params['Params']['Stats_Workflow_2'] = 'OmicScope log2-transformed protein abundances'
    if getattr(self, 'batch_correction', False):
        expression = correct_batch(self, expression)
    choose_variance_model(self, expression)
    # Apply t-test if len(conditions) == 2
    if len(self.Conditions) == 2:
        from .Static_Statistics import ttest
        params = [self.ind_variables, self.ctrl, self.experimental[0], expression, rdata,
                  pvalue, PValue_cutoff]
        data = ttest(self, params=params)
        data = params[4].merge(data, on='Accession')

    # Apply ANOVA if len(conditions) > 2
    elif len(self.Conditions) > 2:
        from .Static_Statistics import anova
        params = [expression, rdata, self.Conditions, pvalue, PValue_cutoff]
        data = anova(self, params=params)
        data = params[1].merge(data, on='Accession')
    data['Test'] = 'intensity'
    if getattr(self, 'partial_expression', None) is not None:
        from .detection import detection_stratum
        data = detection_stratum(self, data)
    data = data.sort_values('pvalue')
    data = data.reset_index(drop=True)
    # # Calculating CV
    data['CV'] = data.iloc[:,data.columns.str.contains('.',regex=False)]\
                .apply(lambda x: np.std(x, ddof=1) / np.mean(x), axis=1)
    # Filtering contaminants
    if self.ExcludeContaminants is True:
        self.Params['Params']['Stats_Warning'] = 'Drop protein contaminants based on Frankenfield, 2022'
        path = os.path.dirname(os.path.abspath(__file__))
        contaminants = list(pd.read_csv(path+'/contaminants.csv')['Accession'])
        data['Accession'] = data.Accession.str.replace(',', ';')
        data = data[~data['Accession'].str.split(';').apply(lambda x: any(gene in contaminants for gene in x))]
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
    expression = expression.replace(np.nan, 0)
    log = copy(self.logTransform)
    rdata = copy(self.rdata)
    pdata = copy(self.pdata)
    pvalue = copy(self.pvalue)
    ctrl = copy(self.ctrl)
    PValue_cutoff = copy(self.PValue_cutoff)
    self.Params['Params']['Stats_Workflow_1'] = 'OmicScope performed Longitudinal Workflow'
    # Log-normalize data
    if log is True:
        expression = expression.replace(0, np.nan)
        expression = np.log2(expression)
        self.Params['Params']['Stats_Workflow_2'] = 'OmicScope log2-transformed protein abundances'
    from .Longitudinal_Stat import Longitudinal_Stats
    data = Longitudinal_Stats(self, assay=expression, pdata=pdata.drop(columns=['technical']),
                              degrees_of_freedom=degrees_of_freedom, pvalue=pvalue, ctrl=ctrl,
                              PValue_cutoff=PValue_cutoff)
    data = rdata.merge(data, on='Accession')
    # longitudinal modules
    data = data.sort_values('pvalue')
    data = data.reset_index(drop=True)
    # # Calculating CV
    data['CV'] = data.iloc[:,data.columns.str.contains('.',regex=False)]\
                .apply(lambda x: np.std(x, ddof=1) / np.mean(x), axis=1)
    # # Filtering Contaminants
    if self.ExcludeContaminants is True:
        self.Params['Params']['Stats_Warning'] = 'Drop protein contaminants based on Frankenfield, 2022'
        path = os.path.dirname(os.path.abspath(__file__))
        contaminants = list(pd.read_csv(path+'/contaminants.csv')['Accession'])
        data['Accession'] = data.Accession.str.replace(',', ';')
        data = data[~data['Accession'].str.split(';').apply(lambda x: any(gene in contaminants for gene in x))]
    data = data.reset_index(drop=True)
    data = copy(data)
    data.iloc[:, data.columns.str.contains(
        '.', regex=False)] = data.iloc[:, data.columns.str.contains('.', regex=False)].replace(np.nan, 0)
    return (data)
