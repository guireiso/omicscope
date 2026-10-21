"""Presence/absence ("on/off") layer of the static workflow.

Entities detected in some conditions but missing (or nearly) in others cannot be tested on intensity,
yet they are often the most interesting biologically. They are tested on detection frequency instead
(exact Fisher / Freeman-Halton test on detected biological replicates, see fast_stats.detection_test),
corrected for multiple testing in their own stratum, labelled in 'Regulation', and given a censored
log2 fold change: absent values are replaced by the limit of detection of each sample (its lowest
detected value), so the reported change is a lower bound (column 'fc_censored').
"""
import numpy as np
import pandas as pd

from . import fast_stats


def _label(row, conditions):
    present = [c for c in conditions if row[c] > 0]
    absent = [c for c in conditions if row[c] == 0]
    if absent and len(present) == 1:
        return f'exclusive to {present[0]}'
    if absent:
        return 'absent in ' + ', '.join(absent)
    return 'detection differs'


def detection_stratum(self, data):
    conditions = list(self.Conditions)
    expr = self.partial_expression.copy()                     # entities x biological samples (linear, 0 = missing)
    counts = self.detection_counts.loc[expr.index, conditions]
    sizes = self.group_sizes.reindex(conditions).to_numpy()
    pvals = fast_stats.detection_test(counts.to_numpy(), sizes)

    # own FDR stratum (the intensity stratum keeps its correction untouched)
    from .Static_Statistics import adjust_pvalues
    message = self.Params['Params'].get('Stats_Workflow_5')
    padj = adjust_pvalues(self, pvals, self.PValue_cutoff)
    if message is not None:
        self.Params['Params']['Stats_Workflow_5'] = message
    self.Params['Params']['Stats_Workflow_7'] = (f'{len(expr)} entities detected only in part of the conditions were tested '
                                                 'on detection frequency (exact test), with their own FDR correction')

    # censored means: missing values sit at the limit of detection of their sample
    log = np.log2(expr.replace(0, np.nan)) if self.logTransform else expr.replace(0, np.nan)
    reference = self.expression.replace(0, np.nan)
    reference = np.log2(reference) if self.logTransform else reference
    lod = pd.concat([reference, log]).min(axis=0)             # per sample
    censored = log.fillna(lod)
    sample_condition = expr.columns.str.split('.').str[-1]
    cmean = pd.DataFrame({c: censored.loc[:, sample_condition == c].mean(axis=1) for c in conditions})
    observed_linear = expr.replace(0, np.nan)

    out = pd.DataFrame(index=expr.index)
    out['pvalue'] = pvals
    out['pAdjusted'] = padj
    out['Test'] = 'detection'
    out['Regulation'] = counts.apply(_label, axis=1, conditions=conditions)
    out['fc_censored'] = True
    out['TotalMean'] = observed_linear.mean(axis=1)
    if len(conditions) == 2:
        ctrl, exp = self.ctrl, self.experimental[0]
        out['mean ' + ctrl] = 2 ** cmean[ctrl] if self.logTransform else cmean[ctrl]
        out['mean ' + exp] = 2 ** cmean[exp] if self.logTransform else cmean[exp]
        out['log2(fc)'] = cmean[exp] - cmean[ctrl]
        out['fc'] = 2 ** out['log2(fc)']
    else:
        out['HighestGroup'], out['LowestGroup'] = cmean.idxmax(axis=1), cmean.idxmin(axis=1)
        out['HighestMean'], out['LowestMean'] = cmean.max(axis=1), cmean.min(axis=1)
        out['log2(fc)'] = out['HighestMean'] - out['LowestMean']
        out['Comparison'] = [[h, lo] for h, lo in zip(out['HighestGroup'], out['LowestGroup'])]
        if self.pvalue == 'pTukey':
            out['pTukey'] = out['pAdjusted']
    out[f'-log10({self.pvalue})'] = -np.log10(out[self.pvalue])
    out = pd.concat([out, observed_linear], axis=1).reset_index().rename(columns={'index': 'Accession'})
    if 'Accession' not in out.columns:
        out = out.rename(columns={out.columns[0]: 'Accession'})
    out = self.partial_rdata.merge(out, on='Accession', how='right')

    data = data.copy()
    data['fc_censored'] = False
    data['Regulation'] = pd.NA
    result = pd.concat([data, out], ignore_index=True)
    for c in conditions:
        result[f'n_detected_{c}'] = result['Accession'].map(self.detection_counts[c]).astype('Int64')
    return result
