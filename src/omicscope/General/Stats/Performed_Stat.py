"""  Statistic Module comprises all statistical tests performed
by OmicScope for differential expression analysis.

    Returns:
        DataFrame: DataFrame with recommended statistical
        analysis performed
"""
import os
from copy import copy

import numpy as np
import pandas as pd


def imported_stat(self, statistics):
    """  If user have already imported statistics from other
    softwares, this function
    1) joins pdata, rdata and assay;
    2) calculates mean of each protein per condition;
    3) calculates the fold-change in relation to ControlGroup;
    4) calculates log2(fc);
    5) calculates -log10(p-value);
    6) If ExcludeContaminants = True, they are excluded.

        Returns:
        DataFrame: Quantitative data
    """
    
    self.Params['Params']['Stats'] = 'Imported from user data'
    for i in statistics:
        if i in self.rdata.columns:
            pdata = copy(self.pdata)
            rdata = copy(self.rdata)
            assay = copy(self.assay)
            statistical_dictionary = {
                'Anova (p)': 'pvalue',
                'pvalue': 'pvalue',
                'p-value': 'pvalue',
                'q-value': 'pAdjusted',
                'q Value': 'pAdjusted',
                'qvalue': 'pAdjusted',
                'padj': 'pAdjusted',
                'p-adjusted': 'pAdjusted',
                'p-adj': 'pAdjusted',
                'pAdjusted': 'pAdjusted'
            }

            sample_annotation = pdata.Sample + '.' + pdata.Condition
            sample_dictionary = dict(zip(pdata.Sample, sample_annotation))
            quant_data = assay
            quant_data = quant_data.rename(columns=sample_dictionary)
            quant_data['TotalMean'] = quant_data.mean(axis=1)
            quant_data = rdata.merge(quant_data, on='Accession')
            quant_data = quant_data.sort_values(i)
            for g in self.Conditions:
                quant_data['mean ' + g] = quant_data.loc[:,
                                                         quant_data.columns.str.endswith(g)].mean(axis=1)
            if len(self.experimental) == 1:
                self.Params['Params']['Stats_Warning1'] = 'Only 1 experimental group identified'
                exp = self.experimental
                quant_data['fc'] = quant_data['mean ' + exp[0]] / quant_data['mean ' + self.ctrl]
                quant_data['log2(fc)'] = np.log2(quant_data['fc'])
                
            elif len(self.experimental) > 1:
                self.Params['Params']['Stats_Warning1'] = 'More than 1 experimental group identified'
                means = quant_data.iloc[:, quant_data.columns.str.startswith('mean ')]
                means = means.replace(0, 0.01)
                fc = means.apply(lambda x: max(x)/min(x), axis='columns')
                quant_data['fc'] = fc
                quant_data['log2(fc)'] = np.log2(quant_data['fc'])
                Comparison = means.apply(lambda x: [x.sort_values().index[-1].split(' ')[-1],
                                                    x.sort_values().index[0].split(' ')[-1]], axis=1)
                quant_data['Comparison'] = Comparison

            quant_data = quant_data.rename(columns=statistical_dictionary)
            quant_data[f'-log10({self.pvalue})'] = -np.log10(quant_data[self.pvalue])
            quant_data['CV'] = quant_data.iloc[:,quant_data.columns.str.contains('.',regex=False)]\
                .apply(lambda x: np.std(x, ddof=1) / np.mean(x), axis=1)

            if self.ExcludeContaminants is True:
                self.Params['Params']['Stats_Warning_2'] = 'Drop protein contaminants based on Frankenfield, 2022'
                path = os.path.dirname(os.path.abspath(__file__))
                contaminants = list(pd.read_csv(path+'/contaminants.csv')['Accession'])
                quant_data['Accession'] = quant_data.Accession.str.replace(',', ';')
                quant_data = quant_data[~quant_data['Accession'].str.split(';').apply(lambda x: any(gene in contaminants for gene in x))]
            break
    return (quant_data)
