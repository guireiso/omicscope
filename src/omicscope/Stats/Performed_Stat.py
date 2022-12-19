"""  Statistic Module comprises all statistical tests performed
by OmicScope for differential expression analysis.

    Returns:
        DataFrame: DataFrame with recommended statistical
        analysis performed
"""

from copy import copy
import numpy as np

def imported_stat(self,statistics):
    """  If user have already imported statistics from other
    softwares, this function 
    1) joins pdata, rdata and assay;
    2) calculates mean of each protein per condition;
    3) calculates the fold-change in relation to ControlGroup;
    4) calculates log2(fc);
    5) calculates -log10(p-value);
    6) If ExcludeKeratins = True, they are excluded.

        Returns:
        DataFrame: Quantitative data
    """
    for i in statistics:
        if i in self.rdata.columns:
            pdata = copy(self.pdata)
            rdata = copy(self.rdata)
            assay = copy(self.assay)

            sample_annotation = pdata.Sample + '.' +pdata.Condition
            sample_dictionary = dict(zip(pdata.Sample,sample_annotation))
            quant_data = assay
            quant_data = quant_data.rename(columns = sample_dictionary)
            quant_data['TotalMean'] = quant_data.mean(axis=1)
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
            quant_data['-log10(pvalue)'] = -np.log10(quant_data['pvalue'])
            if self.ExcludeKeratins is True:
                quant_data = quant_data[~quant_data['Description'].str.contains('Krt|KRT|krt')]
            break
    return(quant_data)