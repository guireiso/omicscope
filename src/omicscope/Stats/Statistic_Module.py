"""  Statistic Module comprises all statistical tests performed
by OmicScope for differential expression analysis.

    Returns:
        DataFrame: DataFrame with recommended statistical
        analysis performed
"""

def perform_stat(self):
    """log2 transformation for expression
    """
    from copy import copy
    import numpy as np
    expression = copy(self.expression)
    log = copy(self.logTransformed)
    rdata = copy(self.rdata)
    pdata = copy(self.pdata)
    # Log-normalize data if it was not
    if log is False:
        expression = expression.replace(0, 0.01)
        expression = np.log2(expression)
    #Perform Longitudinal Analysis or Static Analysis?
    # To perform longitudinal analysis, user must add
    # TimeCourse columns to pdata dataframe
    if True in pdata.columns.str.contains('TimeCourse'):
        pdata = pdata
    else:
        # Apply t-test if len(conditions) == 2
        if len(self.Conditions) == 2:
            from .Static_Statistics import ttest
            # Are variables independent?
            # If true (default) independent t-test is applied.
            # If false related t-test is applied
            
            self.ind_variables = (input('Are variables independent? [Default: True]: ') or 'True')
            params = [self.ind_variables, self.ctrl, self.experimental[0], expression, rdata]
            data = ttest(params=params)
            data = params[4].merge(data, on='Accession')
        # Apply ANOVA if len(conditions) > 2
        elif len(self.Conditions) > 2:     
            from .Static_Statistics import anova
            params = [expression, rdata, self.Conditions]
            data = anova(params=params)
            data = params[1].merge(data, on='Accession')
    data = data.sort_values('pvalue')
    data = data.reset_index(drop=True)
    
    # Filtering Keratin
    data['Description'] = data['Description'].astype(str)
    # Filter Keratin proteins
    data = data[~data['Description'].str.contains('Krt|KRT|krt')]
    data = data.reset_index(drop=True)
    return(data)
