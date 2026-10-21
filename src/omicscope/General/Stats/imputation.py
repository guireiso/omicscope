import sklearn.impute
import pandas as pd
import numpy as np
from copy import copy

def value_imputation(self, expression):
    data = copy(expression)
    if self.logTransform is True:
        data = np.log2(data)
        data = data.replace([0, np.inf, -np.inf], np.nan)
    # remove proteins presenting only NaN
    data = data.T
    if self.imputation_method is None:
        imputated_values = data
    elif self.imputation_method == 'mean':
        imp = sklearn.impute.SimpleImputer(missing_values=np.nan, strategy="mean")
        imputated_values = pd.DataFrame(imp.fit_transform(data.values))
    elif self.imputation_method == 'median':
        imp = sklearn.impute.SimpleImputer(missing_values=np.nan, strategy="median")
        imputated_values = pd.DataFrame(imp.fit_transform(data.values))
    elif self.imputation_method == 'minprob':
        # Perseus-like left-censored imputation: draws from a narrow normal below each sample's
        # detected distribution (down-shift 1.8 sd, width 0.3 sd); reproducible seed
        rng = np.random.default_rng(getattr(self, 'imputation_seed', 0))
        imputated_values = data.copy()
        for sample in data.index:
            row = data.loc[sample]
            mu, sd = row.mean(skipna=True), row.std(skipna=True)
            missing = row.isna()
            imputated_values.loc[sample, missing] = rng.normal(mu - 1.8 * sd, 0.3 * sd, int(missing.sum()))
        self.imputed_mask = data.isna().T
        self.Params['Params']['DataHandling_Imputation'] = 'MinProb left-censored imputation (shift 1.8 sd, width 0.3 sd)'
    elif self.imputation_method == 'knn':
        imp = sklearn.impute.KNNImputer(n_neighbors=2, weights="uniform")
        imputated_values = pd.DataFrame(imp.fit_transform(data.values))
    else:
        raise ValueError("Please, select None or one of following imputation methods: 'mean', 'median', 'knn', 'minprob'.")
    imputated_values.index = data.index
    imputated_values.columns = data.columns
    imputated_values = imputated_values.T
    if self.logTransform is True:
        imputated_values = np.exp2(imputated_values)
        imputated_values = imputated_values.replace(np.nan, 0)
    return imputated_values