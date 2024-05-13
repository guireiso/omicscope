import pandas as pd
import numpy as np
from copy import copy


def quantile_normalization(assay):
    # sort values according to abbundance
    sorted_values = pd.DataFrame(np.sort(assay.values, axis=0))
    # average of protein intesities
    sorted_mean = sorted_values.mean(axis=1)
    # create a dictionary ranking proteins according to abundance
    rank_replace = {}
    rank_replace = {rank + 1: abundance for rank, abundance in sorted_mean.iteritems()}
    # 
    data = assay.rank()
    normalized_data = data.apply(pd.Series.map, arg=rank_replace)
    return normalized_data

def median_normalization(assay):
    # Must be on log-scale
    # consider null as Nan and calculate sample median
    sample_median = assay.copy()
    sample_median = sample_median.replace([0, np.inf, -np.inf], np.nan)
    sample_median = sample_median.median()
    # calculate median_mean
    global_mean = sample_median.mean()
    # calculate global factor
    global_factor = sample_median-global_mean
    # apply normalization factor
    normalized_data = assay - global_factor
    return normalized_data

def average_normalization(assay):
    # Must be on log-scale
    # consider null as Nan and calculate sample mean
    sample_mean = assay.copy()
    sample_mean = sample_mean.replace([0, np.inf, -np.inf], np.nan)
    sample_mean = sample_mean.mean()
    # calculate global average_mean
    global_mean = sample_mean.mean()
    # calculate global factor
    global_factor = sample_mean-global_mean
    # apply normalization factor
    normalized_data = assay - global_factor
    return normalized_data


def normalization(self, expression):
    data = copy(expression)
    if self.logTransform is True:
        data = np.log2(data)
        data = data.replace([0, np.inf, -np.inf], np.nan)
    # apply normalization
    if self.normalization_method is None:
        normalized_data = data
    elif self.normalization_method=='average':
        normalized_data = average_normalization(data)
    elif self.normalization_method=='median':
        normalized_data = median_normalization(data)
    elif self.normalization_method=='quantile':
        normalized_data = quantile_normalization(data)
    else:
        raise ValueError("Please, select None or one of following normalization methods: 'average', 'median', 'quantile'.")
    if self.logTransform is True:
        normalized_data = np.exp2(normalized_data)
        normalized_data = normalized_data.replace(np.nan, 0)
    return normalized_data
