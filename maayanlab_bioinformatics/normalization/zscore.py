import logging
import numpy as np
import pandas as pd
from functools import singledispatch
from scipy.stats import zscore


@singledispatch
def zscore_normalize(mat, ddof=0):
  ''' Compute the z score of each value in the sample, relative to the sample mean and standard deviation.
  In the case of a pd.DataFrame, preserve the index on the output frame.
  '''
  logging.warn('Unrecognized type: ' + type(mat).__name__)
  return zscore_normalize_np(mat, ddof=ddof)

@zscore_normalize.register
def zscore_normalize_np(mat: np.ndarray, ddof=0):
  return zscore(mat, axis=0, ddof=ddof)

@zscore_normalize.register
def zscore_normalize_pd(mat: pd.DataFrame, ddof=0):
  return pd.DataFrame(zscore_normalize_np(mat, ddof=ddof), index=mat.index, columns=mat.columns)
