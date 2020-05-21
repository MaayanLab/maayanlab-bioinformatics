import logging
import numpy as np
import pandas as pd
from functools import singledispatch


@singledispatch
def log2Normalize(mat, offset=1.):
  ''' Compute log normalization of matrix
  Simple `log2(x + offset)`, offset usually set to 1. because log(0) is undefined.
  '''
  logging.warn('Unrecognized type: ' + type(mat).__name__)
  return log2Normalize_np(mat, offset=offset)

@log2Normalize.register
def log2Normalize_np(mat: np.ndarray, offset=1.):
  return np.log2(mat + offset)

@log2Normalize.register
def log2Normalize_pd(mat: pd.DataFrame, offset=1.):
  return np.log2(mat + offset)
