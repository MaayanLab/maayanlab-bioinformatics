import logging
import numpy as np
import pandas as pd
from functools import singledispatch


@singledispatch
def log2_normalize(mat, offset=1.):
  ''' Compute log normalization of matrix
  Simple `log2(x + offset)`, offset usually set to 1. because log(0) is undefined.
  '''
  logging.warn('Unrecognized type: ' + type(mat).__name__)
  return log2_normalize_np(mat, offset=offset)

@log2_normalize.register
def log2_normalize_np(mat: np.ndarray, offset=1.):
  return np.log2(mat + offset)

@log2_normalize.register
def log2_normalize_pd(mat: pd.DataFrame, offset=1.):
  return np.log2(mat + offset)


@singledispatch
def log10_normalize(mat, offset=1.):
  ''' Compute log normalization of matrix
  Simple `log10(x + offset)`, offset usually set to 1. because log(0) is undefined.
  '''
  logging.warn('Unrecognized type: ' + type(mat).__name__)
  return log10_normalize_np(mat, offset=offset)

@log10_normalize.register
def log10_normalize_np(mat: np.ndarray, offset=1.):
  return np.log10(mat + offset)

@log10_normalize.register
def log10_normalize_pd(mat: pd.DataFrame, offset=1.):
  return np.log10(mat + offset)
