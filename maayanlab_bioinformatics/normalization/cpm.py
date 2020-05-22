import logging
import numpy as np
import pandas as pd
from functools import singledispatch


@singledispatch
def cpm_normalize(mat):
  ''' Compute counts-per-million value of counts
  Simple division of each column by the total sum of its counts and multiplying it by 10^6
  '''
  logging.warn('Unrecognized type: ' + type(mat).__name__)
  return cpm_normalize_np(mat)

@cpm_normalize.register
def cpm_normalize_np(mat: np.ndarray):
  return (mat / mat.sum(axis=0)) * 1e6

@cpm_normalize.register
def cpm_normalize_pd(mat: pd.DataFrame):
  return (mat / mat.sum(axis=0)) * 1e6
