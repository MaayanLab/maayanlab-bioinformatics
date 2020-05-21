import logging
import numpy as np
import pandas as pd
from functools import singledispatch


@singledispatch
def cpmNormalize(mat):
  ''' Compute counts-per-million value of counts
  Simple division of each column by the total sum of its counts and multiplying it by 10^6
  '''
  logging.warn('Unrecognized type: ' + type(mat).__name__)
  return cpmNormalize_np(mat)

@cpmNormalize.register
def cpmNormalize_np(mat: np.ndarray):
  return (mat / mat.sum(axis=0)) * 1e6

@cpmNormalize.register
def cpmNormalize_pd(mat: pd.DataFrame):
  return (mat / mat.sum(axis=0)) * 1e6
