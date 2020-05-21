import logging
import numpy as np
import pandas as pd
from functools import singledispatch


@singledispatch
def logNormalize(mat, offset=1.):
  ''' Compute log normalization of matrix
  Simple `log10(x + offset)`, offset usually set to 1. because log(0) is undefined.
  '''
  logging.warn('Unrecognized type: ' + type(mat).__name__)
  return logNormalize_np(mat, offset=offset)

@logNormalize.register
def logNormalize_np(mat: np.ndarray, offset=1.):
  return np.log10(mat + offset)

@logNormalize.register
def logNormalize_pd(mat: pd.DataFrame, offset=1.):
  return np.log10(mat + offset)
