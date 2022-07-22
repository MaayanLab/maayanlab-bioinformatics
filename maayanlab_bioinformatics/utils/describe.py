''' Descriptive statistics on things that aren't pandas data frames.
This can often be a lot more efficient.
'''
import numpy as np

def np_describe(x, axis=0, *, percentiles=[25, 50, 75]):
  ''' Like pandas Series.describe() but operating on numpy arrays / matrices.
  This can be a lot faster especially when working with h5py or sparse data frames.

  :params x: The numpy array to describe
  :params axis: The axis for which to perform describe against
  :returns: dict[str, np.array] A dictionary mapping metric name to results
  '''
  results = {
    'count': (~np.isnan(x)).sum(axis=axis),
    'mean': x.mean(axis=axis),
    'std': x.std(axis=axis),
    'min': x.min(axis=axis),
    'max': x.max(axis=axis),
  }
  if percentiles:
    percentile = np.percentile(x, percentiles, axis=axis)
    results.update({
      f"{p}%": percentile[i]
      for i, p in enumerate(percentiles)
    })
  return results
