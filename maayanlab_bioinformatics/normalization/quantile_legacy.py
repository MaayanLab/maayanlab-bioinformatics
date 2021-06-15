import logging
import numpy as np
import pandas as pd
from functools import singledispatch


@singledispatch
def quantile_normalize(mat):
  ''' Perform quantile normalization on the values of a matrix
  In the case of a pd.DataFrame, preserve the index on the output frame.
  See: https://en.wikipedia.org/wiki/Quantile_normalization
  '''
  logging.warn('Unrecognized type: ' + type(mat).__name__)
  return quantile_normalize_np(mat)

@quantile_normalize.register
def quantile_normalize_np(mat: np.ndarray):
  # sort vector in np (reuse in np)
  sorted_vec = np.sort(mat, axis=0)
  # rank vector in np (no dict necessary)
  rank = sorted_vec.mean(axis=1)
  # construct quantile normalized matrix
  return np.array([
    [
      rank[i]
      for i in np.searchsorted(sorted_vec[:, c], mat[:, c])
    ] for c in range(mat.shape[1])
  ]).T

@quantile_normalize.register
def quantile_normalize_pd(mat: pd.DataFrame):
  return pd.DataFrame(
    quantile_normalize_np(mat.values),
    index=mat.index,
    columns=mat.columns,
  )

def quantile_normalize_h5(in_mat, out_mat, tmp=None):
  import os, tempfile, h5py
  '''
  Maximum memory required (3 * in_mat.shape[1] * sizeof(dtype))
  Storage required 4 * in_mat.size
    - input matrix
    - transposed copy
    - sorted & transposed copy
    - output matrix
  '''
  assert isinstance(in_mat, h5py.Dataset)
  assert isinstance(out_mat, h5py.Dataset)
  assert in_mat.shape == out_mat.shape
  # transpose + sort
  tmp_f = tempfile.mktemp() if tmp is None else tmp
  tmp_h5 = h5py.File(tmp_f, 'w')
  tmp_T_mat = tmp_h5.create_dataset('tmp_T', shape=(in_mat.shape[1], in_mat.shape[0]), dtype=in_mat.dtype)
  tmp_T_sorted_mat = tmp_h5.create_dataset('tmp_T_sorted', shape=(in_mat.shape[1], in_mat.shape[0]), dtype=in_mat.dtype)
  sorted_col_vec_agg_rank = np.zeros(in_mat.shape[0])
  for col in range(in_mat.shape[1]):
    # this single read is potentially expensive but the two writes are cheap
    col_vec = in_mat[:, col]
    tmp_T_mat[col, :] = col_vec
    sorted_col_vec = np.sort(col_vec)
    tmp_T_sorted_mat[col, :] = sorted_col_vec
    sorted_col_vec_agg_rank += sorted_col_vec
  # setup rank matrix
  sorted_col_vec_agg_rank /= in_mat.shape[1]
  # construct output matrix
  for c in range(in_mat.shape[1]):
    # this write is potentially expensive but the reads are cheap
    out_mat[:, c] = [
      sorted_col_vec_agg_rank[i]
      for i in np.searchsorted(tmp_T_sorted_mat[c, :], tmp_T_mat[c, :])
    ]
  # close and remove tmp file
  tmp_h5.close()
  os.remove(tmp_f)
  return out_mat
