# Test based on information here: https://en.wikipedia.org/wiki/Quantile_normalization

import numpy as np
from maayanlab_bioinformatics.normalization.quantile import quantile_normalize, quantile_normalize_h5


def test_quantile_normalization():
  given = np.array([
    [5, 4, 3],
    [2, 1, 4],
    [3, 4, 6],
    [4, 2, 8],
  ])
  expectation = np.array([
    [5.67, 5.17, 2.00],
    [2.00, 2.00, 3.00],
    [3.00, 5.17, 4.67],
    [4.67, 3.00, 5.67],
  ])
  assert np.allclose(quantile_normalize(given), expectation, atol=1e-2)

def test_quantile_normalization_h5():
  import os
  import h5py
  import tempfile
  fname = tempfile.mktemp('test_quantile.h5')
  f = h5py.File(fname, 'w')
  given = f.create_dataset('given', data=np.array([
      [5, 4, 3],
      [2, 1, 4],
      [3, 4, 6],
      [4, 2, 8],
  ]))
  norm = f.create_dataset('norm', shape=given.shape, dtype='float64')
  quantile_normalize_h5(given, norm)
  expectation = np.array([
    [5.67, 4.67, 2.00],
    [2.00, 2.00, 3.00],
    [3.00, 4.67, 4.67],
    [4.67, 3.00, 5.67],
  ])
  assert np.allclose(norm[:], expectation, atol=1e-2)
  os.unlink(fname)
