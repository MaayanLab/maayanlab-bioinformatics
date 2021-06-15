# Test based on information here: https://en.wikipedia.org/wiki/Quantile_normalization

import numpy as np
from maayanlab_bioinformatics.normalization.quantile import quantile_normalize


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
