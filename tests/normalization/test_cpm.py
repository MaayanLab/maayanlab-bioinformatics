import numpy as np
from maayanlab_bioinformatics.normalization.cpm import cpm_normalize

def test_cpm_normalization():
  given = np.array([
    [5, 4, 3],
    [2, 1, 4],
    [3, 4, 6],
    [4, 2, 8],
  ])
  # expectation from R: edgeR::cpm(given)
  expectation = np.array([
    [357142.86, 363636.36, 142857.14],
    [142857.14,  90909.09, 190476.19],
    [214285.71, 363636.36, 285714.29],
    [285714.29, 181818.18, 380952.38],
  ])
  assert np.allclose(cpm_normalize(given), expectation, atol=1e-2)
