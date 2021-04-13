
import numpy as np
from maayanlab_bioinformatics.normalization.cpm import cpm_normalize

def test_cpm_normalization():
  given = np.array([
    [5, 4, 3],
    [2, 1, 4],
    [3, 4, 6],
    [4, 2, 8],
  ])
  from rpy2.robjects import numpy2ri
  from rpy2.robjects.packages import importr
  edgeR = importr('edgeR')
  expectation = numpy2ri.rpy2py(edgeR.cpm(numpy2ri.py2rpy(given)))
  assert np.allclose(cpm_normalize(given), expectation, atol=1e-2)
