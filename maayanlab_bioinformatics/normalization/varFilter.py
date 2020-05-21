import logging
import numpy as np
import pandas as pd
from functools import singledispatch


@singledispatch
def selectHighVariance(mat, top_n=2500, axis=1):
  ''' Select rows with the most variable expression accross all samples.
  Takes a dataframe and returns a filtered dataframe in the same orientation.
  e.g.
        |condition_1|condition_2|
  gene_1|    1      |     1     |
  gene_2|    0      |    10     |

  gene_1 here is *not* variable at all,
  gene_2 here is *very* variable.

  gene_1 will be dropped, while gene_2 is kept.
  '''
  logging.warn('Unrecognized type, note that pd.DataFrame is expected to keep track of the index: ' + type(mat).__name__)
  return selectHighVariance_pd(mat, top_n=top_n)

@selectHighVariance.register
def selectHighVariance_pd(mat: pd.DataFrame, top_n=2500):
  return mat.loc[mat.var(axis=1).sort_values(ascending=False).index[:top_n], :]

