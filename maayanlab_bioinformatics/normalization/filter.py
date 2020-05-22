import logging
import numpy as np
import pandas as pd
from functools import singledispatch
from maayanlab_bioinformatics.normalization.cpm import cpm_normalize


def filter_by_var(mat: pd.DataFrame, top_n=2500, axis=1):
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
  return mat.loc[mat.var(axis=1).sort_values(ascending=False).index[:top_n], :]


def filter_by_expr(mat, design: pd.Series=None, group: pd.DataFrame=None, min_count=10, min_total_count=15, large_n=10, min_prop=0.7, tol=1e-14):
  ''' Ported from R https://rdrr.io/bioc/edgeR/src/R/filterByExpr.R
  '''
  lib_size = mat.sum(axis=0)
  #  Minimum effect sample sample size for any of the coefficients
  if group is None:
    if design is None:
      logging.warn('No group or design set. Assuming all samples belong to one group.')
      min_sample_size = mat.shape[1]
    else:
      min_sample_size = 1 / design.max()
  else:
    min_sample_size = group[group > 0].min()
  #
  if min_sample_size > large_n:
    min_sample_size = large_n + (min_sample_size - large_n) * min_prop
  # CPM cutoff
  median_lib_size = lib_size.median()
  cpm_cutoff = (min_count / median_lib_size) * 1e6
  cpm = cpm_normalize(mat)
  keep_cpm = (cpm >= cpm_cutoff).sum(axis=1) >= (min_sample_size - tol)
  # Total count cutoff
  keep_total_count = mat.sum(axis=1) >= min_total_count - tol
  #
  return mat.loc[keep_cpm & keep_total_count, :]
