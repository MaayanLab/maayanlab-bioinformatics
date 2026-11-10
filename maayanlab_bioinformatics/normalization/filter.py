import logging
import numpy as np
import pandas as pd


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


def filter_by_expr_mask(counts, design=None, group=None, lib_size=None, min_count=10, min_total_count=15, large_n=10, min_prop=0.7, tol=1e-14) -> np.ndarray:
  ''' Port of edgeR's filterByExpr.default https://rdrr.io/bioc/edgeR/src/R/filterByExpr.R

  :param counts: (np.ndarray | pd.DataFrame) raw counts, genes as rows and samples as columns
  :param design: (np.ndarray | pd.DataFrame) design matrix, samples as rows
  :param group: (array-like) group label for each sample, takes precedence over design
  :param lib_size: (array-like) library size for each sample, defaults to column sums of counts
  :return: (np.ndarray) boolean mask, True for each gene (row) to keep
  '''
  y = np.asarray(counts, dtype=np.float64)
  lib_size = y.sum(axis=0) if lib_size is None else np.asarray(lib_size, dtype=np.float64)
  # Minimum effective sample size for any of the coefficients
  if group is None:
    if design is None:
      logging.warning('No group or design set. Assuming all samples belong to one group.')
      min_sample_size = y.shape[1]
    else:
      # R's hat(design) adds an intercept column; the leverage of the most
      #  influential sample gives the smallest group size
      X = np.column_stack([np.ones(y.shape[1]), np.asarray(design, dtype=np.float64)])
      h = np.einsum('ij,ji->i', X, np.linalg.pinv(X))
      min_sample_size = 1 / np.max(h)
  else:
    # size of the smallest non-empty group
    min_sample_size = pd.Series(np.asarray(group)).value_counts().min()
  #
  if min_sample_size > large_n:
    min_sample_size = large_n + (min_sample_size - large_n) * min_prop
  # CPM cutoff
  cpm_cutoff = min_count / np.median(lib_size) * 1e6
  cpm = y / lib_size * 1e6
  keep_cpm = (cpm >= cpm_cutoff).sum(axis=1) >= (min_sample_size - tol)
  # Total count cutoff
  keep_total_count = y.sum(axis=1) >= (min_total_count - tol)
  #
  return keep_cpm & keep_total_count


def filter_by_expr(mat, design=None, group=None, min_count=10, min_total_count=15, large_n=10, min_prop=0.7, tol=1e-14, lib_size=None):
  ''' Port of edgeR's filterByExpr, keeping only genes (rows) with sufficiently large counts.
  Takes a matrix of raw counts (genes as rows, samples as columns) and returns it filtered,
  see `filter_by_expr_mask` for the parameters.
  '''
  keep = filter_by_expr_mask(
    mat, design=design, group=group, lib_size=lib_size,
    min_count=min_count, min_total_count=min_total_count,
    large_n=large_n, min_prop=min_prop, tol=tol,
  )
  if isinstance(mat, pd.DataFrame):
    return mat.loc[keep, :]
  return mat[keep]
