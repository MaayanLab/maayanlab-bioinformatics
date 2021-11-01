import numpy as np
import pandas as pd

def GSEA2003(geneset_membership: pd.Series, gene_difference_metric: pd.Series):
  '''
  Implementation of algorithm described here:
  https://pubmed.ncbi.nlm.nih.gov/12808457/

  :param geneset_membership: (pd.Series) True if in set, False if not, index: all genes
  :param gene_difference_metric: (pd.Series) Difference metric between two classes, e.g. SNR difference
  :return (Tuple[np.array, np.array]) x and y arrays ready to be plotted. ES = y.max()
  '''
  R_i = gene_difference_metric.sort_values(ascending=False) # R_1, ... R_N ordered by difference metric
  S = geneset_membership[R_i.index] # S containing gene_membership members
  G = geneset_membership.sum()
  N = geneset_membership.count()
  X = (
        S  * np.sqrt((N - G) / G) # X_i when member of S
    - (~S) * np.sqrt(G / (N - G)) # X_i when not member of S
  )
  # 0 added to beginning for plotting, doesn't affect sum
  x = np.arange(N + 1)
  y = np.concatenate([[0],np.cumsum(X)])
  return x, y
