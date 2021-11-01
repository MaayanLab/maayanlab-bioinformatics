import numpy as np
import pandas as pd

def GSEA2005(geneset_membership: pd.Series, correlations: pd.Series):
  '''
  Implementation of algorithm described here:
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1239896/

  :param geneset_membership: (pd.Series) True if in set, False if not, index: all genes
  :param correlations: (pd.Series) Correlation of a given gene
  :return (Tuple[np.array, np.array]) x and y arrays ready to be plotted. ES = y.max()
  '''
  r_j = correlations.abs().sort_values(ascending=False) # r_j: correlation of gene_j in ranked order
  S = geneset_membership[correlations.index] # S: geneset mask aligned with r_j
  N = S.count()             # N: number of genes
  N_H = S.sum()             # N_H: number of hits
  N_R = r_j[S].sum()        # N_R: sum of r_j for g_j \in S
  P_hit =    S  * r_j/N_R   # P_hit: fraction of hits weighted by r_j
  P_miss = (~S) * 1/(N-N_H) # P_hit: fraction of misses up to position i
  # 0 added to beginning for plotting, doesn't affect sum
  x = np.arange(N + 1)
  y = np.concatenate([[0],np.cumsum(P_hit - P_miss)])
  return x, y
