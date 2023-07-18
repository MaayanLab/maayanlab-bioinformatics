import pandas as pd
from maayanlab_bioinformatics.normalization import log2_normalize

def logfc_differential_expression(controls_mat: pd.DataFrame, cases_mat: pd.DataFrame):
  ''' NOT RECOMMENDED. Given two separate dataframes (controls, cases) with a shared index (genes),
  we compute the logFC differential expression for all genes, also the average expression.

  :param controls_mat: (pd.DataFrame) the control samples (samples as columns and genes as rows)
  :param cases_mat: (pd.DataFrame) the case samples (samples as columns and genes as rows)
  :return: A data frame with the results
  '''
  assert (controls_mat.index == cases_mat.index).all(), 'Index between controls and cases must be the same'
  df_results = pd.DataFrame({
    'LogFC': log2_normalize(cases_mat.mean(axis=1)) - log2_normalize(controls_mat.mean(axis=1)),
  }, index=controls_mat.index)
  df_results.sort_values('LogFC', key=lambda logfc: logfc.abs(), ascending=False, inplace=True)
  return df_results
