import pandas as pd
import scipy.stats
from maayanlab_bioinformatics.normalization import log2_normalize

def ttest_differential_expression(controls_mat: pd.DataFrame, cases_mat: pd.DataFrame, equal_var=False, alternative='two-sided', log2norm=True):
  ''' Given two separate dataframes (controls, cases) with a shared index (genes),
  we compute the ttest differential expression for all genes. Benjamini-Hochberg Adjusted p-value.

  :param controls_mat: (pd.DataFrame) the control samples (samples as columns and genes as rows)
  :param cases_mat: (pd.DataFrame) the case samples (samples as columns and genes as rows)
  :param equal_var: (bool) Should t-test assume equal variance (default: False)
  :param alternative: (str) Alternative hypothesis (see scipy.stats.ttest_ind) (default: two-sided)
  :param log2norm: (bool) Apply log2norm, typically keep with raw counts but disable if you have normalized data (default: True)
  :return: A data frame with the results
  '''
  assert (controls_mat.index == cases_mat.index).all(), 'Index between controls and cases must be the same'
  if log2norm:
    cases_mat = log2_normalize(cases_mat)
    controls_mat = log2_normalize(controls_mat)
  results = scipy.stats.ttest_ind(cases_mat.T, controls_mat.T, equal_var=equal_var, alternative=alternative)
  df_results = pd.DataFrame({
    'Statistic': results.statistic,
    'Pval': results.pvalue,
  }, index=controls_mat.index)
  df_results['AdjPval'] = scipy.stats.false_discovery_control(df_results['Pval'].fillna(1.), method='bh')
  df_results.sort_values('AdjPval', inplace=True)
  return df_results
