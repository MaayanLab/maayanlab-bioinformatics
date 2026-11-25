import numpy as np
import pandas as pd
from typing import Optional, Tuple
from maayanlab_bioinformatics.normalization.filter import filter_by_expr_mask

# pylimma reports topTable columns in snake_case, map back to R's limma names
_TOP_TABLE_COLUMNS = {
  'log_fc': 'logFC',
  'ave_expr': 'AveExpr',
  'p_value': 'P.Value',
  'adj_p_value': 'adj.P.Val',
  'b': 'B',
}

def _calc_norm_factors_tmm(counts: np.ndarray, lib_size: np.ndarray, logratio_trim=0.3, sum_trim=0.05, a_cutoff=-1e10):
  ''' Port of edgeR's calcNormFactors(method="TMM") with default arguments.
  '''
  x = counts[(counts > 0).any(axis=1)]
  # reference column: the sample whose upper quartile is closest to the mean upper quartile
  f75 = np.quantile(x / lib_size, 0.75, axis=0)
  if np.median(f75) < 1e-20:
    ref = int(np.argmax(np.sqrt(x).sum(axis=0)))
  else:
    ref = int(np.argmin(np.abs(f75 - f75.mean())))
  f = np.array([
    _calc_factor_tmm(x[:, i], x[:, ref], lib_size[i], lib_size[ref], logratio_trim, sum_trim, a_cutoff)
    for i in range(x.shape[1])
  ])
  # factors multiply to one
  return f / np.exp(np.mean(np.log(f)))

def _calc_factor_tmm(obs, ref, n_o, n_r, logratio_trim, sum_trim, a_cutoff):
  import scipy.stats
  with np.errstate(divide='ignore', invalid='ignore'):
    log_r = np.log2((obs / n_o) / (ref / n_r))
    abs_e = (np.log2(obs / n_o) + np.log2(ref / n_r)) / 2
    v = (n_o - obs) / n_o / obs + (n_r - ref) / n_r / ref
  fin = np.isfinite(log_r) & np.isfinite(abs_e) & (abs_e > a_cutoff)
  log_r, abs_e, v = log_r[fin], abs_e[fin], v[fin]
  if np.max(np.abs(log_r)) < 1e-6:
    return 1.0
  n = len(log_r)
  lo_l = np.floor(n * logratio_trim) + 1
  hi_l = n + 1 - lo_l
  lo_s = np.floor(n * sum_trim) + 1
  hi_s = n + 1 - lo_s
  rank_r = scipy.stats.rankdata(log_r)
  rank_e = scipy.stats.rankdata(abs_e)
  keep = (rank_r >= lo_l) & (rank_r <= hi_l) & (rank_e >= lo_s) & (rank_e <= hi_s)
  f = np.nansum(log_r[keep] / v[keep]) / np.nansum(1 / v[keep])
  if np.isnan(f):
    f = 0.0
  return 2 ** f

def _diff_expression(expression: pd.DataFrame, design_dataframe: pd.DataFrame, filter_genes=False, voom_design=False, adjust='BH'):
  ''' Python equivalent of the edgeR/limma pipeline:
  DGEList -> filterByExpr -> calcNormFactors -> voom -> lmFit -> contrasts.fit(B-A) -> eBayes -> topTable
  '''
  import pylimma
  counts = expression.to_numpy(dtype=np.float64)
  genes = expression.index
  design = design_dataframe.to_numpy(dtype=np.float64)
  # like DGEList, library sizes are computed before filtering & kept after
  lib_size = counts.sum(axis=0)
  if filter_genes:
    keep = filter_by_expr_mask(counts, design, lib_size=lib_size)
    counts, genes = counts[keep], genes[keep]
  # TMM normalization factors scale the effective library sizes used by voom
  lib_size = lib_size * _calc_norm_factors_tmm(counts, lib_size)
  # R's voom uses a fixed span by default
  v = pylimma.voom(counts, design if voom_design else None, lib_size=lib_size, adaptive_span=False)
  fit = pylimma.lm_fit(v, design)
  cont_matrix = pylimma.make_contrasts(de='B-A', levels=list(design_dataframe.columns))
  fit = pylimma.contrasts_fit(fit, cont_matrix.to_numpy())
  fit = pylimma.e_bayes(fit)
  limma_dataframe = pylimma.top_table(
    fit, coef=0, number=len(genes), adjust_method=adjust,
    genelist=pd.DataFrame({'gene_symbol': genes}),
  ).rename(columns=_TOP_TABLE_COLUMNS)
  return limma_dataframe[['logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B', 'gene_symbol']]

def limma_voom_differential_expression_design(
  expression: pd.DataFrame,
  design: pd.DataFrame,
  de: Tuple[str, str],
  filter_genes: bool = False,
  voom_design: bool = False,
):
  ''' Use voom and limma (via pylimma) for differential expression.

  Note that this function expects the original raw gene counts.

  alex version: voom_design=True, filter_genes=False
  biojupies version: voom_design=False, filter_genes=True

  :param expression: (pd.DataFrame) the samples
  :param design: (pd.DataFrame) the design dataframe
  :param de: (Tuple[str, str]) ('control_column_name', 'case_column_name')
  :param filter_genes: (bool) Whether to perform edgeR's `filterByExpr` during normalization
  :param voom_design: (bool) Whether to give the voom function the design matrix (supervised)
  :return: A data frame with the results
  '''
  expression = expression.copy()
  design = design.copy().loc[expression.columns]
  expression.columns = design.index = ['s'+str(i) for i, _ in enumerate(expression.columns)]
  design.columns = [
    {de[0]: 'A', de[1]: 'B'}.get(col, 'c'+str(i))
    for i, col in enumerate(design.columns)
  ]
  assert 'A' in design.columns and 'B' in design.columns
  return _diff_expression(
    expression,
    design,
    filter_genes=filter_genes,
    voom_design=voom_design,
  ).sort_values('t', ascending=False).set_index('gene_symbol')

def make_design_matrix(expression_df, controls, cases):
  expression_df = expression_df.copy()
  expression_df.index.name = 'index'
  expression_df = expression_df.reset_index().groupby('index').sum()
  design_df = pd.DataFrame([{'index': 's'+str(i), 'A': int(x in controls), 'B': int(x in cases)} for i, x in enumerate(expression_df.columns)]).set_index('index')
  expression_df.columns = ['s'+str(i) for i, _ in enumerate(expression_df.columns)]
  return expression_df, design_df

def limma_voom_differential_expression(
  controls_mat: pd.DataFrame,
  cases_mat: pd.DataFrame,
  all_data_mat: Optional[pd.DataFrame] = None,
  filter_genes: bool = False,
  voom_design: bool = False,
):
  ''' Use voom and limma (via pylimma) for differential expression.

  Note that this function expects the original raw gene counts.

  alex version: voom_design=True, filter_genes=False
  biojupies version: voom_design=False, filter_genes=True

  :param controls_mat: (pd.DataFrame) the control samples
  :param cases_mat: (pd.DataFrame) the case samples
  :param all_data_mat: (pd.DataFrame) *all* samples (for full experiment normalization)
  :param filter_genes: (bool) Whether to perform edgeR's `filterByExpr` during normalization
  :param voom_design: (bool) Whether to give the voom function the design matrix (supervised)
  :return: A data frame with the results
  '''
  if all_data_mat is None:
    all_data_mat = pd.concat([controls_mat, cases_mat], axis=1)
  all_data_mat.columns = all_data_mat.columns.to_flat_index()
  # transform input into expression/design ready for the limma pipeline
  expression, design = make_design_matrix(all_data_mat, set(controls_mat.columns.to_flat_index()), set(cases_mat.columns.to_flat_index()))
  return _diff_expression(
    expression,
    design,
    filter_genes=filter_genes,
    voom_design=voom_design,
  ).sort_values('t', ascending=False).set_index('gene_symbol')

def up_down_from_limma_voom(expr: pd.DataFrame, top_n: int = 600):
  ''' Given limma_voom_differential_expression output, produce a discrete up/down geneset

  :param top_n: (int) the number of genes in total to produce
  :return: UpDownGeneset a type with `.up` and `.down` methods corresponding to the list of genes.
  '''
  most_significant_expr = expr.sort_values('P.Value').iloc[:top_n]
  return type('UpDownGeneset', tuple(), dict(
    up=list(most_significant_expr[most_significant_expr['logFC'] > 0].dropna().index),
    down=list(most_significant_expr[most_significant_expr['logFC'] < 0].dropna().index),
  ))
