import numpy as np
import pandas as pd
from typing import Optional, Tuple

R = None

def _lazy_load():
  global R
  if R is not None:
    return R
  #
  from rpy2.robjects import r
  #
  r('''
  suppressMessages(library("R.utils"))
  suppressMessages(library("RCurl"))
  suppressMessages(library("DESeq2"))
  suppressMessages(library("limma"))
  suppressMessages(library("statmod"))
  suppressMessages(library("edgeR"))
  diffExpression <- function(expression, design_dataframe, filter_genes=FALSE, voom_design=FALSE, adjust="BH") {
    design <- as.matrix(design_dataframe)
    # turn count matrix into a expression object compatible with edgeR
    dge <- DGEList(counts=expression)
    # filter genes
    if (isTRUE(filter_genes)) {
      keep <- filterByExpr(dge, design)
      dge <- dge[keep,]
    }
    # calculate normalization factors, here the different library sizes per sample are accounted for
    # the normalization factors are appended to the dge object and used in later steps
    dge <- calcNormFactors(dge)
    # to be honest I am not sure what exactly happens here. Limma was developed for affymetrix chips that have
    # values that follow different distributions than read counts. It will apply some sort of log transformation
    # and make it compatible with lmFit
    if (isTRUE(voom_design)) {
      v <- voom(dge, design, plot=FALSE)
    } else {
      v <- voom(dge, plot=FALSE)
    }
    # this is basically just applying a linear fit. The test will calculate how much the differentiation between controls and samples
    # improves the fit over the simplest possible model (average)
    fit <- lmFit(v, design)
    # this is what makes it differential expression from B - A, B and A are set in the design matrix
    cont.matrix <- makeContrasts(de=B-A, levels=design)
    fit <- contrasts.fit(fit, cont.matrix)
    # this will calculate moderated t-statistics using empirical bayes moderation of the standard errors towards a common value
    fit <- eBayes(fit)
    # Get results
    limma_dataframe <- topTable(fit, adjust=adjust, number=nrow(expression))
    limma_dataframe$gene_symbol <- rownames(limma_dataframe)
    #
    return (limma_dataframe)
  }
  ''')
  R = r
  return R

def limma_voom_differential_expression_design(
  expression: pd.DataFrame,
  design: pd.DataFrame,
  de: Tuple[str, str],
  filter_genes: bool = False,
  voom_design: bool = False,
):
  ''' Use R's voom and limma for differential expression.

  Note that this function expects the original raw gene counts.

  alex version: voom_design=True, filter_genes=False
  biojupies version: voom_design=False, filter_genes=True

  :param expression: (pd.DataFrame) the samples
  :param design: (pd.DataFrame) the design dataframe
  :param de: (Tuple[str, str]) ('control_column_name', 'case_column_name')
  :param filter_genes: (bool) Whether to perform R's `filterByExpr` during normalization
  :param voom_design: (bool) Whether to give R's voom function the design matrix (supervised)
  :return: A data frame with the results
  '''
  design = design.copy().loc[expression.columns]
  design.columns = [
    {de[0]: 'A', de[1]: 'B'}.get(col, str(i))
    for i, col in enumerate(design.columns)
  ]
  assert 'A' in design.columns and 'B' in design.columns
  import rpy2.robjects as ro
  from rpy2.robjects import pandas2ri
  from rpy2.robjects.conversion import localconverter
  r = _lazy_load()
  with localconverter(ro.default_converter + pandas2ri.converter):
    return r.diffExpression(
      expression,
      design,
      filter_genes=filter_genes,
      voom_design=voom_design,
    ).sort_values('t', ascending=False).set_index('gene_symbol')

def make_design_matrix(expression_df, controls, cases):
  expression_df = expression_df.copy()
  expression_df.index.name = 'index'
  expression_df = expression_df.reset_index().groupby('index').sum()
  design_df = pd.DataFrame([{'index': str(i), 'A': int(x in controls), 'B': int(x in cases)} for i, x in enumerate(expression_df.columns)]).set_index('index')
  expression_df.columns = [str(i) for i, _ in enumerate(expression_df.columns)]
  return expression_df, design_df

def limma_voom_differential_expression(
  controls_mat: pd.DataFrame,
  cases_mat: pd.DataFrame,
  all_data_mat: Optional[pd.DataFrame] = None,
  filter_genes: bool = False,
  voom_design: bool = False,
):
  ''' Use R's voom and limma for differential expression.

  Note that this function expects the original raw gene counts.

  alex version: voom_design=True, filter_genes=False
  biojupies version: voom_design=False, filter_genes=True

  :param controls_mat: (pd.DataFrame) the control samples
  :param cases_mat: (pd.DataFrame) the case samples
  :param all_data_mat: (pd.DataFrame) *all* samples (for full experiment normalization)
  :param filter_genes: (bool) Whether to perform R's `filterByExpr` during normalization
  :param voom_design: (bool) Whether to give R's voom function the design matrix (supervised)
  :return: A data frame with the results
  '''
  if all_data_mat is None:
    all_data_mat = pd.concat([controls_mat, cases_mat], axis=1)
  all_data_mat.columns = all_data_mat.columns.to_flat_index()
  # transform input into expression/design ready for R functions
  expression, design = make_design_matrix(all_data_mat, set(controls_mat.columns.to_flat_index()), set(cases_mat.columns.to_flat_index()))
  import rpy2.robjects as ro
  from rpy2.robjects import pandas2ri
  from rpy2.robjects.conversion import localconverter
  r = _lazy_load()
  with localconverter(ro.default_converter + pandas2ri.converter):
    return r.diffExpression(
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
