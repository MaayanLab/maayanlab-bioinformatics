import numpy as np
import pandas as pd
from rpy2.robjects import r

# Seemless conversion between numpy and R types
from rpy2.robjects import r, numpy2ri, pandas2ri
numpy2ri.activate()
pandas2ri.activate()

# Load all libraries and create diffExpression function
r('''
library("R.utils")
library("RCurl")
library("DESeq2")
library("limma")
library("statmod")
library("edgeR")
diffExpression <- function(controls_mat, samples_mat, filter_genes=FALSE, contrast=FALSE, adjust="BH") {
  controls <- as.data.frame(controls_mat)
  samples <- as.data.frame(samples_mat)
  # build a single expression matrix from controls and samples Nx(M+V)
  expression = do.call(cbind, append(controls, samples))
  # design matrix is used in linear regression model to define the controls and samples
  # They can be more complicated in when there is other variables other than treatment and control
  # here it will create an intercept and then also a 1 -1 vector to differentiate the samples
  design = c(rep(1,length(controls)), rep(-1, length(samples)))
  dm = model.matrix(~as.factor(design))
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
  v <- voom(dge, dm, plot=FALSE)
  # this is basically just applying a linear fit. The test will calculate how much the differentiation between controls and samples
  # improves the fit over the simplest possible model (average)
  fit <- lmFit(v, dm)
  if (isTRUE(contrast)) {
    cont.matrix <- makeContrasts(de=B-A, levels=design)
    fit <- contrasts.fit(fit, cont.matrix)
  }
  # this will calculate moderated t-statistics using empirical bayes moderation of the standard errors towards a common value
  fit <- eBayes(fit)
  # Get results
  limma_dataframe <- topTable(fit, adjust=adjust, number=nrow(expression))
  limma_dataframe$gene_symbol <- rownames(limma_dataframe)
  #
  return (limma_dataframe)
}
''')
def limmaVoomDifferentialExpression(controls_mat, cases_mat, filter_genes=False, contrast=False):
  assert (controls_mat.index == cases_mat.index).all(), 'Index between controls and cases must be the same'
  # genes on rows, samples on cols
  results = pd.DataFrame(r['diffExpression'](controls_mat, cases_mat, filter_genes=filter_genes, contrast=contrast))
  # recover gene symbols
  ind_lookup = {str(k): v for k, v in enumerate(controls_mat.index, start=1)}
  results.index = results['gene_symbol'].apply(ind_lookup.get)
  # return resulting frame
  return results.drop('gene_symbol', axis=1)

def upDownFromLimmaVoom(expr, top_n=600):
  most_significant_expr = expr.sort_values('P.Value').iloc[:top_n]
  return type('UpDownGeneset', tuple(), dict(
    up=list(most_significant_expr[most_significant_expr['t'] > 0].dropna().index),
    down=list(most_significant_expr[most_significant_expr['t'] < 0].dropna().index),
  ))
