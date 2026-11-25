import pathlib
import numpy as np
import pandas as pd
from maayanlab_bioinformatics.normalization import filter_by_expr

root = pathlib.Path(__file__).parent.parent

def _load():
  df = pd.read_csv(root/'test_example_matrix.txt', sep='\t', index_col=0)
  # the R limma pipeline results contain exactly the genes kept by edgeR's filterByExpr(dge, design)
  expected = pd.read_csv(root/'test_example_matrix_dge_results.txt', sep='\t', index_col=0).index
  design = pd.DataFrame({'A': [1, 1, 1, 0, 0, 0], 'B': [0, 0, 0, 1, 1, 1]}, index=df.columns)
  return df, design, expected

def test_filter_by_expr_design():
  df, design, expected = _load()
  result = filter_by_expr(df, design=design)
  assert set(result.index) == set(expected)
  # order and orientation preserved
  assert list(result.columns) == list(df.columns)
  assert list(result.index) == [g for g in df.index if g in set(expected)]

def test_filter_by_expr_group():
  df, _, expected = _load()
  result = filter_by_expr(df, group=['control']*3 + ['case']*3)
  assert set(result.index) == set(expected)

def test_filter_by_expr_numpy():
  df, design, expected = _load()
  result = filter_by_expr(df.to_numpy(), design=design.to_numpy())
  assert result.shape == (len(expected), df.shape[1])

def test_filter_by_expr_large_n():
  # 12 samples per group exceeds large_n=10: min sample size = 10 + (12-10)*0.7 = 11.4,
  #  so a gene needs 12 samples above the cpm cutoff
  counts = np.full((3, 24), 100)
  counts[1, :11] = 0 # 13 expressed samples -> kept
  counts[2, :13] = 0 # 11 expressed samples -> dropped
  result = filter_by_expr(counts, group=[0]*12 + [1]*12)
  assert (result == counts[:2]).all()
