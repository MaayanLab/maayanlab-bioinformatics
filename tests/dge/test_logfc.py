import pathlib
import pandas as pd
import numpy as np
from maayanlab_bioinformatics.dge.logfc import logfc_differential_expression

def test_logfc():
  df = pd.read_csv(pathlib.Path(__file__).parent.parent/'test_example_matrix.txt', sep='\t', index_col=0)
  df_expected = pd.read_csv(pathlib.Path(__file__).parent.parent/'test_example_matrix_dge_results.txt', sep='\t', index_col=0)
  df_results = logfc_differential_expression(
    df.iloc[:, :3],
    df.iloc[:, 3:],
  )
  df_cmp = pd.concat({'expected': df_expected['logFC'], 'computed': df_results.loc[df_expected.index, 'LogFC']}, axis=1)
  df_cmp_close = np.isclose(df_cmp['expected'], df_cmp['computed'], atol=1.)
  print(df_cmp)
  print(df_cmp[~df_cmp_close])
  # most logfc values are pretty close to those reported by limma
  #  they won't be the same since limma applies a normalization but they should be similar orders of magnitude
  #  so we ensure 95% are relatively close
  assert df_cmp_close.sum() > (0.95*df_expected.shape[0])
