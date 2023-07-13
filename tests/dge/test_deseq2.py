import pathlib
import pandas as pd
from maayanlab_bioinformatics.dge.deseq2 import deseq2_differential_expression

def test_deseq2():
  df = pd.read_csv(pathlib.Path(__file__).parent.parent/'test_example_matrix.txt', sep='\t', index_col=0)
  df_results = deseq2_differential_expression(
    df.iloc[:, :3],
    df.iloc[:, 3:],
  )
  print(df_results)
  assert (df_results['pvalue'] < 0.05).any()
