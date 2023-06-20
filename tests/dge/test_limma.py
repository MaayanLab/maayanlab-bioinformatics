import pathlib
import numpy as np
import pandas as pd
from maayanlab_bioinformatics.dge import limma_voom_differential_expression

def test_limma():
  df = pd.read_csv(pathlib.Path(__file__).parent.parent/'test_example_matrix.txt', sep='\t', index_col=0)
  df_expected = pd.read_csv(pathlib.Path(__file__).parent.parent/'test_example_matrix_dge_results.txt', sep='\t', index_col=0)
  df_results = limma_voom_differential_expression(
    df.iloc[:, :3],
    df.iloc[:, 3:],
    filter_genes=True,
  )
  print(df_expected)
  print(df_results)
  isclose = pd.DataFrame(np.isclose(df_expected, df_results), index=df_results.index, columns=df_results.columns)
  print(isclose.value_counts())
  assert isclose.all().all()

def test_limma_shuffled():
  df = pd.read_csv(pathlib.Path(__file__).parent.parent/'test_example_matrix.txt', sep='\t', index_col=0)
  df_expected = pd.read_csv(pathlib.Path(__file__).parent.parent/'test_example_matrix_dge_results.txt', sep='\t', index_col=0)
  index = df.index.values.copy()
  np.random.shuffle(index)
  controls, cases = df.iloc[:, :3], df.iloc[:, 3:]
  controls_columns = controls.columns.values.copy()
  np.random.shuffle(controls_columns)
  cases_columns = cases.columns.values.copy()
  np.random.shuffle(cases_columns)
  df_results = limma_voom_differential_expression(
    controls.loc[index, controls_columns], cases.loc[index, cases_columns],
    filter_genes=True,
  )
  print(df_expected)
  print(df_results)
  isclose = pd.DataFrame(np.isclose(df_expected, df_results), index=df_results.index, columns=df_results.columns)
  print(isclose.value_counts())
  assert isclose.all().all()

def test_limma_determanism():
  df = pd.read_csv(pathlib.Path(__file__).parent.parent/'test_example_matrix.txt', sep='\t', index_col=0)
  index = df.index.values.copy()
  np.random.shuffle(index)
  controls, cases = df.iloc[:, :3], df.iloc[:, 3:]
  controls_columns = controls.columns.values.copy()
  np.random.shuffle(controls_columns)
  cases_columns = cases.columns.values.copy()
  np.random.shuffle(cases_columns)
  df_expected = limma_voom_differential_expression(
    controls, cases,
    voom_design=True,
  )
  df_results = limma_voom_differential_expression(
    controls.loc[index, controls_columns], cases.loc[index, cases_columns],
    voom_design=True,
  )
  print(df_expected)
  print(df_results)
  isclose = pd.DataFrame(np.isclose(df_expected, df_results), index=df_results.index, columns=df_results.columns)
  print(isclose.value_counts())
  assert isclose.all().all()
