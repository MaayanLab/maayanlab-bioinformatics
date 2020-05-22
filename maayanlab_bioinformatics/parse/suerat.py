import os
import pandas as pd
import scipy.sparse as sp_sparse
from maayanlab_bioinformatics.utils import merge

def suerat_load(base_dir):
  ''' Files prepared for suerat are quite common, this function will load them
  given the directory that contains `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.tsv.gz`.
  '''
  df_barcodes = pd.read_csv(
    os.path.join(base_dir, 'barcodes.tsv.gz'),
    index_col=0,
    header=None,
    sep='\t',
  )
  df_features = pd.read_csv(
    os.path.join(base_dir, 'features.tsv.gz'),
    header=None,
    names=['symbol', 'type'],
    index_col=0,
    sep='\t',
  )
  matrix = pd.read_csv(
    os.path.join(base_dir, 'matrix.mtx.gz'),
    header=None,
    names=['indices', 'indptr', 'data'],
    skiprows=2,
    sep=' ',
  )
  csc_matrix = sp_sparse.csc_matrix(
    (
      matrix['data'].values,
      (
        matrix['indices'].values - 1, # 0 based indexing
        matrix['indptr'].values - 1,  # 0 based indexing
      )
    ),
  )
  df_expression = pd.DataFrame(csc_matrix.todense())
  df_expression.index = df_features.index
  df_expression.columns = df_barcodes.index
  return df_features, df_barcodes, df_expression

def suerat_load_multiple(base_dirs):
  ''' Sets of suerat directories that are meant to be analyzed together are quite common,
  providing all those directories to this function (much like load_suerat_files) will load
  each individually and return a merged version that captures the filename in the barcodes.
  '''
  all_df_features = []
  all_df_barcodes = []
  all_df_expression = []
  #
  for ind, base_dir in enumerate(base_dirs):
    df_features, df_barcodes, df_expression = suerat_load(base_dir)
    df_barcodes['barcode'] = df_barcodes.index
    df_barcodes['file'] = f'File {ind}'
    df_barcodes.index = df_barcodes.index.map(lambda s, ind=ind: f'{ind}:{s}')
    df_expression.columns = df_barcodes.index
    all_df_features.append(df_features)
    all_df_barcodes.append(df_barcodes)
    all_df_expression.append(df_expression)
  #
  df_features = merge(*all_df_features, how='left', suffixes=('', '_')).drop(['symbol_', 'type_'], axis=1)
  df_barcodes = pd.concat(all_df_barcodes)
  df_expression = merge(*all_df_expression)
  #
  return df_features, df_barcodes, df_expression
