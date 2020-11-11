import pandas as pd
from typing import Dict

from maayanlab_bioinformatics.harmonization.ncbi_genes import ncbi_genes_lookup
from maayanlab_bioinformatics.utils.merge import merge

def transcripts_to_genes(df_expression: pd.DataFrame, df_features: pd.DataFrame=None, organism='Mammalia/Homo_sapiens', uppercasegenes=True, lookup_dict: Dict[str, str]=ncbi_genes_lookup):
  ''' Map genes to transcripts given a matrix with transcripts using `ncbi_genes_lookup`
  We take a matrix with genes on the rows.
  In the case of multiple transcript to gene mappings, we keep the one with the highest variance.
  If df_features is provided, we will use 'symbol' column as the transcript names,
   otherwise we will use the df_expression index column.
  The resulting matrix will naturally have fewer samples, corresponding to gene symbols in the
   `lookup_dict` which defaults to official ncbi_gene symbols for homo sapiens.
  '''
  # resolve lookup_dict if necessary
  if callable(lookup_dict):
    lookup_dict = lookup_dict()
  # construct df_features if not provided
  if df_features is None:
    df_features = pd.Series(df_expression.index).to_frame('symbol')
    df_features.index = df_expression.index
  # uppercase genes if necessary
  if uppercasegenes:
    df_features['symbol'] = df_features['symbol'].apply(str.upper)
  # get df_expression but only the highest variance transcript that
  #  corresponds to the same set of genes
  df_transcript_genes = merge(
    df_expression.var(axis=1).to_frame('var'),
    df_features[['symbol']].applymap(lambda s: lookup_dict(s))
  ).groupby('symbol')['var'].idxmax().reset_index()
  df_transcript_genes.index = df_transcript_genes['var']
  df_transcript_genes = df_transcript_genes.drop('var', axis=1)
  # perform the actual mapping
  df_gene_expression = df_expression.loc[df_transcript_genes.index]
  df_gene_expression.index = df_transcript_genes['symbol']
  return df_gene_expression
