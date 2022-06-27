import pandas as pd
from typing import Dict, Optional

from maayanlab_bioinformatics.harmonization.ncbi_genes import ncbi_genes_lookup
from maayanlab_bioinformatics.utils.merge import merge

def transcripts_to_genes(
  df_expression: pd.DataFrame,
  df_features: pd.DataFrame=None,
  strategy='var',
  uppercasegenes=False,
  lookup_dict: Optional[Dict[str, str]]=None,
  organism='Mammalia/Homo_sapiens',
):
  ''' Map gene alternative ids/transcripts to gene symbols using `ncbi_genes_lookup`
  We take a matrix with genes/transcripts on the rows and samples on the columns.
  In the case of multiple gene/transcript to symbol mappings, we adopt the collision strategy specified.
  If df_features is provided, we will use 'symbol' column as the transcript names,
   otherwise we will use the df_expression index column.
  The resulting matrix will naturally have fewer samples, corresponding to gene symbols in the
   `lookup_dict` which defaults to official ncbi_gene symbols for homo sapiens.
  
  :param strategy: ('var'|'sum') collision strategy (select one with highest variance, or sum counts)
  '''
  # resolve lookup_dict if necessary
  if lookup_dict is None:
    lookup_dict = ncbi_genes_lookup(organism=organism)
  elif callable(lookup_dict):
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
  if strategy == 'var':
    df_transcript_genes = merge(
      df_expression.var(axis=1).to_frame('var'),
      df_features[['symbol']].applymap(lambda s: lookup_dict(s))
    ).groupby('symbol')['var'].idxmax().reset_index()
    df_transcript_genes.index = df_transcript_genes['var']
    df_transcript_genes = df_transcript_genes.drop('var', axis=1)
    # perform the actual mapping
    df_gene_expression = df_expression.loc[df_transcript_genes.index]
    df_gene_expression.index = df_transcript_genes['symbol']
  elif strategy == 'sum':
    df_gene_expression = merge(
      df_expression,
      df_features[['symbol']].applymap(lambda s: lookup_dict(s))
    ).groupby('symbol').sum()
  else:
    raise NotImplementedError
  return df_gene_expression
