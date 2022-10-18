import pandas as pd

from maayanlab_bioinformatics.utils.fetch_save_read import fetch_save_read

def mouse_human_homologs(uppercase=False):
  ''' Returns a dataframe with mouse/human gene mappings based on MGI.
  See: http://www.informatics.jax.org/homology.shtml
  
  @param uppercase: bool should mappings be uppercase (i.e. for case insensitive mapping)
  @returns pd.DataFrame
  ```
  |mouse|human|
  |-----|-----|
  |sp140|SP140|
  ```
  '''
  mouse_human_sequence = fetch_save_read(
    'http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt',
    'HOM_MouseHumanSequence.rpt',
    sep='\t',
  )
  mouse_human_sequence_simplified = pd.DataFrame([
    dict(
      mouse=d.loc[d['Common Organism Name'].str.contains('mouse'), 'Symbol'].values,
      human=d.loc[d['Common Organism Name'].str.contains('human'), 'Symbol'].values,
    )
    for _, d in mouse_human_sequence.groupby('DB Class Key')
  ]).explode('mouse').explode('human').dropna()
  if uppercase:
    mouse_human_sequence_simplified['mouse'] = mouse_human_sequence_simplified['mouse'].str.upper()
    mouse_human_sequence_simplified['human'] = mouse_human_sequence_simplified['human'].str.upper()
  return mouse_human_sequence_simplified

def human_expression_to_mouse(human_expression, strategy='sum', uppercase=False):
  ''' Given a human expression matrix, produce a mouse-compatible expression matrix by mapping
  homologs.

  @param human_expression: pd.DataFrame(columns=samples, index=human_genes, values=counts)
  @param strategy: 'sum' -- the strategy to use when aggregating duplicates
  @returns pd.DataFrame(columns=samples, index=mouse_genes, values=counts)
  '''
  if strategy == 'sum':
    mouse_expression = pd.merge(
      left=human_expression.set_index(human_expression.index.str.upper()), left_index=True,
      right=mouse_human_homologs(uppercase=uppercase), right_on='human'
    ).groupby('mouse').sum()
  else:
    raise NotImplementedError
  return mouse_expression

def mouse_expression_to_human(mouse_expression, strategy='sum', uppercase=False):
  ''' Given a mouse expression matrix, produce a human-compatible expression matrix by mapping
  homologs.
  
  @param mouse_expression: pd.DataFrame(columns=samples, index=mouse_genes, values=counts)
  @param strategy: 'sum' -- the strategy to use when aggregating duplicates
  @returns pd.DataFrame(columns=samples, index=human_genes, values=counts)
  '''
  if strategy == 'sum':
    human_expression = pd.merge(
      left=mouse_expression.set_index(mouse_expression.index.str.upper()), left_index=True,
      right=mouse_human_homologs(uppercase=uppercase), right_on='mouse'
    ).groupby('human').sum()
  else:
    raise NotImplementedError
  return human_expression
