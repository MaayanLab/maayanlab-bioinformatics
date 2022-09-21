import pandas as pd
from functools import lru_cache
from maayanlab_bioinformatics.utils import fetch_save_read, merge

@lru_cache()
def ncbi_genes_fetch(organism='Mammalia/Homo_sapiens', filters=lambda ncbi: ncbi['type_of_gene']=='protein-coding'):
  ''' Fetch the current NCBI Human Gene Info database.
  See ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/ for the directory/file of the organism of interest.
  '''
  def maybe_split(record):
    ''' NCBI Stores Nulls as '-' and lists '|' delimited
    '''
    if record in {'', '-'}:
      return set()
    return set(record.split('|'))
  #
  def supplement_dbXref_prefix_omitted(ids):
    ''' NCBI Stores external IDS with Foreign:ID while most datasets just use the ID
    '''
    for id in ids:
      # add original id
      yield id
      # also add id *without* prefix
      if ':' in id:
        yield id.split(':', maxsplit=1)[1]
  #
  ncbi = fetch_save_read(
    'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/{}.gene_info.gz'.format(organism),
    '{}.gene_info.tsv'.format(organism),
    sep='\t',
  )
  if filters and callable(filters):
    ncbi = ncbi[filters(ncbi)]
  #
  ncbi['All_synonyms'] = [
    set.union(
      maybe_split(gene_info['Symbol']),
      maybe_split(gene_info['Symbol_from_nomenclature_authority']),
      maybe_split(str(gene_info['GeneID'])),
      maybe_split(gene_info['Synonyms']),
      maybe_split(gene_info['Other_designations']),
      maybe_split(gene_info['LocusTag']),
      set(supplement_dbXref_prefix_omitted(maybe_split(gene_info['dbXrefs']))),
    )
    for _, gene_info in ncbi.iterrows()
  ]
  return ncbi

@lru_cache()
def ncbi_genes_lookup(organism='Mammalia/Homo_sapiens', filters=lambda ncbi: ncbi['type_of_gene']=='protein-coding'):
  ''' Return a lookup dictionary with synonyms as the keys, and official symbols as the values
  Usage:
  ```python
  ncbi_lookup = ncbi_genes_lookup('Mammalia/Homo_sapiens')
  print(ncbi_lookup('STAT3')) # any alias will get converted into the official symbol
  ```
  '''
  ncbi_genes = ncbi_genes_fetch(organism=organism)
  synonyms, symbols = zip(*{
    (synonym, gene_info['Symbol'])
    for _, gene_info in ncbi_genes.iterrows()
    for synonym in gene_info['All_synonyms']
  })
  ncbi_lookup = pd.Series(symbols, index=synonyms)
  index_values = ncbi_lookup.index.value_counts()
  ambiguous = index_values[index_values > 1].index
  ncbi_lookup_disambiguated = ncbi_lookup[(
    (ncbi_lookup.index == ncbi_lookup) | (~ncbi_lookup.index.isin(ambiguous))
  )]
  return ncbi_lookup_disambiguated.to_dict().get
