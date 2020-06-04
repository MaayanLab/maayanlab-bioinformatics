import pandas as pd
from functools import lru_cache
from maayanlab_bioinformatics.utils import fetch_save_read, merge

@lru_cache()
def ncbi_genes_fetch(organism='Mammalia/Homo_sapiens'):
  ''' Fetch the current NCBI Human Gene Info database.
  See ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/ for the directory/file of the organism of interest.
  '''
  ncbi = fetch_save_read(
    'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/{}.gene_info.gz'.format(organism),
    '{}.gene_info.tsv'.format(organism),
    sep='\t',
  )
  # Ensure nulls are treated as such
  ncbi = ncbi.applymap(lambda v: float('nan') if type(v) == str and v == '-' else v)
  # Break up lists
  split_list = lambda v: v.split('|') if type(v) == str else []
  ncbi['dbXrefs'] = ncbi['dbXrefs'].apply(split_list)
  ncbi['Synonyms'] = ncbi['Synonyms'].apply(split_list)
  ncbi['LocusTag'] = ncbi['LocusTag'].apply(split_list)
  ncbi['Other_designations'] = ncbi['Other_designations'].apply(split_list)
  return ncbi

@lru_cache()
def ncbi_genes_lookup(organism='Mammalia/Homo_sapiens'):
  ''' Return a lookup dictionary with synonyms as the keys, and official symbols as the values
  Usage:
  ```python
  ncbi_lookup = ncbi_genes_lookup('Mammalia/Homo_sapiens')
  print(ncbi_lookup('STAT3')) # any alias will get converted into the official symbol
  ```
  '''
  ncbi_genes = ncbi_genes_fetch(organism=organism)
  ncbi_lookup = {
    sym: row['Symbol']
    for _, row in ncbi_genes.iterrows()
    for sym in [row['Symbol']] + row['Synonyms']
  }
  return ncbi_lookup.get
