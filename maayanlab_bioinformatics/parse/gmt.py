import re
import io
import math as m
import pandas as pd
import contextlib
import logging
import pathlib
import typing
from dataclasses import dataclass
from maayanlab_bioinformatics.utils.maybe_tqdm import maybe_tqdm

def _try_load_number(s):
  try:
    return int(s)
  except ValueError:
    pass
  try:
    return float(s)
  except ValueError:
    pass
  return s

@contextlib.contextmanager
def _ensure_fp(fp, mode):
  if type(fp) == str or isinstance(fp, pathlib.Path):
    with open(fp, mode) as fh:
      yield fh
  else:
    yield fp

def parse_gene_weight(gene):
  ''' A helper to parse the gmt potentially with numeric weights
  '''
  gene, *_weight = re.split(r'[,:;]', gene.strip(), maxsplit=1)
  if _weight:
    _weight, = _weight
    _weight = _try_load_number(_weight)
    if type(_weight) == str:
      gene += _weight
      weight = 1
    else:
      weight = _weight
  else:
    weight = 1
  return gene.strip(), weight

def parse_gene_unweighted(gene):
  ''' A helper to parse the gmt unweighted
  '''
  return gene.strip(), 1

def gmt_read_iter(fh, parse_gene=parse_gene_weight):
  with _ensure_fp(fh, 'r') as fh:
    for n, line in enumerate(fh):
      try:
        term1, term2, genes_str = line.strip().split('\t', maxsplit=2)
      except ValueError:
        logging.warn('Ignoring line {}:{} because it seems empty'.format(n, line))
        continue
      term = '\t'.join(filter(None, map(str.strip, (term1, term2))))
      geneset = {
        k: v
        for k, v in map(parse_gene, genes_str.split('\t'))
        if k
      }
      yield term, geneset

def gmt_read_dict(fh, parse_gene=parse_gene_weight):
  ''' Read .gmt files into a dictionary of the form:
  {
    'term_1\tterm_2': {
      gene_1: weight or 1,
      ...
    },
    ...
  }

  If your genes are encoded in a weird way you can also provide your own `parse_gene` function,
   the current one supports just gene names or gene names with weights separated by non-word/numeric characters.
  '''
  gmt = {}
  for n, (term, geneset) in enumerate(gmt_read_iter(fh, parse_gene=parse_gene)):
    if term in gmt:
      logging.warn('Duplicate term: {}:{}, merging'.format(n, term))
    else:
      gmt[term] = {}
    gmt[term].update(**geneset)
  return gmt

def gmt_read_pd(fh, parse_gene=parse_gene_weight):
  ''' Read .gmt files directly into a data frame.
  '''
  return pd.DataFrame(gmt_read_dict(fh, parse_gene=parse_gene))


def _serialize_gene_weight_pair(gene, weight):
  if weight == 1 or m.isclose(weight, 1.): return gene
  elif m.isclose(weight, 0.) or m.isnan(weight): return None
  else: return '{},{}'.format(gene, weight)

def _ensure_weight(gs):
  if isinstance(gs, dict):
    return gs.items()
  else:
    return ((g, 1) for g in gs)

def gmt_write_dict(gmt, fh, serialize_gene_weight_pair=_serialize_gene_weight_pair):
  ''' Opposite of gmt_read_dict, write a dictionary to a file pointer
  serialize_gene_weight_pair can be used to customize serialization when dealing with weights.
    - it should return the serialized gene,weight pair or None if it should be removed
  By default, 0/nans are dropped, 1s result in a gene (crisp), and everything else uses gene,weight.
  '''
  with _ensure_fp(fh, 'w') as fh:
    for term, geneset in gmt.items():
      if '\t' not in term: serialized_term = term + '\t'
      else: serialized_term = term
      serialized_geneset = '\t'.join(filter(None, (
        serialize_gene_weight_pair(gene, weight)
        for gene, weight in _ensure_weight(geneset)
      )))
      if not serialized_geneset:
        logging.warn('Ignoring term {} because its geneset seems empty'.format(term))
        continue
      print(serialized_term, serialized_geneset, sep='\t', file=fh)

def gmt_write_pd(df, fh, serialize_gene_weight_pair=_serialize_gene_weight_pair):
  ''' Write a pandas dataframe as a gmt, where rows are genes and columns are terms.
  See gmt_write_dict for more information.
  '''
  gmt_write_dict(df.to_dict(), fh, serialize_gene_weight_pair=serialize_gene_weight_pair)

@dataclass
class GMT:
  ''' a data structure for GMTs in memory
  '''
  # the unique set of genes across all gene lists
  background: list[str]
  # first two columns of the GMT
  terms: list[(str, str)]
  # variable gene lists of the GMT
  gene_lists: list[list[str]]

  @staticmethod
  def reader(gmtfile: io.TextIOBase | str | pathlib.Path) -> typing.Iterator[tuple[tuple[str, str], list[str]]]:
    ''' read the .gmt format, a tab separated file with variable columns
    '''
    gene_expr = re.compile(r'^([^:;,]+?)([:;,].+)?$')
    with _ensure_fp(gmtfile, 'r') as fr:
      for line in fr:
        line_split = [cell.strip() for cell in line.strip().split('\t')]
        if len(line_split) < 3: continue
        term, desc, *genes = line_split
        genes = [
          m.group(1)
          for gene in genes
          if gene
          for m in (gene_expr.match(gene),)
          if m
        ]
        yield (term, desc), genes

  @staticmethod
  def from_iter(it: typing.Iterator[tuple[tuple[str, str], list[str]]]):
    ''' initialize a GMT from Iterator[(term, desc), gene_list] (i.e. read_gmt)
    '''
    background = set()
    terms = []
    gene_lists = []
    for (term, desc), genes in maybe_tqdm(it, desc='Reading gmt...'):
      background.update(genes)
      terms.append((term, desc))
      gene_lists.append(genes)
    return GMT(list(background), terms, gene_lists)

  @staticmethod
  def concat(*gmts):
    background = set()
    terms = []
    gene_lists = []
    for gmt in gmts:
      background.update(gmt.background)
      terms += gmt.terms
      gene_lists += gmt.gene_lists
    return GMT(list(background), terms, gene_lists)

  @staticmethod
  def from_file(gmtfile: io.TextIOBase | str | pathlib.Path):
    ''' initialze a GMT from a file
    '''
    return GMT.from_iter(GMT.reader(gmtfile))

  def to_spmatrix(self):
    ''' create a sparse matrix from this GMT
    '''
    import scipy.sparse
    import numpy as np
    spmatrix = scipy.sparse.dok_matrix((len(self.gene_lists), len(self.background)), dtype=np.int8)
    gene_index = { gene: index for index, gene in enumerate(self.background) }
    for i, gene_list in enumerate(maybe_tqdm(self.gene_lists, desc='Building spmatrix...')):
      spmatrix[i, [gene_index[g] for g in gene_list]] = 1
    return spmatrix

  def to_df(self):
    ''' create a sparse pandas dataframe from this GMT
    '''
    import pandas as pd
    return pd.DataFrame.sparse.from_spmatrix(
      self.to_spmatrix(),
      columns=self.background,
      index=self.terms,
    )

  def dedupe(self):
    ''' de-duplicate gene sets in a GMT
    '''
    deduped_terms = []
    deduped_gene_lists = []
    gene_list_hashes = {}
    for term, gene_list in maybe_tqdm(zip(self.terms, self.gene_lists), desc='De-duping...'):
      gene_list_hash = hash(frozenset(gene_list))
      if gene_list_hash not in gene_list_hashes:
        gene_list_hashes[gene_list_hash] = len(deduped_terms)-1
        deduped_terms.append(term)
        deduped_gene_lists.append(gene_list)
      else:
        gene_list_index = gene_list_hashes[gene_list_hash]
        deduped_terms[gene_list_index] = (*deduped_terms[gene_list_index], *term)
    return GMT(self.background, deduped_terms, deduped_gene_lists)
