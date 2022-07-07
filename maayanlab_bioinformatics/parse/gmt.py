import re
import json
import math as m
import pandas as pd
import logging

def _try_json_loads(s):
  try:
    return json.loads(s)
  except:
    return s

def _ensure_fp(fp, mode):
  if type(fp) == str: return open(fp, mode)
  else: return fp

_gene_parser = re.compile(r'^(.+?)([,:;](.+))?$')
def _parse_gene(gene):
  match = _gene_parser.match(gene.strip())
  if match and match.group(3):
    return match.group(1), _try_json_loads(match.group(3))
  else:
    return gene.strip(), 1

def gmt_read_iter(fh, parse_gene=_parse_gene):
  fh = _ensure_fp(fh, 'r')
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

def gmt_read_dict(fh, parse_gene=_parse_gene):
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

def gmt_read_pd(fh, parse_gene=_parse_gene):
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
  fh = _ensure_fp(fh, 'w')
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
