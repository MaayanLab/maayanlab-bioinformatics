import re
import json
import pandas as pd
import logging

def _try_json_loads(s):
  try:
    return json.loads(s)
  except:
    return s

_gene_parser = re.compile(r'^([\w\d]+)[^\w\d]+(.+)$')
def _parse_gene(gene):
  match = _gene_parser.match(gene.strip())
  if match:
    return {match.group(1): _try_json_loads(match.group(2)) }
  else:
    return { gene.strip(): 1 }

def read_gmt(fh, parse_gene=_parse_gene):
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
  for n, line in enumerate(fh):
    try:
      term1, term2, genes_str = line.strip().split('\t', maxsplit=2)
    except:
      logging.warn('Ignoring line {}:{} because it seems empty'.format(n, line))
    term = '\t'.join(filter(None, map(str.strip, (term1, term2))))
    if term in gmt:
      logging.warn('Duplicate term: {}, merging'.format(term))
    else:
      gmt[term] = {}
    for parsed_gene in map(parse_gene, genes_str.split('\t')):
      gmt[term].update(parsed_gene)
  return gmt

def pd_read_gmt(fh, parse_gene=_parse_gene):
  ''' Read .gmt files directly into a data frame.
  '''
  return pd.DataFrame(read_gmt(fh, parse_gene=parse_gene))
