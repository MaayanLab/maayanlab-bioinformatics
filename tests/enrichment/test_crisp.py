import os
import pandas as pd
from maayanlab_bioinformatics.enrichment import enrich_crisp
from maayanlab_bioinformatics.parse import gmt_read_iter

def test_enrich_crisp():
  # TODO: compare pvalues & oddsratio
  geneset = [
    line.strip().upper()
    for line in open(os.path.join(os.path.dirname(__file__), '..', 'test_geneset.txt'), 'r')
  ]
  library_iter = gmt_read_iter(os.path.join(os.path.dirname(__file__), '..', 'test_gmt.gmt'))
  expectation = pd.read_csv(os.path.join(os.path.dirname(__file__), '..', 'test_enrichr_results.txt'), sep='\t').sort_values('P-value')
  expectation_terms = expectation['Term'].tolist()
  results = sorted(enrich_crisp(geneset, library_iter, 20000, True), key=lambda r: r[1].pvalue)
  result_terms = [term for term, _ in results]
  assert expectation_terms == result_terms, f"{expectation_terms} != {result_terms}"
