import time
import requests
import pandas as pd

def enrichr_link_from_genes(genes, description='', enrichr_link='https://maayanlab.cloud/Enrichr', sleep=1):
  ''' Functional access to Enrichr API
  '''
  if sleep:
    time.sleep(sleep)
  resp = requests.post(enrichr_link + '/addList', files={
    'list': (None, '\n'.join(genes)),
    'description': (None, description),
  })
  if resp.status_code != 200:
    raise Exception('Enrichr failed with status {}: {}'.format(
      resp.status_code,
      resp.text,
    ))
  result = resp.json()
  return dict(result, link=enrichr_link + '/enrich?dataset=' + resp.json()['shortId'])

def enrichr_get_top_results(userListId, bg, enrichr_link='https://maayanlab.cloud/Enrichr', sleep=1):
  if sleep:
    time.sleep(sleep)
  resp = requests.get(enrichr_link + '/enrich?userListId={}&backgroundType={}'.format(userListId, bg))
  if resp.status_code != 200:
    raise Exception('Enrichr failed with status {}: {}'.format(
      resp.status_code,
      resp.text,
    ))
  return pd.DataFrame(resp.json()[bg], columns=[
    'rank', 'term', 'pvalue', 'zscore', 'combinedscore', 'overlapping_genes', 'adjusted_pvalue', '', ''
  ])

def enrichr_term_genes(bg, terms, enrichr_link='https://maayanlab.cloud/Enrichr', sleep=1):
  if sleep:
    time.sleep(sleep)
  resp = requests.get(enrichr_link + '/geneSetLibrary?mode=json&libraryName={}&term={}'.format(
      bg,
      ';'.join(terms),
  ))
  if resp.status_code != 200:
    raise Exception('Enrichr failed with status {}: {}'.format(
      resp.status_code,
      resp.text,
    ))
  return resp.json()

class EnrichrUserList:
  ''' Object oriented access to Enrichr results.

  Example:
  ```python
  from maayanlab_bioinformatics.api import EnrichrUserList
  mylist = EnrichrUserList([
    'STAT1', 'ACE2', #...
  ], 'mylist')
  print(mylist.link)
  mylist['GO_Biological_Process_2021'] # returns dataframe with enrichment results
  ```
  '''
  def __init__(self, genes, description='', enrichr_link='https://maayanlab.cloud/Enrichr'):
    self._enrichr_link = enrichr_link
    self._genes = genes
    self._description = description
    self._shortId = None
    self._userListId = None
    self._results = {}

  def __repr__(self):
    return f"EnrichrUserList(..., description={repr(self.description)})"
  
  @property
  def genes(self):
    return self._genes

  @property
  def description(self):
    return self._description

  @property
  def link(self):
    return enrichr_link + '/enrich?dataset=' + self.shortId

  @property
  def shortId(self):
    if self._shortId is None: self._addList()
    return self._shortId

  @property
  def userListId(self):
    if self._userListId is None: self._addList()
    return self._userListId
  
  def __getitem__(self, library):
    if library not in self._results: self._enrich(library)
    return self._results[library].copy()

  def _addList(self):
    resp = requests.post(self._enrichr_link + '/addList', files={
      'list': (None, '\n'.join(self.genes)),
      'description': (None, self.description),
    })
    if resp.status_code != 200:
      raise Exception('Enrichr failed with status {}: {}'.format(
        resp.status_code,
        resp.text,
      ))
    results = resp.json()
    # wait a tinybit before returning link (backoff)
    time.sleep(1)
    #
    self._userListId = results['userListId']
    self._shortId = results['shortId']
  
  def _enrich(self, library):
    resp = requests.get(self._enrichr_link + '/enrich', params={
      'userListId': self.userListId,
      'backgroundType': library,
    })
    if resp.status_code != 200:
      raise Exception('Enrichr failed with status {}: {}'.format(
        resp.status_code,
        resp.text,
      ))
    results = resp.json()
    # wait a tinybit before returning link (backoff)
    time.sleep(1)
    self._results[library] = pd.DataFrame([
      dict(
        rank=rank,
        term=term,
        pvalue=pvalue,
        zscore=zscore,
        combinedscore=combinedscore,
        overlapping_genes=overlapping_genes,
        adjusted_pvalue=adjusted_pvalue,
      )
      for rank, term, pvalue, zscore, combinedscore, overlapping_genes, adjusted_pvalue, _, _ in results[library]
    ])
