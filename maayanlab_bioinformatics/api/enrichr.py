import time
import requests
import pandas as pd

def enrichr_link_from_genes(genes, description='', enrichr_link='https://amp.pharm.mssm.edu/Enrichr', sleep=1):
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

def enrichr_get_top_results(userListId, bg, enrichr_link='https://amp.pharm.mssm.edu/Enrichr', sleep=1):
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

def enrichr_term_genes(bg, terms, enrichr_link='https://amp.pharm.mssm.edu/Enrichr', sleep=1):
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
