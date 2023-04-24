import requests
import pandas as pd

def speedenrich(userlist: list[str], libraries: list[str]=None, background: list[str]=None, description='Example gene list', base_url='https://maayanlab.cloud/speedrichr'):
  ''' Perform enrichment analysis using speedrichr.

  :param userlist: A list of genes (e.g. ['PHF14', 'RBM3', 'MSL1', ...])
  :param libraries: One or more Enrichr Libraries (https://maayanlab.cloud/Enrichr/#libraries) (e.g. ["TRANSFAC_and_JASPAR_PWMs"])
  :param background: An optional background geneset for background correction (e.g. ['PHF14', 'RBM3', 'MSL1', ...])
  :param base_url: If using a different enrichr instance than the public one, specify the base prefix
  :return: A pandas dataframe with enrichment results
  '''
  userlist_response = requests.post(
    base_url+'/api/addList',
    files=dict(
      list=(None, '\n'.join(userlist)),
      description=(None, description),
    )
  ).json()
  if background:
    background_response = requests.post(
      base_url+'/api/addbackground',
      data=dict(background='\n'.join(background)),
    ).json()
    results = {}
    for library in set(libraries):
      results.update(
        requests.post(
          base_url+'/api/backgroundenrich',
          data=dict(
            userListId=userlist_response['userListId'],
            backgroundid=background_response['backgroundid'],
            backgroundType=library,
          )
        ).json()
      )
  else:
    results = {}
    for library in set(libraries):
      results.update(
        requests.get(
          base_url+'/api/enrich',
          params=dict(
            userListId=userlist_response['userListId'],
            backgroundType=library,
          ),
        ).json()
      )
  return pd.DataFrame([
    [l, *r]
    for l, result in results.items()
    for r in result
  ], columns=['library', 'rank', 'term', 'pvalue', 'oddsratio', 'combined score', 'overlap', 'adj pvalue', 'legacy_0', 'legacy_1']).sort_values('rank', ascending=True)
