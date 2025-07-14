import os
import pandas as pd

def fetch_save_read(url, file, cache_dir='', reader=pd.read_csv, sep=',', **kwargs):
  ''' Download file from {url}, save it to {file}, and subsequently read it with {reader} using pandas options on {**kwargs}.
  '''
  if cache_dir or cache_dir == '':
    path = file if cache_dir == '' else cache_dir + os.path.sep + file
    if not os.path.exists(path):
      if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
      df = reader(url, sep=sep, index_col=None)
      df.to_csv(path, sep=sep, index=False)
    return reader(path, sep=sep, **kwargs)
  else:
    return reader(url, sep=sep, **kwargs)

