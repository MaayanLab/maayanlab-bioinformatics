import os
import pandas as pd

def fetch_save_read(url, file, reader=pd.read_csv, sep=',', **kwargs):
  ''' Download file from {url}, save it to {file}, and subsequently read it with {reader} using pandas options on {**kwargs}.
  '''
  if not os.path.exists(file):
    if os.path.dirname(file):
      os.makedirs(os.path.dirname(file), exist_ok=True)
    df = reader(url, sep=sep, index_col=None)
    df.to_csv(file, sep=sep, index=False)
  return pd.read_csv(file, sep=sep, **kwargs)
