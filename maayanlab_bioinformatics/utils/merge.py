import pandas as pd

def merge(*dfs, **kwargs):
  ''' Helper function for many trivial (index based) joins
  '''
  if not dfs:
    return pd.DataFrame()
  #
  left, *rights = dfs
  #
  merged = left
  for right in rights:
    merged = pd.merge(left=merged, left_index=True, right=right, right_index=True, **kwargs)
  #
  return merged
