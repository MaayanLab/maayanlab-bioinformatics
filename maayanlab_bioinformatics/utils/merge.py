import pandas as pd

def merge(left, *rights, **kwargs):
  ''' Helper function for many trivial (index based) joins
  '''
  merged = left
  for right in rights:
    merged = pd.merge(left=merged, left_index=True, right=right, right_index=True, **kwargs)
  return merged
