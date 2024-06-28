def maybe_tqdm(iterable, **kwargs):
  ''' Optional tqdm (omitted if tqdm is not installed)
  '''
  try:
    from tqdm.auto import tqdm
    return tqdm(iterable, **kwargs)
  except ImportError:
    return iterable
