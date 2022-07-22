''' Chunked module has useful helper functions for manipulating ndarrays in chunks, this is
especially useful when working with h5py matrices since operations which respect chunk boundaries
avoid excessive disk random access.
'''
import numpy as np
import itertools as it
import logging
logger = logging.getLogger(__name__)

try:
  from tqdm import tqdm
except ImportError:
  logger.warning('install tqdm for progress bar')
  tqdm = lambda it, **kwargs: it

def chunk_slices(shape, chunks, progress=False):
  '''
  Return slices to chunk through an ndarray.

  :param shape:    The shape of the ndarray or size in 1d.
  :param chunks:   The shape of the chunks or size in all dimensions.
  :param progress: Show tqdm progress bar or not

  :returns: Iterator[slice(start, stop) for each dimension in shape]

  Usage:
  N = np.arange(10)
  [N[s] for s in chunk_slices(len(N), 3)]

  I = np.eye(10)
  [I[i, j] for i, j in chunk_slices(I.shape, 3)]
  '''
  if not progress:
    tqdm = lambda it, **kwargs: it
  # ensure shape is a tuple (we'll flatten it when we're done if not)
  if type(shape) == int:
    shape = (shape,)
    flatten = True
  else:
    flatten = False
  # spread chunks to length of shape
  if type(chunks) == int:
    chunks = (chunks,)*len(shape)
  # ensure shape and chunks are the same dimensions
  assert len(shape) == len(chunks)
  for indices in tqdm(
    # we take the iterable cartesian product, equivalent to a nested for loop on each dimension
    it.product(*[
      # for each dimension, shape / chunks + 1 when chunks are not divisible by shape
      range((s//cs)+int(s%cs!=0))
      for s, cs in zip(shape, chunks)
    ]),
    # the total (for progress bar) is just the actual product of steps in each dimension
    total=np.product([(s//cs)+int(s%cs!=0) for s, cs in zip(shape, chunks)]),
  ):
    # indices have the chunk start index, we can then get the slice (start, end) for each dimension
    slices = tuple([slice(i*cs, min(s, (i+1)*cs)) for i, s, cs in zip(indices, shape, chunks)])
    # if not multiple dimensions, we'll flatten the slice
    if flatten:
      slices, = slices
    yield slices

def chunk_infer(x, chunks=None):
  ''' Helper function for interpreting the chunks param with respect to a matrix x

  :param x:      The matrix (ndarray)
  :param chunks: The chunks parameter,
                 if None (default): Try to infer from chunks attribute (h5py)
                 if int: Use a multiple of the inferred chunks attribute, or alternatively that size in each dimension
                 if tuple: Use the explicit chunks provided for slicing

  :returns: tuple chunks parameter
  '''
  if type(chunks) in (tuple, list):
    assert len(chunks) == len(x.shape), f"Matrix has {len(x.shape)} dimensions but chunks has {len(chunks)}"
    return chunks
  inferred_chunks = getattr(x, 'chunks', None)
  if inferred_chunks is not None:
    if chunks is None:
      chunks = inferred_chunks
    elif type(chunks) == int:
      chunks = np.array(inferred_chunks)*chunks
    else:
      raise NotImplementedError('chunks should be int or tuple')
  elif type(chunks) == int:
    chunks = (chunks,)*len(x.shape)
  else:
    raise NotImplementedError('chunks should be int or tuple')
  return tuple(chunks)

def chunk_applymap(func, x, *, out=None, chunks=None, progress=False):
  ''' Apply function to all elements in a matrix in chunks

  :param func:     The function to apply to each chunk
  :param x:        The matrix to apply it to
  :param out:      The matrix to write to (pass variable to out for inplace)
  :param chunks:   The shape of the chunks in each dimension,
                   can be inferred for h5py arrays based on actual chunks on disk,
                   can be a multiple of an integer value of chunks.
  :param progress: Show tqdm progress bar or not

  :returns:        The augmented matrix (or the original matrix, augmented)
  '''
  if out is None: out = np.zeros(x.shape)
  if x.shape != out.shape: logger.warning('x and out shape do not match, is this what you wanted?')
  if x.dtype != out.dtype: logger.warning('x and out dtype do not match, is this what you wanted?')
  for s in chunk_slices(x.shape, chunks=chunk_infer(x, chunks), progress=progress):
    out[s] = func(x[s])
  return out
