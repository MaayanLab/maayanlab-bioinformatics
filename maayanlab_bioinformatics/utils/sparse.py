import scipy.sparse as sp_sparse

def sp_hdf_dump(hdf, sdf, **kwargs):
  ''' Dump Sparse Pandas DataFrame to h5py object.

  Usage:
  ```python
  import h5py
  import pandas as pd
  import scipy.sparse as sp_sparse

  # write
  f = h5py.File('sparse.h5', 'w')
  sdf = pd.DataFrame.sparse.from_spmatrix(sp_sparse.eye(3))
  sp_hdf_dump(f, sdf)
  f.close()
  ```
  '''
  s = sdf.sparse.to_coo()
  hdf.create_dataset('data', data=s.data, **kwargs)
  hdf.create_dataset('row', data=s.row, **kwargs)
  hdf.create_dataset('col', data=s.col, **kwargs)
  hdf.create_dataset('index', data=sdf.index.values, **kwargs)
  hdf.create_dataset('columns', data=sdf.columns.values, **kwargs)
  hdf.attrs['shape'] = s.shape
  return hdf

def sp_hdf_load(hdf):
  ''' Load Sparse Pandas DataFrame from h5py object.

  Usage:
  ```python
  import h5py
  import pandas as pd
  import scipy.sparse as sp_sparse

  f = h5py.File('sparse.h5', 'r')
  sdf = sp_hdf_load(f)
  f.close()
  ```
  '''
  import pandas as pd
  return pd.DataFrame.sparse.from_spmatrix(
    sp_sparse.coo_array((hdf['data'], (hdf['row'], hdf['col'])), shape=hdf.attrs['shape']),
    index=pd.Series(hdf['index']).str.decode('utf8'),
    columns=pd.Series(hdf['columns']).str.decode('utf8'),
  )

def sp_std(X_ij, ddof=1):
  ''' Standard deviation for a matrix compatible with sparse matrices.
  i is the row index, j is the column index.

  \sigma_j = \sqrt{\frac{\sum(x_ij - \mu_j)^2}{N_j - ddof}}}
  '''
  N_j = X_ij.shape[-1]
  mu_j = X_ij.sum(axis=0) / N_j
  num_j = ((X_ij - mu_j)**2).sum(axis=0)
  denom_j = N_j - ddof
  if sp_sparse.isspmatrix(X_ij):
    return (num_j / denom_j).A.squeeze()**(1/2)
  else:
    return (num_j / denom_j)**(1/2)
