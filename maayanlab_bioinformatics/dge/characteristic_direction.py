import numpy as np
import pandas as pd
from scipy.stats import chi2
from scipy.stats.mstats import zscore
from sklearn.decomposition import PCA

# TODO: revamp _chdir
def _chdir(data, sampleclass, genes, gamma=1., sort=True, calculate_sig=False, nnull=10, sig_only=False, norm_vector=True):
  """ repurposed from https://github.com/MaayanLab/geode/blob/master/geode/geode.py#L10-L115

  Calculate the characteristic direction for a gene expression dataset
  
  Input:
    data: numpy.array, is the data matrix of gene expression where rows correspond to genes and columns correspond to samples
    sampleclass: list or numpy.array, labels of the samples, it has to be consist of 0, 1 and 2, with 0 being columns to be excluded, 1 being control and 2 being perturbation
        example: sampleclass = [1,1,1,2,2,2]
    genes: list or numpy.array, row labels for genes 
    gamma: float, regulaized term. A parameter that smooths the covariance matrix and reduces potential noise in the dataset
    sort: bool, whether to sort the output by the absolute value of chdir
    calculate_sig: bool, whether to calculate the significance of characteristic directions
    nnull: int, number of null characteristic directions to calculate for significance
    sig_only: bool, whether to return only significant genes; active only when calculate_sig is True
    norm_vector: bool, whether to return a characteristic direction vector normalized to unit vector
  Output:
    A list of tuples sorted by the absolute value in descending order characteristic directions of genes.
      If calculate_sig is set to True, each tuple contains a third element which is the ratio of characteristic directions to null ChDir
  """
  
  ## check input
  data.astype(float)
  sampleclass = np.array(list(map(int, sampleclass)))
  # masks
  m_non0 = sampleclass != 0
  m1 = sampleclass[m_non0] == 1
  m2 = sampleclass[m_non0] == 2

  if type(gamma) not in [float, int]:
    raise ValueError("gamma has to be a numeric number")
  if set(sampleclass) != set([1,2]) and set(sampleclass) != set([0,1,2]):
    raise ValueError("sampleclass has to be a list whose elements are in only 0, 1 or 2")
  # if m1.sum()<2 or m2.sum()<2:
  #   raise ValueError("Too few samples to calculate characteristic directions")
  if len(genes) != data.shape[0]:
    raise ValueError("Number of genes does not match the demension of the expression matrix")

  ## normalize data
  data = data[:, m_non0]
  data = zscore(data) # standardize for each genes across samples

  ## start to compute
  n1 = m1.sum() # number of controls
  n2 = m2.sum() # number of experiments

  ## the difference between experiment mean vector and control mean vector.
  meanvec = data[:,m2].mean(axis=1) - data[:,m1].mean(axis=1) 

  ## initialize the pca object
  pca = PCA(n_components=None)
  pca.fit(data.T)

  ## compute the number of PCs to keep
  cumsum = pca.explained_variance_ratio_ # explained variance of each PC
  keepPC = len(cumsum[cumsum > 0.001]) # number of PCs to keep

  v = pca.components_[0:keepPC].T # rotated data 
  r = pca.transform(data.T)[:,0:keepPC] # transformed data

  dd = ( np.dot(r[m1].T,r[m1]) + np.dot(r[m2].T,r[m2]) ) / float(n1+n2-2) # covariance
  sigma = np.mean(np.diag(dd)) # the scalar covariance

  shrunkMats = np.linalg.inv(gamma*dd + sigma*(1-gamma)*np.eye(keepPC))

  b = np.dot(v, np.dot(np.dot(v.T, meanvec), shrunkMats))

  if norm_vector:
    b /= np.linalg.norm(b) # normalize b to unit vector

  grouped = zip([abs(item) for item in b],b,genes)
  if sort:
    grouped = sorted(grouped,key=lambda x: x[0], reverse=True)


  if not calculate_sig: # return sorted b and genes.
    res = [(item[1],item[2]) for item in grouped]
    return res
  else: # generate a null distribution of chdirs
    nu = n1 + n2 - 2
    y1 = np.random.multivariate_normal(np.zeros(keepPC), dd, nnull).T * np.sqrt(nu / chi2.rvs(nu,size=nnull))
    y2 = np.random.multivariate_normal(np.zeros(keepPC), dd, nnull).T * np.sqrt(nu / chi2.rvs(nu,size=nnull))
    y = y2 - y1 ## y is the null of v

    nullchdirs = []
    for col in y.T:
      bn = np.dot(np.dot(np.dot(v,shrunkMats), v.T), np.dot(col,v.T))
      bn /= np.linalg.norm(bn)
      bn = bn ** 2
      bn.sort()
      bn = bn[::-1] ## sort in decending order
      nullchdirs.append(bn)

    nullchdirs = np.array(nullchdirs).T
    nullchdirs = nullchdirs.mean(axis=1)
    b_s = b ** 2 
    b_s.sort()
    b_s = b_s[::-1] # sorted b in decending order
    relerr = b_s / nullchdirs ## relative error
    # ratio_to_null
    ratios = np.cumsum(relerr)/np.sum(relerr)- np.linspace(1./len(meanvec),1,len(meanvec))
    res = [(item[1],item[2], ratio) for item, ratio in zip(grouped, ratios)] 
    # print('Number of significant genes: %s'%(np.argmax(ratios)+1))
    if sig_only:
      return res[0:np.argmax(ratios)+1]
    else:
      return res

def characteristic_direction(controls_mat: pd.DataFrame, cases_mat: pd.DataFrame, gamma=1., nnull=10, norm_vector=True, sort=True, calculate_sig=False):
  ''' Given two separate dataframes (controls, cases) with a shared index (genes), we compute the characteristic direction coefficients for all genes.
  e.g.

  control_mat:
        |control_replicate_1|control_replicate_2|...
  gene_1|        ..         |        ..         |...
  gene_2|        ..         |        ..         |...

  cases_mat:
        |case_replicate_1|case_replicate_2|...
  gene_1|     ..         |     ..         |...
  gene_2|     ..         |     ..         |...
  '''
  assert (controls_mat.index == cases_mat.index).all(), 'Index between controls and cases must be the same'
  n_genes = controls_mat.shape[0]
  # Compute characteristic direction
  results = pd.DataFrame(
    data=_chdir(
      np.concatenate([controls_mat, cases_mat], axis=1),
      np.array(
        [1]*controls_mat.shape[1] + [2]*cases_mat.shape[1]
      ),
      # genes
      np.array(list(controls_mat.index)),
      gamma=gamma,
      nnull=nnull,
      norm_vector=norm_vector,
      sort=False,
      calculate_sig=calculate_sig,
    ),
    columns=['CD-coefficient', 'Index', 'Significance'] if calculate_sig else ['CD-coefficient', 'Index'],
  )
  results.index = results['Index']
  results.index.name = controls_mat.index.name
  if sort:
    results = results.sort_values('CD-coefficient')
  return results.drop('Index', axis=1)

def up_down_from_characteristic_direction(expr: pd.DataFrame, top_n=600):
  ''' Using the output of `characteristic_direction`, we can extract the top n genes
  with the highest absolute characteristic direction coefficients and split them into `up` and `down`.
  '''
  highest_abs_expr = expr.loc[expr.abs().sort_values('CD-coefficient', ascending=False)[:top_n].index]
  return type('UpDownGeneset', tuple(), dict(
    up=list(highest_abs_expr[highest_abs_expr > 0].dropna().index),
    down=list(highest_abs_expr[highest_abs_expr < 0].dropna().index),
  ))
