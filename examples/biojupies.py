#%%[markdown]

# Here we walk through building a biojupies notebook using methods in this library for the purpose of integration testing.

#%%
import sys; sys.path.insert(0, '..')
import os
import numpy as np
import pandas as pd
from IPython.display import display
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import plotly.express as px
from maayanlab_bioinformatics.normalization import filter_by_var, cpm_normalize, zscore_normalize, log10_normalize
from maayanlab_bioinformatics.dge import limma_voom_differential_expression
from maayanlab_bioinformatics.utils import merge
os.environ['R_LIBS_USER'] = '/home/u8sand/.r_libs'

#%%
# biojupies settings
pca_top_genes = 2500
normalization = 'logCPM'
zscore = True
pval_thresh = 0.05
logFC_thresh = 1.5
enrichr_geneset_size = 500
sort_genes_by = 't'

#%%
# ## Load Dataset
df_data = pd.read_csv('../tests/test_example_matrix.txt', sep='\t', index_col=0)
display(df_data.head())

df_metadata = pd.read_csv('../tests/test_example_metadata.txt', sep='\t', index_col=0)
display(df_metadata.head())

#%%
# Biojupies selection
normal = df_data.loc[:, df_metadata['cell type'] == 'normal melanocytes']
perturbation = df_data.loc[:, df_metadata['cell type'] == 'melanoma cell line']

#%%
# ## PCA
# Principal Component Analysis was performed using the PCA function from in the sklearn Python module.
# Prior to performing PCA, the raw gene counts were normalized using the logCPM method
df_data_norm_for_pca = log10_normalize(cpm_normalize(df_data))
# filtered by selecting the 2500 genes with most variable expression
df_data_norm_for_pca = filter_by_var(df_data_norm_for_pca, top_n=2500)
# and finally transformed using the Z-score method.
df_data_norm_for_pca = zscore_normalize(df_data_norm_for_pca)

pca = PCA()
pca.fit(df_data_norm_for_pca.T.values)
df_pca = pd.DataFrame(
  pca.transform(df_data_norm_for_pca.T.values),
  index=df_data_norm_for_pca.columns,
  columns=[
    f"PC{n}({var*100:2.1f}% var. explained)"
    for n, var in enumerate(pca.explained_variance_ratio_)
  ]
)
px.scatter_3d(
  merge(df_pca, df_metadata),
  x=df_pca.columns[0],
  y=df_pca.columns[1],
  z=df_pca.columns[2],
  color='cell type',
  hover_data=[df_pca.index],
)

#%%
# ## clustergrammer

#%%
# ## library size analysis

#%%
# ## Differential Expression Table
dge_table = limma_voom_differential_expression(
  normal,
  perturbation,
  filter_genes=True,
  voom_design=False,
)
dge_table.sort_values('P.Value')

