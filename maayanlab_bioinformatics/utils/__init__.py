'''This module contains general utility functions for convenient analysis
'''

from maayanlab_bioinformatics.utils.describe import np_describe
from maayanlab_bioinformatics.utils.chunked import chunk_slices, chunk_applymap
from maayanlab_bioinformatics.utils.fetch_save_read import fetch_save_read
from maayanlab_bioinformatics.utils.merge import merge
from maayanlab_bioinformatics.utils.sparse import sp_hdf_dump, sp_hdf_load
