''' This module contains functions for differential expression analysis
'''

from maayanlab_bioinformatics.dge.characteristic_direction import characteristic_direction, up_down_from_characteristic_direction
from maayanlab_bioinformatics.dge.limma_voom import limma_voom_differential_expression, up_down_from_limma_voom
