''' This module contains functions for differential expression analysis
'''

import lazy_loader as lazy

__getattr__, __dir__, __all__ = lazy.attach(
    __name__,
    submod_attrs={
        'characteristic_direction': ['characteristic_direction', 'up_down_from_characteristic_direction'],
        'limma_voom': ['limma_voom_differential_expression', 'limma_voom_differential_expression_design', 'up_down_from_limma_voom'],
        'deseq2': ['deseq2_differential_expression'],
        'ttest': ['ttest_differential_expression'],
        'logfc': ['logfc_differential_expression'],
    }
)
