# Ma'ayanlab Bioinformatics

A collection of useful functions for bioinformatics data analysis.

This library contains many functions and methods I use again and again in different analyses including:
- quantile normalization
- other common normalizations
- logcpm, zscore, filter variance
- gmt parser
- single-cell sparse matrix parsing
- transcript to gene conversions
- ...
- various dge including chdir, lima_voom, etc..

## Installation
```
pip install git+https://github.com/Maayanlab/maayanlab-bioinformatics.git

# [OPTIONAL] for some functionality like limma_voom & filter_by_expr
R -e "source('setup.R')"
```
