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
pip3 install --user --upgrade git+https://github.com/Maayanlab/maayanlab-bioinformatics.git
```

If you want to use limma_voom, or some other R only functions, you'll need to install the relevant R dependencies as well.. see setup.R.
