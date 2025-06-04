import rpy2.robjects as ro
ro.r('''
install.packages("R.utils", repos="https://cloud.r-project.org/")
install.packages("RCurl", repos="https://cloud.r-project.org/")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org/")

BiocManager::install("DESeq2")
BiocManager::install("limma")
BiocManager::install("statmod")
BiocManager::install("edgeR")
''')
