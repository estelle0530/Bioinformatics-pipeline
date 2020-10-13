# Useful Functions 

##  Expression Marker Selection  
This project introduces genetic variable selection by a two-fold selection process: 
### 1) Apply principle component analysis to select genes with hihgest loadgins in the top components that explain over 90% of the variance  
### 2) Apply lasso regression for selected features against the outcome of interest (categorical and continuous)

## Parallelized Computation 
This function applies a pre-defined function (regression for pairwise features of interest) to all the cores avaialable in the machine in order to execute non-sequential functions all at once.

## Manage Gene Set Enrichment Leading Edge Genes
This function is written to read in a compiled list of genesets for Gene Set Enrichment Analysis (https://www.gsea-msigdb.org/gsea/index.jsp) and include all member genes in the genesets. This is particularly useful for enrichment analysis on the gene level rather than the gene set level.
