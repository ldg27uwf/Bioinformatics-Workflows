# Read and normalize microarray data
# RMA method

# tidyverse, GEOQuery, affy

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("GEOquery")

# load libraries

library(affy)
library(GEOquery)
library(tidyverse)

# breast cancer dataset SPCA2 knockdown
# .cel files 

getGEOSuppFiles("GSE148537")

# untar files

untar("GSE148537/GSE148537_RAW.tar", exdir = 'data/')

# reading in .cel files
mca_data <- ReadAffy(celfile.path = "data/")

# rna normalization
norm_data <- rma(mca_data)

# ERROR package hgu133plus2cdf not installed
# ‘RSQLite’ 'AnnotationDbi' hgu133plus2cdf' had non zero exit status

BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133plus2cdf")

#get expression estimates
expressdata_norm <- as.data.frame(exprs(norm_data))

# map probe IDs to gene symbols
gse <- getGEO("GSE148537", GSEMatrix = TRUE)
view(gse)

# get feature data to get ID- gene symbol mapping
featureData(gse$GSE148537_series_matrix.txt.gz)

feat.data <- gse$GSE148537_series_matrix.txt.gz@featureData@data

# subset data 
feat.data <- feat.data[c(1,11)]

normal_express <- expressdata %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feat.data, by = 'ID')
