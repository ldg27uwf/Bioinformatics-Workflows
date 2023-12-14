# script to manipulate gene expression data
setwd("~/Desktop/Bioinformagician/RNA-seq_data")

# load libraries

library(dplyr)
library(tidyverse)
library(GEOquery)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

# read the data
dat <- read.csv(file = "../RNA-seq_data/GSE183947_fpkm.csv")

# get metadata 
GSE <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)

metadata <- pData(phenoData(GSE[[1]]))
head(metadata)
colnames(metadata)

metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  #gsub("remove this", "replace with this", in this column)
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis)) 

head(dat) 

#reshape data to merge with metadata

dat.long <- dat %>%
  rename(gene = X) %>%
  #convert wide format to long format
  gather(key = 'samples',value = 'FPKM', -gene) 

#join dataframes 

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

#explore data 

dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  #Store mean FPKM in new column called mean_FPKM
  summarize(mean_FPKM = mean(FPKM), 
            median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM)


