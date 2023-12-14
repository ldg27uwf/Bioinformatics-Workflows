#DESeq2 workflow tutorial 
#Smooth muscle cell seq data
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")
BiocManager::install("DESeq2")

library(DESeq2)
library(airway)
library(tidyverse)

setwd("~/Desktop/Bioinformagician/")

#read in count data
data(airway)
airway

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

#read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)

#read in sample info
colData<- read.csv('sample_info.csv')

#match rownames in colData to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))
#in same order
all(colnames(counts_data) == rownames(colData))

#Construct a DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData =counts_data,
                       colData = colData,
                       design = ~dexamethasone)

#filtering to remove rows with low gene counts
keep <- rowSums(counts(dds))>= 10

dds <- dds[keep,]

#set factor level
dds$dexamethasone <-  relevel(dds$dexamethasone, ref = "untreated")

#collapse technical replicates

#run DESeq
dds <- DESeq(dds)

res <- results(dds)

#explore results
summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

#contrasts
resultsNames(dds)

# compare multiple levels
# results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# MA plot
plotMA(res)
