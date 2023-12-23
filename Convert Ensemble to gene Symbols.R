# Convert Ensemble IDs to gene symbols

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("EnsDb.Hsapiens.v86", force = TRUE)

# annotables not compatible with this version of R

install.packages("devtools")
devtools::install_github("stephenturner/annotables", force = TRUE)
library(biomaRt)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(devtools)
library(annotables)
library(tidyverse)

# input list of Ensembl ID's
ensembl.ids <- read.table("~/Desktop/Bioinformagician/Bioinformatics_Lessons/GENE_ID_LIST", quote="\"", comment.char="", header = F)
View(GENE_ID_LIST)

# method 1: BiomaRt
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
view(datasets)

ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

# query for gene ID and gene name
getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
      filters = "ensembl_gene_id",
      values = ensembl.ids$V1,
      mart = ensembl.con
      )

# method 2: annotables

grch38 %>%
  filter(ensgene %in% ensembl.ids$V1)

# method 3: annotation DBs

keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

mapIds(org.Hs.eg.db,
       keys = ensembl.ids$V1,
       keytype = 'ENSEMBL',
       column = 'SYMBOL')

keytypes(EnsDb.Hsapiens.v86)

columns(EnsDb.Hsapiens.v86)

mapIds(EnsDb.Hsapiens.v86,
       keys = ensembl.ids$V1,
       keytype = 'GENEID',
       column = 'SYMBOL')
