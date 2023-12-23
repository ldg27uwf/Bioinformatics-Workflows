# Convert Ensemble IDs to gene symbols

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("EnsDb.Hsapiens.v86", force = TRUE)

# annotables not compatible with this version of R

install.packages("devtools")
devtools::install_github("stephenturner/annotables")
library(biomaRt)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(devtools)
library(annotables)
library(tidyverse)

# input list of Ensembl ID's
GENE_ID_LIST <- read.table("~/Desktop/Bioinformagician/Bioinformatics_Lessons/GENE_ID_LIST", quote="\"", comment.char="")
View(GENE_ID_LIST)

# method 1: BiomaRt
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
