
# Modify metadata and counts file to match the rownames of metadata and colnames of counts for deseq2 

# Load the libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(GEOquery)

# Read tab-delimited count file
counts <- read.csv("GSE312779_dmd_wt_all_raw_counts.csv.gz")
View(counts)

#  get the metadata using GEOquery package 
gse <- getGEO(GEO = "GSE312779" ,GSEMatrix = TRUE)
length(gse)
gse[1]
metadata <- pData(phenoData(gse[[1]]))

# modify counts matrix
counts$X <- make.unique(as.character(counts$X))
rownames(counts) <- counts[, 1]
counts <- counts[, -1]

# modify metadata
metadata <- metadata %>% 
  select(c("description", "cell line:ch1","genotype:ch1","treatment:ch1")) %>% 
  mutate(description = gsub("Library name: " , "" , description)) 
colnames(metadata) <- gsub(":ch1", "", colnames(metadata))
rownames(metadata) <- metadata[, 1]
metadata <- metadata[, -1]



# match the rownames of metadata to colnames of counts
all(rownames(metadata) %in% colnames(counts))
all(rownames(metadata) == colnames(counts))





