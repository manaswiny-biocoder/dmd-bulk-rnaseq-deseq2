
# This script carries of differential gene expression on the project(GSE312779)



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




           
# match the rownames of metadata to colnames of counts

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


all(rownames(metadata) %in% colnames(counts))
all(rownames(metadata) == colnames(counts))

# Constructing DESeqDataset


library("DESeq2")
metadata$treatment <- factor(metadata$treatment)
metadata$genotype <- factor(metadata$genotype)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ genotype + treatment + genotype:treatment)
View(dds)

resultsNames(dds)

# set factor level
dds$genotype <- relevel(dds$genotype, ref = "WT")
dds$treatment <- relevel(dds$treatment, ref = "Vamorolone")



# filtering genes less then  a threshold read counts

keep <- rowSums(counts(dds)) >= 10
dds <-dds[keep,]
dds
View(dds)







# Run DESeq

dds <- DESeq(dds)
plotDispEsts(dds)


# DMD VS WT 
res1 <- results(dds, contrast = c("genotype", "DMD", "WT"))
summary(res1)



#Log fold change for visualization and ranking
plotMA(res1, ylim = c(-6,6))
reslfc <- lfcShrink(dds, coef = "genotype_DMD_vs_WT" , type = "apeglm")
library(apeglm)
reslfc
plotMA(reslfc, ylim = c(-6,6))


# Visualization of the results

library(pheatmap)
library(EnhancedVolcano)

rld <- rlogTransformation(dds, blind = FALSE)
head(assay(rld))
hist(assay(rld))

PCAA <- plotPCA(rld, intgroup = "genotype")
PCAA

EnhancedVolcano(res1,
                lab = rownames(res1),
                x = 'log2FoldChange',
                y = 'pvalue')

# Turn the result table of DEGS of samples DMD vs WT to dataframe for functional enrichment


res1_data <- as.data.frame(res1)
View(res1_data)

write.table(res1_data, "DMD_WT_res.csv" , sep = "\t", row.names = T)


library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

# GENE SET ENRICHMENT ANALYSIS

deg <- read.table("DMD_WT_res.csv")

deg

deg$gene <- rownames(deg)

deg <- deg[, c("gene" , setdiff(names(deg), "gene"))]

deg

rownames(deg) <- NULL

# rank the genes

deg_order_fc <- deg[order(-deg$log2FoldChange), ]

head(deg_order_fc)

# get logFC
 

logfc <- deg_order_fc$log2FoldChange
logfc

# make a named vector of logFC

names(logfc) <- deg_order_fc$gene
logfc

enrich_go_gsea <- gseGO(geneList = logfc, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, keyType ="SYMBOL", verbose= FALSE)

# save the enrichment table

enrich_go_gsea_df <- enrich_go_gsea@result

enrich_go_gsea@result

dotplot_enrich_go_gsea <- dotplot(enrich_go_gsea, showCategory = 10, orderBy="GeneRatio")

dotplot_enrich_go_gsea
