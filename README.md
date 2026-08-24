DMD vs WT RNA‑seq Differential Expression and GSEA (GSE312779)
This repository contains an R workflow for differential gene expression (DGE) analysis of the GSE312779 RNA‑seq dataset (Duchenne muscular dystrophy vs wild‑type), followed by Gene Ontology (GO) enrichment via GSEA using clusterProfiler.

Overview
Data source: GEO accession GSE312779 (bulk RNA‑seq, DMD vs WT, with treatment information).

Main goals:

Build a DESeq2 model with genotype, treatment, and their interaction.

Identify differentially expressed genes (DEGs) for DMD vs WT.

Perform GO‑GSEA on ranked log2 fold‑changes.

Key packages: DESeq2, apeglm, pheatmap, EnhancedVolcano, clusterProfiler, org.Hs.eg.db, enrichplot, tidyverse, GEOquery.
