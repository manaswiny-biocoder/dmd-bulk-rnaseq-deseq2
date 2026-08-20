# Constructing DESeqDataset and running DESeq2


library("DESeq2")


# set factor level
metadata$treatment <- factor(metadata$treatment)
metadata$genotype <- factor(metadata$genotype)

# set reference level
dds$genotype <- relevel(dds$genotype, ref = "WT")
dds$treatment <- relevel(dds$treatment, ref = "Vamorolone")


# construct DESeqDataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ genotype + treatment + genotype:treatment)
View(dds)
resultsNames(dds)



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
