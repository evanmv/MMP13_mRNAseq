##DESeq2 workflow. Pre-filtering and DGE analysis, w/o independent filtering. PCA analysis. 

##Packages -----
BiocManager::install("regionReport")
BiocManager::install("apeglm")
install.packages("pheatmap")

library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(regionReport)
library(tidyverse)

##DESeq2 -----
#Read in featureCounts data, remove lengths column
fc_counts <- read_tsv("counts.txt") %>%
  mutate(Length = NULL) 
#Change to dataframe, tidyverse didn't want to work to reassign column to rownames which is needed for DESeq2
fc_counts <- as.data.frame(fc_counts)

#Reassign GeneID as rownames and remove column
rownames(fc_counts) <- fc_counts[, 1]
fc_counts <- fc_counts[, -1]
  
#Change col names from feature counts to match study design file
colnames(fc_counts) <- sub("_AlignedReads.sam", "", colnames(fc_counts))
#Read in study design file
col_data <- read.delim("studydesign.txt", row.names = 1)

#Double check for matching row/column names
all(rownames(col_data) == colnames(fc_counts))

dds <- DESeqDataSetFromMatrix(countData = fc_counts,
                              colData = col_data,
                              design = ~ group)
#Pre-filtering 
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >=10) >= smallestGroupSize
dds <- dds[keep, ]

#DGE analysis 
dds_nf <- DESeq(dds)
res_nf <- results(dds_nf, independentFiltering=FALSE)
res_nf

##regionReport -----
#Generate report with regionReport
#Problem with regionReport?
dir.create('DESeq2_nf_out', showWarnings = FALSE, recursive = TRUE)
report <- DESeq2Report(dds_nf, project = 'MMP13 knockdown in DOK', intgroup = 'group', res = 'res.nf', outdir = 'DESeq2.nf_out', output = 'index')
??DESeq2Report

#Shrinkage of effect size with apeglm
resultsNames(dds_nf)
res_lfc <- lfcShrink(dds_nf, 
                    coef = "group_Knockdown_vs_Control",
                    type = "apeglm")
summary(res_lfc)

plotMA(res_lfc, ylim = c(-2,2)) #plot log2 FCs over mean of normalized counts for all samples.
idx <- identify(res$baseMean, res$log2FoldChange) #Id row number of individual genes. Click in plot then click finish
rownames(res)[idx]

#Subset "significant" results
resOrd.nf <- res.nf[order(res.nf$pvalue),]
summary(resOrd.nf)
resSig.nf <- subset(resOrd.nf, padj < 0.05) #Subset of ordered dataset with significance less that 0.05
summary(resSig.nf)

##PCA -----
#Variance stabilizing transformation

vsd.nf <- vst(dds.nf)
head(assay(vsd.nf), 3)

#PCA 
#plotPCA(vsd, intgroup="group") #use returnData=TRUE to save pca data to object
pcaData.nf <- plotPCA(vsd.nf, intgroup="group", returnData=TRUE)
percentVar.nf <- formatC(100 * attr(pcaData.nf, "percentVar"))

#plot w/ GGplot
PCA_nf <- ggplot(pcaData.nf, aes(PC1, PC2, color=group)) +
  geom_point(size=5) +
  xlab(paste0("PC1 (",percentVar.nf[1],"%)")) +
  ylab(paste0("PC2 (",percentVar.nf[2],"%)")) +
  coord_fixed() +
  theme_bw()

ggsave("PCA_nf.png")

#Small multiples plot



