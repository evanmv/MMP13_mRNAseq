##Packages -----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")
BiocManager::install("apeglm")
install.packages("pheatmap")
BiocManager::install("regionReport")
BiocManager::install("goseq")
BiocManager::install("geneLenDataBase") #no luck. Doesn't support hg38
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db") #Needed for GOseq

library(Rsubread)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(regionReport)
library(geneLenDataBase)
library(goseq)
##featureCounts -----
fC <- featureCounts(files = c("alignment_output/B1_AlignedReads.sam","alignment_output/B2_AlignedReads.sam","alignment_output/B3_AlignedReads.sam", "alignment_output/C1_AlignedReads.sam","alignment_output/C2_AlignedReads.sam","alignment_output/C3_AlignedReads.sam"),
              annot.inbuilt = "hg38",
              useMetaFeatures = TRUE,
              isPairedEnd = TRUE,
              nthreads = 8)
View(fC$annotation)
View(fC$stat)
View(fC$counts)

##DESeq2 -----

#Change col names from feature counts to match study design file
colnames(fC$counts) <- sub("_AlignedReads.sam", "", colnames(fC$counts))
#Read in study design file
colData <- read.delim("studydesign.txt", row.names = 1)
?read.delim

#Double check for matching row/column names
all(rownames(colData) == colnames(fC_counts))

dds <- DESeqDataSetFromMatrix(countData = fC_counts,
                              colData = colData,
                              design = ~ group)
#Pre-filtering 
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >=10) >= smallestGroupSize
dds <- dds[keep, ]

#Start DGE analysis
dds <- DESeq(dds)
res <- results(dds)
res

#Shrinkage of effect size
resultsNames(dds)
resLFC <- lfcShrink(dds, 
                    coef = "group_Knockdown_vs_Control",
                    type = "apeglm")

resOrdered <- res[order(res$pvalue),] #Order by smallest p-value
summary(res)

res05 <- results(dds, alpha=0.05)
summary(res05)

resSig <- subset(resOrdered, padj < 0.05) #Subset of ordered dataset with significance less that 0.05

plotMA(resLFC, ylim = c(-2,2)) #plot log2 FCs over mean of normalized counts for all samples.
idx <- identify(res$baseMean, res$log2FoldChange) #Id row number of individual genes. Click in plot then click finish
rownames(res)[idx]

##PCA -----
#Variance stabilizing transformation
?vst

vsd <- vst(dds)
head(assay(vsd), 3)

#PCA 
#plotPCA(vsd, intgroup="group") #use returnData=TRUE to save pca data to object
pcaData <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- formatC(100 * attr(pcaData, "percentVar"))

#plot w/ GGplot
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) +
  coord_fixed() +
  theme_bw()

ggsave("PCA.png")

#Heatmap
select <- order(rowMeans(counts(dds, normalized=TRUE)),
                decreasing = TRUE)[1:20]

df <- colData

ntd <- normTransform(dds)
#Can use various transformations of the data normTransform (ntd) variance stabilizing transformation (vsd)
pheatmap(assay(vsd)[select,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, annotation_col = colData)

?pheatmap

##regionReport -----
#Generate report with regionReport

dir.create('DESeq2_out', showWarnings = FALSE, recursive = TRUE)
report <- DESeq2Report(dds, project = 'MMP13 knockdown in DOK', intgroup = 'group', outdir = 'DESeq2_out', output = 'index')


##GO analysis -----
supportedGeneIDs() #hg38 not supported

DEG <- as.vector(resSig@rownames) #Create character vector of DE genes
ALL <- as.vector(res@rownames) #Create character vector of all genes
gene.vector=as.integer(ALL%in%DEG) #Assign value 1 for genes present in DE genes
names(gene.vector)=ALL #Assign names from ALL
head(gene.vector)

pwf=nullp(gene.vector, "hg38","knownGene")
head(pwf)
GO.wall=goseq(pwf, "hg38","knownGene")
