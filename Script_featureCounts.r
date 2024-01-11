##Packages -----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread")
BiocManager::install("apeglm")
library(Rsubread)
library(DESeq2)

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
all(rownames(colData) == colnames(fC$counts))

dds <- DESeqDataSetFromMatrix(countData = fC$counts,
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

resOrdered <- res[order(res$pvalue),]
summary(res)
