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
BiocManager::install("clusterProfiler") #Alternative for GO/overrepresentation analysis

library(Rsubread)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(regionReport)
library(geneLenDataBase)
library(goseq)
library(GO.db)
library(clusterProfiler)
library(tidyverse)

install.packages('devtools')
require(devtools)
install_version("rvcheck", version = "0.1.8", repos = "http://cran.us.r-project.org")
library(rvcheck)
?patchwork
??rvcheck
install.packages("patchwork")
##featureCounts -----
fC <- featureCounts(files = c("alignment_output/B1_AlignedReads.sam","alignment_output/B2_AlignedReads.sam","alignment_output/B3_AlignedReads.sam", "alignment_output/C1_AlignedReads.sam","alignment_output/C2_AlignedReads.sam","alignment_output/C3_AlignedReads.sam"),
              annot.inbuilt = "hg38",
              useMetaFeatures = TRUE,
              isPairedEnd = TRUE,
              nthreads = 8)
View(fC$annotation)
View(fC$stat)
View(fC$counts)

fC_counts["4322",] #MMP13 counts
fC_counts["22806",]
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
?DESeq
#Shrinkage of effect size
resultsNames(dds)
resLFC <- lfcShrink(dds, 
                    coef = "group_Knockdown_vs_Control",
                    type = "apeglm")

resOrdered <- res[order(res$pvalue),] #Order by smallest p-value
summary(resOrdered)

res05 <- results(dds, alpha=0.05)
summary(res05)

resSig <- subset(resOrdered, padj < 0.05) #Subset of ordered dataset with significance less that 0.05

plotMA(resLFC, ylim = c(-2,2)) #plot log2 FCs over mean of normalized counts for all samples.
idx <- identify(res$baseMean, res$log2FoldChange) #Id row number of individual genes. Click in plot then click finish
rownames(res)[idx]

resSigDown <- subset(resSig, log2FoldChange <= -0.6)
resSigUp <- subset(resSig, log2FoldChange >=0.6)
#Down-regulated genes (significantly) with less than -0.6 Log2 FC
SigDownTidy <- as_tibble(resSigDown) %>%
  mutate(GeneID = resSigDown@rownames, .before = 1)
#Up-regulated genes (sig) with more than 0.6 Log2 FC
sigUpTidy <- as_tibble(resSigUp) %>%
  mutate(GeneID = resSigUp@rownames, .before = 1)
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
  geom_point(size=5) +
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) +
  coord_fixed() +
  theme_bw()

ggsave("PCA_5pt.png")

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


##GO analysis w. goSEQ -----
supportedGeneIDs() #hg38 not supported

#Get data formatted correctly
DEG <- as.vector(resSig@rownames) #Create character vector of DE genes
ALL <- as.vector(res@rownames) #Create character vector of all genes
gene.vector=as.integer(ALL%in%DEG) #Assign value 1 for genes present in DE genes
names(gene.vector)=ALL #Assign names from ALL
head(gene.vector)


DEG2 <- as.vector(c(sigDownTidy$GeneID, sigUpTidy$GeneID))
#GOseq
pwf=nullp(gene.vector, "hg38","knownGene")
head(pwf)
GO.wall=goseq(pwf, "hg38","knownGene")
class(GO.wall)
head(GO.wall)
nrow(GO.wall)

#Filter for significance
enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method='BH')<.05]
head(enriched.GO)

for (go  in enriched.GO[1:10]) {
  print(GOTERM[[go]])
  cat("-----------------------------------------\n")
}

#To save txt file 

capture.output(for (go in enriched.GO[1:100]) {
  print(GOTERM[[go]])
  cat("-----------------------------------------\n")
}
, file = "SigGo.txt")
  
#Molecular function pathways
GO.MF=goseq(pwf, "hg38","knownGene",test.cats=c("GO:MF"))
enriched.GO.mf=GO.MF$category[p.adjust(GO.MF$over_represented_pvalue, method='BH')<0.5]

capture.output(for (go in enriched.GO.mf[1:47]) {
  print(GOTERM[[go]])
  cat("-----------------------------------------\n")
}
, file = "SigGoMF.txt")

#Biological processes pathways
GO.BP=goseq(pwf, "hg38","knownGene",test.cats=c("GO:BP"))
enriched.GO.bp=GO.MF$category[p.adjust(GO.BP$over_represented_pvalue, method='BH')<0.5]

capture.output(for (go in enriched.GO.bp[1:100]) {
  print(GOTERM[[go]])
  cat("-----------------------------------------\n")
}
, file = "SigGoBP.txt")


## ClusterProfiler -----

?clusterProfiler
??clusterProfiler

?enrichGO()
enrich_GO <- enrichGO(
  gene = DEG,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = ALL,
  qvalueCutoff = 0.5,
  readable = TRUE
)
  
enrich_GO_tidy <- enrich_GO %>%
  slot("result") %>%
  tibble::as.tibble()

dotplot(enrich_GO)

ggsave("enrichGO_dotplot.png")

pairwise_termsim(enrich_GO) #?
emapplot(enrich_GO)
DEG

#Try enrichGO with restricted DEG list
#MF
enrich_GO_MF <- enrichGO(
  gene = DEG2,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = ALL,
  qvalueCutoff = 0.5,
  readable = TRUE
)

GO_MF.t <- enrich_GO_MF %>%
  slot("result") %>%
  tibble::as.tibble()

dotplot(enrich_GO_MF,
        title = "GO terms MF")
ggsave("GO.MF_enrichmentDotPlot.png")

MFsimilarityMatrix <- as.tibble(pairwise_termsim(enrich_GO_MF))

emapplot(pairwise_termsim(enrich_GO_MF))
ggsave("GO.MF_enrichmentPlot.png")

#CC

enrich_GO_CC <- enrichGO(
  gene = DEG2,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "CC",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = ALL,
  qvalueCutoff = 0.5,
  readable = TRUE
)

GO_CC.t <- enrich_GO_CC %>%
  slot("result") %>%
  tibble::as.tibble()

dotplot(enrich_GO_CC,
        title = "GO terms CC")
ggsave("GO.CC_enrichmentDotPlot.png")

CCsimilarityMatrix <- as.tibble(pairwise_termsim(enrich_GO_CC))

emapplot(pairwise_termsim(enrich_GO_CC))
ggsave("GO.CC_enrichmentPlot.png")

#BP 

enrich_GO_BP <- enrichGO(
  gene = DEG2,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = ALL,
  qvalueCutoff = 0.5,
  readable = TRUE
)

GO_BP.t <- enrich_GO_BP %>%
  slot("result") %>%
  tibble::as.tibble()

dotplot(enrich_GO_BP,
        title = "GO terms BP")
ggsave("GO.BP_enrichmentDotPlot.png")

BPsimilarityMatrix <- as.tibble(pairwise_termsim(enrich_GO_BP))

emapplot(pairwise_termsim(enrich_GO_BP))
ggsave("GO.BP_enrichmentPlot.png")

#combine 
write.csv(BPsimilarityMatrix, 
          file = "BPsimilarityMatrix.csv")
write.csv(MFsimilarityMatrix,
          file = "MFsimilarityMatrix.csv")
write.csv(CCsimilarityMatrix,
          file = "CCsimilarityMatrix.csv")

write.csv(GO_BP.t,
          file = "GO_BP_results.csv")
write.csv(GO_MF.t,
          file = "GO_MF_results.csv")
write.csv(GO_CC.t,
          file = "GO_CC_results.csv")
##Assign gene symbols from Entrez IDs -----
install.packages("annotate")
library(annotate)
library(org.Hs.eg.db)
?getSYMBOL()
downSymbols <- getSYMBOL(SigDownTidy$GeneID, data = 'org.Hs.eg.db')

sigDownTidy <- 
  SigDownTidy %>%
  mutate(Symbol = downSymbols, .before=1)

upSymbols <- getSYMBOL(sigUpTidy$GeneID, data = 'org.Hs.eg.db')

sigUpTidy <- 
  sigUpTidy %>%
  mutate(Symbol = upSymbols, .before=1)

write.csv(sigUpTidy, file = "upDEGs.csv")
write.csv(sigDownTidy, file = "downDEGs.csv")
library(readxl)

#?
DEGs <- as_tibble(DEG)
allGenes <- as_tibble(ALL)
allGenes <- as_tibble(resOrdered) %>%
  mutate(GeneID = resOrdered@rownames, .before = 1)
fC_counts["4322",]
