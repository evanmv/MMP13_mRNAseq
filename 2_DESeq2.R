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
res_df <- as.data.frame(res_nf)
write_csv(res_df, "res_nf.csv")

##regionReport -----
#Generate report with regionReport
#Problem with regionReport?
dir.create('DESeq2_nf_out', showWarnings = FALSE, recursive = TRUE)
report <- DESeq2Report(dds_nf, project = 'MMP13 knockdown in DOK', intgroup = 'group', res = 'res_nf', outdir = 'DESeq2_nf_out', output = 'index')
??DESeq2Report

#Shrinkage of effect size with apeglm
resultsNames(dds_nf)
res_lfc <- lfcShrink(dds_nf, 
                    coef = "group_Knockdown_vs_Control",
                    type = "apeglm")
summary(res_lfc)

plotMA(res_lfc, ylim = c(-2,2)) #plot log2 FCs over mean of normalized counts for all samples.
idx <- identify(res_nf$baseMean, res_nf$log2FoldChange) #Id row number of individual genes. Click in plot then click finish
rownames(res_nf)[idx] #Helpful after fetching gene names

#Subset "significant" results
res_ord <- res_nf[order(res_nf$pvalue),]
summary(res_ord)
res_sig <- subset(res_ord, padj < 0.05) #Subset of ordered dataset with significance less than 0.05
summary(res_sig)
res_sig_df <- as.data.frame(res_sig)
write_csv(res_sig_df, "res_sig.csv")

##PCA -----
#Variance stabilizing transformation

vsd <- vst(dds_nf)
head(assay(vsd), 3)

#PCA 
#plotPCA(vsd, intgroup="group") #use returnData=TRUE to save pca data to object
pca_data <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percent_var <- formatC(100 * attr(pca_data, "percentVar"))
?plotPCA
#plotPCA#plot w/ GGplot
PCA <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
  geom_point(size=5) +
  xlab(paste0("PC1 (",percent_var[1],"%)")) +
  ylab(paste0("PC2 (",percent_var[2],"%)")) +
  coord_fixed() +
  theme_bw()

ggsave("PCA.png")

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
#Uses top 500 variable genes, same as plotPCA from DESeq2
ntop <- 500 
rv <- rowVars(assay(vsd), useNames = T)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(assay(vsd)[select, ] )

pca_prcomp <- prcomp(mat)
pca <- as.data.frame(pca_prcomp$x)
pca$PC1 #Check that PC1 matches plotPCA
pca_data$PC1
?prcomp

#To get % variance
summ <- summary(pca_prcomp)
percent_var_prcomp <- formatC(100 * summ$importance[2, ])

#Top 4 principal components
pca_res <- pca[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = rownames(col_data),
             group = col_data$group) 
#Add percent variance to colnames
names(pca_res)[names(pca_res) == 'PC1'] <- paste0("PC1 (",percent_var_prcomp[1],"%)")
names(pca_res)[names(pca_res) == 'PC2'] <- paste0("PC2 (",percent_var_prcomp[2],"%)")
names(pca_res)[names(pca_res) == 'PC3'] <- paste0("PC3 (",percent_var_prcomp[3],"%)")
names(pca_res)[names(pca_res) == 'PC4'] <- paste0("PC4 (",percent_var_prcomp[4],"%)")


pca_pivot <- pivot_longer(pca_res, # dataframe to be pivoted
                          cols = 1:4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)
??Magrittr
ggplot(pca_pivot) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot") +
  theme_bw() +
  coord_flip()

ggsave('smallMultiples.png')
?geom_bar
?labs
?facet_wrap
