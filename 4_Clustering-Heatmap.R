##From previous work - histogram/heatmap -----
#Clean up a bit more 
##Packages -----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
install.packages("tidyverse")
BiocManager::install("annotate")
BiocManager::install("DESeq2")
install.packages("gt")
install.packages("gtExtras")
library(tidyverse)
library(limma) 
library(RColorBrewer)
library(gplots) 
library(DESeq2)
library(annotate)
library(org.Hs.eg.db)
library(gt)
library(gtExtras)

## Data -----
#Read in data files 

res_sig <- read.csv("res_sig.csv")
names(res_sig)[names(res_sig) == 'X'] <- 'GeneID'

#Get dds_nf from 2_DESeq2.R lines 8 to 43
vsd <- vst(dds_nf)

vsd_hmap <- assay(vsd) %>%
  as.data.frame()

geneSymbols <- getSYMBOL(rownames(vsd_hmap), data = 'org.Hs.eg.db')

vsd_hmap <- vsd_hmap %>%
  mutate(Symbol = geneSymbols, .before = 1)

vsd_hmap <- na.omit(vsd_hmap) #Elimate NAs from missing gene symbols

deg <- as.vector(res_sig$GeneID)
vsd_hmap_sig <- vsd_hmap[rownames(vsd_hmap) %in% deg,]
#Assign row names from 1st column containing geneSymbols and remove column
row.names(vsd_hmap_sig) <- vsd_hmap_sig[,1] 
vsd_hmap_sig[,1] <- NULL

#Narrow down DEG list to ca. 20 **update - 40 (20240920)

#resSig.nf.Narrow <- subset(resSig.nf, log2FoldChange <= -1 | log2FoldChange >= 1)
#resSig.nf.Narrow <- as_tibble(resSig.nf.Narrow) %>%
 # mutate(GeneID = resSig.nf.Narrow@rownames, .before = 1)

res_sig_08 <- subset(res_sig, log2FoldChange <= -0.8 | log2FoldChange >= 0.8) 
res_sig_08 <- as_tibble(res_sig_08)
geneSymbols_08 <- getSYMBOL(as.character(res_sig_08$GeneID), data = 'org.Hs.eg.db')
res_sig_08 <- res_sig_08 %>% 
  mutate(Symbol = geneSymbols_08, .before = 1) %>%
  mutate(Change = ifelse(log2FoldChange > 0, 'Up', 'Down'))

vsd_hmap_sig_08 <- vsd_hmap[rownames(vsd_hmap) %in% res_sig_08$GeneID,] #match resSig.nf.Narrow

row.names(vsd_hmap_sig_08) <- vsd_hmap_sig_08[, 1] 
vsd_hmap_sig_08[, 1] <- NULL
  
write_csv(res_sig_08, "2024.11.13_DEGs.csv")

#table of DEGs 
#Define myheatcolors first
myheatcolors <- brewer.pal(name="RdBu", n=11)
DEGtable <- res_sig_08 %>% 
  select(Symbol, log2FoldChange, padj, Change) %>%
  gt(groupname_col = "Change") %>%
  tab_header(
    title = md("**DEGs in response to MMP13 knockdown in DOK**"),
    subtitle = "Fold change cutoff +/- 0.8"
  ) %>%
  gt_color_rows(
    columns = log2FoldChange,
    domain = c(-2, 2),
    palette = rev(myheatcolors)
  ) %>%
  tab_row_group(
    label = md("**Up**"),
    rows = res_sig_08$Change == "Up"
  ) %>%
  tab_row_group(
    label = md("**Down**"),
    rows = res_sig_08$Change == "Down"
  ) %>%
  row_group_order(groups = c("**Up**", "**Down**"))
  
gtsave(DEGtable, "tableDEGs2.png") #gtsave doesn't work on fox (no chrome)

## Cluster DEGs ----
#Narrowed non-filtered DEGs (FC 0.8, padj 0.05)
clust_rows <- hclust(as.dist(1-cor(t(vsd_hmap_sig_08), method="pearson")), method="complete") 
clust_columns <- hclust(as.dist(1-cor(vsd_hmap_sig_08, method="spearman")), method="complete") #cluster columns by spearman correlation
module_assign <- cutree(clust_rows, k=2)
module_color <- rainbow(length(unique(module_assign)), start=0.1, end=0.9) 
module_color <- module_color[as.vector(module_assign)] 

## Heatmap -----

myheatcolors <- brewer.pal(name="RdBu", n=11)
#Narrow, non-filtered...
options(bitmapType = 'cairo')
heatmap_nf <- heatmap.2(as.matrix(vsd.nf.hmap.sigNarrow), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row',
          density.info="none", trace="none",  
          cexRow=1, cexCol=2, margins = c(4,6),
          key = TRUE, keysize = 1,
)
 
dev.off()
#Heatmap with narrowed DEG list

heatmap.2(as.matrix(vsd_hmap_sig_08), 
          Rowv=as.dendrogram(clust_rows), 
          Colv=as.dendrogram(clust_columns),
          RowSideColors=module_color,
          col=rev(myheatcolors), scale='row',
          density.info="none", trace="none",  
          cexRow=1, cexCol=2, margins = c(4,8),
          key = TRUE, keysize = 2,
)

#Played around with keysize and cexCol (column label size). This is the best so far. 

#Heatmap with full diff gene set
clustRowsNF <- hclust(as.dist(1-cor(t(vsd.nf.hmap.sig), method="pearson")), method="complete") 
clustColumnsNF <- hclust(as.dist(1-cor(vsd.nf.hmap.sig, method="spearman")), method="complete") #cluster columns by spearman correlation
module.assignNF <- cutree(clustRowsNF, k=2)
module.colorNF <- rainbow(length(unique(module.assignNF)), start=0.1, end=0.9) 
module.colorNF <- module.colorNF[as.vector(module.assignNF)] 

heatmap.2(as.matrix(vsd.nf.hmap.sig), 
          Rowv=as.dendrogram(clustRowsNF), 
          Colv=as.dendrogram(clustColumnsNF),
          RowSideColors=module.colorNF,
          col=rev(myheatcolors), scale='row',
          density.info="none", trace="none",
          labRow = FALSE,
          cexRow=1, cexCol=2, margins = c(4,2),
          key = TRUE, keysize = 1,
)
?heatmap.2
