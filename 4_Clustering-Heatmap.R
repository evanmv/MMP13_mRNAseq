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

vsd.nf.hmap <- assay(vsd.nf) %>%
  as.data.frame()
geneSymbols <- getSYMBOL(vsd.nf.hmap$GeneID, data = 'org.Hs.eg.db')

vsd.nf.hmap <- vsd.nf.hmap %>%
  mutate(Symbol = geneSymbols, .before = 1)

vsd.nf.hmap <- na.omit(vsd.nf.hmap) #Elimate NAs from missing gene symbols

DEG.nf <- as.vector(resSig.nf@rownames)
vsd.nf.hmap.sig <- vsd.nf.hmap[rownames(vsd.nf.hmap) %in% DEG.nf,]
#Assign row names from 1st column containing geneSymbols and remove column
row.names(vsd.nf.hmap.sig) <- vsd.nf.hmap.sig[,1] 
vsd.nf.hmap.sig[,1] <- NULL

#Narrow down DEG list to ca. 20 **update - 40 (20240920)

#resSig.nf.Narrow <- subset(resSig.nf, log2FoldChange <= -1 | log2FoldChange >= 1)
#resSig.nf.Narrow <- as_tibble(resSig.nf.Narrow) %>%
 # mutate(GeneID = resSig.nf.Narrow@rownames, .before = 1)

resSig.nf.Narrow <- subset(resSig.nf, log2FoldChange <= -0.8 | log2FoldChange >= 0.8) 
resSig.nf.Narrow <- as_tibble(resSig.nf.Narrow) %>%
  mutate(GeneID = resSig.nf.Narrow@rownames, .before = 1)
geneSymbols.nf <- getSYMBOL(resSig.nf.Narrow$GeneID, data = 'org.Hs.eg.db')
resSig.nf.Narrow <- resSig.nf.Narrow %>% 
  mutate(Symbol = geneSymbols.nf, .before = 1) %>%
  mutate(Change = ifelse(log2FoldChange > 0, 'Up', 'Down'))

vsd.nf.hmap.sigNarrow <- vsd.nf.hmap[rownames(vsd.nf.hmap) %in% resSig.nf.Narrow$GeneID,] #match resSig.nf.Narrow

row.names(vsd.nf.hmap.sigNarrow) <- vsd.nf.hmap.sigNarrow[,1] 
vsd.nf.hmap.sigNarrow[,1] <- NULL
  
write_csv(resSig.nf.Narrow, "2024.09.20_DEGs.csv")

#table of DEGs 

DEGtable <- resSig.nf.Narrow %>% 
  select(Symbol, log2FoldChange, padj, Change) %>%
  gt(groupname_col = "Change") %>%
  tab_header(
    title = md("**DEGs in response to MMP13 knockdown in DOK**"),
    subtitle = "Fold change cutoff +/- 0.8"
  ) %>%
  gt_color_rows(
    columns = log2FoldChange,
    domain = c(-2, 2),
    palette = myheatcolors
  ) %>%
  tab_row_group(
    label = md("**Up**"),
    rows = resSig.nf.Narrow$Change == "Up"
  ) %>%
  tab_row_group(
    label = md("**Down**"),
    rows = resSig.nf.Narrow$Change == "Down"
  ) %>%
  row_group_order(groups = c("**Up**", "**Down**"))
  
gtsave(DEGtable, "tableDEGs2.png") #gtsave doesn't work on fox (no chrome)

## Cluster DEGs ----
#Narrowed non-filtered DEGs (FC 0.8, padj 0.05)
clustRows <- hclust(as.dist(1-cor(t(vsd.nf.hmap.sigNarrow), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(vsd.nf.hmap.sigNarrow, method="spearman")), method="complete") #cluster columns by spearman correlation
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

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

heatmap.2(as.matrix(vsd.hmap.sigNarrow), 
          Rowv=as.dendrogram(clustRowsN), 
          Colv=as.dendrogram(clustColumnsN),
          RowSideColors=module.colorN,
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
