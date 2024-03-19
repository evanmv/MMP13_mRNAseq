##From previous work - histogram/heatmap -----
##Packages -----
library(tidyverse)
library(limma) 
library(RColorBrewer)
library(gplots) 
library(DESeq2)
library(annotate)
library(org.Hs.eg.db)

## Data -----

vsd.nf.hmap <- assay(vsd.nf) %>%
  as.data.frame()
geneSymbols <- getSYMBOL(rownames(vsd.nf.hmap), data = 'org.Hs.eg.db')

vsd.nf.hmap <- vsd.nf.hmap %>%
  mutate(Symbol = geneSymbols, .before = 1)

vsd.nf.hmap <- na.omit(vsd.nf.hmap) #Elimate NAs from missing gene symbols

DEG.nf <- as.vector(resSig.nf@rownames)
vsd.nf.hmap.sig <- vsd.nf.hmap[rownames(vsd.nf.hmap) %in% DEG.nf,]
#Assign row names from 1st column containing geneSymbols and remove column
row.names(vsd.nf.hmap.sig) <- vsd.nf.hmap.sig[,1] 
vsd.nf.hmap.sig[,1] <- NULL

#Narrow down DEG list to ca. 20

resSig.nf.Narrow <- subset(resSig.nf, log2FoldChange <= -1 | log2FoldChange >= 1)
resSig.nf.Narrow <- as_tibble(resSig.nf.Narrow) %>%
  mutate(GeneID = resSig.nf.Narrow@rownames, .before = 1)

#resSigNarrow <- subset(resSig, log2FoldChange <= -0.8 | log2FoldChange >= 0.8) 
#resSigTidy <- as_tibble(resSigNarrow) %>%
 # mutate(GeneID = resSigNarrow@rownames, .before = 1)

vsd.nf.hmap.sigNarrow <- vsd.nf.hmap[rownames(vsd.nf.hmap) %in% resSig.nf.Narrow$GeneID,]

row.names(vsd.nf.hmap.sigNarrow) <- vsd.nf.hmap.sigNarrow[,1] 
vsd.nf.hmap.sigNarrow[,1] <- NULL
  
?subset

## Cluster DEGs ----
#Narrowed non-filtered DEGs (FC 1, padj 0.05)
clustRows <- hclust(as.dist(1-cor(t(vsd.nf.hmap.sigNarrow), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(vsd.nf.hmap.sigNarrow, method="spearman")), method="complete") #cluster columns by spearman correlation
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

#Cluster with narrowed DEG set

clustRowsN <- hclust(as.dist(1-cor(t(vsd.hmap.sigNarrow), method="pearson")), method="complete") 
clustColumnsN <- hclust(as.dist(1-cor(vsd.hmap.sigNarrow, method="spearman")), method="complete") #cluster columns by spearman correlation
module.assignN <- cutree(clustRowsN, k=2)
module.colorN <- rainbow(length(unique(module.assignN)), start=0.1, end=0.9) 
module.colorN <- module.colorN[as.vector(module.assignN)] 


## Heatmap -----

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
#Narrow, non-filtered...
heatmap.2(as.matrix(vsd.nf.hmap.sigNarrow), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors), scale='row',
          density.info="none", trace="none",  
          cexRow=1, cexCol=2, margins = c(4,6),
          key = TRUE, keysize = 2,
)

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
