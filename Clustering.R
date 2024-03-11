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

vsd.hmap <- assay(vsd) %>%
  as.data.frame()
geneSymbols <- getSYMBOL(rownames(vsd.hmap), data = 'org.Hs.eg.db')

vsd.hmap <- vsd.hmap %>%
  mutate(Symbol = geneSymbols, .before = 1)

vsd.hmap <- na.omit(vsd.hmap) #Elimate NAs from missing gene symbols

vsd.hmap.sig <- vsd.hmap[rownames(vsd.hmap) %in% DEG2,]
#Assign row names from 1st column containing geneSymbols and remove column
row.names(vsd.hmap.sig) <- vsd.hmap.sig[,1] 
vsd.hmap.sig[,1] <- NULL

head(vsd.hmap.sig)

#Narrow down DEG list to ca. 20

resSigNarrow <- subset(resSig, log2FoldChange <= -0.8 | log2FoldChange >= 0.8) 
resSigTidy <- as_tibble(resSigNarrow) %>%
  mutate(GeneID = resSigNarrow@rownames, .before = 1)

vsd.hmap.sigNarrow <- vsd.hmap[rownames(vsd.hmap) %in% resSigTidy$GeneID,]

row.names(vsd.hmap.sigNarrow) <- vsd.hmap.sigNarrow[,1] 
vsd.hmap.sigNarrow[,1] <- NULL
  
?subset

## Cluster DEGs ----

clustRows <- hclust(as.dist(1-cor(t(vsd.hmap.sig), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(vsd.hmap.sig, method="spearman")), method="complete") #cluster columns by spearman correlation
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

heatmap.2(as.matrix(vsd.hmap.sig), 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=2, margins = c(4,4),
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
