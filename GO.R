##Packages -----
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
devtools::install_github("YuLab-SMU/enrichplot")
##GO analysis with clusterProfiler

#Get data formatted correctly
DEG.nf <- as.vector(resSig.nf@rownames) #Create character vector of DE genes
ALL.nf <- as.vector(res.nf@rownames) #Create character vector of all genes
DEG.nf.narrow <- as.vector(resSig.nf.Narrow$GeneID) #Highly differentially expressed
#Relax FC cutoff to 0.8 for GO analysis
resSig.nf.Narrow08 <- subset(resSig.nf, log2FoldChange <= -0.8 | log2FoldChange >= 0.8)
resSig.nf.Narrow08 <- as_tibble(resSig.nf.Narrow08) %>%
  mutate(GeneID = resSig.nf.Narrow08@rownames, .before = 1)
DEG.nf.narrow08 <- as.vector(resSig.nf.Narrow08$GeneID)
#GO_MF
enrich_GO_MFnf <- enrichGO(
  gene = DEG.nf,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = ALL.nf,
  qvalueCutoff = 0.5,
  readable = TRUE
)
head(enrich_GO_MFnf)

GO_MFnf.t <- enrich_GO_MFnf %>%
  slot("result") %>%
  tibble::as.tibble()
?dotplot
dotplot_GO.MF.nf <- dotplot(enrich_GO_MFnf,
        title = "GO terms MF")
ggsave("GO.MFnf_enrichmentDotPlot.png")

MFnfsimilarityMatrix <- as.tibble(pairwise_termsim(enrich_GO_MFnf))

emapplot(pairwise_termsim(enrich_GO_MFnf))
ggsave("GO.MFnf_enrichmentPlot.png")

#GO_CC
enrich_GO_CCnf <- enrichGO(
  gene = DEG.nf,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "CC",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = ALL.nf,
  qvalueCutoff = 0.5,
  readable = TRUE
)
head(enrich_GO_CCnf)

GO_CCnf.t <- enrich_GO_CCnf %>%
  slot("result") %>%
  tibble::as.tibble()
?dotplot
dotplot_GO.CC.nf <- dotplot(enrich_GO_CCnf,
        title = "GO terms CC")
ggsave("GO.CCnf_enrichmentDotPlot.png")

CCnfsimilarityMatrix <- as.tibble(pairwise_termsim(enrich_GO_CCnf))

emapplot(pairwise_termsim(enrich_GO_CCnf))
ggsave("GO.CCnf_enrichmentPlot.png")

#GO_BP
enrich_GO_BPnf <- enrichGO(
  gene = DEG.nf,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = ALL.nf,
  qvalueCutoff = 0.5,
  readable = TRUE
)
head(enrich_GO_BPnf)

GO_BPnf.t <- enrich_GO_BPnf %>%
  slot("result") %>%
  tibble::as.tibble()
?dotplot
dotplot_GO.BP.nf <- dotplot(enrich_GO_BPnf,
        title = "GO terms BP")
ggsave("GO.BPnf_enrichmentDotPlot.png")

BPnfsimilarityMatrix <- as.tibble(pairwise_termsim(enrich_GO_BPnf))

emapplot(pairwise_termsim(enrich_GO_BPnf))
ggsave("GO.BPnf_enrichmentPlot.png")


##Narrow -----

#MF
enrich_GO_MFnf08 <- enrichGO(
  gene = DEG.nf.narrow08,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = ALL.nf,
  qvalueCutoff = 0.5,
  readable = TRUE
)
head(enrich_GO_MFnf08)

GO_MFnf.t <- enrich_GO_MFnf08 %>%
  slot("result") %>%
  tibble::as.tibble()
?dotplot
dotplot(enrich_GO_MFnf08,
        title = "GO terms MF")
ggsave("GO.MFnf08_enrichmentDotPlot.png")

MFnfsimilarityMatrix <- as.tibble(pairwise_termsim(enrich_GO_MFnf08))

emapplot(pairwise_termsim(enrich_GO_MFnf08))
ggsave("GO.MFnf_enrichmentPlot.png")

#Figure with patchwork
#Requires "package "library(patchwork)" and ggplot2

design.GO <- "AABBCC"

dotplot_GO.MF.nf + dotplot_GO.CC.nf + dotplot_GO.BP.nf +
  plot_layout(design = design.GO) +
  plot_annotation(tag_levels = "A")