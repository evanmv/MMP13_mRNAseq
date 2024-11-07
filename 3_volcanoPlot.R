##Load packages -----
library(tidyverse)
library(annotate)
library(org.Hs.eg.db)
library(plotly)
install.packages("patchwork")
library(patchwork)
##Data -----
res_nf <- read.csv("res_nf.csv")

names(res_nf)[names(res_nf) == 'X'] <- 'GeneID'

res_nf_vplot <- as_tibble(res_nf) 

res_nf_symbols <- getSYMBOL(as.character(res_nf_vplot$GeneID), data = 'org.Hs.eg.db')

res_nf_vplot <- mutate(res_nf_vplot, Symbol = res_nf_symbols, .before=1)

##Volcano plot -----

vplot <- ggplot(res_nf_vplot) +
  aes(y=-log10(padj), x=log2FoldChange, text = paste("Symbol:", Symbol)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0.8, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -0.8, linetype="longdash", colour="#2C467A", size=1) +
  ylim(NA, 75) +
  labs(title="Volcano plot",
       subtitle = "MMP13 kd",
      ) +
  theme_bw()

ggplotly(vplot)
ggsave("vplotYlim.nf.png")
## With p-value (not adjusted) to get an idea without all the NAs
#ggplot(res.vplot) +
#  aes(y=-log10(pvalue), x=log2FoldChange, text = paste("Symbol:", Symbol)) +
 # geom_point(size=1) +
  #geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 0.8, linetype="longdash", colour="#2C467A", size=1) +
  #geom_vline(xintercept = -0.8, linetype="longdash", colour="#BE684D", size=1) +
  
  #ylim(NA, 75) +
  #labs(title="Volcano plot",
    #   subtitle = "MMP13 kd",
     #  caption=paste0("produced on ", Sys.time())) +
  #theme_bw()

#ggsave("vplot.png")
?ggplot

#DOWN
vplotdown <- ggplot(res_nf_vplot) +
  aes(y=-log10(padj), x=log2FoldChange, text = paste("Symbol:", Symbol)) +
  geom_point(size=1.5, colour = "#2C467A") +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = -0.8, linetype="longdash", colour="#2C467A", size=1) +
  ylim(1.30, 40) +
  xlim(-2, -0.8) +
  labs(title = "Downregulated genes") +
  theme_bw()

ggplotly(vplotdown) 

#UP
vplotup <- ggplot(res_nf_vplot) +
  aes(y=-log10(padj), x=log2FoldChange, text = paste("Symbol:", Symbol)) +
  geom_point(size=1.5, colour = "#BE684D") +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0.8, linetype="longdash", colour="#BE684D", size=1) +
  ylim(1.30, 10) +
  xlim(0.8, 2) +
  labs(title = "Upregulated genes") +
  theme_bw()

ggplotly(vplotup)
?patchwork
design <- "AAAAAA
AAAAAA
BBBCCC"

vplot + vplotdown + vplotup +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A")
ggsave("vplotCombined.png")


