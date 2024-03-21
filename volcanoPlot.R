##Load packages -----
library(tidyverse)
library(annotate)
library(org.Hs.eg.db)
library(plotly)

##Data -----
res.nf.vplot <- as_tibble(res.nf) %>%
  mutate(GeneID = res.nf@rownames, .before = 1) 

res.nf.Symbols <- getSYMBOL(res.nf.vplot$GeneID, data = 'org.Hs.eg.db')

res.nf.vplot <- mutate(res.nf.vplot, Symbol = res.nf.Symbols, .before=1)

##Volcano plot -----

vplot <- ggplot(res.nf.vplot) +
  aes(y=-log10(padj), x=log2FoldChange, text = paste("Symbol:", Symbol)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1.0, linetype="longdash", colour="#2C467A", size=1) +
  geom_vline(xintercept = -1.0, linetype="longdash", colour="#BE684D", size=1) +
  ylim(NA, 75) +
  labs(title="Volcano plot",
       subtitle = "MMP13 kd",
       caption=paste0("produced on ", Sys.time())) +
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
