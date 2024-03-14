##Load packages -----
library(tidyverse)
library(annotate)
library(org.Hs.eg.db)
library(plotly)

##Data -----
resSymbols <- getSYMBOL(res.vplot$GeneID, data = 'org.Hs.eg.db')
res.vplot <- as_tibble(resOrdered) %>%
  mutate(GeneID = resOrdered@rownames, .before = 1) %>%
  mutate(Symbol = resSymbols, .before=1)

#res.vplot <- na.omit(res.vplot) ? Get rid of independent filtering?

##Volcano plot -----

vplot <- ggplot(res.vplot) +
  aes(y=-log10(padj), x=log2FoldChange, text = paste("Symbol:", Symbol)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0.8, linetype="longdash", colour="#2C467A", size=1) +
  geom_vline(xintercept = -0.8, linetype="longdash", colour="#BE684D", size=1) +
  labs(title="Volcano plot",
       subtitle = "MMP13 kd",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

## With p-value (not adjusted) to get an idea without all the NAs
ggplot(res.vplot) +
  aes(y=-log10(pvalue), x=log2FoldChange, text = paste("Symbol:", Symbol)) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0.8, linetype="longdash", colour="#2C467A", size=1) +
  geom_vline(xintercept = -0.8, linetype="longdash", colour="#BE684D", size=1) +
  
  ylim(NA, 75) +
  labs(title="Volcano plot",
       subtitle = "MMP13 kd",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggsave("vplot.png")
?ggplot
