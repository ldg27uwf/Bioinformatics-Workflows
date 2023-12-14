# script to visualize gene expression data (GSE183947)
setwd("~/Desktop/Bioinformagician/RNA-seq_data")

library(ggplot2)
library(tidyverse)

# basic syntax ggplot

#1. Barplot
dat.long %>%
  filter(gene == 'BRCA1')%>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) +
  geom_col()

#2. Density
dat.long %>%
  filter(gene == 'BRCA1' )%>%
  ggplot(.,aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = .3)

#3 Scatterplot
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x= BRCA1, y = BRCA2, color = tissue)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)

#4. Heatmap  
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

pdf("heatmap_save2.png", width =10, height =8)
#p<- dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill= FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')
#ggsave(p, file = 'heatmap_save1.png', width = 10, height = 8 )

dev.of()
