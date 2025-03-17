####FLC vs Other####
library(dplyr)
library(Seurat)
library(ggplot2)

genes <- c("COL1A1","COL1A2",
           "EPCAM","CDH1","CDH2","CDH3","ICAM1","CLDN3",
           "SNAI1","SNAI2","SNAI3","TWIST1","TWIST2","ZEB1","ZEB2","VIM",
           "CLDN18","SFTPC","AGER","SCGB1A1")

Idents(S) <- S$orig.ident
VlnPlot(S, features = genes,stack=T,pt.size=0,flip = T,add.noise = T,split.by = 'celltype',split.plot = T,
        cols = c("#F8766D", "#729ECE"),)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'top')


