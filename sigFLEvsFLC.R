####FLE vs FLC####
library(tidyverse)
library(Seurat)
S <- readRDS("D:/LAB/RProject/FLC/GSE135893/EPI.RDS")
S@meta.data$orig.ident <- "PF_GSE135893"
counts <- readRDS("counts.RDS")
S <- merge(counts, y =c(S))
S <- S[,S@meta.data$celltype %in% c("FLC","FLE")]

genes <- c("CLDN18","SFTPC","AGER","SCGB1A1")

Idents(S) <- S$orig.ident
VlnPlot(S,features = genes,stack=T,pt.size=0,flip = T,add.noise = T)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')



