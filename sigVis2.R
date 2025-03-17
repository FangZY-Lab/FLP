## FLC signature violin plot
library(dplyr)
library(Seurat)
library(ggplot2)

genes <- c("COL1A1","DCBLD2","EGFR","FAM3C","HSF1","YY1","PTGS2",
           "PLAT","PLAU","PLAUR","CRIP2","ATOX1","OSMR","STAT3",
           "VOPP1","VSTM2L","SRD5A3","UBALD2")

Idents(S) <- S$orig.ident
VlnPlot(S, features = genes,stack=T,pt.size=0,flip = T,add.noise = T,split.by = 'celltype',split.plot = T,
        cols = c("#F8766D", "#729ECE"),)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'top')



