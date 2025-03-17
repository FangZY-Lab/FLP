####other dataset####
library(tidyverse)
library(Seurat)

S1 <- readRDS("D:/LAB/RProject/FLC/GSE138709/EPI.RDS")
S2 <- readRDS("D:/LAB/RProject/FLC/GSE184362/EPI.RDS")
S3 <- readRDS("D:/LAB/RProject/FLC/GSE135893/EPI.RDS")

S1@meta.data$orig.ident <- "ICC"
S2@meta.data$orig.ident <- "PTC"
S3@meta.data$orig.ident <- "PF"

S1$ICC <- S1$celltype
S2$PTC <- S2$celltype
S3$PF <- S3$celltype

DimPlot(S3, reduction = "umap",label = F,group.by = "PF",
        cols = c("#729ECE","#F8766D"))+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),)

FeaturePlot(S3, features = c("COL1A1","COL1A2"),reduction = "umap",ncol = 1,raster=F,pt.size=0.1,
            cols = c("lightgray","red"))

S <- merge(S1, y =c(S2,S3),add.cell.ids = c("A","B","C"))

S$newtype <- "Other"
S@meta.data[S@meta.data$celltype%in%c("FLC","FLE"),]$newtype <- "FLC/FLE"

gene <- c("COL1A1","COL1A2",
          "CDH1","ICAM1","CLDN3",
          "SNAI1","SNAI2","SNAI3","TWIST1","TWIST2","ZEB1","ZEB2","VIM",
          "PLAT","PLAU","PLAUR","CRIP2","ATOX1","OSMR","STAT3")

Idents(S) <- S$orig.ident
VlnPlot(S, features = gene,stack=T,pt.size=0,flip = T,add.noise = T,split.by = 'newtype',split.plot = T,
        cols = c("#F8766D", "#729ECE"),)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'top')

library(tidyverse)
library(Seurat)

S1 <- readRDS("D:/LAB/RProject/FLC/GSE138709/EPI.RDS")
S2 <- readRDS("D:/LAB/RProject/FLC/GSE184362/EPI.RDS")
S1@meta.data$orig.ident <- "ICC_GSE138709"
S2@meta.data$orig.ident <- "PTC_GSE184362"
S <- merge(S1, y =c(S2),add.cell.ids = c("A","B"))

counts <- readRDS("counts.RDS")
counts <- merge(counts, y =c(S))
counts$newtype <- "Other"
counts@meta.data[counts@meta.data$celltype%in%c("FLC","FLE"),]$newtype <- "FLC/FLE"

gene <- c("APOE","CAMK2N1","CLDN3","COL1A1","CRIP1","CRIP2","CTHRC1","CXCL1","CXCL8","DCBLD2","FAM3C","FAM83A",
          "GPR87","HS3ST1","IGFBP3","IRS1","ITGB8","KDR","LINC00511","MUC16","OSMR","PLAT","PLAU","PTGS2",
          "PXDN","RAC1","SLC2A1","SOX11","SPATS2L","SPINT2","SRD5A3","STK24","TNFAIP2","TNNC2","UBALD2",
          "VOPP1","VSTM2L")

Idents(counts) <- counts$orig.ident
VlnPlot(counts, features = gene,stack=T,pt.size=0,flip = T,add.noise = T,split.by = 'newtype',split.plot = T,
        cols = c("#F8766D", "#729ECE"),)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'top')

library(Seurat)
library(viridis)

gene <- c("COL1A1","PLAT","PLAU","CRIP1","CRIP2","DCBLD2","FAM3C","PTGS2")

markers <- list()
markers$FLC <- gene

sce <- AddModuleScore(S,features = markers,name = 
                        'score')
colnames(sce@meta.data)

library(ggpubr)
source("GeomSplitViolin.R")
plot.df <- sce@meta.data[,c("newtype","orig.ident","score1")]
colnames(plot.df)[2] <- "Dataset"
colnames(plot.df)[3] <- "Signature_Score"

summary_data <- plot.df %>%
  group_by(Dataset,newtype) %>%
  summarise(
    mean = mean(Signature_Score),
    min = mean(Signature_Score) - qnorm(0.975)*sd(Signature_Score)/sqrt(n()), #0.95置信区间
    max = mean(Signature_Score) + qnorm(0.975)*sd(Signature_Score)/sqrt(n())
  )

plot.df$Dataset=factor(plot.df$Dataset, levels=c("ICC", "PTC", "PF"))
plot.df %>% ggplot(aes(Dataset, Signature_Score, fill = newtype)) +
  geom_violin(scale = "width", color = NA, trim = F, alpha = 0.5) +
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9)) +
  stat_compare_means(method = "wilcox.test",aes(label = ..p.signif..)) +
  scale_x_discrete("Datasets")+
  scale_y_continuous("Signature Score")+
  scale_fill_manual(values = c("FLC/FLE"="#F8766D","Other"="#729ECE"),breaks = c("FLC/FLE","Other"),labels=c("FLC/FLE","Other"))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.25,"cm"),
    legend.title = element_blank()
  )


S1 <- readRDS("D:/LAB/RProject/FLC/GSE138709/EPI.RDS")
S2 <- readRDS("D:/LAB/RProject/FLC/GSE184362/EPI.RDS")
S3 <- readRDS("D:/LAB/RProject/FLC/GSE135893/EPI.RDS")

gene <- c("COL1A1","PLAT","PLAU","CRIP1","CRIP2","DCBLD2","FAM3C","PTGS2")

markers <- list()
markers$FLC <- gene

sce <- AddModuleScore(S3,features = markers,name = 
                        'score')
colnames(sce@meta.data)
colnames(sce@meta.data)[9] <- "PF FLE score"
FeaturePlot(sce,'PF FLE score', cols = c("#546de5","#ff4757","red"))


