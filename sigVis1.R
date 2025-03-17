##FLC signatures in data
library(dplyr)
library(Seurat)
library(ggplot2)
library(scRNAtoolVis)
gene <- c("APOE","CAMK2N1","CLDN3","COL1A1","CRIP1","CRIP2","CTHRC1","CXCL1","CXCL8","DCBLD2","FAM3C","FAM83A",
          "GPR87","HS3ST1","IGFBP3","IRS1","ITGB8","KDR","LINC00511","MUC16","OSMR","PLAT","PLAU","PTGS2",
          "PXDN","RAC1","SLC2A1","SOX11","SPATS2L","SPINT2","SRD5A3","STK24","TNFAIP2","TNNC2","UBALD2",
          "VOPP1","VSTM2L")
all.genes <- rownames(S)
S <- ScaleData(S, features = all.genes)

DoHeatmap(S, features =gene ,group.by = "celltype",angle = 0,group.colors=c("#F8766D","#729ECE")) + NoLegend() + 
  scale_fill_gradientn(colors=c( '#0099CC', 'white', '#CC0033'))


library(ggpubr)
plot.df <- S@meta.data[,c("celltype","orig.ident","score")]
colnames(plot.df)[2] <- "Dataset"
colnames(plot.df)[3] <- "Signature_Score"

summary_data <- plot.df %>%
  group_by(Dataset,celltype) %>%
  summarise(
    mean = mean(Signature_Score),
    min = mean(Signature_Score) - qnorm(0.975)*sd(Signature_Score)/sqrt(n()), #0.95置信区间
    max = mean(Signature_Score) + qnorm(0.975)*sd(Signature_Score)/sqrt(n())
  )

plot.df %>% ggplot(aes(Dataset, Signature_Score, fill = celltype)) +
  geom_boxplot(width=0.5,outlier.shape = NA,position = position_dodge(0.9)) +
  stat_compare_means(method = "wilcox.test",aes(label = ..p.signif..)) +
  scale_x_discrete("Datasets")+
  scale_y_continuous("FLC Signature Score")+
  scale_fill_manual(values = c("FLC"="#F8766D","Other"="#729ECE"),breaks = c("FLC","Other"),labels=c("FLC","Other"))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.25,"cm"),
    legend.title = element_blank()
  )

features <- c("EPCAM","KRT19","KRT18","CDH1","DCN","ACTA2","COL1A1","COL1A2")
DotPlot(S,features = features,cols = c("lightgray","red"),group.by = "origtype") + RotatedAxis()#24*5


