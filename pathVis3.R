####relation####
library(Seurat)
library(msigdbr)
library(dplyr)
S <- readRDS("./counts.RDS")

markers <- list()
markers$KEGG <- c("APOE","CAMK2N1","CLDN3","COL1A1","CRIP1","CRIP2","CTHRC1","CXCL1","CXCL8","DCBLD2","FAM3C","FAM83A",
                  "GPR87","HS3ST1","IGFBP3","IRS1","ITGB8","KDR","LINC00511","MUC16","OSMR","PLAT","PLAU","PTGS2",
                  "PXDN","RAC1","SLC2A1","SOX11","SPATS2L","SPINT2","SRD5A3","STK24","TNFAIP2","TNNC2","UBALD2",
                  "VOPP1","VSTM2L")

S <- AddModuleScore(S,features = markers,name = 
                      'score')
S$FLC_SIGNATURE <- S$score1
colnames(S@meta.data)

mdb_H <- msigdbr(species = "Homo sapiens", category = "C5")
fgsea_sets<- mdb_H %>% split(x = .$gene_symbol, f = .$gs_name)

markers <- list()
markers$KEGG <- fgsea_sets$GOBP_CELL_CELL_ADHESION

S <- AddModuleScore(S,features = markers,name = 
                      'score')
S$CELL_ADHESION <- S$score1
colnames(S@meta.data)

mdb_H <- msigdbr(species = "Homo sapiens", category = "C5")
mdb_H <- mdb_H[mdb_H$gs_subcat%in%c("GO:BP"),]
fgsea_sets<- mdb_H %>% split(x = .$gene_symbol, f = .$gs_name)

markers <- list()
markers$KEGG <- fgsea_sets$GOBP_CELL_MIGRATION

S <- AddModuleScore(S,features = markers,name = 
                      'score')
S$CELL_MIGRATION <- S$score1
colnames(S@meta.data)


mdb_H <- msigdbr(species = "Homo sapiens", category = "C2")

fgsea_sets<- mdb_H %>% split(x = .$gene_symbol, f = .$gs_name)

markers <- list()
markers$KEGG <- fgsea_sets$SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP

S <- AddModuleScore(S,features = markers,name = 
                      'score')

library(ggpubr)
source("GeomSplitViolin.R")
plot.df <- S@meta.data[,c("celltype","orig.ident","score1")]
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
  geom_violin(scale = "width", color = NA, trim = F, alpha = 0.5) +
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9)) +
  stat_compare_means(method = "wilcox.test",aes(label = ..p.signif..)) +
  scale_x_discrete("Datasets")+
  scale_y_continuous("Signature Score")+
  scale_fill_manual(values = c("FLC"="#F8766D","Other"="#729ECE"),breaks = c("FLC","Other"),labels=c("FLC","Other"))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.25,"cm"),
    legend.title = element_blank()
  )

S$EMT_UP <- S$score1
colnames(S@meta.data)
