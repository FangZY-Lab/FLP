####AddModuleScore####
setwd("D:/LAB/FLC/fang")
library(Seurat)
library(msigdbr)
S <- readRDS("./counts.RDS")

mdb_H <- msigdbr(species = "Homo sapiens", category = "C2")
mdb_H <- mdb_H[mdb_H$gs_subcat%in%c("CP:KEGG"),]
mdb_H$gs_name <- gsub("KEGG_", "", mdb_H$gs_name)
fgsea_sets<- mdb_H %>% split(x = .$gene_symbol, f = .$gs_name)

markers <- list()
markers$KEGG <- fgsea_sets$ADHERENS_JUNCTION

sce <- AddModuleScore(S,features = markers,name = 
                        'score')
colnames(sce@meta.data)

library(ggpubr)
source("GeomSplitViolin.R")
plot.df <- sce@meta.data[,c("celltype","orig.ident","score1")]
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



