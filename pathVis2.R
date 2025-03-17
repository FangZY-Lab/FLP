####progeny ###
library(Seurat)
library(progeny)
library(tidyr)
library(tibble)
library(dplyr)
library(pheatmap)
library(viridis)

#导入数据
S <- readRDS("counts.RDS")
S1 <- S[,S@meta.data$orig.ident %in% c("E-MTAB-6149")]
S2 <- S[,S@meta.data$orig.ident %in% c("GSE131907")]
S3 <- S[,S@meta.data$orig.ident %in% c("GSE143423")]
S4 <- S[,S@meta.data$orig.ident %in% c("GSE189357")]
S5 <- S[,S@meta.data$orig.ident %in% c("PMID34663877")]

S <- S5
Idents(S) <- S$celltype

CellsClusters <- data.frame(cell = names(Idents(S)),
                            celltype = as.character(Idents(S)),
                            stringsAsFactors = F)

S <- progeny(S, scale = F, organism = "Human", top = 500, perm = 1, return_assay = T)
S@assays[["progeny"]]@data[1:6,1:4]

S <- Seurat::ScaleData(S, assay = "progeny")

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(S, slot = "scale.data", assay = "progeny"))) %>%
  rownames_to_column("cell") %>%
  gather(Pathway, Activity, -cell)
dim(progeny_scores_df)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, celltype) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
dim(summarized_progeny_scores)

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = F, stringsAsFactors = F)

paletteLength = 100
mycolor = colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out = ceiling(paletteLength/2) +1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out = floor(paletteLength/2)))


summarized_progeny_scores_df
df_reordered <- summarized_progeny_scores_df[c(11,7,2,6,4,10,5,13,14,1,8,3,12,9)]
pheatmap(t(df_reordered), fontsize = 12,
         fontsize_row = 10,
         cluster_rows = F,
         color = mycolor, breaks = progenyBreaks,
         main = "PMID34663877", angle_col = 45,
         treeheight_col = 0, border_color = NA,legend = F)


