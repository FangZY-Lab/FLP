####merge####
library(dplyr)
library(Seurat)
library(ggplot2)

S1 <- readRDS("D:/LAB/RProject/FLC/E-MTAB-6149/EPI.RDS")
S2 <- readRDS("D:/LAB/RProject/FLC/GSE131907/EPI.RDS")
S3 <- readRDS("D:/LAB/RProject/FLC/GSE143423/EPI.RDS")
S4 <- readRDS("D:/LAB/RProject/FLC/GSE189357/ALL/EPI.RDS")
S5 <- readRDS("D:/LAB/RProject/FLC/PMID 34663877/EPI.RDS")
S1@meta.data$orig.ident <- "E-MTAB-6149"
S2@meta.data$orig.ident <- "GSE131907"
S3@meta.data$orig.ident <- "GSE143423"
S4@meta.data$orig.ident <- "GSE189357"
S5@meta.data$orig.ident <- "PMID34663877"

S <- merge(S1, y =c(S2,S3,S4,S5))
S = S[,S@meta.data$CNV %in% c("Malignant")]
S = S[,S@meta.data$Doublet %in% c("Singlet")]
table(S@meta.data$CNV)
table(S@meta.data$Doublet)

S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S,features=rownames(S))
S <- RunPCA(S, features = VariableFeatures(object = S))
ElbowPlot(S,ndims = 50)
S <- FindNeighbors(S, dims = 1:30)
S <- FindClusters(S, resolution = 0.5)
S <- RunUMAP(S, dims = 1:30)
DimPlot(S, reduction = "umap",label = T)


