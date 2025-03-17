####E-MTAB-6149####
library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
S <- readRDS("counts.RDS")

S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
S <- RunHarmony(S, group.by.vars = "orig.ident")
table(S@meta.data$orig.ident)
S@meta.data$patient <- "null"
S@meta.data[S@meta.data$orig.ident%in%c("1a","1b","1c"),]$patient <- "P1"
S@meta.data[S@meta.data$orig.ident%in%c("2a","2b","2c","2d"),]$patient <- "P2"
S@meta.data[S@meta.data$orig.ident%in%c("3a","3b","3c","3d"),]$patient <- "p3"
S@meta.data[S@meta.data$orig.ident%in%c("4a","4b","4c","4d"),]$patient <- "p4"
S@meta.data[S@meta.data$orig.ident%in%c("5a","5b","5c","5d"),]$patient <- "p5"
table(S@meta.data$patient)
S@meta.data$group <- "tumor"
S@meta.data[S@meta.data$orig.ident%in%c("2d","3d","4d","5d"),]$group <- "non-malignant"
table(S@meta.data$group)
ElbowPlot(S,ndims = 50)
S <- FindNeighbors(S,reduction = "harmony",dims = 1:30)
S <- FindClusters(S,resolution = 0.5)
S <- RunUMAP(S,reduction = "harmony",dims = 1:30)
DimPlot(S,reduction = "umap",label = F,raster=F,pt.size=0.5,)

VlnPlot(S,features = c("CD3D","LYZ","CD79A","MZB1","CLDN5","DCN","TPSAB1","EPCAM","COL1A1"))

new.cluster.ids <- c("T/NK","T/NK","Myeloid","Epi","T/NK","Myeloid","B","Epi","T/NK",
                     "T/NK","B","Endo","Fibroblast","Epi","Mast","Epi","Epi","Epi")
names(new.cluster.ids) <- levels(S)
S <- RenameIdents(S, new.cluster.ids)
S$celltype <- S@active.ident
saveRDS(S,file="counts.RDS")
Endo = S[,S@meta.data$celltype %in% c("Endo")]
saveRDS(Endo,file="Endo.RDS")
Epi = S[,S@meta.data$seurat_clusters %in% c(3,7,13,15,16,17)]
saveRDS(Epi,file="Epi.RDS")

library(dplyr)
library(Seurat)
library(ggplot2)
Epi <- readRDS("Epi.RDS")
Epi <- NormalizeData(Epi, normalization.method = "LogNormalize", scale.factor = 10000)
Epi <- FindVariableFeatures(Epi, selection.method = "vst", nfeatures = 2000)
Epi <- ScaleData(Epi)
Epi <- RunPCA(Epi, features = VariableFeatures(object = Epi))
ElbowPlot(Epi,ndims = 50)
Epi <- FindNeighbors(Epi, dims = 1:20)
Epi <- FindClusters(Epi, resolution = 1)
Epi <- RunUMAP(Epi, dims = 1:20)
DimPlot(Epi,reduction = "umap",label = T,raster=F,pt.size=1)
DimPlot(Epi,reduction = "umap",label = T,raster=F,pt.size=1,group.by = "patient")
DimPlot(Epi,reduction = "umap",label = T,raster=F,pt.size=1,group.by = "group")
VlnPlot(Epi, features = c("EPCAM","DCN","COL1A1","COL1A2"),ncol = 2)
FeaturePlot(Epi, features = c("EPCAM","DCN","COL1A1","COL1A2"),reduction = "umap")
new.cluster.ids <- c("Other","Other","Other","Other","Other","Other","Other","Other","FLC",
                     "Other","Other","Other","Other","Other","Other","Other","Other","Other")
names(new.cluster.ids) <- levels(Epi)
Epi <- RenameIdents(Epi, new.cluster.ids)
Epi$celltype <- Epi@active.ident
DimPlot(Epi,reduction = "umap",label = T,raster=F,pt.size=1)
FeaturePlot(Epi, features = c("COL1A1","COL1A2"),reduction = "umap",ncol = 2,raster=F,
            cols = c("lightgray","tomato"))
saveRDS(Epi,file="EPI.RDS")

Idents(Epi) <- Epi$seurat_clusters
VlnPlot(Epi, features = c("CLDN18","SFTPC","CAPS","COL1A1","COL1A2","CD79A","LYZ","CD3D","EPCAM"))
DimPlot(Epi,reduction = "umap",label = T,raster=F,pt.size=1)
new.cluster.ids <- c("Malignant","Malignant","Malignant","Malignant","Malignant","Non-Malignant",
                     "Non-Malignant","Malignant","Malignant","Malignant","Malignant","Malignant",
                     "Non-Malignant","Malignant","Non-Malignant","Non-Malignant","Malignant",
                     "Non-Malignant")
names(new.cluster.ids) <- levels(Epi)
Epi <- RenameIdents(Epi, new.cluster.ids)
Epi$CNV <- Epi@active.ident
DimPlot(Epi,reduction = "umap",label = F,raster=F,group.by = "CNV")
saveRDS(Epi,file="EPI.RDS")

FLC.markers <- FindMarkers(Epi, ident.1 = "FLC", min.pct = 0.25)
FLC.markers %>% top_n(n = 1000, wt = avg_log2FC) -> M1
M1 <- M1[rownames(M1)%in%c("COL1A1","PLAU","PGM2L1","CRIP2","CRIP1","PTGS2","CLDN3","EGFR","PLAT","ACTN1",
                           "NCOA7","AGRN","LAMP1","PSAP","P4HB","HES1"),]
saveRDS(FLC.markers,file="FLC.markers.RDS")

library(DoubletFinder)
EPI <- readRDS("Epi.RDS")
sweep.res.list <- paramSweep(EPI, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
DoubletRate = ncol(EPI)*8*1e-6 
homotypic.prop <- modelHomotypic(EPI$seurat_clusters) 
nExp_poi <- round(DoubletRate*nrow(EPI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

EPI <- doubletFinder(EPI, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
EPI@meta.data[,"Doublet"] <- EPI@meta.data$DF.classifications_0.25_0.02_520
table(EPI@meta.data$Doublet)
DimPlot(EPI, reduction = "umap", group.by ="Doublet")
saveRDS(EPI,file="EPI.RDS")

####GSE131907####
library(dplyr)
library(Seurat)
library(ggplot2)
S.data <- readRDS("tLung.RDS")
S =CreateSeuratObject(counts = S.data,project = 'SP',min.cells = 0, min.features = 0)
remove(S.data)

S[["percent.mt"]] <- PercentageFeatureSet(S, pattern = "^MT-")
VlnPlot(S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
S <- subset(S, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
ElbowPlot(S,ndims = 40)
S <- FindNeighbors(S, dims = 1:20)
S <- FindClusters(S, resolution = 0.5)
S <- RunUMAP(S, dims = 1:20)
DimPlot(S, reduction = "umap",label = T)

saveRDS(S,file="counts.RDS")
S <- readRDS("counts.RDS")
S@meta.data$patient <- substr(rownames(S@meta.data),23,25)
table(S@meta.data$patient)

VlnPlot(S,features = c("CD3D","LYZ","CD79A","MZB1","CLDN5","DCN","TPSAB1","EPCAM","COL1A1"),group.by = "seurat_clusters")
Idents(S) <- S$seurat_clusters
new.cluster.ids <- c("T/NK","Myeloid","T/NK","Epi","B","T/NK","Mast","T/NK","B","Fibroblast",
                     "Myeloid","Myeloid","Myeloid","Endo","Epi","Epi","Fibroblast","Epi","Myeloid",
                     "Epi","B","Epi","Epi")
names(new.cluster.ids) <- levels(S)
S <- RenameIdents(S, new.cluster.ids)
S$celltype <- S@active.ident
saveRDS(S,file="counts.RDS")
Endo = S[,S@meta.data$celltype %in% c("Endo")]
saveRDS(Endo,file="Endo.RDS")

features <- c("CD3D","CD3E","TRBC2",
              "KLRD1","NKG7","GNLY",
              "CD79A","MS4A1","BANK1",
              "LYZ","AIF1","MS4A7",
              "TPSAB1","TPSB2","MS4A2",
              "EPCAM","SCGB3A1","SFTPB",
              "CLDN5","VWF","RAMP2",
              "DCN","COL1A2","COL1A1",
              "MZB1","JCHAIN","IGHA1")
DotPlot(S,features = features,cols = c("lightgray","coral2")) + RotatedAxis()#24*5



EPI = S[,S@meta.data$celltype %in% c("Epithelial")]
saveRDS(EPI,file="EPI.RDS")
EPI <- readRDS("EPI.RDS")
EPI <- NormalizeData(EPI, normalization.method = "LogNormalize", scale.factor = 10000)
EPI <- FindVariableFeatures(EPI, selection.method = "vst", nfeatures = 2000)
EPI <- ScaleData(EPI)
EPI <- RunPCA(EPI, features = VariableFeatures(object = EPI))
ElbowPlot(EPI,ndims = 50)
EPI <- FindNeighbors(EPI, dims = 1:30)
EPI <- FindClusters(EPI, resolution = 0.5)
EPI <- RunUMAP(EPI, dims = 1:30)
DimPlot(EPI, reduction = "umap",label = T)
DimPlot(EPI, reduction = "umap",label = T,group.by = "patient")
VlnPlot(EPI, features = c("FAM183A","AGER","SFTPC","SCGB1A1","COL1A1","EPCAM","PTPRC","LYZ"))
new.cluster.ids <- c("Other","Other","Other","Other","Other","Other","Other","FLC","Other","Other","Other","Other",
                     "FLC","FLC","Other","Other","Other","Other","Other","Other")
names(new.cluster.ids) <- levels(EPI)
EPI <- RenameIdents(EPI, new.cluster.ids)
EPI$celltype <- EPI@active.ident
DimPlot(EPI, reduction = "umap",label = T)

saveRDS(EPI,file="EPI.RDS")
FeaturePlot(EPI, features = c("COL1A1","COL1A2"),reduction = "umap",ncol = 2,raster=F,
            cols = c("lightgray","tomato"))
FeaturePlot(EPI, features = c("EPCAM","DCN","COL1A1","COL1A2"),reduction = "umap")

Idents(EPI) <- EPI$seurat_clusters
VlnPlot(EPI, features = c("CLDN18","SFTPC","CAPS","COL1A1","COL1A2","CD79A","LYZ","CD3D","EPCAM"))
DimPlot(EPI,reduction = "umap",label = T,raster=F,pt.size=1)
new.cluster.ids <- c("Malignant","Malignant","Non-Malignant","Malignant","Malignant","Malignant","Malignant",
                     "Malignant","Malignant","Malignant","Malignant","Malignant","Malignant",
                     "Malignant","Non-Malignant","Non-Malignant","Non-Malignant","Malignant",
                     "Non-Malignant","Non-Malignant")
names(new.cluster.ids) <- levels(EPI)
EPI <- RenameIdents(EPI, new.cluster.ids)
EPI$CNV <- EPI@active.ident
DimPlot(EPI,reduction = "umap",label = F,raster=F,group.by = "CNV")
saveRDS(EPI,file="EPI.RDS")

library(DoubletFinder)
EPI <- readRDS("EPI.RDS")
sweep.res.list <- paramSweep(EPI, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
DoubletRate = ncol(EPI)*8*1e-6 
homotypic.prop <- modelHomotypic(EPI$seurat_clusters) 
nExp_poi <- round(DoubletRate*nrow(EPI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

EPI <- doubletFinder(EPI, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
EPI@meta.data[,"Doublet"] <- EPI@meta.data$DF.classifications_0.25_0.26_219
table(EPI@meta.data$Doublet)
DimPlot(EPI, reduction = "umap", group.by ="Doublet")
saveRDS(EPI,file="EPI.RDS")


####GSE143423####
library(dplyr)
library(Seurat)
S.data <- read.csv(file="GSE143423_lbm_scRNAseq_gene_expression_counts.csv.gz", header=T, sep=",")
rownames(S.data) <- S.data[,1]
S.data <- S.data[,-1]
S =CreateSeuratObject(counts = S.data,project = 'SP',min.cells = 0, min.features = 0)
S[["percent.mt"]] <- PercentageFeatureSet(S, pattern = "^MT-") #MT-  MT.
VlnPlot(S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
S <- subset(S, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 20)
S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
ElbowPlot(S,ndims = 50)
S <- FindNeighbors(S, dims = 1:20)
S <- FindClusters(S, resolution = 1)
S <- RunUMAP(S, dims = 1:20)
DimPlot(S, reduction = "umap",label = T)

saveRDS(S,file="counts.RDS")
S <- readRDS("counts.RDS")
S@meta.data$patient <- substr(rownames(S@meta.data),1,4)
table(S@meta.data$patient)

VlnPlot(S,features = c("CD3D","LYZ","CD79A","MZB1","CLDN5","DCN","TPSAB1","EPCAM","COL1A1"),group.by = "seurat_clusters")

Idents(S) <- S$seurat_clusters
new.cluster.ids <- c("Epi","Epi","Myeloid","Epi","Epi","Epi","Epi","Epi","Epi","Fibroblast",
                     "Endo","Fibroblast","Myeloid","T/NK","Myeloid","Epi")
names(new.cluster.ids) <- levels(S)
S <- RenameIdents(S, new.cluster.ids)
S$celltype <- S@active.ident
saveRDS(S,file="counts.RDS")
Endo = S[,S@meta.data$celltype %in% c("Endo")]
saveRDS(Endo,file="Endo.RDS")
#tumor
Malignant = S[,S@meta.data$seurat_clusters %in% c(0,1,2,3,4,5,7,8,15)]
Malignant <- NormalizeData(Malignant, normalization.method = "LogNormalize", scale.factor = 10000)
Malignant <- FindVariableFeatures(Malignant, selection.method = "vst", nfeatures = 2000)
Malignant <- ScaleData(Malignant)
Malignant <- RunPCA(Malignant, features = VariableFeatures(object = Malignant))
ElbowPlot(Malignant,ndims = 50)
Malignant <- FindNeighbors(Malignant, dims = 1:20)
Malignant <- FindClusters(Malignant, resolution = 0.5)
Malignant <- RunUMAP(Malignant, dims = 1:20)
DimPlot(Malignant, reduction = "umap",label = T)
DimPlot(Malignant, reduction = "umap",label = T,group.by = "patient")
VlnPlot(Malignant, features = c("CD3D","FAM183A","AGER","SFTPC","SCGB1A1","S100A11","COL1A1","CRABP2","EPCAM"))

new.cluster.ids <- c("Other","FLC","Other","FLC","Other","FLC","Other","Other")
names(new.cluster.ids) <- levels(Malignant)
Malignant <- RenameIdents(Malignant, new.cluster.ids)
Malignant$celltype <- Malignant@active.ident
DimPlot(Malignant, reduction = "umap",label = T)
saveRDS(Malignant,file="EPI.RDS")
EPI <- readRDS("EPI.RDS")

FeaturePlot(Malignant, features = c("COL1A1","COL1A2"),reduction = "umap",ncol = 2,raster=F,
            cols = c("lightgray","tomato"))
FeaturePlot(Malignant, features = c("EPCAM","DCN","COL1A1","COL1A2"),reduction = "umap")


Idents(EPI) <- EPI$seurat_clusters
DimPlot(EPI,reduction = "umap",label = T,raster=F,pt.size=1)
new.cluster.ids <- c("Malignant","Malignant","Malignant","Malignant","Malignant","Malignant","Malignant",
                     "Malignant")
names(new.cluster.ids) <- levels(EPI)
EPI <- RenameIdents(EPI, new.cluster.ids)
EPI$CNV <- EPI@active.ident
DimPlot(EPI,reduction = "umap",label = F,raster=F,group.by = "CNV")
saveRDS(EPI,file="EPI.RDS")


library(DoubletFinder)
EPI <- readRDS("EPI.RDS")
sweep.res.list <- paramSweep(EPI, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
DoubletRate = ncol(EPI)*8*1e-6 
homotypic.prop <- modelHomotypic(EPI$seurat_clusters) 
nExp_poi <- round(DoubletRate*nrow(EPI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

EPI <- doubletFinder(EPI, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
EPI@meta.data[,"Doublet"] <- EPI@meta.data$DF.classifications_0.25_0.3_689
table(EPI@meta.data$Doublet)
DimPlot(EPI, reduction = "umap", group.by ="Doublet")
saveRDS(EPI,file="EPI.RDS")

####GSE189357####
library(dplyr)
library(Seurat)
TD1 <- readRDS("TD1EPI.RDS")
TD2 <- readRDS("TD2EPI.RDS")
TD3 <- readRDS("TD3EPI.RDS")
TD4 <- readRDS("TD4EPI.RDS")
TD6 <- readRDS("TD6EPI.RDS")
TD8 <- readRDS("TD8EPI.RDS")
TD9 <- readRDS("TD9EPI.RDS")
EPI <- merge(TD1, y =c(TD2,TD3,TD4,TD6,TD8,TD9), 
             add.cell.ids = c("TD1","TD2","TD3","TD4","TD6","TD8","TD9"))
remove(TD1,TD2,TD3,TD4,TD6,TD8,TD9)
EPI <- NormalizeData(EPI, normalization.method = "LogNormalize", scale.factor = 10000)
EPI <- FindVariableFeatures(EPI, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(EPI)
EPI <- ScaleData(EPI, features = all.genes)
EPI <- RunPCA(EPI, features = VariableFeatures(object = EPI))
ElbowPlot(EPI,ndims = 50)
EPI <- FindNeighbors(EPI, dims = 1:20)
EPI <- FindClusters(EPI, resolution = 0.5)
EPI <- RunUMAP(EPI, dims = 1:20)
DimPlot(EPI, reduction = "umap",label = T)
VlnPlot(EPI, features = c("EPCAM","DCN","COL1A1","COL1A2"))
VlnPlot(EPI, features = c("COL1A1","PLAU","PGM2L1","CRIP2","CRIP1","PTGS2","CLDN3","EGFR","PLAT","ACTN1",
                          "NCOA7","AGRN","LAMP1","PSAP","P4HB","HES1"))
new.cluster.ids <- c("other","other","other","FLC","other","other","other","other","other","other","other",
                     "other","other","other","other","other","FLC","other","other")
names(new.cluster.ids) <- levels(EPI)
EPI <- RenameIdents(EPI, new.cluster.ids)
EPI$celltype <- EPI@active.ident
DimPlot(EPI, reduction = "umap",label = T)
saveRDS(EPI,file="EPI.RDS")
FLC.markers <- FindMarkers(EPI, ident.1 = "FLC", min.pct = 0.25)
FLC.markers %>% top_n(n = 500, wt = avg_log2FC) -> M1
M1 <- M1[rownames(M1)%in%c("COL1A1","PLAU","PGM2L1","CRIP2","CRIP1","PTGS2","CLDN3","EGFR","PLAT","ACTN1",
                           "NCOA7","AGRN","LAMP1","PSAP","P4HB","HES1"),]

####PMID34663877####
library(dplyr)
library(Seurat)
library(tidyr)
S <- readRDS("counts.RDS")
S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
S <- ScaleData(S)
S <- RunPCA(S, features = VariableFeatures(object = S))
S <- RunHarmony(S, group.by.vars = "orig.ident")
ElbowPlot(S,ndims = 50)
S <- FindNeighbors(S,reduction = "harmony",dims = 1:30)
S <- FindClusters(S,resolution = 0.5)
S <- RunUMAP(S,reduction = "harmony",dims = 1:30)
DimPlot(S,reduction = "umap",label = T,raster=F,pt.size=0.5,group.by = "seurat_clusters")
VlnPlot(S,features = c("CD3D","LYZ","CD79A","MZB1","CLDN5","DCN","TPSAB1","EPCAM","COL1A1"),group.by = "seurat_clusters")

Idents(S) <- S$seurat_clusters
new.cluster.ids <- c("T/NK","Epi","T/NK","Myeloid","Myeloid","B","Myeloid","B","T/NK","T/NK",
                     "Myeloid","Mast","Epi","Epi","Fibroblast","Endo","Other","Other","Other",
                     "Other","Other","Myeloid","Other","Other","B","Other")
names(new.cluster.ids) <- levels(S)
S <- RenameIdents(S, new.cluster.ids)
S$celltype <- S@active.ident
S = S[,S@meta.data$celltype %in% c("Epi","T/NK","Myeloid","Endo","Fibroblast","B","Mast")]

saveRDS(S,file="counts.RDS")
Endo = S[,S@meta.data$celltype %in% c("Endo")]
saveRDS(Endo,file="Endo.RDS")

Malignant = S[,S@meta.data$seurat_clusters %in% c(6,10,11,12,22,23)]
Malignant <- NormalizeData(Malignant, normalization.method = "LogNormalize", scale.factor = 10000)
Malignant <- FindVariableFeatures(Malignant, selection.method = "vst", nfeatures = 2000)
Malignant <- ScaleData(Malignant)
Malignant <- RunPCA(Malignant, features = VariableFeatures(object = Malignant))
ElbowPlot(Malignant,ndims = 50)
Malignant <- FindNeighbors(Malignant, dims = 1:30)
Malignant <- FindClusters(Malignant, resolution = 1)
Malignant <- RunUMAP(Malignant, dims = 1:30)
DimPlot(Malignant, reduction = "umap",label = T)
VlnPlot(Malignant, features = c("EPCAM","DCN","COL1A1","COL1A2"))
VlnPlot(Malignant, features = c("COL1A1","PLAU","PGM2L1","CRIP2","CRIP1","PTGS2","CLDN3","EGFR","PLAT","ACTN1",
                                "NCOA7","AGRN","LAMP1","PSAP","P4HB","HES1"))
new.cluster.ids <- c("Other","Other","FLC","Other","Other","Other","FLC","FLC","Other",
                     "Other","Other","Other","Other","FLC","Other","FLC","Other","Other","Other","Other","Other")
names(new.cluster.ids) <- levels(Malignant)
Malignant <- RenameIdents(Malignant, new.cluster.ids)
Malignant$celltype <- Malignant@active.ident
DimPlot(Malignant, reduction = "umap",label = T)
DimPlot(Malignant, reduction = "umap",label = T,group.by = "orig.ident")
FeaturePlot(Malignant, features = c("EPCAM","DCN","COL1A1","COL1A2"),reduction = "umap")

saveRDS(Malignant,file="EPI.RDS")
Malignant <- readRDS("EPI.RDS")
FeaturePlot(Malignant, features = c("COL1A1","COL1A2"),reduction = "umap",ncol = 2,raster=F,
            cols = c("lightgray","tomato"))

FLC <- Malignant[["RNA"]]@counts
FLC <- FLC[rownames(FLC)%in%c("EPCAM","COL1A1","COL1A2","DCN"),]
FLC <- as.data.frame(t(as.matrix(FLC)))
FLC <- FLC[FLC$COL1A2>1,]
Malignant@meta.data$del <- "save"
Malignant@meta.data[rownames(Malignant@meta.data)%in%rownames(FLC),]$del <- "del"
Malignant= Malignant[,Malignant@meta.data$del %in% c("save")]
DimPlot(Malignant, reduction = "umap",label = F,group.by = "celltype")

Idents(EPI) <- EPI$seurat_clusters
VlnPlot(EPI, features = c("CLDN18","SFTPC","CAPS","COL1A1","COL1A2","CD79A","LYZ","CD3D","EPCAM"))
DimPlot(EPI,reduction = "umap",label = T,raster=F,pt.size=1)
new.cluster.ids <- c("Malignant","Malignant","Malignant","Malignant","Malignant","Malignant","Malignant",
                     "Malignant","Non-Malignant","Non-Malignant","Malignant","Malignant","Non-Malignant",
                     "Malignant","Malignant","Malignant","Non-Malignant","Malignant","Malignant","Malignant",
                     "Non-Malignant")
names(new.cluster.ids) <- levels(EPI)
EPI <- RenameIdents(EPI, new.cluster.ids)
EPI$CNV <- EPI@active.ident
DimPlot(EPI,reduction = "umap",label = F,raster=F,group.by = "CNV")
saveRDS(EPI,file="EPI.RDS")


library(DoubletFinder)
EPI <- readRDS("EPI.RDS")
sweep.res.list <- paramSweep(EPI, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
DoubletRate = ncol(EPI)*8*1e-6 
homotypic.prop <- modelHomotypic(EPI$seurat_clusters) 
nExp_poi <- round(DoubletRate*nrow(EPI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

EPI <- doubletFinder(EPI, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)

EPI@meta.data[,"Doublet"] <- EPI@meta.data$DF.classifications_0.25_0.005_549
table(EPI@meta.data$Doublet)
DimPlot(EPI, reduction = "umap", group.by ="Doublet")
saveRDS(EPI,file="EPI.RDS")

FLC.markers <- FindMarkers(Malignant, ident.1 = "FLC", min.pct = 0.25)
FLC.markers %>% top_n(n = 500, wt = avg_log2FC) -> M1
M1 <- M1[rownames(M1)%in%c("COL1A1","PLAU","PGM2L1","CRIP2","CRIP1","PTGS2","CLDN3","EGFR","PLAT","ACTN1",
                           "NCOA7","AGRN","LAMP1","PSAP","P4HB","HES1"),]

