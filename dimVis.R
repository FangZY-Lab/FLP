##visualization on 2D
data_min100S.top1000Gavg.supervised_umap_v2 = readRDS('D:/LAB/FLC/fang/pipeline/data_min100S.top1000Gavg.supervised_umap_v2_s10.rds')
Y <- readRDS('D:/LAB/FLC/fang/pipeline/data_min100S.meta.rds')
S <- readRDS("D:/LAB/FLC/counts.RDS")

df = data.frame(U=data_min100S.top1000Gavg.supervised_umap_v2[[1]],Y[,-(1:2)])
ggplot(df,aes(U.1,U.2))+geom_point(aes(colour=celltype))

meta <- df

S$cellid <- rownames(S@meta.data)
S = S[,S@meta.data$cellid %in% rownames(meta)]

meta = meta[S$cellid,]
umap1 <- meta$U.1
names(umap1) <- rownames(meta)
S@reductions[["umap"]]@cell.embeddings[,1] <- umap1
umap2 <- meta$U.2
names(umap2) <- rownames(meta)
S@reductions[["umap"]]@cell.embeddings[,2] <- umap2
S@meta.data$celltype <- meta$celltype
Idents(S)='celltype'
DimPlot(S, reduction = 'umap', label = F, repel = T)
FeaturePlot(S, features = c("EPCAM","DCN","COL1A1","COL1A2"),reduction = "umap",ncol = 4,raster=F,pt.size=0.1,
            cols = c("lightgray","red"))
features <- c("EPCAM","KRT19","KRT18","CDH1","DCN","ACTA2","COL1A1","COL1A2")
DotPlot(S,features = features,cols = c("lightgray","red"),group.by = "origtype") + RotatedAxis()#24*5



