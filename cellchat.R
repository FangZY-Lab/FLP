####cellchat####
library(dplyr)
library(Seurat)
library(ggplot2)

S <- readRDS("counts.RDS")

S1 <- S[,S@meta.data$orig.ident %in% c("E-MTAB-6149")]
counts1 <- readRDS("D:/LAB/RProject/FLC/E-MTAB-6149/counts.RDS")

table(counts1$celltype)
counts1 <- counts1[,counts1@meta.data$celltype %in% c("T/NK","Myeloid","B","Endo","Fibroblast","Mast")]

FLC <- S1[,S1@meta.data$celltype %in% c("FLC")]
FLC$sample <- substr(rownames(FLC@meta.data),1,3)
table(FLC$sample)
table(counts1$orig.ident)
counts1 = counts1[,counts1@meta.data$orig.ident %in% c("p019t","p023t","p024t","p031t")]

S <- merge(counts1,S1,add.cell.ids = c("A","B"))

Idents(S) <- S$celltype

remove(counts1,S1,FLC)

library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)

cellchat <- createCellChat(S@assays$RNA@data)
meta <- data.frame(cellType = S$celltype, row.names =  Cells(S))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "celltype")
cellchat <- setIdent(cellchat, ident.use = "celltype")
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
groupSize

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pathways.show <- unique(df.net$pathway_name)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_scatter(cellchat)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

mat <- cellchat@net$weight
par(mfrow = c(2,4), xpd=TRUE)
par(mar=c(2,2,1,1))
for (i in 1:nrow(mat)) {#
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])}


netVisual_bubble(cellchat, sources.use = c("FLC"), 
                 targets.use = c("T/NK","B","Myeloid","Endo","Mast","Fibroblast"), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = c("FLC"), 
                 targets.use = c("T/NK","B","Myeloid","Endo","Mast","Fibroblast"), 
                 signaling = c("VEGF","VISFATIN","MIF","MK","KIT","GDF","EDN","ANGPTL"), 
                 remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = c("T/NK","B","Myeloid","Endo","Mast","Fibroblast","FLC"), 
                 targets.use = c("FLC"), 
                 signaling = c("MK","ANGPTL","EGF","TWEAK","PARs","PTN","PERIOSTIN","GRN","OSM"), 
                 remove.isolate = FALSE)


