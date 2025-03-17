####FLC vs Fib####
Fib1 <- readRDS("D:/LAB/RProject/FLC/E-MTAB-6149/Fib.RDS")
Fib2 <- readRDS("D:/LAB/RProject/FLC/GSE131907/Fib.RDS")
Fib3 <- readRDS("D:/LAB/RProject/FLC/GSE143423/Fib.RDS")
Fib4 <- readRDS("D:/LAB/RProject/FLC/GSE189357/ALL/Fib.RDS")
Fib5 <- readRDS("D:/LAB/RProject/FLC/PMID 34663877/Fib.RDS")

Fib1@meta.data$orig.ident <- "E-MTAB-6149"
Fib2@meta.data$orig.ident <- "GSE131907"
Fib3@meta.data$orig.ident <- "GSE143423"
Fib4@meta.data$orig.ident <- "GSE189357"
Fib5@meta.data$orig.ident <- "PMID34663877"

Fib <- merge(Fib1, y =c(Fib2,Fib3,Fib4,Fib5))

Epi <- readRDS("counts.RDS")
S <- merge(Epi, y =c(Fib),add.cell.ids = c("Epi","Fib"))
table(S$celltype)
S = S[,S@meta.data$celltype %in% c("FLC","Fibroblast")]

genes <- c("DCN","ACTA2","FGF7","S100A4","FAP","PDGFRA","PDGFRB","THY1","PDPN",
           "CAV1","ITGB1","SNAI1","TWIST1","ZEB1","VIM","COL1A1","COL1A2","COL3A1","COL4A1","COL5A1")

celllevel <- c("FLC","Fibroblast")
S$celltype <- factor(S$celltype, levels=celllevel)

Idents(S) <- S$orig.ident

VlnPlot(S, features = genes,stack=T,pt.size=0,flip = T,add.noise = T,split.by = 'celltype',split.plot = T,
        cols = c("#F8766D","#66CC66"),)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'top')



