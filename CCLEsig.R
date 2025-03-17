####CCLE####
library(data.table)
library(ggplot2)

dat = data.table::fread("OmicsExpressionProteinCodingGenesTPMLogp1.csv",data.table = F)
dat[1:4,1:4]
rownames(dat) <- dat[,1]
dat <- dat[,-1]
exp <- t(dat)
exp[1:4,1:4]

clinical = data.table::fread("Model.csv",data.table = F)
clinical[1:5,c(1:5,11:15)]

ModelID = intersect(colnames(exp),clinical$ModelID)
exp = exp[,ModelID]
clinical = clinical[match(ModelID,clinical$ModelID),]
table(clinical$ModelID %in% colnames(exp))
identical(clinical$ModelID,colnames(exp))

rownames(exp) <- gsub("\\(.*?\\)","",rownames(exp))
rownames(exp) = trimws(rownames(exp))
expr <- exp

table(clinical$OncotreePrimaryDisease)

clinical <- clinical[clinical$OncotreePrimaryDisease%in%c("Bladder Squamous Cell Carcinoma",
                                                          "Bladder Urothelial Carcinoma",
                                                          "Colorectal Adenocarcinoma",
                                                          "Endometrial Carcinoma",
                                                          "Esophageal Squamous Cell Carcinoma",
                                                          "Esophagogastric Adenocarcinoma",
                                                          "Head and Neck Squamous Cell Carcinoma",
                                                          "Hepatocellular Carcinoma",
                                                          "Invasive Breast Carcinoma",
                                                          "Non-Small Cell Lung Cancer",
                                                          "Ovarian Epithelial Tumor",
                                                          "Pancreatic Adenocarcinoma",
                                                          "Prostate Adenocarcinoma"),]
expr <- exp[,colnames(exp)%in%clinical$ModelID]


COL1A1 <- data.frame(COL1A1 = as.numeric(expr["COL1A1",]))
COL1A2 <- data.frame(COL1A2 = as.numeric(expr["COL1A2",]))
DCN <- data.frame(DCN = as.numeric(expr["DCN",]))
CRIP2 <- data.frame(CRIP2 = as.numeric(expr["CRIP2",]))
PLAU <- data.frame(PLAU = as.numeric(expr["PLAU",]))
PLAT <- data.frame(PLAT = as.numeric(expr["PLAT",]))
EGFR <- data.frame(EGFR = as.numeric(expr["EGFR",]))
OSMR <- data.frame(OSMR = as.numeric(expr["OSMR",]))
PTGS2 <- data.frame(PTGS2 = as.numeric(expr["PTGS2",]))
IGFBP3 <- data.frame(IGFBP3 = as.numeric(expr["IGFBP3",]))

pdata <- cbind(COL1A1,COL1A2,DCN,CRIP2,PLAU,PLAT,EGFR,OSMR,PTGS2,IGFBP3)

pdata$row <- colnames(expr)
pdata$group <- "null"
pdata[pdata$COL1A1<1&pdata$COL1A2<1,]$group <- "COL1A1-COL1A2-"
pdata[pdata$COL1A1>1&pdata$COL1A2<1,]$group <- "COL1A1+COL1A2-"
pdata[pdata$COL1A1>1&pdata$COL1A2>1,]$group <- "COL1A1+COL1A2+"

table(pdata$group)
pdata <- pdata[pdata$group%in%c("COL1A1-COL1A2-","COL1A1+COL1A2-"),]

ggplot(pdata,aes(x = group,y = COL1A1)) +
  theme(legend.position = "none")+
  geom_boxplot(aes(fill = group),alpha = 0.7)+
  geom_jitter(aes(color = group))+
  scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  scale_color_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  theme(panel.grid = element_blank())


pdata$OncotreePrimaryDisease <- clinical[match(pdata$row,clinical$ModelID),]$OncotreePrimaryDisease
A <- pdata
expr <- exp[,colnames(exp)%in%A$row]

library(tidyverse)
library(limma)
group <- A[,c(11,12)]
table(group$group)
group$group <- factor(group$group, levels = c("COL1A1+COL1A2-","COL1A1-COL1A2-")) 

design <- cbind(control = ifelse(group$group == "COL1A1-COL1A2-", 1, 0), 
                FLC = ifelse(group$group == "COL1A1+COL1A2-", 1, 0))
contrast.matrix <- makeContrasts(contrasts = 'FLC-control', levels = design)


fit <- lmFit(expr, design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = nrow(fit2), adjust="fdr")

saveRDS(tempOutput,"CCLE FLC marker.RDS")

gene <- c("APOE","CAMK2N1","CLDN3","COL1A1","CRIP1","CRIP2","CTHRC1","CXCL1","CXCL8","DCBLD2","FAM3C","FAM83A",
          "GPR87","HS3ST1","IGFBP3","IRS1","ITGB8","KDR","LINC00511","MUC16","OSMR","PLAT","PLAU","PTGS2",
          "PXDN","RAC1","SLC2A1","SOX11","SPATS2L","SPINT2","SRD5A3","STK24","TNFAIP2","TNNC2","UBALD2",
          "VOPP1","VSTM2L")

output <- tempOutput[tempOutput$logFC>1,]
com <- Reduce(intersect, list(rownames(output),gene))

output <- tempOutput
output$change = ifelse(output$adj.P.Val < 0.05 & abs(output$logFC) >= 1, 
                       ifelse(output$logFC> 1 ,'Up','Down'),
                       'Stable')
geneList0 <- com
geneList <- output[geneList0,]

library('ggplot2')
ggplot(
  output, aes(x = logFC, y = -log10(adj.P.Val), colour=change)) +
  geom_point(alpha=0.5, size=3.5) +
  scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  geom_point(data=geneList,aes(x = logFC, y = -log10(adj.P.Val)),colour="yellow",size=3.5)+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="bottom", 
        legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.key.size = unit(1, 'cm'),
        legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
        axis.title.x =element_text(size=18), 
        axis.title.y=element_text(size=18),
        axis.text=element_text(size=14,face = "bold"))

geneList1 <- output[rownames(output) %in% geneList0,]
geneList1 <- subset(geneList1, select = -change)
geneList1$label <- rownames(geneList1)

library(ggrepel)
p + geom_label_repel(data = geneList1, 
                     aes(x = logFC, y = -log10(adj.P.Val), label = label),
                     size = 3,color="black",
                     box.padding = unit(0.4, "lines"), 
                     segment.color = "black",
                     segment.size = 0.4,
)



