##cancer cell lines
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
identical(clinical$ModelID,colnames(exp))## [1] TRUE
rownames(exp) <- gsub("\\(.*?\\)","",rownames(exp))
rownames(exp) = trimws(rownames(exp))
expr <- exp

COL1A1 <- data.frame(Expression = as.numeric(exp[grep("COL1A1",rownames(exp)),]),clinical)
COL1A2 <- data.frame(Expression = as.numeric(exp[grep("COL1A2",rownames(exp)),]),clinical)
EPCAM <- data.frame(Expression = as.numeric(exp[grep("EPCAM",rownames(exp)),]),clinical)
DCN <- data.frame(Expression = as.numeric(exp[grep("DCN",rownames(exp)),]),clinical)

COL1A1$gene <- "COL1A1"
COL1A2$gene <- "COL1A2"

A <- rbind(COL1A1,COL1A2)
A <- A[,c("StrippedCellLineName","Expression","gene")]
colnames(A) <- c("CellLineName","Expression","gene")

B <- A[A$CellLineName%in%c("A427","A549","ABC1","CALU6","LUDLU1","NCIH1395"),]
C <- A[A$CellLineName%in%c("DM3","HLFA","HS229T","HS618T","RS5","TIG3TD"),]

ggplot(B,aes(CellLineName,Expression))+
  geom_col(aes(fill=CellLineName))+
  scale_fill_manual(values=c("#00c16e","#729ECE","#FF9E4A","#CDCC5D","#ED665D","#AD8BC9"))+
  facet_grid(~gene)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab("")+
  theme(panel.grid = element_blank())


