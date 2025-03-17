####TCGA####
setwd("D:/LAB/LUAD TCGA")
library(tidyverse)

fpkm_01A <- read.table("LUAD_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
fpkm_11A <- read.table("LUAD_fpkm_mRNA_11A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

gene <- c("APOE","CAMK2N1","CLDN3","COL1A1","CRIP1","CRIP2","CTHRC1","CXCL1","CXCL8","DCBLD2","FAM3C","FAM83A",
          "GPR87","HS3ST1","IGFBP3","IRS1","ITGB8","KDR","LINC00511","MUC16","OSMR","PLAT","PLAU","PTGS2",
          "PXDN","RAC1","SLC2A1","SOX11","SPATS2L","SPINT2","SRD5A3","STK24","TNFAIP2","TNNC2","UBALD2",
          "VOPP1","VSTM2L")
a <- fpkm_01A[gene,]
b <- fpkm_11A[gene,]
a <- na.omit(a)
b <- na.omit(b)
a <- t(a)
b <- t(b)
class(a)
a <- as.data.frame(a)
b <- as.data.frame(b)

a$Sample <- rownames(a)
a$Group <- "Tumor"
b$Sample <- rownames(b)
b$Group <- "Normal"
TCGA <- rbind(a,b)
pdat <- gather(TCGA,key=Gene,value = Expression,-c(Group,Sample))

ggboxplot(pdat, x = "Gene", y = "Expression",outlier.shape = 20,
          fill = "Group", palette = "lancet")+
  stat_compare_means(aes(group = Group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1))


