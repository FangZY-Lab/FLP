####ssGSEA####
library(GSVA)
library(tidyverse)
setwd("D:/LAB/LUAD TCGA")
expr <- data.table::fread("LUAD_fpkm_mRNA_all.txt",data.table = F) #读取表达文件
expr <- as.matrix(expr)
markers <- list()
markers$FLC <- c("APOE","CAMK2N1","CLDN3","COL1A1","CRIP1","CRIP2","CTHRC1","CXCL1","CXCL8","DCBLD2","FAM3C","FAM83A",
                 "GPR87","HS3ST1","IGFBP3","IRS1","ITGB8","KDR","LINC00511","MUC16","OSMR","PLAT","PLAU","PTGS2",
                 "PXDN","RAC1","SLC2A1","SOX11","SPATS2L","SPINT2","SRD5A3","STK24","TNFAIP2","TNNC2","UBALD2",
                 "VOPP1","VSTM2L")
gsva_data <- gsva(expr,markers, method = "ssgsea")

a <- gsva_data %>% t() %>% as.data.frame()
a$group =  ifelse(str_sub(rownames(a),14,15)<10,"tumor","normal")
colnames(a)[2] <- "LUAD"
a <- gather(a,key=ssGSEA,value = Expression,-c(group,sample))


expr <- data.table::fread("LUSC_fpkm_mRNA_all.txt",data.table = F)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.matrix(expr)
gsva_data <- gsva(expr,markers, method = "ssgsea")
b <- gsva_data %>% t() %>% as.data.frame()
write.table(b, file = "LUSC ssGSEA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
b$group =  ifelse(str_sub(rownames(b),14,15)<10,"tumor","normal")
b <- b %>% rownames_to_column("sample")
colnames(b)[2] <- "LUSC"
b <- gather(b,key=ssGSEA,value = Expression,-c(group,sample))
com <- rbind(a,b)

ggboxplot(com, x = "ssGSEA", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1))


