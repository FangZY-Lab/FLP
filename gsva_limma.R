
####GSVA####
library(msigdbr) #MSigdb gene sets
library(GSVA)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(limma)

genesets <- msigdbr(species = "Homo sapiens", category = "C2")
genesets <- genesets[genesets$gs_subcat%in%c("CP:KEGG"),]
genesets$gs_name <- gsub("KEGG_", "", genesets$gs_name)
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

Idents(S) <- S$origtype
expr <- AverageExpression(S,assays = "RNA",slot = "data")[[1]]
expr <- as.matrix(expr)
gsva.res <- gsva(expr, genesets,method = "gsva")

#group info
group_list <- c("Other","FLC","Other","FLC","Other","FLC","Other","FLC","FLC","Other")

#limma
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(expr)

contrast.matrix<-makeContrasts(contrasts = "FLC-Other",levels = design)
fit <- lmFit(gsva.res,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

diff$group <- ifelse( diff$logFC > 0 & diff$P.Value < 0.05 ,"up" ,
                      ifelse(diff$logFC < 0 & diff$P.Value < 0.05 ,"down","noSig"))
diff2 <- diff %>% 
  mutate(hjust2 = ifelse(t>0,1,0)) %>% 
  mutate(nudge_y = ifelse(t>0,-0.1,0.1)) %>% 
  filter(group != "noSig") %>%
  arrange(t) %>% 
  rownames_to_column("ID")
diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))

ggplot(diff2, aes(ID, t,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","up"), #设置颜色
                    values = c("#008020","#08519C"))+
  geom_text(data = diff2, aes(label = diff2$ID, #添加通路标签
                              y = diff2$nudge_y),
            nudge_x=0,nudge_y=0,hjust =diff2$hjust,
            size=3)+
  labs(x="",
       y=paste0("t value of GSVA score\n","FLC vs Other"),
       title = "GSVA")+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(), #主题微调
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
  )


