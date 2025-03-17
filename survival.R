####survival####
library(tidyverse)

surv = read.table(file = 'TCGA-LUAD.survival.tsv', sep = '\t', header = TRUE)

surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]

expr <- read.table("LUAD ssGSEA.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(rownames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[comgene,,drop=F]
surv <- surv[comgene,]

exp_sur <- cbind(expr,surv)

exp_sur$OS.time <- exp_sur$OS.time/365
exp_sur_01A <- exp_sur[substr(rownames(exp_sur),14,16) == "01A",]

surv <- exp_sur_01A
surv$OS.time <- surv$OS.time*12

median(surv$`FLC`)
surv$group <- ifelse(surv$`FLC` > median(surv$`FLC`),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)


library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)

library(survminer)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           palette = "jco",
           legend.labs = c("FLC Low", "FLC High"),
           size = 1,
           xlim = c(0,120),
           break.time.by = 20,
           legend.title = "",
           surv.median.line = "hv",
           ylab = "Survival probability (%)",
           xlab = "Time (Months)",
           ncensor.plot = F,
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)


