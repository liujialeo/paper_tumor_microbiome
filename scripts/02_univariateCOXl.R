rm(list = ls())
#install.packages('survival')
library(survival)
library(tidyr)
library(tidyverse)
library(dplyr)

load(file = "output/exprSet.Rdata")
load(file = "output/metadata.Rdata")

exprSet_change <- as.data.frame(t(exprSet*100))
exprSet_change <- log(exprSet_change+0.0001)
range(exprSet_change)
exprSet <- cbind(rownames(exprSet_change),exprSet_change)
colnames(exprSet)[1] <- "Aliquot_barcode_RNAseq"

# ################################################
# metadata <- metadata[metadata$Sample_type=="cancer",]
# ################################################

clinical <- metadata[,c(4,11,10)]
rt <- inner_join(clinical,exprSet, by = "Aliquot_barcode_RNAseq")
rt <- rt %>% 
  column_to_rownames("Aliquot_barcode_RNAseq")
colnames(rt)[1:2] <- c("futime","fustat")

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="output/diffgene_uniCoxResult.txt",sep="\t",row.names=F,quote=F)

sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<0.05,] #P值的筛选阈值可设为0.10
# 删除含有无穷大值的行

write.table(sigTab,file="output/diffgene_uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)

sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="output/diffgene_uniSigExp.txt",sep="\t",row.names=F,quote=F)
