rm(list = ls())
library(survival)
library(dplyr)
rt=read.table("output/diffgene_lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
colnames(rt) <- gsub(" ", ".", colnames(rt))
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
###################################
outTab <- as.data.frame(outTab)
outTab=outTab[as.numeric(as.vector(outTab$pvalue))<0.01,] #P值的筛选阈值可设为0.10
#####################################
write.table(outTab,file="output/diffgene_multiCox.xls",sep="\t",row.names=F,quote=F)

selected <- c("futime", "fustat", rownames(outTab))

rt <- rt %>%
  select(selected)

#####################################################################################
###重新计算多因素方差
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
###################################


#####################################
write.table(outTab,file="output/diffgene_multiCox.xls",sep="\t",row.names=F,quote=F)

#####################################################################################

#计算出每个患者的风险评分，并按中位数分为高低危
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="output/risk.txt",
            sep="\t",
            quote=F,
            row.names=F)
