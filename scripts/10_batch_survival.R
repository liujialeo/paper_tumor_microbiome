### 本节任务: 批量生存分析
#############################################################################=
rm(list = ls())
rt=read.table("output/risk.txt",header=T,sep="\t")
selected <- colnames(rt)[2:19] ###更改列数

load("output/rt_plot.Rdata")
test <- rt_plot[1:10,1:10]
colnames(rt_plot) <- gsub(" ", ".", colnames(rt_plot))
rt_plot <- rt_plot%>%
  select(selected)

library(dplyr)
coxdata <- rt_plot %>% 
  # 去掉小于30天的
  filter(futime >= 30) %>% 
  mutate(futime = futime/365)
test <- coxdata[,1:10]
## 单因素cox分析，正常情况下这里所有的因素都会小于0.05
## gsub
colnames(coxdata) <- gsub("-","_",colnames(coxdata))
colnames(coxdata) <- gsub(":","_",colnames(coxdata))
genes <- colnames(coxdata)[-c(1:2)][1:16] ###更改列数
library(survival)
res <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(futime, fustat)~', genes[i]))
  x = coxph(surv, data = coxdata)
  x = summary(x)
  p.value=signif(x$wald["pvalue"], digits=2)
  HR =signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = HR
  res[i,3] = CI
  res[i,4] = p.value
}
names(res) <- c("ID","HR","95% CI","p.value")
# 删除包含 NA 和 Inf 值的行
res <- na.omit(res)

save(res,file = "output/res_HR.Rdata")

### 批量生存分析log-rank
rm(list = ls())
library(survival)
library(survminer)
library(dplyr)
rt=read.table("output/risk.txt",header=T,sep="\t")
selected <- colnames(rt)[2:19] ###更改列数

load("output/rt_plot.Rdata")
test <- rt_plot[1:10,1:10]
colnames(rt_plot) <- gsub(" ", ".", colnames(rt_plot))
rt <- rt_plot%>%
  select(selected)

res.cut <- surv_cutpoint(rt, 
                         time = "futime", 
                         event = "fustat", 
                         variables = names(rt)[3:ncol(rt)], 
                         minprop = F) 
res.cat <- surv_categorize(res.cut)

test <- res.cat[1:10,1:10]
### 也可以使用median来批量做
colnames(rt) <- gsub("-","_",colnames(rt))
colnames(rt) <- gsub(":","_",colnames(rt))
genes <- colnames(rt)[-c(1:2)][1:16] ###更改列数

res2 <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(futime, fustat)~', genes[i]))
  x = survdiff(surv, data = res.cat)
  pValue=1-pchisq(x$chisq,df=1)
  res2[i,1] = genes[i]
  res2[i,2] = pValue
}
names(res2) <- c("ID","pValue_log")
load(file = "output/res_HR.Rdata")
res1 <- res
res <- merge(res1,res2,by="ID")

### 联合筛选
library(dplyr)
res_filter <- res %>% 
  filter(p.value < 0.01) %>% 
  filter(pValue_log < 0.05)
## 最终得到少量位点，进行下一步筛选。
save(res_filter,file = "output/res_filter.Rdata")

######################################################
### 我们做批量生存
rm(list = ls())
library(survival)
library(survminer)
library(dplyr)
rt=read.table("output/risk.txt",header=T,sep="\t")
selected <- colnames(rt)[2:19] ###更改列数

load("output/rt_plot.Rdata")
test <- rt_plot[1:10,1:10]
colnames(rt_plot) <- gsub(" ", ".", colnames(rt_plot))
rt <- rt_plot%>%
  select(selected)

rt <- rt %>% 
  # 去掉小于30天的
  filter(futime >= 30) %>% 
  mutate(futime = futime/365)

for (i in 3:ncol(rt)) {
  geneindex <- colnames(rt)[i]
  print(colnames(rt)[i])
  rt_batch <- rt[,c("futime","fustat",geneindex)]
  #save(rt,file = "output/rt.Rdata")
  #ifelse 联合median快速二分类
  rt_batch$risk <- ifelse(rt_batch[,geneindex] > median(rt_batch[,geneindex]),"High","Low")

  ### logrank的方法
  ### 首先Surv函数用于创建生存数据对象
  library(survival)
  surv_object = Surv(rt_batch$futime, rt_batch$fustat)
  ## 生存数据拟合survfit
  fit1 <- survfit(surv_object ~ risk, data = rt_batch)
  summary(fit1)
  library(survminer)
  p <- ggsurvplot(fit1,
                  data=rt_batch,
                  surv.median.line = "hv",  # 添加中位生存时间线
                  legend = c(0.8,0.85), # 指定图例位置
                  conf.int=TRUE,
                  pval=TRUE,
                  pval.size=6,
                  risk.table=F,
                  legend.labs=c("High", "Low"),
                  legend.title=paste0(geneindex," Abundance"),
                  xlab="Time(years)",
                  break.time.by = 1,
                  ggtheme = theme_light(),
                  risk.table.y.text.col = T,
                  risk.table.height = 0.2,
                  risk.table.y.text = F,
                  ncensor.plot = F,
                  ncensor.plot.height = 0.2,
                  conf.int.style = "ribbon")
  p

  library(export)
  #graph2pdf(file = paste0("./output/","survival_",geneindex,".pdf"),width = 8,height = 8)
  ggsave(file = paste0("./output/","survival_",geneindex,".pdf"),width = 8,height = 6)
}


# ##################################################################################
# ##突发状况如何处理???
# ### 批量生存分析log-rank
# rm(list = ls())
# if(T){
#   load("resources/dd1.Rda")
#   library(survival)
#   library(survminer)
#   library(dplyr)
#   rt <- dd1 %>% 
#     filter(futime >= 30) %>% # 去掉小于30天的
#     mutate(futime = futime/365)
#   test <- rt[1:10,1:10]
#   #colnames(rt) <- gsub("-","_",colnames(rt))
#   #colnames(rt) <- gsub(":","_",colnames(rt))
# }
# 
# genes <- colnames(rt)[-c(1:2)][1:2000]
# res2 <- data.frame()
# for (i in 1:length(genes)) {
#   print(i)
#   surv =as.formula(paste('Surv(futime, fustat)~', "group"))
#   group = ifelse(rt[,genes[i]]>median(rt[,genes[i]]),"high","low")
#   #if(length(table(group))==1) next
#   data = cbind(rt[,1:2],group)
#   x = survdiff(surv, data = data)
#   pValue=1-pchisq(x$chisq,df=1) 
#   res2[i,1] = genes[i]
#   res2[i,2] = pValue
# }
# names(res2) <- c("ID","pValue_log")
# 
# #######################
# ### 去掉表达量均一的基因
# dd <- rt[,-c(1,2)]
# dd1 <- dd[,apply(dd,2,var)!=0]
# 
# rt <- cbind(rt[,c(1,2)],dd1)
# 
# ## 8秒完成2万个基因的生存分析，人人都可以！
# ## https://mp.weixin.qq.com/s/o4e1HzG4zPIQoGT6-7D0ug