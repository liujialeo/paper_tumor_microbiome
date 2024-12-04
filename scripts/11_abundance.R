#############不同肿瘤类型
rm(list = ls())
library(survival)
library(survminer)
library(dplyr)
rt=read.table("output/risk.txt",header=T,sep="\t")
selected <- colnames(rt)[4:19]

load("output/rt_plot.Rdata")
test <- rt_plot[1:10,1:10]
colnames(rt_plot) <- gsub(" ", ".", colnames(rt_plot))
rt <- rt_plot%>%
  select(selected)

load(file = "output/metadata.Rdata")
metadata <- metadata[,c(4,1)]


rt <- cbind(rownames(rt),rt)
colnames(rt)[1] <- "Aliquot_barcode_RNAseq"
dd <- inner_join(metadata,rt, by = "Aliquot_barcode_RNAseq")


dd1 <- dd %>% 
  pivot_longer(cols=3:18,
               names_to= "genus",
               values_to = "Abundance")
dd1 <- dd1[,-1]
colnames(dd1)[1] <- "sample"

library(ggplot2)
library(ggpubr)


### 箱线图
ggplot(data =dd1, aes(x = genus, y = Abundance))+
  geom_boxplot(aes(fill = sample),outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")

### 小提琴
ggplot(data =dd1, aes(x = genus, y = Abundance))+
  geom_violin(aes(fill = sample),position = position_dodge(1),scale = "width")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")

### 混合叠加
ggplot(data =dd1, aes(x = genus, y = Abundance))+
  geom_boxplot(aes(fill = sample),position = position_dodge(1),width=.3,outlier.shape = NA)+
  geom_violin(aes(colour = sample),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), 
                     #method = "t.test",
                     method = "anova",
                     label = "p.signif")+
  scale_fill_manual(values = c("#fee0d2","#fc9272","#ef3e2f","#a50f15","#08519c","#9ecae1")) +  # 修改填充颜色
  scale_color_manual(values = c("#fee0d2","#fc9272","#ef3e2f","#a50f15","#08519c","#9ecae1"))    # 修改边框颜色




#############肿瘤vs正常
rm(list = ls())
library(survival)
library(survminer)
library(dplyr)
rt=read.table("output/risk.txt",header=T,sep="\t")
selected <- colnames(rt)[4:19]

load("output/rt_plot.Rdata")
test <- rt_plot[1:10,1:10]
colnames(rt_plot) <- gsub(" ", ".", colnames(rt_plot))
rt <- rt_plot%>%
  select(selected)

load(file = "output/metadata.Rdata")
metadata <- metadata[,c(4,5)]


rt <- cbind(rownames(rt),rt)
colnames(rt)[1] <- "Aliquot_barcode_RNAseq"
dd <- inner_join(metadata,rt, by = "Aliquot_barcode_RNAseq")


dd1 <- dd %>% 
  pivot_longer(cols=3:18,
               names_to= "genus",
               values_to = "Abundance")
dd1 <- dd1[,-1]
colnames(dd1)[1] <- "sample"

library(ggplot2)
library(ggpubr)


### 箱线图
ggplot(data =dd1, aes(x = genus, y = Abundance))+
  geom_boxplot(aes(fill = sample),outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  #ylim(0,10)
  stat_compare_means(aes(group=sample), label = "p.signif")

### 小提琴
ggplot(data =dd1, aes(x = genus, y = Abundance))+
  geom_violin(aes(fill = sample),position = position_dodge(1),scale = "width")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")

### 混合叠加
ggplot(data =dd1, aes(x = genus, y = Abundance))+
  geom_boxplot(aes(fill = sample),position = position_dodge(1),width=.3,outlier.shape = NA)+
  geom_violin(aes(colour = sample),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), 
                     #method = "t.test",
                     method = "anova",
                     label = "p.signif")+
  scale_fill_manual(values = c("#fee0d2","#fc9272","#ef3e2f","#a50f15","#08519c","#9ecae1")) +  # 修改填充颜色
  scale_color_manual(values = c("#fee0d2","#fc9272","#ef3e2f","#a50f15","#08519c","#9ecae1"))    # 修改边框颜色



####heatdata
rm(list = ls())
library(dplyr)
library(pheatmap)
load(file = "output/exprSet.Rdata")
load(file = "output/metadata.Rdata")

exprSet_change <- as.data.frame(t(exprSet*100))
exprSet_change <- log(exprSet_change+0.0001)
range(exprSet_change)
exprSet <- cbind(rownames(exprSet_change),exprSet_change)
colnames(exprSet)[1] <- "Aliquot_barcode_RNAseq"

metadata <- metadata[,c(4,1,5)]
heatdata <- inner_join(metadata,exprSet, by = "Aliquot_barcode_RNAseq")


heatdata <- as.data.frame(heatdata)
rownames(heatdata) <- heatdata[,1]
metadata <- heatdata[,1:3]
heatdata <- heatdata[,-c(1:3)]
#####只做特定的菌
rt=read.table("output/risk.txt",header=T,sep="\t")
selected <- colnames(rt)[4:19]
colnames(heatdata) <- gsub(" ", ".", colnames(heatdata))
heatdata <- heatdata%>%
  select(selected)

heatdata <- as.data.frame(t(heatdata))

annotation_col <- metadata
annotation_col <- as.data.frame(annotation_col)
rownames(annotation_col) <- annotation_col[,1]
annotation_col <- annotation_col[,-1]
annotation_col$Sample_type <- factor(annotation_col$Sample_type, levels = c("normal","cancer"))
### 如果注释出界, 可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("#08519c", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 0.2, cellheight = 10,# 格子比例
         fontsize = 10)

