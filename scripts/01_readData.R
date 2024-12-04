rm(list = ls())
library(limma)
library(dplyr)
library(tidyr)
library(tidyverse)
exprSet <- data.table::fread(file = "data/BIC genus expression matrix.txt")

class(exprSet)
exprSet <- as.data.frame(exprSet)
colnames(exprSet)[1] <- "gene_id"
### 查看分组
### https://dwz.cn/WVgQUqfw
### 样本名称

rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
### 我们发现了7个转移的样本，本次分析，我们关注的是癌症和癌旁，先把转移的样本去掉
### 原发和转移的对比作为家庭作业
# 
# TCGA_id <- TCGA_id[substring(TCGA_id,14,15)!="06"]
# TCGA_id <- TCGA_id[substring(TCGA_id,14,15)!="02"]
# 
# exprSet <- cbind(exprSet$gene_id,exprSet[,TCGA_id])
# TCGA_id <- colnames(exprSet)[-1]
# table(substring(TCGA_id,14,15))

# ### 创建metadata
# sample <- ifelse(substring(TCGA_id,14,15)=="01","cancer","normal")
# sample <- factor(sample,levels = c("normal","cancer"),ordered = F)
# metadata <- data.frame(TCGA_id,sample) 
# save(metadata,file = "./output/COAD_metadata.Rdata")
# colnames(exprSet)[1] <- "gene_id"

###将样本名miRNA标注改为RNA标注，方便与转录组数据合并分析
metadata <- data.table::fread(file = "data/BIC metadata matrix.txt")
metadata <- metadata[,c(2,3,4,5,6,9,11,15:19)]
metadata$Sample_type <- ifelse( metadata$Sample_type =="Primary Tumor","cancer","normal")
###筛选消化道肿瘤
target_projects <- c("COAD", "STAD", "LIHC", "ESCA", "PAAD", "CHOL")

# 使用filter函数筛选出Project_name列中是目标项目名称的行
metadata <- metadata %>% filter(Project_name %in% target_projects)

save(metadata,file = "output/metadata.Rdata")
exprSet <- as.data.frame(t(exprSet))
exprSet_change <- cbind(rownames(exprSet),exprSet)
colnames(exprSet_change)[1] <- "Aliquot_barcode_miRNAseq"
metadata_change <- metadata[,c(3,4)]
exprSet_change <- as.data.frame(inner_join(metadata_change,exprSet_change, by="Aliquot_barcode_miRNAseq")) 
exprSet_change <- exprSet_change %>%
  select(-Aliquot_barcode_miRNAseq) %>%
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  ## 把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  ## 去重，Aliquot_barcode_miRNAseq留下第一个
  distinct(Aliquot_barcode_RNAseq,.keep_all = T) %>% 
  ## 反向选择去除rowMean这一列
  select(-rowMean) %>% 
  ## 列名转行名
  column_to_rownames("Aliquot_barcode_RNAseq")
exprSet <- as.data.frame(t(exprSet_change))

TCGA_id <- colnames(exprSet)
table(substring(TCGA_id,14,15))
group <- ifelse(substring(TCGA_id,14,15)=="01","cancer","normal")
group <- factor(group,levels = c("normal","cancer"),ordered = F)


##删除行全部是0的行
zero_rows <- apply(exprSet, 1, function(x) all(x == 0))
table(zero_rows)
exprSet <- exprSet[!zero_rows, ]
save(exprSet,file = "output/exprSet.Rdata")

# res.pca <- prcomp(t(exprSet), scale = TRUE)
# library(factoextra)
# fviz_pca_ind(res.pca,col.ind = group)

### 构建比较矩阵
design <- model.matrix(~group)
### 比较矩阵命名
colnames(design) <- levels(group)
design

### 2.线性模型拟合
fit <- lmFit(exprSet,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
### 4.输出差异分析结果,其中coef的数目不能操过design的列数
### 此处的2代表的是design中第二列和第一列的比较
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
### 这个数据很重要需要保存一下
save(allDiff,file = "output/allDiff.Rdata")

###################################################################################
### 定义差异基因：差异倍数2倍，矫正后的p值小于0.05
load("output/allDiff.Rdata")
library(dplyr)
diffgene <- allDiff %>% 
  filter(P.Value < 0.05) %>% 
  filter(abs(logFC) >0)

### 如果出现行名丢失的情况，需要先把行名变成列，处理完毕后再把列变成行名
### 这个工作是由tibble这个包里面的rownames_to_column()和column_to_rownames()完成的
library(tibble)
diffgene <- allDiff %>% 
  rownames_to_column() %>% 
  filter(P.Value < 0.05) %>% 
  filter(abs(logFC) >0) %>% 
  column_to_rownames()

### 可选方案:使用subset直接获取,&是and的意思
diffgene <- subset(allDiff,abs(logFC) >0 & P.Value < 0.05)
test <- allDiff[allDiff$P.Value < 0.05 & abs(allDiff$logFC)>0,]
### 该数据也需要保存，此处一次性保存两个数据，如果是多个，一次写入变量名称即可。
save(diffgene,group,file = "output/diffgene.Rdata")
### 到此差异基因的分析就结束了
####################################################################################
####################################################################################
## 作图环节
## 1.把现在数据调整成可以作图的格式
### 这个技能是data wrangling部分重点掌握的技能
### 复习一下流程：输入数据是表达量，经过三步
### 1.探针ID转换，2.行列转置，3，添加分组信息。最终获得的是数据框

### 行列转置
exprSet <- as.data.frame(t(exprSet))
### 添加分组信息
exprSet_change <- cbind(rownames(exprSet),exprSet)
colnames(exprSet_change)[1] <- "Aliquot_barcode_RNAseq"
dd <- inner_join(metadata,exprSet_change, by ="Aliquot_barcode_RNAseq")
#dd <- cbind(group=group,exprSet)
### 截取部分展示,这就是清洁数据
test = dd[,1:10]

## 2.作图展示
library(ggplot2)
ggplot(data = dd,aes(x=Project_name,y=Dorea,fill=Project_name))+
  geom_boxplot()+
  geom_point()+
  theme_bw()

## 3.steal plot
my_comparisons <- list(
  c("normal", "cancer")
)
library(ggpubr)
ggboxplot(
  dd, x = "group", y = "Dorea",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

## 改写成函数
diffplot <- function(gene){
  my_comparisons <- list(
    c("normal", "cancer")
  )
  library(ggpubr)
  ggboxplot(
    dd, x = "group", y = gene,
    color = "group", palette = c("#00AFBB", "#E7B800"),
    add = "jitter"
  )+
    #ylim(-0.5,1)+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

diffplot("Dorea")
diffplot("Acinetobacter")

