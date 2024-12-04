#####################################################################
##########################################################一、数据清洗

rm(list = ls())


### 1.加载数据
load("data/TCGA-LIHC.Rdata")    #####修改肿瘤缩写
expr_df <- as.data.frame(countsdata)


### 行名变成第一列，cbind十分好用
expr_df <- cbind(gene_id= rownames(expr_df),expr_df)

### 2.基因名称转换
load("data/gtfdata.Rdata")
gtfdata$gene_id <- gsub("\\..*", "",gtfdata$gene_id)

library(dplyr)
exprSet <- gtfdata %>% 
  ## 和表达量的数据交叉合并，等同于merge
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  ## 去掉多余列
  dplyr::select(-gene_id) %>% 
  ## 以下是为了删除重复的行(有些基因名称相同)
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean) 

save(exprSet,file = "output/countsdata_TCGA-LIHC.Rdata")
### 3.准备分类信息
#############################
#制作metadata，不要管这个单词，这一步就是区别肿瘤和正常组
#要对TCGA的id有一点了解，其中第14和15位的数字很重要
#其中01-09是tumor，也就是癌症样本；其中10-29是normal，也就是癌旁

TCGA_id <- colnames(exprSet)[-1]
## 使用table来统计频次么
table(substring(TCGA_id,14,15))

## 分类信息
## 讲解ifelse
sample <- ifelse(substring(TCGA_id,14,15)=="01","cancer","normal")
sample <- factor(sample,levels = c("normal","cancer"))
metadata <- data.frame(TCGA_id,sample) 

### 添加分组信息
colnames(metadata) <- c("sample","group")


#########################################################
### 3.核心环节，构建dds对象，前面的操作都是铺垫
### 要记住四个参数
### 要有数据countData，这里就是exprSet
### 要有分组信息，在colData中，这里是metadata
### design部分是分组信息，格式是~group
### 第一列如果是基因名称，需要自动处理，设置参数tidy=TRUE
### 对象是个复合体
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                             colData=metadata, 
                             design=~group,
                             tidy=TRUE)
nrow(dds)
rownames(dds)
### 筛选样本，counts函数提取表达量数据
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

#########################################################
### 4.数据质量判断
### vst标准化处理
vsd <- vst(dds, blind = FALSE)
### 内置函数plotPCA进行主成分分析画图
plotPCA(vsd, "group")


### 用内置函数plotCounts来进行快速，简易作图
### 找到阳性基因,此处ESR1
### dds来自于上一步
### gene 输入的是ensemble ID
### intgroup 输入的是metadata中分组信息
# plotCounts(dds, gene = "ENSMUSG00000064370", intgroup=c("group"))

### 导出标准化后的表达数据
### assay函数提取vst标准化后的数据，保存数据用于热图
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:5,1:5]
### 保存数据,用于表达量作图，比如差异分析，热图
save(exprSet_vst,file = "output/exprSet_vst_LIHC.Rdata")  #####修改肿瘤缩写

#####################################################################
##########################################################二、GSVA

rm(list = ls())
library(dplyr)
### 1.加载marker
load(file = "data/cellMarker_ssGSEA.Rdata")
### 2.加载表达数据，这是vst标准化后的数据
load(file = "output/exprSet_vst_LIHC.Rdata") #####修改肿瘤缩写
expr <- as.matrix(exprSet_vst)
library(GSVA)
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
##转置
tcga_gsva <- as.data.frame(t(gsva_data))

load("output/metadata.Rdata")
clin <- metadata[,c(4,1,5)]
colnames(clin) <- c("sample","type","group")

tcga_gsva <- cbind(rownames(tcga_gsva),tcga_gsva)
colnames(tcga_gsva)[1] <- "sample"

tcga_gsva <- inner_join(clin,tcga_gsva, by= "sample")
save(tcga_gsva,file = "output/tcga_gsva.Rdata")
### 对照组与实验组的差异

save(tcga_gsva,file = "output/tcga_gsva_LIHC.Rdata") #####修改肿瘤缩写

#####################################################################
###########################################################三、Immune
rm(list = ls())
##加载所有癌症的免疫浸润结果
load("output/tcga_gsva_CHOL.Rdata")
CHOL <- tcga_gsva

load("output/tcga_gsva_COAD.Rdata")
COAD <- tcga_gsva

load("output/tcga_gsva_ESCA.Rdata")
ESCA <- tcga_gsva

load("output/tcga_gsva_STAD.Rdata")
STAD <- tcga_gsva

load("output/tcga_gsva_PAAD.Rdata")
PAAD <- tcga_gsva

load("output/tcga_gsva_LIHC.Rdata")
LIHC <- tcga_gsva

##比较两个数据框列名和列顺序是否完全一致
identical(colnames(CHOL), colnames(COAD)) && all.equal(colnames(CHOL), colnames(COAD))
identical(colnames(COAD), colnames(ESCA)) && all.equal(colnames(COAD), colnames(ESCA))
identical(colnames(ESCA), colnames(STAD)) && all.equal(colnames(ESCA), colnames(STAD))
identical(colnames(STAD), colnames(PAAD)) && all.equal(colnames(STAD), colnames(PAAD))
identical(colnames(PAAD), colnames(LIHC)) && all.equal(colnames(PAAD), colnames(LIHC))

##行合并
gsva <- rbind(CHOL,COAD,ESCA,STAD,PAAD,LIHC)
save(gsva, file = "output/gsva.Rdata")

#####################################################################
###########################################################四、association of Immune and bacteria
rm(list = ls())
library(dplyr)
library(pheatmap)

load(file = "output/exprSet.Rdata")
load(file = "output/metadata.Rdata")
load(file = "output/gsva.Rdata")

exprSet_change <- as.data.frame(t(exprSet*100))
exprSet_change <- log(exprSet_change+0.0001)
range(exprSet_change)
exprSet <- cbind(rownames(exprSet_change),exprSet_change)
colnames(exprSet)[1] <- "Aliquot_barcode_RNAseq"

metadata <- metadata[,c(4,1,5)]
exprSet <- inner_join(metadata,exprSet, by = "Aliquot_barcode_RNAseq")
exprSet <- exprSet[,-c(2:3)]
colnames(exprSet)[1] <- "sample"
#####只做特定的菌
rt=read.table("output/risk.txt",header=T,sep="\t")
selected <- c("sample",colnames(rt)[4:19])
colnames(exprSet) <- gsub(" ", ".", colnames(exprSet))
exprSet <- exprSet%>%
  select(selected)


exprSet <- inner_join(gsva,exprSet, by = "sample")
colnames(exprSet) <- gsub(" ", ".", colnames(exprSet))
exprSet <- as.data.frame(exprSet)
class(exprSet)
save(exprSet, file = "output/exprSet_immune_gsva_genus.Rdata")
library(ggstatsplot)
### data是数据，x是横坐标，y是纵坐标
abc <- as.data.frame(colnames(exprSet))
ggbetweenstats(
  data = exprSet,
  x = group,
  y = Dorea
)

### x变成亚组信息
ggbetweenstats(
  data = exprSet,
  x = group,
  y = Eosinophil
)
### 2.相关性分析
### 两基因相关性作图
### 还是使用ggstatsplot
library(ggstatsplot)
ggscatterstats(
  data = exprSet,
  x = Dorea,
  y = Activated.CD4.T.cell,
  marginal.type = "densigram"
)

### marginal.type可以控制形状
ggscatterstats(
  data = exprSet,
  x = FOXA1,
  y = ESR1,
  marginal.type = "boxplot"
)

gene = colnames(exprSet)[4:47]
### 选取数据相关性分析
M<-cor(exprSet[,gene])
### 作图
library(corrplot)
library(ggplot2)
library(export)

### method是展示形式，order是是否聚类，type限定上下三角
corrplot(M, method="circle")
corrplot(M, method="circle", order="hclust")
corrplot(M, method="pie", order="hclust")
corrplot(M, method="color", order="hclust",col = COL2(n=20))
graph2pdf(file = "output/correlation_immune_genus.pdf",width = 16,height = 16)

corrplot(M, method="number", order="hclust")
graph2pdf(file = "output/correlation_immune_genus_number.pdf",width = 36,height = 36)

corrplot(M, type="upper", order="hclust")
corrplot(M, type="lower", order="hclust")
corrplot(M, method="pie", type="upper", order="hclust")


############################################################################################
######################展示某个菌GSEA-HALLMARK信号通路的相关性，利用
#####################首先计算GSEA信号通路的GSVA结果
rm(list = ls())
library(GSVA)
library(ggplot2)
library(ggcor)
library(data.table)
library(dplyr)
gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}

### 1.加载marker
gset <- gmt2list("data/h.all.v2024.1.Hs.symbols.gmt")
gset.list <- gset[names(gset)]
# 例如提取带有`METABOLISM`字样的通路
#gset.list <- gset[names(gset) %like% "METABOLISM"]
### 2.加载表达数据，这是vst标准化后的数据
load(file = "output/exprSet_vst_LIHC.Rdata") #####修改肿瘤缩写
expr <- as.matrix(exprSet_vst)
gsva_data <- gsva(expr,gset.list, method = "ssgsea")
tcga_gsva <- as.data.frame(gsva_data)

save(tcga_gsva,file = "output/tcga_gsva_hallmark_LIHC.Rdata") #####修改肿瘤缩写


rm(list = ls())
##加载所有癌症的免疫浸润结果
load("output/tcga_gsva_hallmark_CHOL.Rdata")
CHOL <- tcga_gsva


load("output/tcga_gsva_hallmark_COAD.Rdata")
COAD <- tcga_gsva


load("output/tcga_gsva_hallmark_ESCA.Rdata")
ESCA <- tcga_gsva


load("output/tcga_gsva_hallmark_STAD.Rdata")
STAD <- tcga_gsva


load("output/tcga_gsva_hallmark_PAAD.Rdata")
PAAD <- tcga_gsva


load("output/tcga_gsva_hallmark_LIHC.Rdata")
LIHC <- tcga_gsva


##比较两个数据框列名和列顺序是否完全一致
identical(rownames(CHOL), rownames(COAD)) && all.equal(rownames(CHOL), rownames(COAD))
identical(rownames(COAD), rownames(ESCA)) && all.equal(rownames(COAD), rownames(ESCA))
identical(rownames(ESCA), rownames(STAD)) && all.equal(rownames(ESCA), rownames(STAD))
identical(rownames(STAD), rownames(PAAD)) && all.equal(rownames(STAD), rownames(PAAD))
identical(rownames(PAAD), rownames(LIHC)) && all.equal(rownames(PAAD), rownames(LIHC))

##行合并
gsva <- cbind(CHOL,COAD,ESCA,STAD,PAAD,LIHC)
save(gsva, file = "output/gsva_hallmark.Rdata")

#######################
rm(list = ls())
library(dplyr)
library(pheatmap)

load(file = "output/exprSet.Rdata")
load(file = "output/metadata.Rdata")
load(file = "output/gsva_hallmark.Rdata")

exprSet_change <- as.data.frame(t(exprSet*100))
exprSet_change <- log(exprSet_change+0.0001)
range(exprSet_change)
exprSet <- cbind(rownames(exprSet_change),exprSet_change)
colnames(exprSet)[1] <- "Aliquot_barcode_RNAseq"

metadata <- metadata[,c(4,1,5)]

exprSet <- inner_join(metadata,exprSet, by = "Aliquot_barcode_RNAseq")
#exprSet <- exprSet[,-c(2:3)]
colnames(exprSet)[1:3] <-  c("sample","type","group")
#####只做特定的菌
rt=read.table("output/risk.txt",header=T,sep="\t")
selected <- c("sample","type","group",colnames(rt)[4:19])
colnames(exprSet) <- gsub(" ", ".", colnames(exprSet))
exprSet <- exprSet%>%
  select(selected)

gsva <- as.data.frame(t(gsva))
gsva <- cbind(rownames(gsva),gsva)
colnames(gsva)[1] <- "sample"

exprSet <- inner_join(exprSet, gsva, by = "sample")
colnames(exprSet) <- gsub(" ", ".", colnames(exprSet))
exprSet <- as.data.frame(exprSet)
class(exprSet)
save(exprSet, file = "output/exprSet_hallmark_gsva_genus.Rdata")
library(ggstatsplot)
### data是数据，x是横坐标，y是纵坐标
abc <- as.data.frame(colnames(exprSet)[4:19])

ggbetweenstats(
  data = exprSet,
  x = group,
  y = Oceanithermus
)
ggsave("output/relative_abundance_Oceanithermus_group.pdf",width = 5, height = 5)


### x变成亚组信息
ggbetweenstats(
  data = exprSet,
  x = type,
  y = Oceanithermus
)
ggsave("output/relative_abundance_Oceanithermus_cancertype.pdf",width = 8, height = 5)
### 2.相关性分析
### 两基因相关性作图
### 还是使用ggstatsplot
library(ggstatsplot)
ggscatterstats(
  data = exprSet,
  x = Dorea,
  y = HALLMARK_ANGIOGENESIS,
  marginal.type = "densigram"
)

### marginal.type可以控制形状
ggscatterstats(
  data = exprSet,
  x = Dorea,
  y = HALLMARK_E2F_TARGETS,
  marginal.type = "boxplot"
)

gene = colnames(exprSet)[4:69]
### 选取数据相关性分析
M<-cor(exprSet[,gene])
### 作图
library(corrplot)
library(ggplot2)
library(export)

### method是展示形式，order是是否聚类，type限定上下三角
corrplot(M, method="circle")
corrplot(M, method="circle", order="hclust")
corrplot(M, method="pie", order="hclust")
corrplot(M, method="color", order="hclust",col = COL2(n=20))
graph2pdf(file = "output/correlation_hallmark_genus.pdf",width = 25,height = 25)

corrplot(M, method="number", order="hclust")
graph2pdf(file = "output/correlation_hallmark_genus_number.pdf",width = 36,height = 36)

corrplot(M, type="upper", order="hclust")
corrplot(M, type="lower", order="hclust")
corrplot(M, method="pie", type="upper", order="hclust")

###展示数据来源的sankey图
library(ggalluvial)
df <- exprSet[,c(2:3)]
#定义足够多的颜色，后面从这里选颜色
mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)

#格式转换
UCB_lodes <- to_lodes_form(df[,1:ncol(df)],
                           axes = 1:ncol(df),
                           id = "Cohort")
dim(UCB_lodes)
ggplot(UCB_lodes,
       aes(x = x, stratum = stratum, alluvium = Cohort,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/8) + #线跟方块间空隙的宽窄
  geom_stratum(alpha = .9,width = 1/10) + #方块的透明度、宽度
  geom_text(stat = "stratum", size = 3,color="black") + #文字大小、颜色
  
  #不喜欢默认的配色方案，用前面自己写的配色方案
  scale_fill_manual(values = mycol) +
  
  xlab("") + ylab("") +
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #去掉坐标轴
  ggtitle("")+
  guides(fill = FALSE) 
ggsave(filename = "output/sankey.pdf",height = 8,width = 6)


