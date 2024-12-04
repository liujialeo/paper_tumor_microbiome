
#############1.加载微生物组数据
#######################
rm(list = ls())
library(dplyr)
library(pheatmap)

library(devtools)
library(ggplot2)
library(plyr)
library(ggord)

load(file = "output/exprSet.Rdata")
load(file = "output/metadata.Rdata")
load(file = "output/gsva_hallmark.Rdata")

exprSet_change <- as.data.frame(t(exprSet*100))
#exprSet_change <- log(exprSet_change+0.0001)
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


exprSet$sum <- rowSums(exprSet[,c(4:19)])
# 计算sum列的平均值
mean_sum <- median(exprSet$sum)
# 根据条件给score列赋值
exprSet$score <- ifelse(exprSet$sum > mean_sum, "high", "low")
exprSet$sample <- substr(exprSet$sample,1,15)
exprSet$SAMPLE_BARCODE <- substr(exprSet$sample,1,12)
exprSet <- exprSet%>%
  distinct(SAMPLE_BARCODE,.keep_all = T)%>%
  distinct(sample,.keep_all = T)
exprSet <- exprSet[exprSet$group == "cancer"]
exprSet <- exprSet[,c(22,1,2,21)]
colnames(exprSet) <- c("PATIENT_BARCODE","SAMPLE_BARCODE","DISEASE","SUBTYPE")



# 让SAMPLE_BARCODE作为rowname
setIndex <- function(df, index = 'SAMPLE_BARCODE'){
  df <- as.data.frame(df)
  rownames(df) <- df[,index]
  df[,index] <- NULL
  return(df)
}

#生成颜色
multiColorMapping <- function(x, type, colors = c("blue", "red"), log = T){
  if(!is.list(colors) && !is.factor(type)){
    type <- as.factor(type)
  }
  if(log){
    x <- log(x+1)
  }
  x_color <- rep("white", length(x))
  for(i in 1:length(x)){
    x_color[i] <- rgb(colorRamp(colors = c("white", 
                                           colors[type[i]]))(x[i]/max(x[type == type[i]])), 
                      maxColorValue = 255)
  }
  return(x_color)
}
# sheet1 is a matrix of patient * alteration (that is type.gene)
data_alteration <- setIndex(readxl::read_excel("NIHMS957693-supplement-8.xlsx", skip = 2))
data_alteration <- data_alteration[rownames(data_alteration) %in% exprSet$SAMPLE_BARCODE,]
write.csv(data_alteration, file = "data/data_alteration.csv", row.names = TRUE)

# sheet2 is a matrix of patient * gene，包含187个基因
data_gene <- setIndex(readxl::read_excel("NIHMS957693-supplement-8.xlsx", skip = 0, sheet = 2))
data_gene <- data_gene[rownames(data_gene) %in% exprSet$SAMPLE_BARCODE,]
write.csv(data_gene, file = "data/data_gene.csv", row.names = TRUE)

# sheet3 is a matrix of patient * pathway，包含10个pathway
data_pathway <- setIndex(readxl::read_excel("NIHMS957693-supplement-8.xlsx", skip = 0, sheet = 3, na = "NA"))
data_pathway <- data_pathway[rownames(data_pathway) %in% exprSet$SAMPLE_BARCODE,]
write.csv(data_pathway, file = "data/data_pathway.csv", row.names = TRUE)


exprSet <- exprSet[exprSet$SAMPLE_BARCODE %in% rownames(data_alteration),]
write.csv(exprSet, file = "data/NIHMS957693-supplement-5.csv", row.names = FALSE)
