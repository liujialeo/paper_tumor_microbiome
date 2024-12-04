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

meta_df <- exprSet[,1:3]
rownames(meta_df) <- meta_df[,1]
meta_df <- meta_df[,-1]

rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-c(1:3)]
res.pca <- prcomp(exprSet, scale = TRUE)

#定义足够多的颜色，用于展示分组

mycol <- c("orangered","deepskyblue")

#用ggord画基本PCA图和置信区间背景色
ggord(res.pca, grp_in = meta_df$group, repel=TRUE,
      alpha = 0.8,#点和置信区间背景设为半透明，以凸显箭头和文字
      #或者单独修改置信区间背景的透明度
      #alpha_el = 0.3,
      #    ellipse = F,
      #ellipse_pro = 0.95,#置信区间
      legend.position = 'top',#将图例放在图的上面
      #poly = F, #去掉圈内填充颜色
      #ellipse = FALSE, hull = TRUE, #改变不同类群的背景形状
      #obslab = TRUE, #将点标注出来
      #ptslab = TRUE,
      size = 2.5, #点的大小
      xlims = c(-12,13),
      ylims = c(-12,10),
      #    hull = T,
      #如果想用默认的颜色，就在下面这行前面加个#
      cols = mycol,
      alpha_el = 0.15, ##置信区间颜色透明度
      #facet = T,
      arrow=NULL, #箭头的头的大小
      #vec_ext = 5,#箭头尾巴长短
      #veccol="brown",#箭头颜色
      txt=NULL) + #箭头指向的基因名的字体大小
  # theme(panel.grid =element_blank()) + 
  theme(legend.position = 'top') #+
#用yyplot继续添加虚线的置信区间
#  geom_ord_ellipse(ellipse_pro = .95, #先画个.95的圆圈
#                   #color='darkgrey', #圈圈的颜色
#                   size=0.5, 
#                   lty= 2 ) + #画成虚线，可以用1-6的数字设置为其他线型
#   geom_ord_ellipse(ellipse_pro = .98, #再画个.98的圆圈
#                   #color='grey', #把这行注释掉，就是跟点一样的颜色
#                   size=0.5, lty=2 ) 

ggsave(filename = "output/PCA_genus.pdf",width = 10,height = 10)

mycol <- c("orangered","deepskyblue","#984ea3","#4daf4a","#08519c","#a50f15")
meta_df$type <- factor(meta_df$type,levels = c("COAD","CHOL","ESCA","LIHC","PAAD","STAD"))

ggord(res.pca, grp_in = meta_df$type, repel=TRUE,
      alpha = 0.8,#点和置信区间背景设为半透明，以凸显箭头和文字
      #或者单独修改置信区间背景的透明度
      #alpha_el = 0.3,
      #    ellipse = F,
      #ellipse_pro = 0.95,#置信区间
      legend.position = 'top',#将图例放在图的上面
      #poly = F, #去掉圈内填充颜色
      #ellipse = FALSE, hull = TRUE, #改变不同类群的背景形状
      #obslab = TRUE, #将点标注出来
      #ptslab = TRUE,
      size = 2.5, #点的大小
      xlims = c(-12,13),
      ylims = c(-12,10),
      #    hull = T,
      #如果想用默认的颜色，就在下面这行前面加个#
      cols = mycol,
      alpha_el = 0.15, ##置信区间颜色透明度
      #facet = T,
      arrow=NULL, #箭头的头的大小
      #vec_ext = 5,#箭头尾巴长短
      #veccol="brown",#箭头颜色
      txt=NULL) + #箭头指向的基因名的字体大小
  # theme(panel.grid =element_blank()) + 
  theme(legend.position = 'top') #+
#用yyplot继续添加虚线的置信区间
#  geom_ord_ellipse(ellipse_pro = .95, #先画个.95的圆圈
#                   #color='darkgrey', #圈圈的颜色
#                   size=0.5, 
#                   lty= 2 ) + #画成虚线，可以用1-6的数字设置为其他线型
#   geom_ord_ellipse(ellipse_pro = .98, #再画个.98的圆圈
#                   #color='grey', #把这行注释掉，就是跟点一样的颜色
#                   size=0.5, lty=2 ) 

ggsave(filename = "output/PCA_cancer.pdf",width = 10,height = 10)
