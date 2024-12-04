
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
mean_sum <- mean(exprSet$sum)
# 根据条件给score列赋值
exprSet$score <- ifelse(exprSet$sum > mean_sum, "high", "low")
meta <- exprSet[,c(1:3,21)]


############2.加载转录组数据
load(file = "output/countsdata_TCGA-CHOL.Rdata")
CHOL <- exprSet

load(file = "output/countsdata_TCGA-COAD.Rdata")
COAD <- exprSet

load(file = "output/countsdata_TCGA-ESCA.Rdata")
ESCA <- exprSet

load(file = "output/countsdata_TCGA-STAD.Rdata")
STAD <- exprSet

load(file = "output/countsdata_TCGA-PAAD.Rdata")
PAAD <- exprSet

load(file = "output/countsdata_TCGA-LIHC.Rdata")
LIHC <- exprSet

CHOL_COAD <- inner_join(CHOL,COAD, by = "gene_name")
CHOL_COAD_ESCA <- inner_join(CHOL_COAD,ESCA, by = "gene_name")
CHOL_COAD_ESCA_STAD <- inner_join(CHOL_COAD_ESCA,STAD, by = "gene_name")
CHOL_COAD_ESCA_STAD_PAAD <- inner_join(CHOL_COAD_ESCA_STAD,PAAD, by = "gene_name")
CHOL_COAD_ESCA_STAD_PAAD_LIHC <- inner_join(CHOL_COAD_ESCA_STAD_PAAD,LIHC, by = "gene_name")
exprSet_rna <- CHOL_COAD_ESCA_STAD_PAAD_LIHC

############3.合并score和转录组数据
rownames(exprSet_rna) <- exprSet_rna[,1]
exprSet_rna <- exprSet_rna[,-1]
exprSet_rna <- as.data.frame(t(exprSet_rna))
exprSet_rna <- cbind(rownames(exprSet_rna),exprSet_rna)
colnames(exprSet_rna)[1] <- "sample"
exprSet_rna <- inner_join(meta,exprSet_rna, by = "sample")
exprSet_rna <- as.data.frame(exprSet_rna)

genus_score <- exprSet_rna[,c(1,4)]
colnames(genus_score) <- c("sample","group")

rownames(exprSet_rna) <- exprSet_rna[,1] 
exprSet_rna <- exprSet_rna[,-c(1:4)]
exprSet_rna <- as.data.frame(t(exprSet_rna))
exprSet_rna <- cbind(rownames(exprSet_rna),exprSet_rna)
colnames(exprSet_rna)[1] <- "gene_name"

##############4.根据genus_Score计算差异
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=exprSet_rna, 
                             colData=genus_score, 
                             design=~group,
                             tidy=TRUE)
nrow(dds)
rownames(dds)
### 筛选样本，counts函数提取表达量数据
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

### vst标准化处理
vsd <- vst(dds, blind = FALSE)
### 内置函数plotPCA进行主成分分析画图
plotPCA(vsd, "group")


### 用内置函数plotCounts来进行快速，简易作图
### 找到阳性基因,此处ESR1
### dds来自于上一步
### gene 输入的是ensemble ID
### intgroup 输入的是metadata中分组信息
plotCounts(dds, gene = "HMGCR", intgroup=c("group"))

### 导出标准化后的表达数据
### assay函数提取vst标准化后的数据，保存数据用于热图
exprSet_vst <- as.data.frame(assay(vsd))
dds <- DESeq(dds)

####logFC矫正，RNAseq很重要的一步

### contrast参数设置
### 依次是，1.分组信息(metadata中的列) 2.处理组，3.对照组
contrast=c("group", "high", "low")
### results函数获取差异分析的结果
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
### 内置函数plotMA作图
plotMA(dd1, ylim=c(-5,5))
### logFC矫正
dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
plotMA(dd2, ylim=c(-5,5))

### 导出差异分析的结果
library(dplyr)
library(tibble)
library(tidyr)
res <- dd2 %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_name") 
colnames(res) <- c("gene","baseMean","logFC","lfcSE","P.Value","adj.P.Val")
save(res,file = "output/res_genus_score.Rdata")


####################5.GSEA分析
rm(list = ls())

##############################################################
### GSEA 分析
load(file = "output/res_genus_score.Rdata")

library(dplyr)
gene_df <- res %>%
  dplyr::select(logFC,gene) %>%
  ## 去掉NA
  filter(gene!="") %>%
  ## 去掉重复
  distinct(gene,.keep_all = T)

### 1.获取基因logFC
geneList <- gene_df$logFC
### 2.命名
names(geneList) = gene_df$gene
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)
library(clusterProfiler)


### 准备基因集
library(msigdbr)
msigdbr_species()
msigdbr_collections()
library(dplyr)
## Hallmarks gene set
h_df <- msigdbr(species = "Homo sapiens") %>% 
  filter(gs_cat == "H") %>% 
  dplyr::select(gs_name,gene_symbol)

y <- GSEA(geneList,TERM2GENE =h_df)
yd <- as.data.frame(y)
library(ggplot2)
dotplot(y,showCategory=15,split=".sign")+facet_grid(~.sign)
ggsave("output/genus_GSEA_dotploot.pdf", width = 9, height = 10)

### 自定义画图
ggplot(y, showCategory = 20, aes(NES, forcats::fct_reorder(Description, NES))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  #scale_color_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() + 
  xlab("Normalized Enrichment Score") +
  ylab(NULL)
ggsave("output/genus_GSEA_dotploot_2.pdf", width = 8, height = 5)
### 火山图展示
hallmarks <- h_df
y <- GSEA(geneList,TERM2GENE =hallmarks,pvalueCutoff = 1)
yd <- as.data.frame(y)
### 看整体分布
dotplot(y,showCategory=30,split=".sign")+facet_grid(~.sign)
y1 <- filter(y,p.adjust < 0.05)
dotplot(y1,showCategory=30,split=".sign")+facet_grid(~.sign)
yd1 <- as.data.frame(y1)

library(ggrepel)
yd$sig = ifelse(yd$p.adjust<0.05,ifelse(yd$NES>0,"pos","neg"),"none")
yd$sig = factor(yd$sig,levels = c("pos","neg","none"),
                labels = c("pos","neg","none"))

ggplot(yd,aes(-log10(p.adjust),NES)) + 
  geom_point(aes(color=sig,size = NES))+
  geom_vline(xintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_hline(yintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
  ## 主题
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  #labs(x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')+
  ## 标签
  geom_text_repel(data=subset(yd, p.adjust < 0.05), 
                  aes(label=ID),col="black",
                  force        = 1.2,
                  nudge_x      = 0.3,
                  direction    = "y",
                  hjust        = 0,
                  segment.size = 0.3)
ggsave("output/genus_GSEA_dotploot_volcano.pdf", width = 15, height = 12)
#######单独画特定通路的GSEA图

library(enrichplot)
pathway.id = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
ggsave("output/genus_GSEA_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf", width = 12, height = 8)


###批量GSEA作图
### 新建文件夹，
dir.create("output/GSEA_out")
for (i in 1:nrow(yd1) ) {
  print(i)
  pathway.id = yd1$ID[i]
  p = gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
  ggsave(p, filename = paste0("output/GSEA_out/GSEA_",pathway.id,".pdf"), width = 12, height = 8)
  ggsave(p, filename = paste0("output/GSEA_out/GSEA_",pathway.id,".tiff"), width = 12, height = 8)
}
