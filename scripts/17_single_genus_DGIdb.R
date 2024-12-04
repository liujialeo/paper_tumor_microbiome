

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
df <- data.table::fread(file = "data/DGIdb-interactions.tsv")
df <- df[,c(3,1)]
# #将每个药物与其对应的基因一一对应
# drug_gene_list <- lapply(strsplit(df$TARGET, ", "), trimws)  # trimws函数用于去除基因名称中可能存在的空格
# drug_gene_df <- data.frame(
#   DRUG_NAME = rep(df$DRUG_NAME, sapply(drug_gene_list, length)),
#   GENE = unlist(drug_gene_list)
# )
colnames(df) <- c("trem","gene")

y <- GSEA(geneList,TERM2GENE =df,pvalueCutoff = 0.08)
yd <- as.data.frame(y)
library(ggplot2)
dotplot(y,showCategory=15,split=".sign")+facet_grid(~.sign)

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
 ggsave("output/genus_GSEA_DGIdb_dotploot_2.pdf", width = 13, height = 6)
# ### 火山图展示
# hallmarks <- df
# y <- GSEA(geneList,TERM2GENE =hallmarks,pvalueCutoff = 1)
# yd <- as.data.frame(y)
# ### 看整体分布
# dotplot(y,showCategory=30,split=".sign")+facet_grid(~.sign)
# y1 <- filter(y,p.adjust < 0.05)
# dotplot(y1,showCategory=30,split=".sign")+facet_grid(~.sign)
# yd1 <- as.data.frame(y1)
# 
# library(ggrepel)
# yd$sig = ifelse(yd$p.adjust<0.05,ifelse(yd$NES>0,"pos","neg"),"none")
# yd$sig = factor(yd$sig,levels = c("pos","neg","none"),
#                 labels = c("pos","neg","none"))
# 
# ggplot(yd,aes(-log10(p.adjust),NES)) + 
#   geom_point(aes(color=sig,size = NES))+
#   geom_vline(xintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
#   geom_hline(yintercept = 0,lty=4,lwd=0.6,alpha=0.8)+
#   ## 主题
#   theme_bw()+
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),   
#         axis.line = element_line(colour = "black"))+
#   #labs(x="log2 (fold change)",y="-log10 (q-value)")+
#   theme(plot.title = element_text(hjust = 0.5))+
#   theme(legend.position='none')+
#   ## 标签
#   geom_text_repel(data=subset(yd, p.adjust < 0.05), 
#                   aes(label=ID),col="black",
#                   force        = 1.2,
#                   nudge_x      = 0.3,
#                   direction    = "y",
#                   hjust        = 0,
#                   segment.size = 0.3)
# ggsave("output/genus_GSEA_dotploot_volcano.pdf", width = 15, height = 12)
#######单独画特定通路的GSEA图

library(enrichplot)
pathway.id = "TICLOPIDINE"
gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
ggsave("output/genus_DGIdb_TICLOPIDINE.pdf", width = 12, height = 8)


###批量GSEA作图
### 新建文件夹，
dir.create("output/GSEA_out")
for (i in 1:nrow(yd) ) {
  print(i)
  pathway.id = yd$ID[i]
  p = gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
  ggsave(p, filename = paste0("output/GSEA_out/DGI_",pathway.id,".pdf"), width = 12, height = 8)
  ggsave(p, filename = paste0("output/GSEA_out/DGI_",pathway.id,".tiff"), width = 12, height = 8)
}
