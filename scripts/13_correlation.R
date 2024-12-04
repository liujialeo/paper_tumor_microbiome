rm(list = ls())
source("data/GTBA_functions.R")
load(file = "output/exprSet_immune_gsva_genus.Rdata")
### 准备数据
tcga_splitdata <- split(exprSet,exprSet$type)

### 准备菌/免疫细胞
gene1 ="Granulicella"
gene2 = "Natural.killer.T.cell"
tcga_plotdf <- getpancordata(gene1,gene2,data=tcga_splitdata,method = "spearman")
pancorplot(tcga_plotdf)

### 单组织作图
sincorplot("Granulicella","Natural.killer.T.cell","COAD",data=exprSet,source = "TCGA")
library(ggstatsplot)
ggscatterstats(
  data = exprSet,
  x = Granulicella,
  y = Natural.killer.T.cell,
  marginal.type = "densigram"
)


gene1 ="Granulicella"
gene2 = "Mast.cell"
tcga_plotdf <- getpancordata(gene1,gene2,data=tcga_splitdata,method = "spearman")
pancorplot(tcga_plotdf)

### 单组织作图
sincorplot("Granulicella","Mast.cell","COAD",data=exprSet,source = "TCGA")
library(ggstatsplot)
ggscatterstats(
  data = exprSet,
  x = Granulicella,
  y = Mast.cell,
  marginal.type = "densigram"
)


gene1 ="Granulicella"
gene2 = "Granulicella"
tcga_plotdf <- getpancordata(gene1,gene2,data=tcga_splitdata,method = "spearman")
pancorplot(tcga_plotdf)

### 单组织作图
sincorplot("Granulicella","Mast.cell","COAD",data=exprSet,source = "TCGA")
library(ggstatsplot)
ggscatterstats(
  data = exprSet,
  x = Granulicella,
  y = Mast.cell,
  marginal.type = "densigram"
)



genus <- as.data.frame(colnames(exprSet)[32:47]) 
immune <- as.data.frame(colnames(exprSet)[4:31]) 

for (i in 1:(nrow(genus) - 1)) {
  for (j in (i + 1):nrow(immune)) {
    gene1 <- genus[i,]
    gene2 <- immune[j,]
    
    tcga_plotdf <- getpancordata(gene1, gene2, data = tcga_splitdata, method = "spearman")
    p <- pancorplot(tcga_plotdf)
    
    filename <- paste0("./output/immune_genus_correlation/", paste0(gene1,"_", gene2), ".pdf")
    ggsave(filename, plot = p, width = 8, height = 6)
  }
}


####不同菌分开画图,菌相对丰度，组合到一起
library(tidyr)
dd <- exprSet[,c(3,32:47)]
data <- dd %>%
  pivot_longer(cols=-c(1),
               names_to= "genus",
               values_to = "Abundance")

## 作图
ggplot(data = data,aes(x=genus,y=Abundance,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  stat_compare_means(aes(group=group), label = "p.signif", method = "t.test")

## 尝试更清晰的展示
my_comparisons <- list(
  c("normal", "cancer")
)
ggplot(data = data,aes(x=genus,y=Abundance,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(.~genus)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
