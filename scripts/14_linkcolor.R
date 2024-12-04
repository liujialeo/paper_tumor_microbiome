######################展示某个菌和immune cell的相关性，利用GSVA结果
rm(list = ls())
library(data.table)
library(GSVA)
library(ggplot2)
library(ggcor)
selectedgenus="Oceanithermus"  ###修改菌名

load(file = "output/exprSet_immune_gsva_genus.Rdata")
rownames(exprSet) <- exprSet[,1]

exprSet_immune <- exprSet[,4:31]
exprSet_immune <- as.data.frame(t(exprSet_immune)) 


exprSet_genus <- as.data.frame(exprSet[,colnames(exprSet) == selectedgenus, drop = FALSE])
exprSet_genus <- as.data.frame(t(exprSet_genus))

immPath.score <- rbind(exprSet_immune,exprSet_genus)

immCorSiglec15 <- NULL
for (i in rownames(immPath.score)) {
  cr <- cor.test(as.numeric(immPath.score[i,]),
                 as.numeric(exprSet_genus),
                 method = "pearson")
  immCorSiglec15 <- rbind.data.frame(immCorSiglec15,
                                     data.frame(gene = selectedgenus,  ###修改菌名
                                                path = i,
                                                r = cr$estimate,
                                                p = cr$p.value,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
}
immCorSiglec15$sign <- ifelse(immCorSiglec15$r > 0,"pos","neg")
immCorSiglec15$absR <- abs(immCorSiglec15$r)
immCorSiglec15$rSeg <- as.character(cut(immCorSiglec15$absR,c(0,0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"),include.lowest = T))
immCorSiglec15$pSeg <- as.character(cut(immCorSiglec15$p,c(0,0.001,0.01,0.05,1),labels = c("<0.001","<0.01","<0.05","ns"),include.lowest = T))
immCorSiglec15[nrow(immCorSiglec15),"pSeg"] <- "Not Applicable"

immCorSiglec15$rSeg <- factor(immCorSiglec15$rSeg, levels = c("0.25","0.50","0.75","1.00"))
immCorSiglec15$pSeg <- factor(immCorSiglec15$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
immCorSiglec15$sign <- factor(immCorSiglec15$sign, levels = c("pos","neg"))

p1 <- quickcor(t(immPath.score), 
               type = "lower",
               show.diag = TRUE) + 
  geom_colour() +
  add_link(df = immCorSiglec15, 
           mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
           spec.key = "gene",
           env.key = "path",
           diag.label = FALSE) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  #remove_axis("x") +
  scale_color_manual(values = c("#19A078","#a50f15","#08519c","#7570B4","#fccde5")) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) 
p1



ggsave(filename = paste0("output/ggcor_immune_",selectedgenus,".pdf"), width = 10,height = 8)




#########################
rm(list = ls())
library(data.table)
library(GSVA)
library(ggplot2)
library(ggcor)
library(dplyr)
selectedgenus="Oceanithermus"  ###修改菌名

load(file = "output/exprSet_immune_gsva_genus.Rdata")
load(file = "output/gsva_hallmark.Rdata")

rownames(exprSet) <- exprSet[,1]
exprSet_genus <- as.data.frame(exprSet[,colnames(exprSet) == selectedgenus, drop = FALSE]) 
exprSet_genus <- as.data.frame(t(exprSet_genus))

exprSet_genus <- as.data.frame(t(exprSet_genus))
exprSet_genus <- cbind(rownames(exprSet_genus),exprSet_genus)
colnames(exprSet_genus) <- "id"

gsva <- as.data.frame(t(gsva))
gsva <- cbind(rownames(gsva),gsva)
colnames(gsva)[1] <- "id"

gsva <- inner_join(exprSet_genus,gsva, by = "id")
rownames(gsva) <- gsva[,1]
gsva <- gsva[,-c(1:2)]
gsva <- as.data.frame(t(gsva))

exprSet_genus <- as.data.frame(exprSet[,colnames(exprSet) == selectedgenus, drop = FALSE]) 
exprSet_genus <- as.data.frame(t(exprSet_genus))

immPath.score <- rbind(gsva,exprSet_genus)

immCorSiglec15 <- NULL
for (i in rownames(immPath.score)) {
  cr <- cor.test(as.numeric(immPath.score[i,]),
                 as.numeric(exprSet_genus),
                 method = "pearson")
  immCorSiglec15 <- rbind.data.frame(immCorSiglec15,
                                     data.frame(gene = selectedgenus, 
                                                path = i,
                                                r = cr$estimate,
                                                p = cr$p.value,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
}
immCorSiglec15$sign <- ifelse(immCorSiglec15$r > 0,"pos","neg")
immCorSiglec15$absR <- abs(immCorSiglec15$r)
immCorSiglec15$rSeg <- as.character(cut(immCorSiglec15$absR,c(0,0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"),include.lowest = T))
immCorSiglec15$pSeg <- as.character(cut(immCorSiglec15$p,c(0,0.001,0.01,0.05,1),labels = c("<0.001","<0.01","<0.05","ns"),include.lowest = T))
immCorSiglec15[nrow(immCorSiglec15),"pSeg"] <- "Not Applicable"

immCorSiglec15$rSeg <- factor(immCorSiglec15$rSeg, levels = c("0.25","0.50","0.75","1.00"))
immCorSiglec15$pSeg <- factor(immCorSiglec15$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
immCorSiglec15$sign <- factor(immCorSiglec15$sign, levels = c("pos","neg"))

p1 <- quickcor(t(immPath.score), 
               type = "upper",
               show.diag = TRUE) + 
  geom_colour() +
  add_link(df = immCorSiglec15, 
           mapping = aes(colour = pSeg, size = rSeg, linetype = sign),
           spec.key = "gene",
           env.key = "path",
           diag.label = FALSE) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  #remove_axis("x") +
  scale_color_manual(values = c("#19A078","#a50f15","#08519c","#7570B4","#fccde5")) +
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) 
p1

ggsave(filename = paste0("output/ggcor_hallmark_",selectedgenus,".pdf"), width = 16,height = 11)
