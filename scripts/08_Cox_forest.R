rm(list = ls())
library(ggplot2)
#install.packages("survival")
#install.packages("survminer")
library(tidyr)
library(dplyr)
library(survival)
library("survminer")

rt=read.table("output/risk.txt",header=T,sep="\t")
colnames(rt)[1] <- "TCGA_id"
load(file = "output/metadata.Rdata")
metadata <- metadata[,c(4,12)]
colnames(metadata) <- c("TCGA_id","age")
rt <- inner_join(metadata,rt, by = "TCGA_id")


s = as.formula(paste("Surv(futime, fustat) ~ ",paste(colnames(rt)[c(2,5:20)],collapse = "+"))) ###这里需要修改
model <- coxph(s, data = rt )


options(scipen=1)
p <- ggforest(model, data =rt, 
              main = "Hazard ratio", 
              cpositions = c(0.10, 0.22, 0.4), 
              fontsize = 1.8, 
              refLabel = "1", noDigits = 4)
p
ggsave(p,file = "output/multivariateCOX_forest.pdf",height = 13,width = 25)

