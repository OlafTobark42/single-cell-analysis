
##==## Part 10: Scenic 转录因子调控网络 


library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
rm(list=ls())

install.packages("arrow")
packageVersion("AUCell")
packageVersion("RcisTarget")
packageVersion("GENIE3")

##==分析准备==##
dir.create("SCENIC")
dir.create("SCENIC/int")
scRNA <- readRDS("scRNA_annoted.rds")
setwd("./SCENIC") 
##准备细胞meta信息
cellInfo <- data.frame(scRNA@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")
##准备表达矩阵
#为了节省计算资源，随机抽取1000个细胞的数据子集
subcell <- sample(colnames(scRNA),3000)
scRNAsub <- scRNA[,subcell]
saveRDS(scRNAsub, "scRNAsub.rds")
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
##设置分析环境
mydbDIR <- "~/projects/cistarget"
mydbs <- c("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
           "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
names(mydbs) <- c("10kb","500bp")
scenicOptions <- initializeScenic(org="mgi", 
                                  nCores=8,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNSCC")
saveRDS(scenicOptions, "int/scenicOptions.rds")
