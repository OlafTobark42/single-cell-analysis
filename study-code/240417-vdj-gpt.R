# vdj-chatGPT
library(tidyverse)
library(Seurat)
# 选择一个r d s
seu <- readRDS("rds-old/scRNA_annoted.rds")
# 读取VDJ重排信息
vdj_data <- read.csv("vdj/sb/airr_rearrangement.tsv", sep = "\t")

# 读取克隆型数据
clonotypes <- read.csv("vdj/sb/clonotypes.csv")

# 假设Seurat对象的条形码在元数据中的列名为"orig.ident"
# 需要确保VDJ数据框的条形码列（这里假设列名为"barcode"）格式与Seurat对象一致
# vdj_data$barcode <- gsub(pattern = "-1", replacement = "", x = vdj_data$barcode)
# seu的行名像这样MCAO_MCAO_AAACCTGTCCAAAGTC-1 
# 所以vdj_data新增一列barcode，改下名
vdj_data$barcode <- paste0("Sham_Sham_", vdj_data$cell_id)
# 检查cell_id列是否有重复的条形码
table(vdj_data$barcode)
# 有重复的 from yunlai zhou
vdj_data <- vdj_data[!duplicated(vdj_data$barcode), ]
# another method of GPT
if (F) {library(dplyr)
  
  # 修改名称并去重
  vdj_data$barcode <- paste0("Sham_Sham_", vdj_data$cell_id)
  vdj_data <- vdj_data %>%
    distinct(barcode, .keep_all = TRUE)  # 去除重复的行
  
  # 将处理后的vdj_data添加到Seurat对象
  rownames(vdj_data) <- vdj_data$barcode
  seu <- AddMetaData(object = seu, metadata = vdj_data)
  head(seu$barcode)
}
# 以克隆型信息为例，添加到Seurat对象的元数据
# 假设Seurat对象名为seu
# 确保vdj_data的行名是条形码
rownames(vdj_data) <- vdj_data$barcode
seu <- AddMetaData(object = seu, metadata = vdj_data)

# 因为单转细胞数很多，而bcr细胞数较少
table(rownames(seu@meta.data) %in% rownames(vdj_data))





# 对细胞进行聚类
seu <- FindClusters(seu, resolution = 0.5)

# 绘制UMAP，颜色按克隆类型
DimPlot(seu, group.by = "clone_id")

FeaturePlot(seu,features="frequency")
