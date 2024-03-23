rm(list = ls())
# 使用一个命令安装和加载多个包
packages <- c("plyr", "ggplot2", "readxl", "ggpubr", "ggsci", "tidyverse", "dplyr", "harmony",
              "rstatix", "Seurat", "reshape2", "RColorBrewer", "stringr", "writexl", "tidyr")

# 检查并安装未安装的包
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
# 加载所有包
lapply(packages, library, character.only = TRUE)

# 设置文件路径
file_path <- "./rawdata"

# 获取所有CSV文件名
file_names <- list.files(file_path, pattern = "*.csv", full.names = TRUE)

# 初始化一个空的Seurat对象列表
seurat_objects <- list()

# 遍历文件名，导入数据并创建Seurat对象
for (file_name in file_names) {
  # 读取表达矩阵
  expr_matrix <- read.csv(file_name, row.names = 1)
  
  # 从文件名中提取信息
  parts <- strsplit(basename(file_name), "_")[[1]]
  wt_or_ko <- parts[3]
  time_after_stroke <- parts[2]
  ipsil_or_contra <- gsub("_counts\\.csv", "", parts[4])
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)
  
  # 将提取的信息添加到Seurat对象的元数据中
  seurat_obj$wt_or_ko <- wt_or_ko
  seurat_obj$time_after_stroke <- time_after_stroke
  seurat_obj$ipsil_or_contra <- ipsil_or_contra
  
  # 将Seurat对象添加到列表中
  seurat_objects[[file_name]] <- seurat_obj
}

# 合并所有的Seurat对象成为一个
combined_seurat_object <- Reduce(function(x, y) merge(x, y), seurat_objects)
dim(combined_seurat_object)  # genes cells
table(combined_seurat_object$orig.ident)
# 现在combined_seurat_object包含了所有的数据和元数据
# 可以进行进一步的分析

# 假设已经创建了combined_seurat_object

# 1. 质量控制
# 选择细胞和基因的QC阈值通常取决于数据的特定情况
# combined_seurat_object <- subset(combined_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 2. 数据规范化
combined_seurat_object <- NormalizeData(combined_seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# 3. 批次效应去除
# 假设文件名中包含的批次信息已经添加到元数据中
combined_seurat_object <- FindVariableFeatures(combined_seurat_object, selection.method = "vst", nfeatures = 2000)
combined_seurat_object <- ScaleData(combined_seurat_object)
combined_seurat_object <- RunPCA(combined_seurat_object, features = VariableFeatures(object = combined_seurat_object))

# 使用Harmony进行批次校正
combined_seurat_object <- RunHarmony(combined_seurat_object, group.by.vars = "wt_or_ko", project.dim = FALSE)

# 4. 降维
combined_seurat_object <- RunUMAP(combined_seurat_object, dims = 1:20)

# 5. 聚类
combined_seurat_object <- FindNeighbors(combined_seurat_object, dims = 1:20)
combined_seurat_object <- FindClusters(combined_seurat_object, resolution = 0.5)

# 6. 分群注释
# 通常需要先了解你的数据来适当标注每个群体
# 加载所需的库
library(SingleR)
library(celldex)
# 准备用于SingleR注释的数据
# 假设你已经有了一个名为combined_seurat_object的Seurat对象
# Seurat V4
# data_for_annotation <- as.data.frame(GetAssayData(combined_seurat_object, slot = "data"))
# Seurat V5
# 直接从Seurat对象的RNA测定中获取经过归一化的表达数据
data_for_annotation <- combined_seurat_object[["RNA"]]$scale.data

# 载入参考数据集，这里以MouseRNAseqData为例
# ref_data <- celldex::MouseRNAseqData()
load("~/projects/celldex/ImmGenData.Rdata")
ref_data <- ref
# 运行SingleR进行自动注释
singleR_results <- SingleR(test = data_for_annotation, ref = ref_data, 
                           labels = ref_data$label.main)

# 查看SingleR的注释结果
table(singleR_results$labels)

# 将SingleR的注释结果添加到Seurat对象的元数据中
combined_seurat_object <- AddMetaData(combined_seurat_object, metadata = singleR_results$labels, col.name = "SingleR_annotation")

# 可以使用Seurat的可视化功能查看注释结果
DimPlot(combined_seurat_object, reduction = "umap", group.by = "SingleR_annotation")

# 可视化
DimPlot(combined_seurat_object, reduction = "umap", group.by = "seurat_clusters")

# 至此，你有了经过预处理和初步分析的Seurat对象
# 首先保存UMAP图
pdf(file = "./visual/umap_plot.pdf")
DimPlot(combined_seurat_object, reduction = "umap", group.by = "SingleR_annotation")
dev.off()

# 计算每个样本中每种细胞类型的比例
# 首先根据wt_or_ko, ipsil_or_contra, time_after_stroke分组
cell_counts <- table(combined_seurat_object$SingleR_annotation, combined_seurat_object$wt_or_ko, combined_seurat_object$ipsil_or_contra, combined_seurat_object$time_after_stroke)

# 将表格转换为数据框，以便于绘图
cell_counts_df <- as.data.frame(cell_counts)

# 为每个样本类型（wt和ko）绘制细胞比例柱状图
for (type in c("wt", "ko")) {
  # 筛选对应类型的数据
  subset_df <- subset(cell_counts_df, Var2 == type)
  
  # 计算每个组的总细胞数
  group_sums <- aggregate(Freq ~ Var2 + Var3 + Var4, data = subset_df, sum)
  
  # 合并总数回原数据框以计算比例
  subset_df <- merge(subset_df, group_sums, by = c("Var2", "Var3", "Var4"))
  
  # 计算比例
  subset_df$Freq <- with(subset_df, Freq.x / Freq.y)
  
  # 查看计算后的数据
  head(subset_df, 30)

  # 绘制柱状图
  prop_plot <- ggplot(subset_df, aes(x = Var4, y = Freq, fill = Var3)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = sprintf("%s Cell Type Proportions", type), x = "", y = "Proportion") +
    facet_wrap(~ Var1, scales = "free_y") +
    theme_minimal()
  ggsave(filename = sprintf("./visual/%s_cell_proportions.pdf", type),
         prop_plot, width = 8, height = 8)
}

save.image("rds/240321-singleR.rda")
