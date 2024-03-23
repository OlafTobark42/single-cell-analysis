rm(list = ls())
# 使用一个命令安装和加载多个包
packages <- c("plyr", "ggplot2", "readxl", "ggpubr", "ggsci", "tidyverse", "dplyr", "harmony",
              "rstatix", "Seurat", "reshape2", "RColorBrewer", "stringr", "writexl", "tidyr")

# 检查并安装未安装的包
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
# 加载所有包
lapply(packages, library, character.only = TRUE)

load("~/r/d1d2-197731/rds/240321-singleR.rda")

# 已有Seurat对象 combined_seurat_object，筛选wt和ko的h48 ipsil样本
wt_ipsil <- subset(combined_seurat_object, 
                   subset = wt_or_ko == 'wt' & 
                     ipsil_or_contra == 'ip' & time_after_stroke == "h48")
ko_ipsil <- subset(combined_seurat_object, 
                   subset = wt_or_ko == 'ko' & 
                     ipsil_or_contra == 'ip' & time_after_stroke == "h48")

# 合并wt和ko以进行差异表达分析
ipsil_combined <- merge(wt_ipsil, ko_ipsil)
# 合并数据层以进行差异表达分析
ipsil_combined <- JoinLayers(ipsil_combined)

# 人工调整注释 ------------------------------------------------------------------

# 重新分类前，最好先备份原始注释
ipsil_combined$orig_SingleR_annotation <- ipsil_combined$SingleR_annotation

# 修改分类
ipsil_combined$SingleR_annotation <- ifelse(ipsil_combined$SingleR_annotation == "B cells, pro", "B cells", ipsil_combined$SingleR_annotation)
ipsil_combined$SingleR_annotation <- ifelse(ipsil_combined$SingleR_annotation == "Tgd", "T cells", ipsil_combined$SingleR_annotation)
ipsil_combined$SingleR_annotation <- ifelse(ipsil_combined$SingleR_annotation == "Basophils", "Neutrophils", ipsil_combined$SingleR_annotation)

# 移除"Epithelial cells"分类的细胞
ipsil_combined <- subset(ipsil_combined, SingleR_annotation != "Epithelial cells")

# 检查修改后的细胞分类统计
table(ipsil_combined$SingleR_annotation, ipsil_combined$wt_or_ko)

# 设置标识符
Idents(ipsil_combined) <- ipsil_combined$wt_or_ko

# 合并数据层以进行差异表达分析
ipsil_combined <- JoinLayers(ipsil_combined)

# 计算所有细胞类型的差异表达基因
de_genes_all <- FindMarkers(ipsil_combined, ident.1 = 'ko', ident.2 = 'wt', min.pct = 0.25, logfc.threshold = 0.25)

# 创建文件夹结构
dir.create("ipsil-deg", showWarnings = FALSE)
dir.create("ipsil-deg/visual", showWarnings = FALSE)
dir.create("ipsil-deg/result", showWarnings = FALSE)

# 导出差异表达基因结果到CSV
write.csv(de_genes_all, "ipsil-deg/result/de_genes_h48_wt_vs_ko.csv")

# 对差异表达基因进行排序，分别获取高表达和低表达的前15个基因
top_genes_up <- rownames(head(de_genes_all[order(de_genes_all$avg_log2FC, decreasing = TRUE), ], 15))
top_genes_down <- rownames(head(de_genes_all[order(de_genes_all$avg_log2FC, decreasing = FALSE), ], 15))

# 合并高表达和低表达基因为绘图所用的基因列表
top_genes <- c(top_genes_up, top_genes_down)
top_genes

# 整体热图 --------------------------------------------------------------------
# 计算平均表达量
avg_expr <- AverageExpression(ipsil_combined, return.seurat = TRUE)
# 获取平均表达矩阵
avg_expr_data <- GetAssayData(avg_expr, assay = "RNA", slot = "data")
# avg_expr_data 是包含平均表达值的矩阵
# 首先计算每个基因在所有条件下的平均表达值
genes_mean <- rowMeans(avg_expr_data)

# 然后计算每个基因在每个条件下与其平均表达值的差
diff_expr_data <- sweep(avg_expr_data, 1, genes_mean, FUN = "-")
# 现在 diff_expr_data 包含的是每个基因在每个条件下相对于其平均表达值的差异
# 确保 top_genes 中的基因都存在于 avg_expr_data 中
valid_genes <- top_genes[top_genes %in% rownames(avg_expr_data)]

# 如果你希望查看哪些基因是有效的，可以打印 valid_genes
print(valid_genes)

# 然后使用 valid_genes 索引 diff_expr_data
diff_expr_data <- diff_expr_data[valid_genes, ]

# 使用pheatmap绘制热图
# 定义颜色：红色代表上调，蓝色代表下调
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(10)

# 绘制热图
pheatmap(diff_expr_data, color = colors, border_color = NA, display_numbers = FALSE, 
         cluster_rows = TRUE, cluster_cols = F)

#virdis方案
library(viridis)
library(pheatmap)
# 获取viridis颜色映射
colors <- viridis(10)

# 使用pheatmap绘制热图，这次使用viridis颜色映射
p <- pheatmap(diff_expr_data, color = colors, border_color = NA, display_numbers = FALSE, 
         cluster_rows = TRUE, cluster_cols = FALSE)
print(p)
ggsave("ipsil-deg/visual/热图配色选择/virdis.pdf", p, width = 8, height = 10)

# 配色方案选择 ------------------------------------------------------------------
library(RColorBrewer)
# 假设ipsil_combined和top_genes已经定义
# 创建存储热图的目录
dir.create("ipsil-deg/visual/热图配色选择", showWarnings = FALSE)

# 定义序列和发散调色板名称
sequential_palettes <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", 
                         "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")
diverging_palettes <- c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")

# 生成并保存热图
for (palette in c(sequential_palettes, diverging_palettes)) {
  colors <- brewer.pal(7, palette)
  p <- # 绘制热图
    pheatmap(diff_expr_data, color = colors, border_color = NA, display_numbers = FALSE, cluster_rows = TRUE, cluster_cols = F)
  pdf(file = paste0("ipsil-deg/visual/热图配色选择/", palette, ".pdf"))
  print(p)
  dev.off()
}

# 各细胞类型Top-DEG和热图 ---------------------------------------------------------
# 对每个组合进行差异表达基因分析并生成热图
# 按细胞类型分别计算差异基因，并导出结果
Idents(ipsil_combined) <- ipsil_combined$wt_or_ko
cell_types <- unique(ipsil_combined$SingleR_annotation)
table(ipsil_combined$SingleR_annotation, ipsil_combined$wt_or_ko)
for (cell_type in cell_types) {
  # 分别计算每种细胞类型的差异表达基因
  # sub_data <- subset(ipsil_combined, subset = cell_type == ipsil_combined$SingleR_annotation)
  # 假设你已经定义了一个变量cell_type，它包含了你想要筛选的细胞类型的名称
  # 正确的调用方法应该像这样：
  sub_data <- subset(ipsil_combined, subset = SingleR_annotation == cell_type)
  
  de_genes <- FindMarkers(sub_data, ident.1 = 'ko', ident.2 = 'wt', min.pct = 0.25, logfc.threshold = 0.25)
  
  # 保存差异基因数据
  write.csv(de_genes, sprintf("ipsil-deg/result/de_genes_%s.csv", cell_type))
  
  # 对差异表达基因进行排序，分别获取高表达和低表达的前15个基因
  top_genes_up <- rownames(head(de_genes_all[order(de_genes_all$avg_log2FC, decreasing = TRUE), ], 15))
  top_genes_down <- rownames(head(de_genes_all[order(de_genes_all$avg_log2FC, decreasing = FALSE), ], 15))
  
  # 合并高表达和低表达基因为绘图所用的基因列表
  top_genes <- c(top_genes_up, top_genes_down)
  
  # 导出pdf
  pdf(sprintf("ipsil-deg/visual/heatmap_%s.pdf", cell_type))
  # 获取viridis颜色映射
  colors <- viridis(10)
  # 使用pheatmap绘制热图，这次使用viridis颜色映射
  pheatmap(diff_expr_data, color = colors, border_color = NA, display_numbers = FALSE, 
                cluster_rows = TRUE, cluster_cols = FALSE)
  dev.off()
  
  # 循环中提示
  print(paste0(cell_type, " heatmap produced"))
}

# 富集分析 --------------------------------------------------------------------

# 富集分析及可视化
library(clusterProfiler)

# 总体的富集分析
de_genes_list <- rownames(de_genes_all[de_genes_all$p_val_adj < 0.05, ])
go_enrichment_all <- enrichGO(gene = de_genes_list, OrgDb = 'org.Mm.eg.db', 
                              keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", 
                              qvalueCutoff = 0.05)

# 导出富集结果到CSV
write.csv(as.data.frame(go_enrichment_all), "ipsil-deg/result/go_enrichment_all.csv")

# 富集结果可视化并导出为PDF
pdf("ipsil-deg/visual/go_enrichment_all.pdf")
barplot(go_enrichment_all, showCategory = 20)
dev.off()

prefix <- "ipsil-deg/"
# 按细胞类型进行富集分析，并导出结果
Idents(ipsil_combined) <- ipsil_combined$wt_or_ko
cell_types <- unique(ipsil_combined$SingleR_annotation)
table(ipsil_combined$SingleR_annotation, ipsil_combined$wt_or_ko)
for (cell_type in cell_types) {
  # 分别计算每种细胞类型的差异表达基因
  sub_data <- subset(ipsil_combined, subset = SingleR_annotation == cell_type)
  de_genes_by_type <- FindMarkers(sub_data, ident.1 = 'ko', ident.2 = 'wt', min.pct = 0.25, logfc.threshold = 0.25)
  
  de_genes_list <- rownames(de_genes_by_type[de_genes_by_type$p_val < 0.05,])
  go_enrichment <- enrichGO(gene = de_genes_list, OrgDb = 'org.Mm.eg.db', keyType = "SYMBOL", 
                            ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  # 导出结果
  write.csv(as.data.frame(go_enrichment), paste0(prefix, "result/",
    "go_enrichment_", cell_type, ".csv"))
  
  # 可视化
  pdf(paste0(prefix, "visual/",
    "go_enrichment_", cell_type, ".pdf"))
  barplot(go_enrichment, showCategory = 20)
  dev.off()
  
  # 循环中提示
  print(paste0(cell_type, " GO enriched"))
}
