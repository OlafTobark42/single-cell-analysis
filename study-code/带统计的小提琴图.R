library(Seurat)
library(ggplot2)
library(ggsci)
library(ggpubr)
rm(list = ls())
seu <- readRDS("data/rawdata/mouse_blood_day2.rds")
# 定义基因和分组顺序
show_genes <- c("Mcemp1", "Clec4d", "Cacna1e")
# seu$orig.ident <- factor(seu$orig.ident, levels = c("sham", "1h", "3h"))
seu$orig.ident <- factor(seu$treatment, levels = c("Sham","MCAO"))
# seu <- subset(seu, orig.ident != "3h")
# seu$orig.ident <- factor(seu$orig.ident, levels = c("sham", "1h"))
# 创建小提琴图，并在中间添加箱线图和统计信息
VlnPlot(seu, features = show_genes, group.by = "celltype", pt.size = 0, 
               split.by = "orig.ident") &
  # ylim(0.1,3) &
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) & 
  stat_compare_means(aes(group = "orig.ident"), method = "t.test", label = "p.signif") &
  ggsci::scale_fill_d3() & 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# 比较两个不同orig.ident
vln <- VlnPlot(seu, features = "Mcemp1", group.by = "celltype", pt.size = 0,split.by = "orig.ident") +
  stat_compare_means(method = "wilcox.test", label = "p.format") &
  ggsci::scale_fill_d3() & 
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) & 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

vln

# 三组比较 --------------------------------------------------------------------
for (gene in show_genes) {
  vln <- VlnPlot(seu, features = gene, group.by = "celltype", pt.size = 0,split.by = "orig.ident") +
    stat_compare_means(method = "wilcox.test", label = "p.format") &
    ggsci::scale_fill_d3() & 
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) & 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  ggsave(filename = paste0("Vln_", gene, ".pdf"),
         vln, width = 10, height = 6)
}


# 两组比较 --------------------------------------------------------------------
for (gene in show_genes) {
  vln <- VlnPlot(seu, features = gene, group.by = "celltype", pt.size = 0,split.by = "orig.ident") +
    stat_compare_means(method = "kruskal.test", label = "p.format") &
    ggsci::scale_fill_d3() & 
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) & 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  ggsave(filename = paste0("Vln_", gene, ".pdf"),
         vln, width = 8, height = 6)
}



# 保存图像
# ggsave("Cacna1e.pdf", vln, width = 8, height = 6)



# 手动提取后绘图 -----------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr); library(tidyverse)

# 定义基因和分组顺序
show_genes <- c("Mcemp1", "Clec4d", "Cacna1e")
seu$orig.ident <- factor(seu$orig.ident, levels = c("sham", "1h", "3h"))

# 提取数据用于 ggplot
plot_data <- FetchData(seu, vars = c(show_genes, "celltype", "orig.ident"))
plot_data <- plot_data %>% 
  pivot_longer(cols = all_of(show_genes), names_to = "gene", values_to = "expression")

# 绘制小提琴图并添加箱线图和统计信息
vln <- ggplot(plot_data, aes(x = celltype, y = expression, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  facet_wrap(~gene, scales = "free_y") +
  stat_compare_means(aes(group = orig.ident), method = "t.test", label = "p.signif",
                     hide.ns = TRUE, size = 3, position = position_dodge(width = 0.9)) +
  ggsci::scale_fill_d3() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
vln
# 保存图像
ggsave("features_vln_with_stats.pdf", vln, width = 18, height = 8)

