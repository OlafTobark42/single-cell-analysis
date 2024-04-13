##==## Part 8: Monocle 拟时序分析 以选取T cell亚群为例

#monocle order function error so have to load from changed monocle file
devtools::load_all("/home/data/t220329/r/x86_64-pc-linux-gnu-library/4.2/monocle")
library(tidyverse)
library(patchwork)
#rm(list=ls())

head(scRNAsub)
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=40, relative_expr = TRUE)
#mycds <- detectGenes(mycds, min_expr = 2)  #很多教程不用
saveRDS(mycds,"mycds.rds")
mycds <- readRDS("mycds.rds")

#选择发育差异表达基因,选择clusters差异表达基因,选择离散程度高的基因,自定义发育marker基因

##使用clusters差异表达基因
diff.genes <- read.csv('cell_identify/all_markers_diff_genes_wilcox.csv')
diff.genes <- subset(diff.genes,p_val_adj<0.01)$gene
mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)
##使用seurat选择的高变基因
var.genes <- VariableFeatures(scRNAsub)
mycds <- setOrderingFilter(mycds, var.genes)
p2 <- plot_ordering_genes(mycds)
##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
##结果对比
p1|p2|p3

dir.create("./pseudotime")
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)

saveRDS(mycds,"ordered_reduced.rds")
# mycds <- readRDS("mycds_ImmGen.rds")

dir.create("./pseudotime/Tcell")
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("pseudotime/Tcell/State.pdf", plot = plot1, width = 6, height = 5)
ggsave("pseudotime/Tcell/State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("pseudotime/Tcell/Cluster.pdf", plot = plot2, width = 6, height = 5)
ggsave("pseudotime/Tcell/Cluster.png", plot = plot2, width = 6, height = 5)
##Cluster轨迹分布图
plot4 <- plot_cell_trajectory(mycds, color_by = "celltype")
ggsave("pseudotime/Tcell/Celltype_imm.pdf", plot = plot4, width = 6, height = 5)
ggsave("pseudotime/Tcell/Celltype_imm.png", plot = plot4, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("pseudotime/Tcell/Pseudotime.pdf", plot = plot3, width = 6, height = 5)
ggsave("pseudotime/Tcell/Pseudotime.png", plot = plot3, width = 6, height = 5)
##合并作图
plotc <- plot1|plot2|plot3
ggsave("pseudotime/Tcell/Combination.pdf", plot = plotc, width = 10, height = 3.5)
ggsave("pseudotime/Tcell/Combination.png", plot = plotc, width = 10, height = 3.5)
##保存结果
write.csv(pData(mycds), "pseudotime/Tcell/pseudotime.csv")

s.genes <- c("Itga4","Ccr7","Rpl12","Rps4x","Vim")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("pseudotime/Tcell/plot_genes_jitter.pdf", plot = p1, width = 10, height = 3.5)
ggsave("pseudotime/Tcell/plot_genes_violin.pdf", plot = p2, width = 10, height = 3.5)
ggsave("pseudotime/Tcell/plot_genes_in_pseudotime.pdf", plot = p3, width = 10, height = 3.5)
ggsave("pseudotime/Tcell/genes_visual.png", plot = plotc, width = 8, height = 4.5)

