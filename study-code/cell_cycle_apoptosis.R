library(Seurat);library(ggplot2)
seu <- readRDS('/Users/zacchan/Library/Mobile Documents/com~apple~CloudDocs/古熊组/RNA水平/单细胞项目/spleen-scRNA-seq-240418jimmy备份/rds/annoted-done.rds')

head(seu)

colnames(seu@meta.data)

seu <- subset(seu, singler.mouseRNA.fine == "B cells")


# 加载细胞周期基因集
cc.genes <- Seurat::cc.genes.updated.2019

# 转换 cc.genes 为小鼠基因名格式
mouse_s_genes <- stringr::str_to_title(Seurat::cc.genes$s.genes)  # S期基因
mouse_g2m_genes <- stringr::str_to_title(Seurat::cc.genes$g2m.genes)  # G2/M期基因


# 查看数据集中的基因名
gene_names <- rownames(seu)

# 检查 S 期和 G2/M 期基因是否在数据集中
mouse_s_genes_in_data <- intersect(mouse_s_genes, gene_names)
mouse_g2m_genes_in_data <- intersect(mouse_g2m_genes, gene_names)

# 打印匹配到的数据
cat("匹配到的 S 期基因数：", length(mouse_s_genes_in_data), "\n")
cat("匹配到的 G2/M 期基因数：", length(mouse_g2m_genes_in_data), "\n")

# 运行 CellCycleScoring
seu <- CellCycleScoring(
  object = seu,
  s.features = mouse_s_genes_in_data,
  g2m.features = mouse_g2m_genes_in_data,
  set.ident = TRUE
)



# 设置分组因子
seu$Group <- factor(seu$orig.ident, levels = c("Sham", "MCAO"))

# 查看细胞周期得分的分布
VlnPlot(seu, features = c("S.Score", "G2M.Score"), group.by = "Group", pt.size = 0) + 
  # ggtitle("Cell Cycle Scores in B cells") + 
  theme_minimal()

# 进行统计分析
s_score_test <- t.test(seu$S.Score ~ seu$Group)
g2m_score_test <- t.test(seu$G2M.Score ~ seu$Group)

# 输出结果
print(s_score_test)
print(g2m_score_test)


# 统计细胞周期阶段分布
table_phase <- table(seu$Phase, seu$Group)

# 转换为百分比
prop_phase <- prop.table(table_phase, margin = 2)

# 可视化
barplot(prop_phase, beside = TRUE, col = c("skyblue", "orange", "pink"),
         legend = rownames(prop_phase),
        main = "Cell Cycle Phase Proportion in B cells",
        ylab = "Proportion")



# 凋亡 ----------------------------------------------------------------------

library(Seurat)

# 1. 筛选 B 细胞
# seu <- subset(seu, subset = celltype == "B cell")

# 2. 准备凋亡相关基因集（示例基因，可以根据研究需求调整）
apoptosis_genes <- c("Bax", "Bak1", "Casp3", "Casp7", "Tp53", "Fas", "Faslg", "Bid", "Bcl2l11", "Casp8")  # 凋亡相关基因

# 3. 检查基因是否存在于数据集中
apoptosis_genes_in_data <- intersect(apoptosis_genes, rownames(seu))

# 4. 计算凋亡评分
seu <- AddModuleScore(seu, features = list(apoptosis_genes_in_data), name = "Apoptosis_Score")

# 5. 可视化不同组（如 Sham vs MCAO）之间的凋亡评分
VlnPlot(seu, features = "Apoptosis_Score1", group.by = "orig.ident", pt.size = 0) + 
  ggtitle("Apoptosis Score in B cells") + 
  theme_minimal()

