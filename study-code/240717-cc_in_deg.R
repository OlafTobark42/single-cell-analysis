# 240717 提取细胞互作差异L_R对，检测在deg中的高低表达
#合并cellchat对象:
object.list <- list(Sham = cc_sham,
                    MCAO = cc_mcao)  # 对照组sham在前，比较组mcao在后，注意顺序
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

# 绘图
# 4.配受体对水平通讯差异
### 4.1 总配受体对概率差异气泡图
levels(cellchat@idents$joint) #查看细胞亚群
cellchat@group
levels(cellchat@idents)
netVisual_bubble(cellchat,
                 sources.use = "B cells",
                 targets.use = c("T cells",  "Dendritic cells"),
                 comparison = c(1,2),
                 angle.x = 45)

### 4.2 区分上下调配受体对

p7 <- netVisual_bubble(cellchat,
                       sources.use = 7,
                       targets.use = c(1:2),
                       comparison = c(1, 2),
                       max.dataset = 2,
                       title.name = "Increased signaling in LS",
                       angle.x = 45,
                       remove.isolate = T) #Increased为比较组通讯概率更强的配受体对信息
p8 <- netVisual_bubble(cellchat,
                       sources.use = 7,
                       targets.use = c(1:2),
                       comparison = c(1, 2),
                       max.dataset = 1,
                       title.name = "Decreased signaling in LS",
                       angle.x = 45,
                       remove.isolate = T) #Decreased为对照组通讯概率更强的配受体对信息
p7 + p8
# 在气泡图中，横坐标为细胞对，颜色区分样本；纵坐标为配受体。
# 气泡的大小表示p值，p值越小气泡越大。颜色表示通讯概率的大小。

# 确定源和目标细胞群组
sources.use <- "T cells"  # T细胞的名称
targets.use <- "B cells"  # B细胞的名称

# 提取Sham组的互作数据
df.sham <- subsetCommunication(cc_sham, sources.use = sources.use, targets.use = targets.use)

# 提取MCAO组的互作数据
df.mcao <- subsetCommunication(cc_mcao, sources.use = sources.use, targets.use = targets.use)

# 将两个组的数据合并在一起
df.sham$group <- "Sham"
df.mcao$group <- "MCAO"
df_combined <- rbind(df.sham, df.mcao)

# 比较Sham组和MCAO组之间的互作差异
df_combined_diff <- merge(df.sham, df.mcao, by = c("ligand", "receptor"), suffixes = c(".sham", ".mcao"))

# 计算互作差异
df_combined_diff$diff <- df_combined_diff$prob.mcao - df_combined_diff$prob.sham

# 确定互作在MCAO组中的通讯是增强还是减弱
df_combined_diff$change <- ifelse(df_combined_diff$diff > 0, "增强", "减弱")

# 查看结果
print(df_combined_diff)


# 分别在DEG中检查它们的表达差异

# 定义配体和受体基因列表
ligands <- df_combined_diff$ligand
receptors <- df_combined_diff$receptor

dir.create("visual/cc_deg")
# 绘制并保存配体基因在T细胞中的表达差异图
for (gene in ligands) {
  p <- VlnPlot(seu, features = gene, group.by = "orig.ident") + 
    ggtitle(paste("Expression of", gene, "in T cells (Sham vs MCAO)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = paste0("visual/cc_deg/","VlnPlot_", gene, "_T_cells.png"), plot = p)
}

# 绘制并保存受体基因在B细胞中的表达差异图
for (gene in receptors) {
  p <- VlnPlot(seu, features = gene, group.by = "orig.ident") + 
    ggtitle(paste("Expression of", gene, "in B cells (Sham vs MCAO)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = paste0("visual/cc_deg/","VlnPlot_", gene, "_B_cells.png"), plot = p)
}

# 安装和加载必要的包
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
library(Seurat)

# 定义配体和受体基因列表
ligands <- df_combined_diff$ligand
receptors <- df_combined_diff$receptor

# 筛选出T细胞和B细胞的子集
seu_T_cells <- subset(seu, idents = "T cells", subset = singler.mouseRNA.fine == "T cells")
seu_B_cells <- subset(seu, idents = "B cells", subset = singler.mouseRNA.fine == "B cells")

# 绘制并保存配体基因在T细胞中的表达差异图
for (gene in ligands) {
  p <- VlnPlot(seu_T_cells, features = gene, group.by = "orig.ident") + 
    ggtitle(paste("Expression of", gene, "in T cells (Sham vs MCAO)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = paste0("visual/cc_deg/","VlnPlot2_", gene, "_T_cells.png"), plot = p, width = 8, height = 6)
}

# 绘制并保存受体基因在B细胞中的表达差异图
for (gene in receptors) {
  p <- VlnPlot(seu_B_cells, features = gene, group.by = "orig.ident") + 
    ggtitle(paste("Expression of", gene, "in B cells (Sham vs MCAO)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = paste0("visual/cc_deg/","VlnPlot2_", gene, "_B_cells.png"), plot = p, width = 8, height = 6)
}


