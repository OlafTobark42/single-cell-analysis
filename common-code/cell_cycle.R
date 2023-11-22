# Cell Cycle 
mouse_cc_genes <- readRDS("rds/mouse_cell_cycle_genes.rds")
str(mouse_cc_genes)
s.genes=mouse_cell_cycle_genes[[1]]
g2m.genes=mouse_cell_cycle_genes[[2]]
scRNA <- CellCycleScoring(scRNA, s.features = s.genes, g2m.features = g2m.genes)
head(scRNA)
scRNA@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()

DimPlot(scRNA,reduction = "umap" ,group.by = "Phase", split.by = "orig.ident")
# DimPlot(PBMC.all2,reduction = "tsne" ,group.by = "Phase")
FeaturePlot(scRNA,features = c("Pcna","Mki67"),reduction = "umap")
FeaturePlot(PBMC.all2,features = c("Pcna","Mki67"),reduction = "tsne")

Idents(scRNA) <- scRNA$singler.mouseRNA
table(scRNA$singler.mouseRNA)
table(scRNA$singler.immgendata.fine)

RidgePlot(scRNA, 
          features = c("Pcna", "Top2a", "Mcm6", "Mki67"), 
          cols = pal_npg("nrc", alpha = 0.7)(10),
          ncol = 2)

RidgePlot(scRNA, "Mki67")
Idents(scRNA) <- scRNA$Phase

table(scRNA$Phase,scRNA$orig.ident ,scRNA$singler.mouseRNA) %>% proportions(1)

t_sub <- subset(scRNA, singler.mouseRNA == "T cells")  # ident 不加引号
table(t_sub$orig.ident, t_sub$Phase) %>% proportions(1)

b_sub <- subset(scRNA, singler.mouseRNA == "B cells")
table(b_sub$orig.ident, b_sub$Phase) %>% proportions(1)
RidgePlot(b_sub, "Pcna")

