# Judge a dataset
head(sce)
table(sce$orig.ident)
table(sce$celltype)
DimPlot(sce, group.by = "celltype",
        repel = T, label = T)

