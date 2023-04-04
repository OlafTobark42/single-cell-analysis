# Judge a dataset
SumSce <- function(sce){
  head(sce)
  print(table(sce$orig.ident))
  print(table(sce$celltype))
  print(table(sce$orig.ident,sce$celltype))
  DimPlot(sce, group.by = "celltype",
          repel = T, label = T)
  
}

