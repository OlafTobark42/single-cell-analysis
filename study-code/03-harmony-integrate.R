library(harmony)
scRNA <- NormalizeData(scRNA) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)
scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1))

library(clustree)
cl_tree <- clustree(scRNA)
ggsave("visual/umap/clustree.pdf", width = 16, height = 9)

Idents(scRNA) <- scRNA$RNA_snn_res.0.5

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:25)

table(scRNA$seurat_clusters)
table(scRNA$orig.ident)
saveRDS(scRNA,"rds/scRNA-harmony-cluster.rds")
