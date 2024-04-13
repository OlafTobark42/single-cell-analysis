# 08-Vlnplot and FeaturePlot of Markers
names(immune_gene)
dir.create("visual/vln-feat-cellt")
for (cellt in names(immune_gene)) {
  vln <- VlnPlot(scRNA,features = immune_gene[[cellt]],ncol = 4)
  feap <- FeaturePlot(scRNA,features = immune_gene[[cellt]],ncol = 4)
  ggsave(paste0("visual/vln-feat-cellt/",cellt,"_vln.pdf"),
         vln,width = 16,height = 9)
  ggsave(paste0("visual/vln-feat-cellt/",cellt,"_feap.pdf"),
         feap,width = 16,height = 9)
}


