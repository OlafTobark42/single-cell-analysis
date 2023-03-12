scRNA$nsfc[scRNA$nsfc=="Granulocytes" & scRNA$orig.ident=="Sham"][1:1100] <- NA


scRNA <- scRNA[,scRNA@meta.data$orig.ident %in% c("MCAO","Sham") &
                 is.na(scRNA$nsfc)==F]


length(scRNA$nsfc[scRNA$nsfc=="Granulocytes" & scRNA$orig.ident=="Sham"])


scRNA$celltype[scRNA$celltype=="Granulocytes"] <- "Astrocytes"
unique(scRNA$celltype)
scRNA$celltype[scRNA$celltype=="NK cells"] <- "T cells"


cd4_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(0,2)]
cd4_sce2 = sce[, Idents(sce) %in% c( "Naive CD4 T" ,  "Memory CD4 T" )]
cd4_sce3 = subset(sce,seurat_clusters %in% c(0,2))
## 较简单，核心原理就是R里取子集的3种策略：逻辑值，坐标，名字

MI.integrated <- FindClusters(object = MI.integrated, resolution = seq(0.1,1,0.1))
clus.tree.out <- clustree(MI.integrated) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set3")+
  coord_flip()
clus.tree.out
