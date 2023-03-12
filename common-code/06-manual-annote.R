## Annotation 

marker_list <- list(
  "DC"=c("H2-Eb1","H2-Ab1","H2-Aa","Cd74","Vim"),
  "Monocytes"=c("Vim","Ly6c2","Plac8","Ifitm3","S100a4"),
  "T cell"=c("Cd3d","Trbc2","Cd3e","Trac","Ms4a1","Ms4a4b","Ccl5"),
  "B cell"=c("Ms4a1","Igkc","Cd79a","Cd79b"),
  "Neutrophil"=c("S100a8","Retnlg","S100a9","Mmp8","Mmp9"),
  "Granulo"=c("Lyz2","Fn1","Hp","S100a8"),
  "Macrophage"=c("Ifit3","Isg15","Cd74","H2-Ab1"),
  "NK"=c("Nkg7","AW112010","Klrk1","Prf1","Gzma")
)
names(marker_list)
dir.create("visual/markers")
for (cellt in names(marker_list)) {
  vlnp <- VlnPlot(scRNA, features = marker_list[[cellt]], 
                  pt.size = 0, group.by = "celltype",ncol = 2)
  ggsave(paste0("visual/markers/vln_",cellt,".pdf"), vlnp,
         width = 16, height = 9)
  featp <- FeaturePlot(scRNA, features = marker_list[[cellt]],
                       ncol = 2)
  ggsave(paste0("visual/markers/featp_",cellt,".pdf"), featp,
         width = 16, height = 9)
}

