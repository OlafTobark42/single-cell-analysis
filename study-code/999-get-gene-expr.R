## gene expression
"S100A8" %in% rownames(sce)
sce@active.ident
sce@meta.data["S100A8",]
scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == 0),'manual'] <- "Microglia"
head(scRNA)
sce@meta.data[which(AverageExpression(s))]
AverageExpression(sce, "S100A8",)

which(FetchData(sce, "S100A8"))
expr <- AverageExpression(sce, assays = "SCT", features = "S100A8", group.by = "copykat.pred")
expr[[1]][1,"aneuploid"]
expr[[1]][1,"diploid"]
sce@active.ident
sce@meta.data[which(FetchData(sce, "S100A8")>expr[[1]][1,"aneuploid"] & 
                      sce$copykat.pred == "aneuploid") , "S100A8"] <- "high"
sce@meta.data[which(FetchData(sce, "S100A8")<expr[[1]][1,"aneuploid"] & 
                      sce$copykat.pred == "aneuploid") , "S100A8"] <- "low"

table(sce$S100A8)
sce$S100A8 <- "NA"
ifelse()
FetchData(sce, vars = c("SCT", "copykat.pred"))
AverageExpression(sce, assays = "SCT", features = "S100A8", group.by = "copykat.pred")
sce[["SCT"]]@counts["S100A8",]>0

sce@meta.data[which(sce[["SCT"]]@counts["S100A8",] > expr[[1]][1,"aneuploid"] & 
                      sce$copykat.pred == "aneuploid"), "S100A8"] <- "High"
sce@meta.data[which(sce[["SCT"]]@counts["S100A8",] < expr[[1]][1,"aneuploid"] & 
                      sce$copykat.pred == "aneuploid"), "S100A8"] <- "Low"
sce@meta.data[which(sce[["SCT"]]@counts["S100A8",] > expr[[1]][1,"diploid"] & 
                      sce$copykat.pred == "diploid"), "S100A8"] <- "High"
sce@meta.data[which(sce[["SCT"]]@counts["S100A8",] < expr[[1]][1,"diploid"] & 
                      sce$copykat.pred == "diploid"), "S100A8"] <- "Low"
table(sce$S100A8)
Idents(sce) <- sce$S100A8
deg <- FindMarkers(sce, ident.1 = "High", ident.2 = "Low")  # ident.2 是 sham ??!!
dir.create("result")
write.csv(deg, "result/s100a8-deg.csv") 
saveRDS(deg, "rds/deg.rds")
