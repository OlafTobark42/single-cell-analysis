##==##       Part 4: Annotation  ##==##
library(SingleR)
#look cellex , first mouse cell 粗筛, then immune data label=label.fine
#ref <- MouseRNAseqData()
load("~/projects/celldex/ImmGenData.Rdata")
#testdata <- GetAssayData(scRNA, slot="data")
#clusters <- scRNA@meta.data$RNA_snn_res.0.5
#clusters
pred.mus <- SingleR(test = GetAssayData(scRNA, slot="data"), 
                    ref = ref, 
                    labels = ref$label.main, 
                    clusters = scRNA@meta.data$RNA_snn_res.0.5, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")
#pred.mus
celltype = data.frame(ClusterID=rownames(pred.mus), celltype=pred.mus$labels, stringsAsFactors = F)
dir.create("./results")
write.csv(celltype,"results/celltype.csv",row.names = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),
                  'celltype'] <- celltype$celltype[i]}


