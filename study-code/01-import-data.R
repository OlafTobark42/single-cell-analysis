# 01-import-data from cellranger "barcode" "matrix "features" to Seurat Object

library(Seurat)
# First Specify the clasee of scRNAlist
scRNAlist <- list()
for (i in 1:length(ds.dir)) {
  scRNAlist[[i]] <- CreateSeuratObject(counts=Read10X(data.dir = ds.dir[i]), 
                                       project=names(ds.dir)[i], min.cells=3, min.features = 200)
  # 给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = names(ds.dir)[i])   
  # 计算线粒体基因比例
  if (T) {    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")   # Hm: MT
  }
  # 计算核糖体基因比例
  if (T) {
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^Rp[sl]")  # Hm: RP[SL]
  }
  #计算红细胞基因比例
  if (T) {
    library(stringr)
    HB.genes <- str_to_title(tolower(c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")))
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
  }
  # this subset step should be in 02-qc
  # scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt < 10) 
}

# Save scRNAlist to RDS file
if (dir.exists("rds")) saveRDS(scRNAlist,"rds/scRNAlist.rds") else {
  dir.create("rds")
  saveRDS(scRNAlist,"rds/scRNAlist.rds")
}

# Merge multi-samples into one object
scRNA <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
# Save file
if (dir.exists("rds")) saveRDS(scRNA,"rds/scRNA-raw-merge.rds") else {
  dir.create("rds")
  saveRDS(scRNA,"rds/scRNA-raw-merge.rds")
}

# Clean Up scRNAlist
rm(scRNAlist)
gc()
