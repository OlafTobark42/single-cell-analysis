# RNA velocity
BiocManager::install("SeuratWrappers")
remotes::install_github("satijalab/seurat-wrappers")
library(devtools)
library(velocyto.R)
library(SeuratWrappers)

sce <- readRDS("annoted.rds")
bm <- RunVelocity(sce, spliced = deltaT = 1, kCells =25,
                  fit.quantile = 0.02)
sce@meta.data
