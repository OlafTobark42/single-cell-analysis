counts <- read.csv("GSM7060815_Brain_GR180716_counts.csv")
metadata <- read.csv("GSM7060815_Brain_GR180716_metadata.csv")
head(counts,4)
str(counts)
head(metadata)

library(Seurat)
sham2d1 <- CreateSeuratObject(counts = c1, meta.data = m1)
# project name == orig.ident

table(sham2d1$parent)

c2 <- read.csv("GSM7060816_Brain_GR181128_counts.csv")
m2 <- read.csv("GSM7060816_Brain_GR181128_metadata.csv")
library(tidyverse)
c2 <- column_to_rownames(c2, var = "X")
m2 <- column_to_rownames(m2, var = "X")

sham2d2 <- CreateSeuratObject(counts = c2, meta.data = m2,
                              project = "sh2d2")
table(sham2d2$parent)
table(sham2d2$orig.ident)
head(sham2d2)

# batch import
getwd()
fs=list.files('./inputdata',full.names=T)
scRNAlist <- list()
for (dir in fs) {
  datas <- list.files(dir, full.names = T)
  # datas[1] == counts; datas[2] == metadata
  count <- read.csv(datas[1]) %>% column_to_rownames(var = "X")
  meta <- read.csv(datas[2]) %>% column_to_rownames(var = "X")
  scRNAlist[[dir]] <- CreateSeuratObject(counts = count, 
                                         meta.data = meta,
                                         project = tail(str_split(dir, "/")[[1]], 1))
}
