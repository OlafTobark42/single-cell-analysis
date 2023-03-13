# DoubletFinder

head(scRNA)
table(scRNA$orig.ident)
library(DoubletFinder)
## 寻找最优pK值
pc.num <- 1:25
sweep.res.list <- paramSweep_v3(scRNA, PCs = pc.num, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.076                     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(scRNA$celltype)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(scRNA)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
scRNA <- doubletFinder_v3(scRNA, PCs = pc.num, pN = 0.28, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
head(scRNA)
## 结果展示，分类结果在scRNA@meta.data中
tail(colnames(scRNA@meta.data),1)
DimPlot(scRNA, reduction = "umap", group.by = tail(colnames(scRNA@meta.data), 1))
#"DF.classifications_0.25_0.3_171"

unique(scRNA$DF.classifications_0.28_0.28_1294)
singlelet <- subset(scRNA, 
                    DF.classifications_0.28_0.28_1294 == "Singlet")

unique(singlelet$DF.classifications_0.28_0.28_1294)


