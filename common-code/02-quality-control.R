#02-quality-control
library(patchwork)
library(ggplot2)
if (!dir.exists("QC")) dir.create("QC")

### 绘制质控小提琴图
# 设置可能用到的主题
theme.set = theme(axis.title.x=element_blank())
# 设置绘图元素
plot.featrures = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb", "percent.HB")
col.num <- length(unique(scRNA$orig.ident))  # How many samples

# Before QC
dim(scRNA)   #查看基因数和细胞总数
table(scRNA@meta.data$orig.ident)  #查看每个样本的细胞数
# 质控前小提琴图
for (i in seq_along(plot.featrures)) {
  # Single
  vln.plot = VlnPlot(scRNA, group.by="orig.ident", pt.size = 0,
                       features = plot.featrures[i]) + theme.set + NoLegend()
  # y="nCount_RNA", x=one feature
  fp.plot <- FeatureScatter(scRNA, feature1 ="nCount_RNA", feature2 = plot.featrures[i])
  
  ggsave(paste0("QC/before_vln_",plot.featrures[i],".pdf"), 
         plot = vln.plot, width = 10, height = 8)
  ggsave(paste0("QC/before_fp_",plot.featrures[i],".pdf"), 
         plot = fp.plot, width = 10, height = 8)
  }


# Set quality control standards
minGene <- 500; maxGene <- 4000; pctMT <- 10; maxUMI <- 15000; pctHB <- 1

### 数据质控并绘制小提琴图
scRNA <- subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)

# After QC
dim(scRNA)   #查看基因数和细胞总数
table(scRNA@meta.data$orig.ident)  #查看每个样本的细胞数
# 质控前小提琴图
for (i in seq_along(plot.featrures)) {
  # Single
  vln.plot = VlnPlot(scRNA, group.by="orig.ident", pt.size = 0,
                     features = plot.featrures[i]) + theme.set + NoLegend()
  # y="nCount_RNA", x=one feature
  fp.plot <- FeatureScatter(scRNA, feature1 ="nCount_RNA", feature2 = plot.featrures[i])
  
  ggsave(paste0("QC/after_vln_",plot.featrures[i],".pdf"), 
         plot = vln.plot, width = 10, height = 8)
  ggsave(paste0("QC/after_fp_",plot.featrures[i],".pdf"), 
         plot = fp.plot, width = 10, height = 8)
  }
