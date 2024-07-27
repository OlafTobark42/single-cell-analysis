# 2023-11-03 Restart
rm(list = ls())
sessionInfo()
set.seed(seed = 1234)


# Step 1: Data Import -----------------------------------------------------

ds.dir=list.files('./inputdata',full.names=T)

# Classic Import
if (F) {
  scRNAlist = lapply(fs,function(folder){ 
    CreateSeuratObject(counts = Read10X(folder), 
                       project = folder )
  })
}
# Import with QC
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
# Save into RDS
if (F) {
  # Save scRNAlist to RDS file
  if (!dir.exists("rds")) dir.create("rds")
  # saveRDS(scRNAlist,"rds/scRNAlist.rds")
  
  # Merge multi-samples into one object
  scRNA <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
  # Save file
  if (!dir.exists("rds")) dir.create("rds")
  saveRDS(scRNA,"rds/scRNA-raw-merge.rds")
}

# Change Order
table(scRNA$orig.ident)
scRNA$orig.ident <- factor(scRNA$orig.ident,
                           levels = c("Naive", "Sham", "MCAO"))

# Step 2: Quality Control ####
if (T){
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
  for (i in seq_along(plot.featrures)) {  # seq_along生成与对象长度相同的序列，可以等同为1:length（）
    # Single
    vln.plot = VlnPlot(scRNA, group.by="orig.ident", pt.size = 0,
                       features = plot.featrures[i]) + theme.set + NoLegend()
    # y="nCount_RNA", x=one feature
    fp.plot <- FeatureScatter(scRNA, feature1 ="nCount_RNA", feature2 = plot.featrures[i])
    
    ggsave(paste0("QC/before_vln_",plot.featrures[i],".pdf"), 
           plot = vln.plot, width = 10, height = 7)
    ggsave(paste0("QC/before_fp_",plot.featrures[i],".pdf"), 
           plot = fp.plot, width = 10, height = 7)
  }
  
  
  # Set quality control standards
  minGene <- 500; maxGene <- 4000; pctMT <- 10; maxUMI <- 15000; pctHB <- 1
  
  ### 数据质控并绘制小提琴图
  scRNA <- subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                    nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)
  
  # After QC
  dim(scRNA)   #查看基因数和细胞总数
  table(scRNA@meta.data$orig.ident)  #查看每个样本的细胞数
  # 质控后小提琴图
  for (i in seq_along(plot.featrures)) {
    # Single
    vln.plot = VlnPlot(scRNA, group.by="orig.ident", pt.size = 0,
                       features = plot.featrures[i]) + theme.set + NoLegend()
    # y="nCount_RNA", x=one feature
    fp.plot <- FeatureScatter(scRNA, feature1 ="nCount_RNA", feature2 = plot.featrures[i])
    
    ggsave(paste0("QC/after_vln_",plot.featrures[i],".pdf"), 
           plot = vln.plot, width = 10, height = 7)
    ggsave(paste0("QC/after_fp_",plot.featrures[i],".pdf"), 
           plot = fp.plot, width = 10, height = 7)
  }
  
}

# Step 3: Harmony Remove Batch Effect ####
if (T) {
  library(harmony)
  scRNA <- NormalizeData(scRNA) %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(verbose=FALSE)
  scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")
  scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:25) %>% 
    FindClusters(resolution = seq(0.1, 1, by = 0.1))
  
  library(clustree)
  cl_tree <- clustree(scRNA)
  if (!dir.exists("visual/umap")) dir.create("visual/umap")
  ggsave("visual/umap/clustree.pdf", width = 10, height = 7)
  
  Idents(scRNA) <- scRNA$RNA_snn_res.0.8  # 0.5 or 0.8? see clustertree
  
  scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:25)
  scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:25)
  
  table(Idents(scRNA))
  table(scRNA$orig.ident, Idents(scRNA))
  saveRDS(scRNA,"rds/scRNA-harmony-cluster.rds")
}

# Step 4: UMAP Preview ####
DimPlot(scRNA, reduction = "tsne", label = T)
DimPlot(scRNA, reduction = "umap")

##作图
library(RColorBrewer)
library(paletteer)
d_palettes<- palettes_d_names
which(d_palettes$length==length(unique(Idents(scRNA))))
alt <- d_palettes[which(d_palettes$length==length(unique(Idents(scRNA)))),]
col_set <- paletteer_d(paste0(alt[1,1],"::",alt[1,2]))
p <- DimPlot(scRNA,reduction = 'tsne',label = T, 
             cols = col_set, repel=T, label.size=5) +
  NoAxes() +  labs(title = "") +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab), 
            fontface = "plain")  # fontface = 'italic'

# Step 5: Auto Annotation -------------------------------------------------
library(SingleR)
#look cellex, first mouse cell 粗筛, then immune data label=label.fine
#ref <- MouseRNAseqData()
load("~/projects/celldex/MouseRNAseqData.Rdata")
#testdata <- GetAssayData(scRNA, slot="data")
#clusters <- scRNA@meta.data$RNA_snn_res.0.5
#clusters
pred.mus <- SingleR(test = GetAssayData(scRNA, slot="data"), 
                    ref = ref, 
                    labels = ref$label.main, 
                    clusters = Idents(scRNA),
                    # clusters = scRNA@meta.data$RNA_snn_res.0.8,  
                    # Warning: 0.5 or 0.8 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")
# 结果转换为表格
celltype = data.frame(ClusterID=rownames(pred.mus), 
                      celltype=pred.mus$labels, 
                      stringsAsFactors = F)
#dir.create("./results")
#write.csv(celltype,"results/celltype.csv",row.names = F)

# 将注释插入进meta data
scRNA@meta.data$singler.mouseRNA = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(Idents(scRNA) == celltype$ClusterID[i]),
                  'singler.mouseRNA'] <- celltype$celltype[i]}

# 检验效果
table(scRNA$orig.ident, scRNA$singler.mouseRNA)

# 进一步label.fine
load("~/projects/celldex/ImmGenData.Rdata")
pred.mus.fine <- SingleR(test = GetAssayData(scRNA, slot="data"), 
                    ref = ref, 
                    labels = ref$label.fine, 
                    clusters = Idents(scRNA),
                    # clusters = scRNA@meta.data$RNA_snn_res.0.8,  
                    # Warning: 0.5 or 0.8 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")
celltype.fine = data.frame(ClusterID=rownames(pred.mus.fine), 
                      celltype=pred.mus.fine$labels, 
                      stringsAsFactors = F)
scRNA@meta.data$singler.immgendata.fine = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(Idents(scRNA) == celltype.fine[i,"ClusterID"]),
                  'singler.immgendata.fine'] <- celltype.fine[i,"celltype"]}
table(scRNA$orig.ident, scRNA$singler.immgendata.fine)
table(scRNA$orig.ident, Idents(scRNA))


# Step 6: UMAP Review -----------------------------------------------------
library(ggsci)
DimPlot(scRNA, reduction = "tsne",
        split.by = "orig.ident", shuffle = T)+
  ggsci::scale_color_npg()  # use ggsci

table(scRNA$orig.ident, scRNA$singler.immgendata.fine)
table(scRNA$orig.ident, scRNA$singler.mouseRNA)

DimPlot(scRNA, reduction = "umap",
        group.by = "singler.immgendata.fine",
        label = T)
library("scales")
a <- show_col(pal_npg()(10))
cols <- pal_npg()(length(unique(scRNA$singler.mouseRNA)))
show_col(cols)

# Step 7: Doublet Finder --------------------------------------------------
# DoubletFinder
head(scRNA)
table(scRNA$orig.ident)
library(DoubletFinder)
## 寻找最优pK值
pc.num <- 1:25
sweep.res.list <- paramSweep_v3(scRNA, PCs = pc.num, sct = F)
# time use 3min
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.076                     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(scRNA$singler.mouseRNA)   # 最好提供celltype
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

unique(scRNA@meta.data[,tail(colnames(scRNA@meta.data),1)])
singlelet <- subset(scRNA, 
                    DF.classifications_0.28_0.27_1275 == "Singlet")

unique(singlelet@meta.data[,tail(colnames(scRNA@meta.data),1)])
saveRDS(scRNA, file = "rds/doubletfinded.rds")



# Step 8: UMAP Review 2 ---------------------------------------------------
if (F){
  singlelet$singler.mouseRNA
  p <- DimPlot(singlelet,reduction = 'umap',
               label = T, group.by = "singler.mouseRNA",
               repel=T, label.size=5) +
    NoAxes() +  
    theme(aspect.ratio = 1)
  p  
  Idents(singlelet) <- singlelet$RNA_snn_res.0.5
  singlelet <- RunUMAP(singlelet, reduction = "harmony",
                       dims = 1:25)
  DimPlot(singlelet, reduction = "umap",
          label =T, group.by = "singler.mouseRNA")
  
  colnames(singlelet@meta.data)
}

celltype[c(13,16),2] <- "T cells"

scRNA <- subset(scRNA, 
                    DF.classifications_0.28_0.27_1275 == "Singlet")
table(scRNA$singler.mouseRNA)
Idents(scRNA) <- scRNA$RNA_snn_res.0.8
UMAPPlot(scRNA, label =T)

test <- scRNA
test$RNA_snn_res.0.8 <- factor(test$RNA_snn_res.0.8,
                               levels = c(1,21,11,12,2,15,10,19,
                                          0,3,4,7,16,
                                          6,9,5,17,
                                          13,26,24,23,8,
                                          22,25,20,14,18))
Idents(test) <- test$RNA_snn_res.0.8
UMAPPlot(test, label =T)

library(RColorBrewer)
library(paletteer)
d_palettes<- palettes_d_names
# actually >= is ok, not must ==
which(d_palettes$length>=length(unique(test$singler.immgendata.fine)))
alt <- d_palettes[which(d_palettes$length>=length(unique(test$singler.immgendata.fine))),]
# export many pics to choose color
if (T) {
  for (i in 1:nrow(alt)) {
    if (!dir.exists("visual/umap/choose_color")) dir.create("visual/umap/choose_color")
    paletter_name <- paste0(alt[i,1],"::",alt[i,2])
    col_set <- paletteer_d(paste0(alt[i,1],"::",alt[i,2]))
    umap.singler <- DimPlot(test, reduction = "umap", group.by='cs',
                            repel=T, label=T, label.size=5,
                            cols = col_set)+labs(title = paste0(alt[i,1],"::",alt[i,2]))
    ggsave(paste0("visual/umap/choose_color/",sprintf("%03d",i),"_",paletter_name,".pdf"),
           umap.singler, height = 6.18, width = 10)
  }
}

# cluster 23 Pre B is not good
test <- subset(test, singler.immgendata.fine != "B cells (preB.FrC)")

# ggthemes::Classic_20 is good
# col_set <- paletteer_d("ggthemes::Classic_20")

# Step 9: Plot UMAP  ------------------------------------------------------
library(paletteer)
col_set <- paletteer_d("ggthemes::Classic_20")
# ggsci::category20c_d3  备选
# extract plot data and merge with celltype
if (T) {
  umap12 <- Embeddings(test, reduction = "umap") %>% as.data.frame() %>%
    cbind(cell_type = test@meta.data$singler.immgendata.fine)
  
  # 构造坐标轴需要的标签和位置信息:
  # get botomn-left coord
  lower <- floor(min(min(umap12$UMAP_1),min(umap12$UMAP_2))) 
  
  # get relative line length
  linelen <- abs(0.3*lower) + lower
  
  # mid point
  mid <- abs(0.3*lower)/2 + lower
  
  # axies data
  axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                     group = c(1,1,2,2),
                     label = rep(c('UMAP_1','UMAP_2'),each = 2))
  
  # axies label
  label <- data.frame(lab = c('UMAP_1','UMAP_2'),angle = c(90,0),
                      x = c(lower - 1,mid),y = c(mid,lower - 1))
  # cell_nums <- format(dim(scRNA)[2], big.mark = ",")
  label[3,] <- data.frame(lab = paste0(format(dim(test)[2], big.mark = ","),
                                       " cells"),
                          angle = 0, x = mid, y = lower + 1)
  
}

p <- DimPlot(test,reduction = 'umap',
             label = T, 
             group.by = "singler.immgendata.fine",
             cols = col_set, 
             repel=T, label.size=5) +
  NoAxes() +  labs(title = "") +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab), 
            fontface = "plain")  # fontface = 'italic'
cells.located <- CellSelector(plot = p)
cells.located
ggsave("visual/umap/col-good-nonaive-3.pdf", p, width = 10, height = 6.18)

saveRDS(test, file = "rds/havenaive.rds")

test <- subset(test,
               orig.ident != "Naive")
test$orig.ident <- factor(test$orig.ident,
                          levels = c("Sham", "MCAO"))
table(test$orig.ident, test$singler.immgendata.fine)
table(test$singler.mouseRNA.fine, test$orig.ident)

test$singler.immgendata.fine <- factor(test$singler.immgendata.fine,
                                       levels = c("B cells (B.Fo)","B cells (B.MZ)","B cells (B.T2)","B cells (B.T3)","B cells (B1a)","DC (DC.8-4-11B-)","DC (DC.PDC.8+)",
                                                  "Macrophages (MF.RP)","Monocytes (MO.6C-II-)","Monocytes (MO.6C+II-)","Basophils (BA)","Neutrophils (GN.ARTH)","Neutrophils (GN)","NK cells (NK.DAP10-)",
                                                  "NKT (NKT.4+)","T cells (T.8Mem)","T cells (T.8Nve)","T cells (T.CD4TESTCJ)","T cells (T.Tregs)"))

Idents(test) <- test$singler.immgendata.fine
sample_split <- DimPlot(test, reduction = "umap",
        split.by = "orig.ident",
        cols = col_set) + 
  NoAxes() +  labs(title = "") +
  theme(legend.position = "none")
ggsave("visual/umap/sample_split-1.pdf", sample_split,
       width = 12, height = 6.18)

a <- levels(test$singler.immgendata.fine)
mode(a)
cat(a,sep = '\",\"')
paste(a,sep = ",")


# Step 10: Proportion ------------------------------------------------------
cellpro <- table(test$orig.ident,
                 test$singler.immgendata.fine) %>% proportions(1) %>% as.data.frame()
fig.prop <- list()
dir.create("visual/fig1")
library(ggrepel)
for (sample.orig in unique(cellpro$Var1)) {
  sample_prop <- cellpro[cellpro$Var1==sample.orig,]
  sample_prop$ymax <- cumsum(sample_prop$Freq)
  sample_prop$ymin <- c(0,head(sample_prop$ymax,n=-1))
  #设置标签
  lab=paste0(sample_prop$Var2,'\n',round(sample_prop$Freq*100,1),"%")
  lab <- paste0(round(sample_prop$Freq*100,1),"%")
  sample_prop$lab=lab
  sample_prop
  
  fig.prop[[sample.orig]] <- ggplot(sample_prop,aes(ymax=ymax,ymin=ymin,
                                                    xmax=4,xmin=2))+
    geom_rect(aes(fill=Var2))+
    theme_void()+
    xlim(0,4)+
    coord_polar(theta="y")+
    scale_fill_manual(values = col_set) +
    geom_text_repel(aes(y= (ymax + ymin)/2,x=3,label=lab),size=4)+
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5, vjust = -6))+
    ggtitle(sample.orig)
  # p1 + geom_text_repel(aes(y= (ymax + ymin)/2,x=4,label=lab),size=4)
  ggsave(paste0("visual/fig1/",sample.orig,".pdf"),
         fig.prop[[sample.orig]],
         width = 16, height = 9)
}

prop.plot <- fig.prop[[1]]

for (i in 2:length(fig.prop)) {
  prop.plot <- prop.plot+fig.prop[[i]]
  
}

ggsave("visual/fig1/patch.pdf", prop.plot, width = 16, height = 9)

### bar plot
## Cell Proportion
scRNA1 <- test
cellprop <- table(scRNA1$orig.ident,scRNA1$singler.immgendata.fine)
cellprop <- proportions(cellprop,1)
as.data.frame(cellprop)
cellprop
proportions(cellprop,1)
colourCount = length(unique(cellprop$Var1))
col_set
ggplot(as.data.frame(cellprop))+
  geom_bar(aes(x =Var1, y= Freq, fill = Var2,),
           stat = "identity",width = 0.7,size = 0.5)+ 
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = col_set)+
  NoLegend()+theme_classic()

# Step 11: FeaturePlot ----------------------------------------------------
feat_gene <- c("Cd3d", "Cd4", "Cd8a", "Foxp3",  # T cells
               "Ms4a1", "Cd79a",  # B cells
               "Cd68", "Spp1","Il1b", "Ly6g",
               "S100a8", "S100a9")

f1 <- FeaturePlot(test, "Cd3d")
f1+scale_color_viridis()
library(viridis)
# Step 12: Find All Markers -----------------------------------------------
markers_seur = FindAllMarkers(scRNA, only.pos = TRUE)
library(tidyverse)
markers_group <- markers_seur %>% group_by(cluster)
markers_top5 <- markers_seur %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
markers_top5 <- markers_top5[!duplicated(markers_top5$gene),]
dotp <- DotPlot(scRNA, features = markers_top5$gene) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("visual/feat.pdf", dotp, width = 20, height = 10)




# Step 13: DEG of each celltypes ------------------------------------------
deg_celltype <- list()
Idents(seu) <- seu$orig.ident
for (t in unique(seu$singler.immgendata.fine)) {
  deg_celltype[[t]] <- FindMarkers(subset(seu,
                                        singler.immgendata.fine == t),
                                 ident.1 = "MCAO",
                                 ident.2 = "Sham")
  
}
deg_celltype[[2]]

# Step 13: Cellchat -------------------------------------------------------
sham <- subset(seu, ident = "Sham")
mcao <- subset(seu, ident = "MCAO")
library(CellChat)
colnames(seu@meta.data)
# sham and mcao
if (T) {
  cellchat <- createCellChat(object = mcao,
                             meta =  mcao@meta.data,
                             group.by = 'singler.mouseRNA.fine',
                             assay = 'RNA')
  str(cellchat)
  
  #showDatabaseCategory()
  CellChatDB <- CellChatDB.mouse
  head(CellChatDB$interaction)
  cellchat@DB <- CellChatDB
  dplyr::glimpse(CellChatDB$interaction)
  
  subsetDB(CellChatDB.mouse, search = 'Secreted Signaling')
  subsetDB(CellChatDB.human, search = '^KEGG', key = 'evidence') #regular expression here!!!
  
  library(future)
  library(parallel)
  detectCores()
  plan(strategy = 'multisession', workers = 100)
  # plan("sequential")
  
  #subset
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat) # 相当于suerat中的FindAllMarkers
  # above time use: 10min
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse) #储存上一步的结果到cellchat@LR$LRsig
  Sys.time()
  cellchat <- computeCommunProb(cellchat)   #core procedure take some time
  Sys.time()  # 100 cores 10min
  # cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  levels(cellchat@idents)            #查看细胞顺序
  # vertex.receiver = c(3, 6)          #指定靶细胞的索引
  cellchat@netP$pathways       
}

cc_sham <- cellchat 
cc_mcao <- cellchat
save(cc_sham,cc_mcao, file = "cellchat-split.rda")
# 查看两个数据集细胞分群，并检查是否一致：
identical(levels(cc_sham@idents),
          levels(cc_mcao@idents))

#合并cellchat对象:
object.list <- list(Sham = cc_sham,
                    MCAO = cc_mcao)  # 对照组sham在前，比较组mcao在后，注意顺序
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

dplyr::glimpse(cellchat)
##合并内容:'data.signaling','images','net','netP','meta','idents','var.features','DB',and 'LR'

# 1.总体比较：通讯数目与通讯强度差异
p1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
p2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1+ ggsci::scale_fill_npg() + p2 + ggsci::scale_fill_npg()

# 2.细胞亚群水平的通讯差异
### 2.1 细胞通讯差异网络图
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# 展示了样本间所有细胞亚群中的配受体对数目差异(左)及通讯概率差异(右)；
# 外周圆不同颜色代表不同细胞亚群，大小表示亚群的配受体对数目，圈越大，细胞间配受体对数目比值越大.
# 蓝线表示对照组通讯较强，红线表示比较组通讯较强。线越粗表示通讯变化程度越强。

### 2.2 细胞通讯差异热图
p3 <- netVisual_heatmap(cellchat)
p4 <- netVisual_heatmap(cellchat, measure = "weight")
p3 + p4

# 3.信号通路水平的通讯差异
### 3.1 组间富集信号通路差异条形图
##基于信息流或互作数对信号通路进行排序
p5 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) #堆叠
p6 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE) #不堆叠
p5 + p6
# 当比较（LS）组与对照组（NL）的通路概率总和的比值＜0.95且差异pval＜0.05（秩和检验）时，
# 则该通路在对照组中的通讯强度显著增加（纵坐标为红色）；当比较组与对照组的通路概率总和的比值＞1.05且差异pval＜0.05时，
# 则该通路在比较组中的通讯强度显著增加（纵坐标为蓝色）；
# 纵坐标为黑色表示该通路在两组间没有差异。左侧为比例图，右侧为实际数值比对图。

### 3.2 传出信号通路水平热图
# complexheatmap is needed
library(CellChat)
library(ComplexHeatmap)
library(patchwork)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways,
                       object.list[[i+1]]@netP$pathways)
pathway.union

ccs <- netAnalysis_computeCentrality(cc_sham)
ccm <- netAnalysis_computeCentrality(cc_mcao)
save(ccs,ccm, file = "rds/cellchat-divided.rda")
object.list <- list(Sham = ccs,
                    MCAO = ccm)
# plot
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "outgoing", #传出
                                        signaling = pathway.union,
                                        title = names(object.list)[i],
                                        width = 5,
                                        height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                        pattern = "outgoing", #传出
                                        signaling = pathway.union,
                                        title = names(object.list)[i+1],
                                        width = 5,
                                        height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

### 3.3 传入信号通路水平热图
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "incoming", #传入
                                        signaling = pathway.union,
                                        title = names(object.list)[i],
                                        width = 5, height = 6,
                                        color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                        pattern = "incoming", #传入
                                        signaling = pathway.union,
                                        title = names(object.list)[i+1],
                                        width = 5, height = 6,
                                        color.heatmap = "GnBu")
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

### 3.4 总体信号通路水平热图
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "all", #总体
                                        signaling = pathway.union,
                                        title = names(object.list)[i],
                                        width = 5, height = 6,
                                        color.heatmap = "OrRd")
ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                        pattern = "all", #总体
                                        signaling = pathway.union,
                                        title = names(object.list)[i+1],
                                        width = 5, height = 6,
                                        color.heatmap = "OrRd")
draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))

# 纵坐标为信号通路，横坐标为细胞亚群，热图颜色代表信号强度，颜色越深，通讯越强。上侧和右侧的柱子是纵轴和横轴强度的累积。
# 左图NL组，右图LS组。

# 4.配受体对水平通讯差异
### 4.1 总配受体对概率差异气泡图
levels(cellchat@idents$joint) #查看细胞亚群
cellchat@group
levels(cellchat@idents)
netVisual_bubble(cellchat,
                 sources.use = "T cells",
                 targets.use = c("B cells",  "Dendritic cells"),
                 comparison = c(1,2),
                 angle.x = 45)

### 4.2 区分上下调配受体对

p7 <- netVisual_bubble(cellchat,
                       sources.use = 7,
                       targets.use = c(1:2),
                       comparison = c(1, 2),
                       max.dataset = 2,
                       title.name = "Increased signaling in LS",
                       angle.x = 45,
                       remove.isolate = T) #Increased为比较组通讯概率更强的配受体对信息
p8 <- netVisual_bubble(cellchat,
                       sources.use = 7,
                       targets.use = c(1:2),
                       comparison = c(1, 2),
                       max.dataset = 1,
                       title.name = "Decreased signaling in LS",
                       angle.x = 45,
                       remove.isolate = T) #Decreased为对照组通讯概率更强的配受体对信息
p7 + p8
# 在气泡图中，横坐标为细胞对，颜色区分样本；纵坐标为配受体。
# 气泡的大小表示p值，p值越小气泡越大。颜色表示通讯概率的大小。

# 5.单个/特定信号通路水平差异可视化
#使用网络图：
pathways.show <- c("CCL") #选择目标信号通路
weight.max <- getMaxWeight(object.list,
                           slot.name = c("netP"),
                           attribute = pathways.show) 
#控制不同数据集的边权重
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]],
                      signaling = pathways.show,
                      layout = "circle",
                      edge.weight.max = weight.max[1],
                      edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

#使用热图：
pathways.show <- c("CCL")
par(mfrow = c(1,2), xpd = TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ht7 <- ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# plot export
if (!dir.exists("visual/cellchat")) dir.create("visual/cellchat")
ls() %>% cat(sep = ",")
for (ht in list(ht1,ht2,ht3,ht4,ht5,ht6,ht7)) {
  pdf(file = paste0("visual/cellchat/",ht,".pdf"),
      width = 10, height = 6.18)
  ht
  dev.off()
}

