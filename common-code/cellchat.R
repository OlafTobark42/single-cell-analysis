##==## Part 9: Cellchat 细胞通讯 以选取T cell亚群为例


library(CellChat)
library(Seurat)
library(SeuratData)

scRNA <- readRDS("scRNA_annoted.rds")
#data <- GetAssayData(object = pbmc3k.final, slot = 'data')
cellchat <- createCellChat(object = scRNA,
                           meta =  scRNA@meta.data,
                           group.by = 'celltype',
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
detectCores()
plan(strategy = 'multisession', workers = 40)
plan("sequential")
#subset
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat) # 相当于suerat中的FindAllMarkers
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse) #储存上一步的结果到cellchat@LR$LRsig

cellchat <- computeCommunProb(cellchat)   #core procedure take some time
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)  ###自己输出一把，很好理解
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))##配体细胞1群2群，受体细胞4群5群
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))##WNT及TGFb相关网络

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
saveRDS(cellchat,"cellchat_done.rds")
##Visulize
dir.create("./cellchat")
net_count <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
net_weight <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
saveRDS(net_count,"cellchat/net_count.rds")
#ggsave("cellchat/net_count.pdf",net_count,dpi=600)
#ggsave("cellchat/net_weight.png",net_weight,dpi=600)
#net_weight

mat <- cellchat@net$weight
#par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])###通過edge.weight.max参数比较不同网络的边缘权重
}


levels(cellchat@idents)            #查看细胞顺序
vertex.receiver = c(3, 6)          #指定靶细胞的索引
cellchat@netP$pathways             #查看富集到的信号通路
pathways.show <- c("CCL","CXCL","ICAM","CD226")            #指定需要展示的通路
# Hierarchy plot
png(filename = "cellchat/sig_pathway_hierarchy.png", width = 1000, height = 650)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = 2)
dev.off()
# Circle plot
png(filename = "cellchat/sig_pathway_cricle.png", width = 650, height = 600)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
dev.off()
# 计算配体-受体对信号网络的贡献度
png(filename = "cellchat/sig_pathway_L-R.png", width = 800, height = 600)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
# 分析细胞在信号网络中角色
cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")
png(filename = "cellchat/sig_pathway_role.png", width = 800, height = 600)
netVisual_signalingRole(cellchat, signaling = pathways.show)
dev.off()

# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


saveRDS(cellchat, file = "cellchat.rds")

plotGeneExpression(cellchat, signaling = "MID")
violin_plot <- plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)
ggsave("cellchat/cscl.svg",violin_plot)
#?aggregateNet
#output
aggregateNet(
  object,
  sources.use = NULL,
  targets.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  remove.isolate = TRUE,
  thresh = 0.05,
  return.object = TRUE
)

cellchat <- readRDS("cellchat_done.rds")
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 6
cellchat <- identifyCommunicationPatterns(object = cellchat,
                                          pattern = "outgoing",
                                          k = nPatterns)
# river plot
riverplot <- netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
dotplot <- netAnalysis_dot(cellchat, pattern = "outgoing")
ggsave("cellchat/riverplot.png",riverplot,height = 20,width = 20,dpi = "retina")
ggsave("cellchat/riverplot.pdf",riverplot,height = 20,width = 20,dpi = "retina")
ggsave("cellchat/dotplot.pdf",dotplot,height = 20,width = 20,dpi = "retina")


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(1,2,4,6)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# 气泡图（全部配体-受体）

levels(cellchat@idents)
# show all the significant interactions (L-R pairs)

# 需要制定受体细胞核配体细胞
p1=netVisual_bubble(cellchat, sources.use = c(3,4,5),
                    targets.use = c(1,7),remove.isolate = FALSE)
ggsave("cellchat/L-R.pdf",p1,height = 10,width = 8,dpi = "retina")

#比如指定CCL和CXCL这两个信号通路
p2 <- netVisual_bubble(cellchat,sources.use =c(3,4,5),targets.use = c(1,7),
                       signaling = pathways.show,remove.isolate = F)
ggsave("cellchat/L-R_CCL.pdf",p2,height = 10,width = 8,dpi = "retina")



