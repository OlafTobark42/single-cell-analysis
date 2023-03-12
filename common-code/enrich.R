## Enrichment
dir.create("enrich")

scRNA <- readRDS("scRNA.rds")
mycds <- readRDS("mycds.rds")
#比较cluster0和cluster1的差异表达基因
dge.cluster <- FindMarkers(scRNA,ident.1 = 8,ident.2 = 11)
dge.celltype <- dge.cluster
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>1)
#比较B_cell和T_cells的差异表达基因
dge.celltype <- FindMarkers(scRNA, ident.1 = "B cells", ident.2 = 'T cells', group.by = 'celltype')
sig_dge.celltype <- subset(dge.celltype, p_val_adj<0.01&abs(avg_log2FC)>1)
#比较拟时State1和State3的差异表达基因
p_data <- subset(pData(mycds),select='State')
scRNAsub <- subset(scRNA, cells=row.names(p_data))
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
dge.State <- FindMarkers(scRNAsub, ident.1 = 1, ident.2 = 3, group.by = 'State')
sig_dge.State <- subset(dge.State, p_val_adj<0.01&abs(avg_log2FC)>1)



#差异基因GO富集分析
ego_ALL <- enrichGO(gene          = row.names(sig_dge.celltype),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Mm.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(ego_all,'enrich/enrichGO.csv')           
ego_CC <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP,showCategory = 8) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 8) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 8) + ggtitle("barplot for Molecular function")
plotc <- p_BP/p_CC/p_MF
ggsave('enrich/enrichGO_celltype.png', plotc, width = 12,height = 10)
ggsave('enrich/enrichGO_celltype.png', plotc, width = 12,height = 10)

library(clusterProfiler)
library(org.Mm.eg.db)
#kegg
genelist <- bitr(row.names(degs), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Mm.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist$ENTREZID)
p1 <- barplot(ekegg, showCategory=10)
p2 <- dotplot(ekegg, showCategory=10)
plotc = p1/p2
ggsave("enrich/enrichKEGG_celltype.png", plot = plotc, width = 12, height = 10)
#save.image("DEG_Enrich.RData")

s.genes <- c("ITGB1","CCR7","KLRB1","GNLY")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("pseudotime/genes_visual.png", plot = plotc, width = 8, height = 4.5)

