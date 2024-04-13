##==## Part 5: Find Markers
## find all markers
# dir.create("cell_markers")
markers_seur = FindAllMarkers(scRNA, only.pos = TRUE)
dir.create("results")
write.csv(markers_seur,"results/all_markers_diff_genes_wilcox.csv",row.names=F)

library(tidyverse)
markers_seur <- markers_seur[markers_seur$p_val<0.05,]
top50 = markers_seur %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top10 = markers_seur %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "results/top10_diff_genes_wilcox.csv", row.names = F)
write.csv(top50, "results/top50_diff_genes_wilcox.csv", row.names = F)
write.csv(markers_seur,"results/all_markers_p-sig.csv",row.names=F)

