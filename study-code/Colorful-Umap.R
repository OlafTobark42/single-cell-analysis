DimPlot(seu, group.by = "celltype")+ scale_color
library(SingleR)
#look cellex , first mouse cell 粗筛, then immune data label=label.fine
#ref <- MouseRNAseqData()
load("~/projects/celldex/MouseRNAseqData.Rdata")
#testdata <- GetAssayData(seu, slot="data")
#clusters <- seu@meta.data$RNA_snn_res.0.5
#clusters
pred.mus <- SingleR(test = GetAssayData(seu, slot="data"), 
                    ref = ref, 
                    labels = ref$label.main, 
                    clusters = seu@meta.data$seurat_clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")
#pred.mus
celltype = data.frame(ClusterID=rownames(pred.mus), celltype=pred.mus$labels, stringsAsFactors = F)
dir.create("./results")
write.csv(celltype,"results/celltype.csv",row.names = F)
seu@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  seu@meta.data[which(seu@meta.data$seurat_clusters == celltype$ClusterID[i]),
                  'celltype'] <- celltype$celltype[i]}

unique(seu$celltype)

# 转化因子顺序 ------------------------------------------------------------------

# 确保 celltype 是一个因子
seu$celltype <- as.character(seu$celltype)

# 修改 "Granulocytes" 为 "Neutrophils"
seu$celltype[seu$celltype == "Granulocytes"] <- "Neutrophils"

# 确定新的顺序，将 "Microglia" 放到最前面
new_levels <- c(c( "Neutrophils","Microglia"),
                setdiff(unique(seu$celltype), c("Microglia", "Neutrophils")))

# 转换 celltype 为 factor，并应用新的顺序
seu$celltype <- factor(seu$celltype, levels = new_levels)

# 验证修改后的 celltype
levels(seu$celltype)

# 多配色方案col.list -------------------------------------------------------------------

source("code/Spatial_MCAO/codebase/supporting_functions_MCAO.R")

# 加载必要的库
library(Seurat)
library(ggplot2)

# 确定输出目录
output_dir <- "visual/umap"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 遍历 col.list，生成每种配色方案的 DimPlot
for (palette_name in names(col.list)) {
  # 提取当前的配色方案
  current_palette <- col.list[[palette_name]]
  
  # # 确保配色数量与 celltype 数量匹配
  # celltypes <- unique(Idents(seu))
  # if (length(current_palette) < length(celltypes)) {
  #   stop(paste("配色方案", palette_name, "中的颜色数量不足以覆盖所有 celltype"))
  # }
  
  # 为 DimPlot 设置自定义配色
  dimplot <- DimPlot(
    seu,
    group.by = "celltype",
    label = TRUE,       # 显示标签
    label.size = 4      # 标签字体大小
  ) +
    scale_color_manual(values = current_palette[1:length(celltypes)]) + # 应用配色方案
    theme_minimal() +  # 使用简单的背景主题
    theme(
      legend.position = "right",         # 图例显示在右侧
      panel.grid = element_blank(),      # 去掉网格线
      axis.text = element_blank(),       # 去掉坐标轴文字
      axis.ticks = element_blank(),      # 去掉坐标轴刻度
      axis.title = element_blank()       # 去掉坐标轴标题
    )
  
  # 保存为 PDF，文件名包含配色方案的名称
  output_file <- file.path(output_dir, paste0("DimPlot_", palette_name, ".pdf"))
  ggsave(
    filename = output_file,
    plot = dimplot,
    width = 8,  # 设置 PDF 宽度
    height = 6  # 设置 PDF 高度
  )
  
  message("已保存: ", output_file)
}


# ggsci和rcolorbrewer的多配色方案 ------------------------------------------------

# 加载必要的库
library(Seurat)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

# 设置输出目录
output_dir <- "visual/umap-colorful"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 获取所有细胞类型
celltypes <- unique(seu$celltype)
n_celltypes <- length(celltypes)

# 获取 ggsci 包中所有可用的配色方案
ggsci_palettes <- list(
  "ggsci_npg" = pal_npg("nrc"),
  "ggsci_aaas" = pal_aaas("default"),
  "ggsci_jco" = pal_jco("default"),
  "ggsci_lancet" = pal_lancet("lanonc"),
  "ggsci_ucscgb" = pal_ucscgb("default"),
  "ggsci_d3" = pal_d3("category20"),
  "ggsci_material_blue" = pal_material("blue"), # 示例颜色主题
  "ggsci_material_red" = pal_material("red"),   # 示例颜色主题
  "ggsci_futurama" = pal_futurama("planetexpress"),
  "ggsci_rickandmorty" = pal_rickandmorty("schwifty")
)
ggsci_palettes[[1]]

# 加载必要的库
library(RColorBrewer)

# 获取 RColorBrewer 中所有 "qual" 和 "div" 类型的配色方案
rcb_palettes <- rownames(brewer.pal.info[brewer.pal.info$category %in% c("qual", "div"), ])

# 获取每个配色方案的最大颜色数
rcb_palette_info <- brewer.pal.info[rownames(brewer.pal.info) %in% rcb_palettes, "maxcolors"]

# 检查结果
rcb_palettes
rcb_palette_info


# 创建补充颜色函数
# 创建补充颜色函数
get_filled_palette <- function(colors, n_needed) {
  # 如果颜色不足，补充随机颜色；否则返回前 n_needed 个颜色
  if (length(colors) < n_needed) {
    extra_colors <- grDevices::colors()[sample(1:length(grDevices::colors()), n_needed - length(colors))]
    c(colors, extra_colors)
  } else {
    colors[1:n_needed]
  }
}

# 遍历 ggsci 配色方案
for (palette_name in names(ggsci_palettes)) {
  palette_function <- ggsci_palettes[[palette_name]]
  
  # 获取当前配色方案颜色，补足颜色数量
  palette_colors <- get_filled_palette(palette_function(n_celltypes), n_celltypes)
  
  # 创建 DimPlot 并应用配色方案
  dimplot <- DimPlot(
    seu,
    group.by = "celltype",
    label = TRUE,       # 显示标签
    label.size = 4      # 标签字体大小
  ) +
    scale_color_manual(values = palette_colors) + # 应用配色方案
    theme_minimal() +  # 使用简单的背景主题
    theme(
      legend.position = "right",         # 图例显示在右侧
      panel.grid = element_blank(),      # 去掉网格线
      axis.text = element_blank(),       # 去掉坐标轴文字
      axis.ticks = element_blank(),      # 去掉坐标轴刻度
      axis.title = element_blank(),      # 去掉坐标轴标题
      plot.title = element_blank()       # 去掉 UMAP 图标题
    )
  
  # 保存为 PDF，文件名包含配色方案的名称
  output_file <- file.path(output_dir, paste0("DimPlot_", palette_name, ".pdf"))
  ggsave(
    filename = output_file,
    plot = dimplot,
    width = 8,  # 设置 PDF 宽度
    height = 6  # 设置 PDF 高度
  )
  
  message("已保存: ", output_file)
}

# 遍历 RColorBrewer 配色方案
for (palette_name in rcb_palettes) {
  max_colors <- brewer.pal.info[palette_name, "maxcolors"]
  
  # 如果配色方案颜色不足，则补充；否则使用完整配色
  if (max_colors < n_celltypes) {
    message("配色方案 ", palette_name, " 颜色不足，补充随机颜色。")
    palette_colors <- get_filled_palette(brewer.pal(max_colors, palette_name), n_celltypes)
  } else {
    palette_colors <- brewer.pal(n_celltypes, palette_name)
  }
  
  # 创建 DimPlot 并应用配色方案
  dimplot <- DimPlot(
    seu,
    group.by = "celltype",
    label = TRUE,       # 显示标签
    label.size = 4      # 标签字体大小
  ) +
    scale_color_manual(values = palette_colors) + # 应用配色方案
    theme_minimal() +  # 使用简单的背景主题
    theme(
      legend.position = "right",         # 图例显示在右侧
      panel.grid = element_blank(),      # 去掉网格线
      axis.text = element_blank(),       # 去掉坐标轴文字
      axis.ticks = element_blank(),      # 去掉坐标轴刻度
      axis.title = element_blank(),      # 去掉坐标轴标题
      plot.title = element_blank()       # 去掉 UMAP 图标题
    )
  
  # 保存为 PDF，文件名包含配色方案的名称
  output_file <- file.path(output_dir, paste0("DimPlot_RColorBrewer_", palette_name, ".pdf"))
  ggsave(
    filename = output_file,
    plot = dimplot,
    width = 8,  # 设置 PDF 宽度
    height = 6  # 设置 PDF 高度
  )
  
  message("已保存: ", output_file)
}


# split.by=GSE ------------------------------------------------------------
names(ggsci_palettes)
palette_function <- ggsci_palettes[["ggsci_npg"]]
# 获取当前配色方案颜色，补足颜色数量
palette_colors <- get_filled_palette(palette_function(n_celltypes), n_celltypes)
palette_colors

dimplot <- DimPlot(
  seu,
  group.by = "celltype", split.by = "geo_accession",
  label = F       # no显示标签
  # label.size = 4      # 标签字体大小
) +
  scale_color_manual(values = palette_colors) + # 应用配色方案
  theme_minimal() +  # 使用简单的背景主题
  theme(
    legend.position = "right",         # 图例显示在右侧
    panel.grid = element_blank(),      # 去掉网格线
    axis.text = element_blank(),       # 去掉坐标轴文字
    axis.ticks = element_blank(),      # 去掉坐标轴刻度
    axis.title = element_blank(),      # 去掉坐标轴标题
    plot.title = element_blank()       # 去掉 UMAP 图标题
  )

# seu$geo_accession
# DimPlot(seu,split.by = "geo_accession", group.by = "celltype")
ggsave("visual/umap-split.pdf", width = 20, height = 8)

