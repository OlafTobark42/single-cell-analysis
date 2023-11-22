sessionInfo()
##作图
library(RColorBrewer)
library(paletteer)
d_palettes<- palettes_d_names
which(d_palettes$length==length(unique(scRNA$manual)))
alt <- d_palettes[which(d_palettes$length==length(unique(scRNA$manual))),]
# export many pics to choose color
if (T) {
  for (i in 1:nrow(alt)) {
    dir.create("visual/fig1/choose_color")
    col_set <- paletteer_d(paste0(alt[i,1],"::",alt[i,2]))
    umap.singler <- DimPlot(scRNA, reduction = "umap", group.by='manual',
                            repel=T, label=T, label.size=5,
                            cols = col_set)+labs(title = "")
    ggsave(paste0("visual/fig1/choose_color/",i,".pdf"),
           umap.singler, height = 8, width = 8)
    
  }
}

# col4all package
if (F) {
  remotes::install_github("mtennekes/cols4all")
  BiocManager::install("colorblindcheck")
  library(cols4all)
  # c4a_gui()
  # can't show colors
  
  palette.colors()
  
}

i <- 1;  col_set <- paletteer_d(paste0(alt[i,1],"::",alt[i,2]))
library(ggplot2)
dir.create("visual/fig1")
library(RColorBrewer)
# col_set <- brewer.pal(n=length(unique(scRNA$celltype)),name = "Set2")
# Fig 1a1
library(ggsci)
# extract plot data and merge with celltype
if (T) {
  umap12 <- Embeddings(scRNA, reduction = "umap") %>% as.data.frame() %>%
    cbind(cell_type = scRNA@meta.data$manual)
  
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
  label[3,] <- data.frame(lab = paste0(format(dim(scRNA)[2], big.mark = ","),
                                       " cells"),
                          angle = 0, x = mid, y = lower + 1)
  
}

# traditional plot method 
if (F) {
  umap.singler <- DimPlot(scRNA, reduction = "umap", group.by='celltype',
                          repel=T, label=T, label.size=5,
                          cols = col_set)+labs(title = "")
}


p <- DimPlot(scRNA,reduction = 'umap',label = T, group.by = "manual",
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

ggsave("visual/fig1/col-good.pdf", p, width = 8, height = 8)


# Fig 1a2
umap.orig <- DimPlot(scRNA, reduction = "umap", group.by='orig.ident',
                     repel=F, label=F, label.size=5)+labs(title = "")
ggsave("visual/fig1/orig.pdf", umap.orig, width = 16, height = 9)
#Fig 1a3
umap.seur <- DimPlot(scRNA, reduction = "umap", group.by='seurat_clusters',
                     repel=T, label=T, label.size=5)+labs(title = "")
ggsave("visual/fig1/seur.pdf", umap.seur, width = 16, height = 9)
#Fig 1a4
umap.singler.split <- DimPlot(scRNA, reduction = "umap", group.by='celltype',
                              split.by = "orig.ident", ncol = 2, 
                              repel=T, label=T, label.size=5)+labs(title = "")
ggsave("visual/fig1/singler-split.pdf", umap.singler.split, width = 16, height = 9)


## Fig1b proportion
## Cell Proportion
cellpro <- table(scRNA$orig.ident,scRNA$manual) %>% proportions(1) %>% as.data.frame()
# cellprop <- proportions(cellprop,1)
# cellpro <- as.data.frame(cellprop)
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
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, vjust = -6))+
    ggtitle(sample.orig)
 
  ggsave(paste0("visual/fig1/",sample.orig,".pdf"),
         fig.prop[[sample.orig]],
         width = 16, height = 9)
  
}

prop.plot <- fig.prop[[1]]

for (i in 2:length(fig.prop)) {
  prop.plot <- prop.plot+fig.prop[[i]]
  
}

ggsave("visual/fig1/patch.pdf", prop.plot, width = 16, height = 9)


