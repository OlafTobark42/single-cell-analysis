library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(cols4all)
new_scRNA$orig.ident <- factor(new_scRNA$orig.ident,levels = c("Sham","Ruptured"))
cellprop <- table(new_scRNA$orig.ident,new_scRNA$celltype)
cellprop <- proportions(cellprop,1)
cellprop <- as.data.frame(cellprop)
df <- cellprop
#转换为因子，指定绘图顺序：
df$Var2 <- factor(df$Var2,levels = unique(df$Var2))
p1 <- ggplot(df, aes(x = Var1, y=Freq, fill = Var2,
                     stratum = Var2, alluvium = Var2)) +
  scale_fill_npg() + 
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  geom_col(width = 0.6,
           color = 'white', size = 0.5) + #同时赋予白色描边，适当加粗
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0,
            color = 'white', size = 0.5) #同时赋予白色描边，适当加粗
pdf(file = "Cell Proportion.pdf",width = 6,height = 6)
p1
dev.off()