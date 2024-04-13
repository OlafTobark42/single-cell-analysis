# 09-make-fig1c featp
unique(scRNA$celltype)

a <- readline(prompt = "input markers you want to show")
input_features <- strsplit(a, ",")[[1]]
#features = c("Cx3cr1","P2ry12","Hexb","S100a8","Fn1","Plac8","Ms4a1","Cd3d")

##Fig1c 

fig1c <- FeaturePlot(scRNA,features = input_features,ncol=4)
#rm(scRNA)
#save.image("rds/makefig.rds")
ggsave("visual/fig1c.pdf",fig1c,width = 16,height = 9)


