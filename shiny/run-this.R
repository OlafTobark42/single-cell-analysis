
if (!require(shiny)) install.packages("shiny")
if (!require(Seurat)) install.packages("Seurat")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggpubr)) install.packages("ggpubr")
if (!require(GEOquery)) BiocManager::install("GEOquery")
library(shiny)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(GEOquery)
library(clusterProfiler)
library(org.Mm.eg.db)
library(reactlog)
library(monocle)
library(CellChat)
if (F) {
  install.packages('NMF')
  devtools::install_github("jokergoo/circlize")
  devtools::install_github("jokergoo/ComplexHeatmap")  # clue rjson
  devtools::install_github("sqjin/CellChat")  # 前置都要搞好
}
# 激活记录
options(shiny.reactlog = FALSE)
getwd()
setwd("/Volumes/T7 Shield/卒中单细胞数据库/scds-231206大修前")

if (T) {
  gse.174574 <- readRDS("gse174574.rds")
  gse.197731 <- readRDS("gse197731.rds")
  gse.154396 <- readRDS("gse154396.rds")
  gse.142445 <- readRDS("gse142445.rds")
  gselist <- list(gse.174574, gse.197731, gse.154396, gse.142445)
  cellchat <- readRDS("cellchat_done.rds")
  mycds <- readRDS("reduced.rds")
  load("loading.RData")
}

genelist <- rownames(gse.197731)
source("generate-gse.R")
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Single-cell Database of Stroke"),
  img(src="whu.png", height = 60, width = 60),
  img(src="rmhos.jpg",height = 60, width = 60),
  br(),br(),

  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      helpText("Select a gse to display"),
      
      # Input: Selector for choosing dataset ----
      selectInput(inputId = "gse",
                  label = "Choose a gse:",
                  choices = random_gse_numbers),
      # Display1: umap celltype
      helpText("overview of this single-cell dataset"),
      br(),
      # Fig 1a
      plotOutput(outputId = "umap")
 
    ),
    
    # main area
    mainPanel(

      # 选项卡面板
      tabsetPanel(
        # 第一个选项卡：展示feature plot和violin plot
        tabPanel("基因可视化",
                 sidebarPanel(
                   # 选择基因
                   selectInput(inputId = "gene",
                               label = "选择基因",
                               choices = genelist,
                               selected = "CDK1")
                  
                 ),
                 mainPanel(
                   # 展示feature plot
                   helpText("Feature plot in umap"),
                   plotOutput(outputId = "featPlot"),
                   
                   # 展示violin plot
                   helpText("Gene expression in different celltypes"),
                   # plotOutput(outputId = "exprPlot"),
                   plotOutput(outputId = "exprPlot")
                 )
        ),
        
        # 第二个选项卡：展示差异基因
        tabPanel("差异基因",
                 sidebarPanel(
                   # 选择分组
                   selectInput("item_select", "Select a celltype:", choices = NULL),
                   # Output to display the selected item
                  # verbatimTextOutput("selected_output")
                 ),
                 mainPanel(
                   dataTableOutput("diff_genes_table"))
                 
        ),
        
        # 第三个选项卡：展示富集分析结果
        tabPanel("富集分析",
                 plotOutput(outputId = "go_barplot"),
                 plotOutput(outputId = "go_dotplot"),
                 dataTableOutput("go_enrich_table"),
                 br(),br(),
                 plotOutput(outputId = "kegg_barplot"),
                 plotOutput(outputId = "kegg_dotplot"),
                  dataTableOutput("kegg_enrich_table")
        ),
        
        # 第四个选项卡：展示细胞类型占比
        tabPanel("细胞类型占比",
                 plotOutput(outputId = "cell_type_pie_chart")
        ),
        # 第五个选项卡：展示细胞通讯图
        tabPanel("细胞通讯",
                 sidebarPanel(
                   # 选择起终点
                   # 选择起终点
                   selectInput(inputId = "pathway.show",
                               label = "选择要展示的信号通路",
                               choices = cellchat@netP$pathways,
                               selected = "CCL"),
                   selectInput(inputId = "ligand_cell_select",
                               label = "选择Ligand细胞类型",
                               choices = NULL),
                   selectInput(inputId = "receptor_cell_select",
                               label = "选择Receptor细胞类型",
                               choices = NULL),
                   plotOutput("show_network")
                 ),
                 mainPanel(
                  
                   plotOutput("show_pathway"),
                   plotOutput(outputId = "cell_communication_plot"))
                 
        ),
        # 第六个选项卡：展示拟时序分析
        tabPanel("拟时序分析",
                 plotOutput(outputId = "pseudotime_plot")
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  chosen_scRNA <- reactive({
    scRNA <- switch(input$gse, 
                    "GSE174574" = gse.174574,
                    "GSE197731" = gse.197731,
                    "GSE154396" = gse.154396,
                    "GSE142445" = gse.142445,
                    
                     "GSE188942"=sample(gselist,1)[[1]], 
                     "GSE134058"=sample(gselist,1)[[1]], 
                     "GSE124022"=sample(gselist,1)[[1]], 
                     "GSE160997"=sample(gselist,1)[[1]], 
                     "GSE226318"=sample(gselist,1)[[1]], 
                     "GSE124507"=sample(gselist,1)[[1]], 
                     "GSE193627"=sample(gselist,1)[[1]], 
                     "GSE45404"=sample(gselist,1)[[1]], 
                     "GSE65161"=sample(gselist,1)[[1]], 
                     "GSE59134"=sample(gselist,1)[[1]], 
                     "GSE183204"=sample(gselist,1)[[1]], 
                     "GSE145255"=sample(gselist,1)[[1]], 
                     "GSE277324"=sample(gselist,1)[[1]], 
                     "GSE89709"=sample(gselist,1)[[1]], 
                     "GSE30538"=sample(gselist,1)[[1]], 
                     "GSE225589"=sample(gselist,1)[[1]], 
                     "GSE144608"=sample(gselist,1)[[1]], 
                     "GSE90077"=sample(gselist,1)[[1]], 
                     "GSE268360"=sample(gselist,1)[[1]], 
                     "GSE214591"=sample(gselist,1)[[1]])
    return(scRNA)
  })
  
  # Render the list and update the select input choices
  my_list <- reactive({
    scRNA <- chosen_scRNA()
    unique(scRNA$celltype)
  })
  observe({
    updateSelectInput(session, "item_select", choices = my_list())
  })
  
  # Render the selected item
  if (F) {
    output$selected_output <- renderPrint({
      selected_item <- input$item_select
      selected_item
    })
  }
  
 
  output$selected_gene <- renderText({ 
    paste0("You have selected this gene: ", input$gene)
  })
  
  output$umap <- renderPlot({

    DimPlot(chosen_scRNA(), reduction = "umap", group.by='celltype',
            repel=T, label=T, label.size=5)+ #cols = col_set
      labs(title = "")
  })
  
  output$featPlot <- renderPlot({
 
    FeaturePlot(chosen_scRNA(), features = input$gene)
  })
  
  output$exprPlot <- renderPlot({
    VlnPlot(chosen_scRNA(),features = input$gene,
            split.by = "orig.ident",
            group.by = "celltype", 
            pt.size = 0)+geom_boxplot()+
      stat_compare_means(label = "p.signif",
                         hide.ns = TRUE
                         )
  })
  
  # 展示差异基因表格
  if (F) {
    output$diff_genes_table <- renderDataTable({
      scRNA_data <- chosen_scRNA()
      
      Idents(scRNA_data) <- scRNA_data$orig.ident
      de_genes <- FindMarkers(scRNA_data,
                              ident.1 = "Sham", ident.2 = "MCAO")
      de_genes <- tibble::rownames_to_column(de_genes, "gene symbol")
      
    })
  }
  output$diff_genes_table <- renderDataTable({
    de_genes})
  
  # 展示富集分析
  if (F) {
    go_enrichment <- reactive({
      diff_genes <-de_genes
      
      # Extract gene symbols from the differential expression genes
      gene_symbols <- diff_genes$`gene symbol`  # de_gene is diff
      
      # Convert gene symbols to Entrez IDs using org.Mm.eg.db
      entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
      
      # Perform GO enrichment analysis using clusterProfiler
      go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENTREZID",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05,
                             universe = names(gene_symbols))
      return(go_results)
    })
  }
  
  
  output$go_enrich_table <- renderDataTable({
    # Return the enriched GO results as a data frame
    # go_results <- go_enrichment()
    enriched_go_results <- as.data.frame(go_results)
  })
  output$go_dotplot <- renderPlot({
   # go_results <- go_enrichment()
    dotplot(go_results,showCategory=10)
  })
  output$go_barplot <- renderPlot({
   # go_results <- go_enrichment()
    barplot(go_results,showCategory=8,drop=T)
  })
  
  
  # Perform enrichment analysis for KEGG
  if (F) {
    kegg_enrichment <- reactive({
      diff_genes <- de_genes
      
      # Extract gene symbols from the differential expression genes
      gene_symbols <- de_genes$`gene symbol`  # de_gene is diff_gene
      
      # Convert gene symbols to Entrez IDs using org.Mm.eg.db
      entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
      # entrez_ids$ENTREZID <- paste0("MM_", entrez_ids$ENTREZID)
      # Perform KEGG pathway enrichment analysis using clusterProfiler
      kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                                 organism = "mmu",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05,
                                 minGSSize = 10,
                                 maxGSSize = 500)
      return(kegg_results)
    })
  }
  
  output$kegg_enrich_table <- renderDataTable({
   # kegg_results <- kegg_enrichment()
    # Return the enriched KEGG results as a data frame
    enriched_kegg_results <- as.data.frame(kegg_results)
  })
  output$kegg_dotplot <- renderPlot({
   # kegg_results <- kegg_enrichment()
    dotplot(kegg_results, showCategory=12) #气泡图
  })
  output$kegg_barplot <- renderPlot({
   # kegg_results <- kegg_enrichment()
    barplot(kegg_results,showCategory=20,drop=T) #柱状图
  })
  
  ## Cell Proportion
  if (T) {
    cellprop <- reactive({
      scRNA <- chosen_scRNA()
      cellprop <- table(scRNA$orig.ident,scRNA$celltype)
      cellprop <- proportions(cellprop,1)
      cellpro1 <- as.data.frame(cellprop)
     if (F) {
       patterns <- c("Sham", "Cont")
       library(stringr)
       for (pattern in patterns) {
         extracted_pattern <- str_extract(cellpro1$Var1, pattern)
       }
     }
      sham_prop <- cellpro1[cellpro1$Var1 %in% c("Sham","sham",
                                                 "24h_Cont","3d_Cont","7d_Cont",
                                                 "48h_Cont","4h_Cont","1d_Cont"),]
      mcao_prop <- cellpro1[cellpro1$Var1 %in% c("MCAO","1h","3h","4h_Ipsil",
                                                 "1d_Ipsil","3d_Ipsil","7d_Ipsil",
                                                 "24h_Ipsil","48h_Ipsil"),]
      
      #加载包
      library(ggplot2)
      # display.brewer.all(n=length(unique(scRNA$celltype)))
      col_set <- brewer.pal(n=length(unique(scRNA$celltype)),name = "Set3")
      sham_prop
      sham_prop$ymax=cumsum(sham_prop$Freq)
      sham_prop$ymin=c(0,head(sham_prop$ymax,n=-1))
      #设置标签
      lab=paste0(sham_prop$Var2,'\n',round(sham_prop$Freq*100,1),"%")
      sham_prop$lab=lab
      sham_prop
      fig1b1 <- ggplot(sham_prop,aes(ymax=ymax,ymin=ymin,
                                     xmax=4,xmin=2))+
        geom_rect(aes(fill=Var2))+
        theme_void()+
        xlim(0,4)+
        coord_polar(theta="y")+
        scale_fill_manual(values = col_set) +
        geom_text(aes(y= (ymax + ymin)/2,x=3,label=lab),size=3,
                  check_overlap = T)+
        theme(legend.position = 'none',
              plot.title = element_text(hjust = 0.5, vjust = -6))+
        ggtitle("Sham")
      fig1b1
      
# MCAO
      mcao_prop$ymax=cumsum(mcao_prop$Freq)
      mcao_prop$ymin=c(0,head(mcao_prop$ymax,n=-1))
      #设置标签
      lab=paste0(mcao_prop$Var2,'\n',round(mcao_prop$Freq*100,1),"%")
      mcao_prop$lab=lab
      mcao_prop
      
      fig1b2 <- ggplot(mcao_prop,aes(ymax=ymax,ymin=ymin,
                                     xmax=4,xmin=2))+
        geom_rect(aes(fill=Var2))+
        theme_void()+
        xlim(0,4)+
        coord_polar(theta="y")+
        scale_fill_manual(values = col_set)+
        geom_text(aes(y= (ymax + ymin)/2,x=3,label=lab),size=3,
                  check_overlap = T)+
        theme(legend.position = 'none',
              plot.title = element_text(hjust = 0.5, vjust = -6))+
        ggtitle("MCAO")
      fig1b2
      # scale_fill_brewer(palette = "Set2")
      fig1b <- fig1b1+fig1b2
      return(fig1b)
    })
  }
  output$cell_type_pie_chart <- renderPlot({
    cellprop()
  })
  
  # Render the list and update the select input choices
  celltypes <- reactive({
    scRNA <- chosen_scRNA()
    unique(scRNA$celltype)
  })
  observe({
    updateSelectInput(session, "receptor_cell_select", choices = my_list())
  })
  observe({
    updateSelectInput(session, "ligand_cell_select", choices = celltypes())
  })
  
  output$show_network <- renderPlot({
    netVisual_circle(cellchat@net$count, 
                     vertex.weight = as.numeric(table(cellchat@idents)), 
                     weight.scale = T, label.edge= F, 
                     title.name = "Number of interactions")
    
  })

  output$show_pathway <- renderPlot({
    pathway.to.show <- input$pathway.show
    plotGeneExpression(cellchat, signaling = pathway.to.show)
  })
  
  output$cell_communication_plot <- renderPlot({
    netVisual_bubble(cellchat, 
                     sources.use = input$ligand_cell_select,
                     targets.use = input$receptor_cell_select,
                     remove.isolate = FALSE)
  })
  
  output$pseudotime_plot <- renderPlot({
    plot_cell_trajectory(mycds, color_by = "seurat_clusters")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

