library(shiny)

# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Single-cell RNA-seq database of ischemic stroke"),
  
  # Sidebar layout with a input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Selector for choosing dataset ----
      selectInput(inputId = "gene",
                  label = "Choose a gene:",
                  choices = rownames(sce)),
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "exprPlot")
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  

  # Show the first "n" observations ----
  output$exprPlot <- renderPlot({
    VlnPlot(sce,features = input$gene,
                             split.by = "orig.ident",
                             group.by = "celltype", 
            pt.size = 0)+geom_boxplot()
    })
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)





