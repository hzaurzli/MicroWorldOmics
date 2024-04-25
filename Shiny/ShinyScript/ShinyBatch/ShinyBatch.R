library(shiny)

ui <- fluidPage(
  titlePanel("PLSDA-batch (correct batch effects in microbiome data)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Step 1: Choose abundance CSV File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      fileInput("file2", "Step 2: Choose metadata CSV File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      hr(),
      selectInput(inputId ="Detection", label ="Step 3: Batch effect detection", 
                  choices = c('PCA','Boxplots','density plots','Heatmap')),
      actionButton("update", "Selected"),  
      hr(),
      sliderInput("ncomp_trt", "Step 4: The number of treatment associated dimensions",
                  min = 0, max = 10,
                  value = 1),
      sliderInput("ncomp_bat", "Step 5: The number of batch associated dimensions",
                  min = 0, max = 10,
                  value = 4),
      sliderInput("keepX_trt", "Step 6: The number of variables to keep in X-loadings, sPLSDA_batch only",
                  min = 0, max = 500,
                  value = 100),
      hr(),
      radioButtons("DownloadType", "Final step: Download type (matrix)",
                   c("PLSDA"="p", "sPLSDA"="s")
      ),
      downloadButton("downloadData", "Download"),
      hr(),
      h5('Developer:'),
      h6(' EvaYiwenWang (method), Small runze (shiny app)'),
      br(),
      h5('Github: '),
      h6('https://github.com/EvaYiwenWang (EvaYiwenWang)'),
      h6('https://github.com/hzaurzli (Small runze)'),
      br(),
      h5('Cition:'),
      h6('PLSDA-batch: a multivariate framework to correct for batch effects in microbiome data')
    ),
    mainPanel(
      h4("Before correcting"),
      br(),
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("detectfig")),
                  tabPanel("Summary (pRDA)", verbatimTextOutput("pRDAtable"))
      ),
      br(),
      hr(),
      br(),
      h4("After correcting"),
      br(),
      tabsetPanel(type = "tab",
                  tabPanel("PLSDA-batch (PCA)", plotOutput("PLSDA")),
                  tabPanel("sPLSDA-batch (PCA)", plotOutput("sPLSDA"))
      )
    )
  )
)



server <- function(input, output, session) {
  filedata <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,row.names = 1)
  })
  
  metadata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,row.names = 1)
  })
  
  library(PLSDAbatch)
  library(pheatmap)
  library(vegan)
  library(gridExtra)
  library(mixOmics)
  library(Biobase)
  
  terms <- reactive({
    # Change when the "update" button is pressed...
    input$update
    # ...but not for anything else
    isolate({
      withProgress({
        setProgress(message = "Processing corpus...")
        input$Detection
      })
    })
  })
  
  
  datasetInput <- reactive({
    data = filedata() 
    metadata = metadata()
    
    if(is.null(data) | is.null(metadata)){
      warning("Please upload files!")
    } 
    else{
      ad.count <- as.matrix(data)
      
      ad.metadata <- data.frame(metadata)
      ad.batch = factor(ad.metadata$sequencing_run_date, 
                        levels = unique(ad.metadata$sequencing_run_date))
      ad.trt = as.factor(ad.metadata$initial_phenol_concentration.regroup)
      names(ad.batch) <- names(ad.trt) <- rownames(ad.metadata)
      
      ad.filter.res <- PreFL(data = ad.count)
      ad.filter <- ad.filter.res$data.filter
      
      ad.filter.res$zero.prob
      sum(ad.filter == 0)/(nrow(ad.filter) * ncol(ad.filter))
      
      ad.clr <- logratio.transfo(X = ad.filter, logratio = 'CLR', offset = 1) 
      class(ad.clr) = 'matrix'
      
      ad.factors.df <- data.frame(trt = ad.trt, batch = ad.batch)
      class(ad.clr) <- 'matrix'
      ad.rda.before <- varpart(ad.clr, ~ trt, ~ batch, 
                               data = ad.factors.df, scale = TRUE)
      ad.rda.before$part$indfract
    }
  })
  

  output$pRDAtable <- renderPrint({
    datasetInput()
  })
  
  
  output$detectfig <- renderPlot({

    v <- terms()
    print(v)
    
    data = filedata() 
    metadata = metadata()
    
    if(is.null(data) | is.null(metadata)){
      warning("Please upload files!")
    } 
    else{
      ad.count <- as.matrix(data)
      
      ad.metadata <- data.frame(metadata)
      ad.batch = factor(ad.metadata$sequencing_run_date, 
                        levels = unique(ad.metadata$sequencing_run_date))
      ad.trt = as.factor(ad.metadata$initial_phenol_concentration.regroup)
      names(ad.batch) <- names(ad.trt) <- rownames(ad.metadata)
      
      ad.filter.res <- PreFL(data = ad.count)
      ad.filter <- ad.filter.res$data.filter
      
      ad.filter.res$zero.prob
      sum(ad.filter == 0)/(nrow(ad.filter) * ncol(ad.filter))
      
      ad.clr <- logratio.transfo(X = ad.filter, logratio = 'CLR', offset = 1) 
      class(ad.clr) = 'matrix'
      
      ad.pca.before <- pca(ad.clr, ncomp = 3, scale = TRUE)
    }
    
    if (v[1] == "PCA") {
      if(!exists("ad.pca.before")){
        warning("Please upload files!")
      } 
      else{
        p1 = Scatter_Density(object = ad.pca.before, batch = ad.batch, trt = ad.trt, 
                             title = 'Before correcting (PCA)', trt.legend.title = 'Phenol conc.')
      }
    } 
    else if (v[1] == "Boxplots") {
      if(!exists("ad.pca.before")){
        warning("Please upload files!")
      } 
      else{
        ad.OTU.name <- selectVar(ad.pca.before, comp = 1)$name[1]
        ad.OTU_batch <- data.frame(value = ad.clr[,ad.OTU.name], batch = ad.batch)
        box_plot(df = ad.OTU_batch, title = paste(ad.OTU.name, '(Before correcting)'), 
                 x.angle = 30)
      }
    }
    else if (v[1] == "density plots") {
      if(!exists("ad.pca.before")){
        warning("Please upload files!")
      } 
      else{
        ad.OTU.name <- selectVar(ad.pca.before, comp = 1)$name[1]
        ad.OTU_batch <- data.frame(value = ad.clr[,ad.OTU.name], batch = ad.batch)
        density_plot(df = ad.OTU_batch, title = paste(ad.OTU.name, '(Before correcting)'))
      }
    }
    else if (v[1] == "Heatmap"){
      if(!exists("ad.clr")){
        warning("Please upload files!")
      } 
      else{
        ad.clr.s <- scale(ad.clr, center = TRUE, scale = TRUE)
        ad.clr.ss <- scale(t(ad.clr.s), center = TRUE, scale = TRUE)
        
        ad.anno_col <- data.frame(Batch = ad.batch, Treatment = ad.trt)
        ad.anno_colors <- list(Batch = color.mixo(seq_len(5)), 
                               Treatment = pb_color(seq_len(2)))
        names(ad.anno_colors$Batch) = levels(ad.batch)
        names(ad.anno_colors$Treatment) = levels(ad.trt)
        
        pheatmap(ad.clr.ss, 
                 cluster_rows = FALSE, 
                 fontsize_row = 4, 
                 fontsize_col = 6,
                 fontsize = 8,
                 clustering_distance_rows = 'euclidean',
                 clustering_method = 'ward.D',
                 treeheight_row = 30,
                 annotation_col = ad.anno_col,
                 annotation_colors = ad.anno_colors,
                 border_color = 'NA',
                 main = 'Before correcting - Scaled')
      }
    }
  })
  
  output$PLSDA <- renderPlot({
    data = filedata() 
    metadata = metadata()
    
    
    if(is.null(data) | is.null(metadata)){
      warning("Please upload files!")
    } 
    else{
      ad.count <- as.matrix(data)
      
      ad.metadata <- data.frame(metadata)
      ad.batch = factor(ad.metadata$sequencing_run_date, 
                        levels = unique(ad.metadata$sequencing_run_date))
      ad.trt = as.factor(ad.metadata$initial_phenol_concentration.regroup)
      names(ad.batch) <- names(ad.trt) <- rownames(ad.metadata)
      
      ad.filter.res <- PreFL(data = ad.count)
      ad.filter <- ad.filter.res$data.filter
      dim(ad.filter)
      
      ad.filter.res$zero.prob
      sum(ad.filter == 0)/(nrow(ad.filter) * ncol(ad.filter))
      
      ad.clr <- logratio.transfo(X = ad.filter, logratio = 'CLR', offset = 1) 
      class(ad.clr) = 'matrix'
      
      ad.PLSDA_batch.res <- PLSDA_batch(X = ad.clr, 
                                        Y.trt = ad.trt, Y.bat = ad.batch,
                                        ncomp.trt = input$ncomp_trt, ncomp.bat = input$ncomp_bat)
      ad.PLSDA_batch <- ad.PLSDA_batch.res$X.nobatch
      ad.pca.PLSDA_batch <- pca(ad.PLSDA_batch, ncomp = 3, scale = TRUE)
      ad.pca.PLSDA_batch.plot <- Scatter_Density(object = ad.pca.PLSDA_batch, 
                                                 batch = ad.batch, 
                                                 trt = ad.trt, 
                                                 title = 'PLSDA-batch (after correcting)')
    }
  })
  
  output$sPLSDA <- renderPlot({
    data = filedata() 
    metadata = metadata()
    
    
    if(is.null(expression) | is.null(metadata)){
      warning("Please upload files!")
    } 
    else{
      ad.count <- as.matrix(data)
      
      ad.metadata <- data.frame(metadata)
      ad.batch = factor(ad.metadata$sequencing_run_date, 
                        levels = unique(ad.metadata$sequencing_run_date))
      ad.trt = as.factor(ad.metadata$initial_phenol_concentration.regroup)
      names(ad.batch) <- names(ad.trt) <- rownames(ad.metadata)
      
      ad.filter.res <- PreFL(data = ad.count)
      ad.filter <- ad.filter.res$data.filter
      dim(ad.filter)
      
      ad.filter.res$zero.prob
      sum(ad.filter == 0)/(nrow(ad.filter) * ncol(ad.filter))
      
      ad.clr <- logratio.transfo(X = ad.filter, logratio = 'CLR', offset = 1) 
      class(ad.clr) = 'matrix'
      
      ad.sPLSDA_batch.res <- PLSDA_batch(X = ad.clr, 
                                         Y.trt = ad.trt, Y.bat = ad.batch,
                                         ncomp.trt = input$ncomp_trt, keepX.trt = input$keepX_trt,
                                         ncomp.bat = input$ncomp_bat)
      ad.sPLSDA_batch <- ad.sPLSDA_batch.res$X.nobatch
      ad.pca.sPLSDA_batch <- pca(ad.sPLSDA_batch, ncomp = 3, scale = TRUE)
      ad.pca.sPLSDA_batch.plot <- Scatter_Density(object = ad.pca.sPLSDA_batch, 
                                                  batch = ad.batch, 
                                                  trt = ad.trt, 
                                                  title = 'sPLSDA-batch (after correcting)')
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      data = filedata() 
      metadata = metadata()
      ad.count <- as.matrix(data)
      
      ad.metadata <- data.frame(metadata)
      ad.batch = factor(ad.metadata$sequencing_run_date, 
                        levels = unique(ad.metadata$sequencing_run_date))
      ad.trt = as.factor(ad.metadata$initial_phenol_concentration.regroup)
      names(ad.batch) <- names(ad.trt) <- rownames(ad.metadata)
      
      ad.filter.res <- PreFL(data = ad.count)
      ad.filter <- ad.filter.res$data.filter
      dim(ad.filter)
      
      ad.filter.res$zero.prob
      sum(ad.filter == 0)/(nrow(ad.filter) * ncol(ad.filter))
      
      ad.clr <- logratio.transfo(X = ad.filter, logratio = 'CLR', offset = 1) 
      class(ad.clr) = 'matrix'
      
      if (input$DownloadType == 'p') {
        ad.PLSDA_batch.res <- PLSDA_batch(X = ad.clr, 
                                          Y.trt = ad.trt, Y.bat = ad.batch,
                                          ncomp.trt = input$ncomp_trt, ncomp.bat = input$ncomp_bat)
        ad.PLSDA_batch <- ad.PLSDA_batch.res$X.nobatch
        
        write.csv(data.frame(ad.PLSDA_batch), file, row.names = TRUE)
      } else if (input$DownloadType == 's') {
        ad.sPLSDA_batch.res <- PLSDA_batch(X = ad.clr, 
                                           Y.trt = ad.trt, Y.bat = ad.batch,
                                           ncomp.trt = input$ncomp_trt, keepX.trt = input$keepX_trt,
                                           ncomp.bat = input$ncomp_bat)
        ad.sPLSDA_batch <- ad.sPLSDA_batch.res$X.nobatch
        
        write.csv(data.frame(ad.sPLSDA_batch), file, row.names = TRUE)
      }
    }
  )

}

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50325)

