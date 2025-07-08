library(shiny)
library(SpiecEasi)
library(igraph)

ui <- tagList(
  fluidPage(
    titlePanel("Network analysis based SpiecEasi"),
    sidebarLayout(
      sidebarPanel(
        # uiOutput 做上传文件的 ui, 对应后面的 output$file1
        uiOutput('file1'),
        
        selectInput(inputId ="method", 
                    label ="Select Method:", 
                    choices = c('MB','Glasso','Normal'),
                    selected = 'MB'),
        
        uiOutput('file2'),
        br(),
        actionButton('reset', 'RESET'),
        actionButton('start', 'START'),
        hr(),
        h5(strong('Download all results:')),
        downloadButton("downloadData", "Download Table"),
        br(),
        br(),
        radioButtons('extPlot', 'Plot output format',
                     choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
        downloadButton("downloadGraph", "Download graph"),
        hr(),
        h5('Developer:'),
        h6('Small runze (shiny app)'),
        br(),
        h5('Github: '),
        h6('https://github.com/hzaurzli (Small runze)'),
      ),
      mainPanel(
        h4("Relationship matrix and network graph"),
        shinycssloaders::withSpinner(
          dataTableOutput("table")
        ),
        hr(),
        plotOutput(outputId = "detectfig")
        
      )
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  otu_dat <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T,row.names = 1)
  })
  
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 1: Choose otu abundance matrices",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  # 通过 input$method 的值来判断应该渲染哪一个 UI
  observeEvent(input$method, {
    print(input$method)
    values$file <- NULL
    if(input$method == 'Normal'){
      output$file2 <- renderUI({
        sliderInput("num", "Corrlation cutoff in Normal pattern",
                    value = 0.3, min = 0, max = 1)
      })
    }
    else{
      output$file2 <- renderUI({
        h5(strong(paste("You have selected", input$method)))
      })
    }
  }, ignoreNULL = F)
  
  
  observeEvent(input$start, {
    
    input_method <<- input$method
    
    otu_dat = otu_dat()
    
    if(is.null(otu_dat)){
      warning("Please upload files!")
    }
    else{
      output$table <- renderDataTable({
        amgut1.filt <<- as.matrix(otu_dat)
        depths <- rowSums(amgut1.filt)
        amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
        amgut1.filt.cs <- round(amgut1.filt.n * min(depths))
        
        d <- ncol(amgut1.filt.cs)
        n <- nrow(amgut1.filt.cs)
        e <- d
        
        
        set.seed(10010)
        graph <- SpiecEasi::make_graph('cluster', d, e)
        Prec <- graph2prec(graph)
        Cor <- cov2cor(prec2cov(Prec))
        
        if (input_method == 'MB') {
          se.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2,
                                    nlambda=20, pulsar.params=list(rep.num=50))
          
          ig.mb <<- adj2igraph(getRefit(se.mb.amgut))
          df = get.data.frame(ig.mb)
          
        }
        else if (input_method == 'Glasso') {
          se.gl.amgut <- spiec.easi(amgut1.filt, method='Glasso', lambda.min.ratio=1e-2,
                                    nlambda=20, pulsar.params=list(rep.num=50))
          
          ig.gl <<- adj2igraph(getRefit(se.gl.amgut))
          df = get.data.frame(ig.gl)
          
        }
        else{
          sparcc.amgut <- sparcc(amgut1.filt)
          ## Define arbitrary threshold for SparCC correlation matrix for the graph
          sparcc.graph <- abs(sparcc.amgut$Cor) >= input$num
          diag(sparcc.graph) <- 0
          library(Matrix)
          sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
          
          ig.sparcc <<- adj2igraph(sparcc.graph)
          
          df = get.data.frame(ig.sparcc)
          
        }
        
        index = data.frame(id = seq(1:ncol(amgut1.filt)),
                           name = colnames(amgut1.filt))
        
        for (i in 1:nrow(df)) {
          df[i,1] = index[which(df[i,1] == index[,1]),2]
        }
        
        for (i in 1:nrow(df)) {
          df[i,2] = index[which(df[i,2] == index[,1]),2]
        }
        
        final <<- df[,c(1,2)]
      }, options = list(pageLength = 10))
    }
    
    output$detectfig <- renderPlot({
      if (input_method == 'MB') {
        library(igraph)
        ## set size of vertex proportional to clr-mean
        vsize <- rowMeans(clr(amgut1.filt, 1))+6
        am.coord <- layout.fruchterman.reingold(ig.mb)
        
        par(mar = c(0.1, 0, 1, 0))
        plot(ig.mb, layout=am.coord, vertex.size=vsize, 
             vertex.label=NA, main="MB")
      }
      else if (input_method == 'Glasso') {
        library(igraph)
        ## set size of vertex proportional to clr-mean
        vsize <- rowMeans(clr(amgut1.filt, 1))+6
        am.coord <- layout.fruchterman.reingold(ig.gl)
        
        par(mar = c(0.1, 0, 1, 0))
        plot(ig.gl, layout=am.coord, vertex.size=vsize, 
             vertex.label=NA, main="Glasso")
      }
      else{
        library(igraph)
        ## set size of vertex proportional to clr-mean
        vsize <- rowMeans(clr(amgut1.filt, 1))+6
        am.coord <- layout.fruchterman.reingold(ig.sparcc)
        
        par(mar = c(0.1, 0, 1, 0))
        plot(ig.sparcc, layout=am.coord, 
             vertex.size=vsize, vertex.label=NA, main="Sparcc")
      } 
      
    })
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final,file,row.names = T,quote = F)
    }
  )
  
  
  output$downloadGraph <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot)
    },
    content = function(file){
      if(input$extPlot == 'pdf'){
        pdf(file)
      }else if(input$extPlot == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      
      if (input_method == 'MB') {
        library(igraph)
        ## set size of vertex proportional to clr-mean
        vsize <- rowMeans(clr(amgut1.filt, 1))+6
        am.coord <- layout.fruchterman.reingold(ig.mb)
        
        par(mar = c(0.1, 0, 1, 0))
        plot(ig.mb, layout=am.coord, vertex.size=vsize, 
             vertex.label=NA, main="MB")
      }
      else if (input_method == 'Glasso') {
        library(igraph)
        ## set size of vertex proportional to clr-mean
        vsize <- rowMeans(clr(amgut1.filt, 1))+6
        am.coord <- layout.fruchterman.reingold(ig.gl)
        
        par(mar = c(0.1, 0, 1, 0))
        plot(ig.gl, layout=am.coord, vertex.size=vsize, 
             vertex.label=NA, main="Glasso")
      }
      else{
        library(igraph)
        ## set size of vertex proportional to clr-mean
        vsize <- rowMeans(clr(amgut1.filt, 1))+6
        am.coord <- layout.fruchterman.reingold(ig.sparcc)
        
        par(mar = c(0.1, 0, 1, 0))
        plot(ig.sparcc, layout=am.coord, 
             vertex.size=vsize, vertex.label=NA, main="Sparcc")
      } 
      
      dev.off()
    }
  )
  
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50736)
