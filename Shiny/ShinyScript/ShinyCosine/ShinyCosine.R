library(shiny)
library(igraph)
library(lsa)


ui <- tagList(
  fluidPage(
    titlePanel("Network analysis based cosine similarity"),
    sidebarLayout(
      sidebarPanel(
        # uiOutput 做上传文件的 ui, 对应后面的 output$file1
        uiOutput('file1'),
        actionButton('reset', 'RESET'),
        actionButton('start', 'START'),
        hr(),
        downloadButton("downloadData", "Download Table"),
        hr(),
        h5('Developer:'),
        h6('Small runze (shiny app)'),
        br(),
        h5('Github: '),
        h6('https://github.com/hzaurzli (Small runze)'),
      ),
      mainPanel(
        h4("Relationship matrix"),
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
  
  
  observeEvent(input$start, {
    
    otu_dat = otu_dat()
    
    if(is.null(otu_dat)){
      warning("Please upload files!")
    }
    else{
      output$table <- renderDataTable({
        final = as.matrix(t(otu_dat))
        
        myAdjacencyMatrix = as.matrix(cosine(final))
        g  <- graph.adjacency(myAdjacencyMatrix,weighted=TRUE)
        df_dat <<- get.data.frame(g)
      },options = list(pageLength = 10))
    }
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(df_dat,file,row.names = T,quote = F)
    }
  )
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50750)

