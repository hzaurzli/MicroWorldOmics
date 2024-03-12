library(shiny)

ui <- tagList(
  fluidPage(
    titlePanel("Time prediction based bacterial abundance"),
    sidebarLayout(
      sidebarPanel(
        # uiOutput 做上传文件的 ui, 对应后面的 output$file1
        uiOutput('file1'),
        # 对应后面的 output$file2
        uiOutput('file2'),
        
        uiOutput('file3'),
        
        numericInput("num", "Number of comps",
                     value = 6, min = 2, max = 15),
        
        actionButton('reset', 'RESET'),
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
        h4("Time prediction table"),
        br(),
        br(),
        shinycssloaders::withSpinner(
          dataTableOutput("table")
        )
      )
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  otu_train <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T,row.names = 1)
  })
  
  
  timer <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T)
  })
  
  
  otu_predict <- reactive({
    infile <- input$file3
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
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file2 <- renderUI({
      fileInput("file2", "Step 2: Choose time table",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  

  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file3 <- renderUI({
      fileInput("file3", "Step 3: Choose prediction otu abundance matrices",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table <- renderDataTable({
    library(pls)
    
    otu_train <- otu_train()
    otu_predict <- otu_predict()
    timer <- timer()
    
    if(is.null(otu_train) | is.null(otu_predict) | is.null(timer)){
      warning("Please upload files!")
    } 
    else{
      cpdata = otu_train
      cptimes = timer
      predict_dat = otu_predict
      
      t_cpdata = as.matrix(t(cpdata))
      t_cpdata_frame = data.frame(cptimes, t_cpdata=I(t_cpdata))
      # pls建模
      
      cptimes = cptimes[,-1]
      mod = plsr(cptimes~t_cpdata,ncomp = input$num, scale = TRUE, 
                 validation = "none")
      summary(mod)
      
      predT = predict(mod,ncomp=2,type='response',
                      newdata = t(predict_dat))
      
      predT_tmp = data.frame(predT)
      predT_table <<- data.frame(Time_point = row.names(predT_tmp), 
                                 Time_val = predT_tmp[,1])
      
    }
  }, options = list(pageLength = 10))
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(predT_table,file,row.names = T,quote = F)
    }
  )
  
  
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50735)

