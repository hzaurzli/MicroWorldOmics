library(shiny)

ui <- fluidPage(
  titlePanel("MCA dimensionality reduction"),
  sidebarLayout(
    sidebarPanel(
      # uiOutput 做上传文件的 ui, 对应后面的 output$file1
      uiOutput('file1'),
      # 对应后面的 output$file2
      uiOutput('file2'),
      
      actionButton('reset', 'RESET'),
      hr(),
      numericInput("num1", "Legend size",
                   value = 15, min = 0, max = 100),
      numericInput("num2", "Coordinate size",
                   value = 15, min = 0, max = 100),
      numericInput("num3", "Text one",
                   value = 20, min = 0, max = 100),
      sliderInput("size", "Points size",
                  min = 0, max = 25,
                  value = 8),
      hr(),
      radioButtons('extPlot', 'Plot output format',
                   choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
      downloadButton("plotDown", "Download figure"),
      hr(),
      downloadButton("downloadData", "Download MCA vector"),
      hr(),
      downloadButton("downloadData_rank", "Download OTU rank"),
      hr(),
      h5('Developer:'),
      h6('Small runze (shiny app)'),
      br(),
      h5('Github: '),
      h6('https://github.com/hzaurzli (Small runze)')
    ),
    mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Plot", shinycssloaders::withSpinner(plotOutput("detectfig", width = "100%", height = '600px'))),
                    tabPanel("MCA vector",shinycssloaders::withSpinner(dataTableOutput("table"))),
                    tabPanel("OTU rank",shinycssloaders::withSpinner(dataTableOutput("table_rank")))
        )
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  filedata <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  metadata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 1: Choose abundance matrix",
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
      fileInput("file2", "Step 2: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$detectfig <- renderPlot({
    library(MicroMCA)
    library(ggplot2)
    
    
    data = filedata()
    group = metadata()
    
    if(is.null(data) | is.null(group)){
      warning("Please upload files!")
    } 
    else{
        mca = MicroMCA::run_mca_dim(data)
        dat = mca$samplesCoordinates[,c(1,2)]
        dat = cbind(dat,group)
        
        p1 <<- ggplot(dat,aes(x = MCA_1, y = MCA_2, 
                              color = Group, group = Group, 
                              fill = Group)) +
                geom_point(size = input$size) + 
                theme_classic() + 
                theme(legend.title=element_text(size=input$num1),
                      axis.text=element_text(size=input$num2),
                      axis.title.x = element_text(size=input$num3),
                      axis.title.y = element_text(size=input$num3),
                      legend.text = element_text(size=input$num2))
                
        
        print(p1)
      }
    })
    
  ## 下载图片写法
  output$plotDown <- downloadHandler(
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
      # 打印全局变量 p1 (ggplot对象)
      print(p1)
      dev.off()
    }
  )
  
  output$table <- renderDataTable({
    library(MicroMCA)
    
    data = filedata()
    group = metadata()
    
    if(is.null(data) | is.null(group)){
      warning("Please upload files!")
    } 
    else{
      mca = MicroMCA::run_mca_dim(data)
      dat = data.frame(mca$samplesCoordinates[,c(1,2,3)])
      dat = data.frame(cbind(dat,group))
      
	  dat_final <<- dat
    }
  }, options = list(pageLength = 10))
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(dat_final,file,row.names = T,quote = F)
    }
  )
  
  
  output$table_rank <- renderDataTable({
    library(MicroMCA)
    
    data = filedata()
    group = metadata()
    
    if(is.null(data) | is.null(group)){
      warning("Please upload files!")
    } 
    else{
      mca = MicroMCA::run_mca_dim(data)
      imp = MicroMCA::calculate_importance(dat = mca,info = group)
      
      result = data.frame()
      
      for (i in 1:length(imp)) {
        result[i,1] = names(imp[i])
        for (j in 1:length(imp[[1]])) {
          result[i,j+1] = paste(names(imp[[i]][j]),imp[[i]][j],sep = ':')
        }
      }
      colnames(result) = c('otu',
                           paste('rank',rep(1:length(imp[[1]])),
                                 sep = '_')
                           )
      final <<- result
    }
  }, options = list(pageLength = 10))
  
  
  output$downloadData_rank <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final,file,row.names = T,quote = F)
    }
  )
}

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50633)

