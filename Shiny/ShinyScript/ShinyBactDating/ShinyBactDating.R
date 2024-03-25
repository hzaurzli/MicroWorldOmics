library(shiny)
library(BactDating)
library(ape)


ui <- fluidPage(
  titlePanel("Molecular dating by BactDating"),
  sidebarLayout(
    sidebarPanel(
      # uiOutput 做上传文件的 ui, 对应后面的 output$file1
      uiOutput('file1'),
      # 对应后面的 output$file2
      uiOutput('file2'),
      radioButtons("inSelect", "if TRUE, the most distant tip from the root is 
                                considered as the origin of the time scale;
                                if FALSE, this is the root node.",
                   c("True","False"), selected = "False"),
      radioButtons("inSelect1", "Show labels",
                   c("True","False"), selected = "False"),
      actionButton('reset', 'RESET'),
      hr(),
      radioButtons('extPlot', 'Plot output format',
                   choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
      downloadButton("plotDown", "Download molecular dating tree"),
      br(),
      br(),
      downloadButton("plotDown1", "Download temporal signal"),
      hr(),
      h5('Developer:'),
      h6('Small runze (shiny app)'),
      br(),
      h5('Github: '),
      h6('https://github.com/hzaurzli (Small runze)')
    ),
    mainPanel(
      h3('Molecular dating results'),
      shinycssloaders::withSpinner(dataTableOutput("table")),
      uiOutput('file3'),
      plotOutput("detectfig"),
      uiOutput('file4'),
      plotOutput("detectfig1")
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  treedata <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.tree(infile$datapath)
  })
  
  datedata <- reactive({
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
      fileInput("file1", "First: Choose tree (newick)",
                accept = c(
                  "newick","nwk",
                  ".nwk,text/plain",
                  ".tree")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file2 <- renderUI({
      fileInput("file2", "Step 2: Choose date file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$file <- NULL
    output$file3 <- renderUI({
      NULL
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$file <- NULL
    output$file4 <- renderUI({
      NULL
    })
  }, ignoreNULL = F)
  
  
  output$table <- renderDataTable({

    data = treedata()
    date = datedata()
    
    if(is.null(data) | is.null(date)){
      warning("Please upload files!")
    } 
    else{
      loc = match(data$tip.label,date[,1])
      d = date[loc,2]
      names(d) = data$tip.label
      res = roottotip(data,d,permTest = F)
      res_t = bactdate(data,d)
      
      d_1 <<- d
      t_1 <<- data
      res_t_1 <<- res_t
      
      
      output$detectfig <- renderPlot({
        roottotip(t_1,d_1)
      })
      
      output$detectfig1 <- renderPlot({
        if(input$inSelect == 'True') {
          if(input$inSelect1 == 'True') {
            plot(res_t_1, show.tip.label = T)
            axisPhylo(backward = T)
          }
          else {
            plot(res_t_1, show.tip.label = F)
            axisPhylo(backward = T)
          }
        }
        else {
          if(input$inSelect1 == 'True') {
            plot(res_t_1, show.tip.label = T)
            axisPhylo(backward = F)
          }
          else {
            plot(res_t_1, show.tip.label = F)
            axisPhylo(backward = F)
          }
        }
      })
      
      output$file3 <- renderUI({
        h5('Molecular temporal signal')
      })
      
      output$file4 <- renderUI({
        h5('Molecular dating tree')
      })
      
      rec = data.frame(Rate = res$rate,Intercept = res$ori, 
                       P_val = res$pvalue,Rootdate = res_t$rootdate)
      row.names(rec) = 'Item'
      final <<- rec
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
      plot(res_t_1, show.tip.label = F)
      axisPhylo(backward = F)
      dev.off()
    }
  )
  
  ## 下载图片写法
  output$plotDown1 <- downloadHandler(
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
      roottotip(t_1,d_1)
      dev.off()
    }
  )

  
  output$downloadData <- downloadHandler(
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
runApp(app, port = 51733)

