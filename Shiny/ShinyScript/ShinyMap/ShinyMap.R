library(shiny)

ui <- fluidPage(
  titlePanel("Shiny Map"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Step 1: Choose world (RData) File",
                accept = c(
                  "RData",
                  "RData/rdata,text/plain",
                  ".RData")
      ),
      fileInput("file2", "Step 2: Choose metadata File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      hr(),
      sliderInput("num1", "Legend size",
                   value = 15, min = 0, max = 25),
      sliderInput("num2", "Legend text size",
                  value = 13, min = 0, max = 25),
      ## 适用于少量文本
      textInput("low", "color low", value = "red"),
      ## 适用于少量文本
      textInput("mid", "color mid", value = "white"),
      ## 适用于少量文本
      textInput("high", "color high", value = "blue"),
      
      numericInput("mid_point", "Mid point value",
                   value = 15, min = 0, max = 100),
      
      radioButtons('extPlot', 'Plot output format',
                   choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
      downloadButton("plotDown", "Download"),
      hr(),
      h5('Developer:'),
      h6('Small runze (shiny app)'),
      br(),
      h5('Github: '),
      h6('https://github.com/hzaurzli (Small runze)')
    ),
    mainPanel(
      plotOutput("detectfig",width = "100%",height = '300px')
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  # 上传文件后自动读取文件
  mapread <- reactive({
    sessionEnvir <- sys.frame()
    if (!is.null(input$file1)) eval(parse(text = load(input$file1$datapath, sessionEnvir)))
  })
  
  
  metadata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
			 fileEncoding = 'utf-8')
  })
  
  
  
  output$detectfig <- renderPlot({
    library(ggplot2)
    theme_set(theme_bw())
    library(sf)

    
    map = mapread()
    meta = metadata()
    
    print(is.null(map) | is.null(meta))
    
    if(is.null(map) | is.null(meta)){
      warning("Please upload files!")
    } 
    else{
      name = colnames(meta)[2] 
      colnames(meta) = c('Country','Data')
      colnames(map)[1] = c('Country')
      
      data = merge(map,meta,by='Country',all = T)
      
      options(scipen=200)
      p = ggplot(data = data) +
            geom_sf(aes(fill = Data)) + 
            scale_fill_gradient2(low = input$low,
                                 mid = input$mid,
                                 high = input$high,
                                 midpoint = input$mid_point) + 
            theme(legend.title=element_text(size=input$num1),
                  legend.text=element_text(size=input$num2)) +
            labs(fill = name)
      
      print(p)
    }
    
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
  })
}

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50329)
