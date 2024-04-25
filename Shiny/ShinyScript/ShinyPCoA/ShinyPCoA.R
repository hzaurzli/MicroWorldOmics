library(shiny)

ui <- fluidPage(
  titlePanel("PCoA"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Step 1: Choose data File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      fileInput("file2", "Step 2: Choose metadata File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      hr(),
      radioButtons("Type", "Confidence ellipse",
                   c("Yes"="y", "No"="n")
      ),
      hr(),
      numericInput("num1", "Legend size",
                   value = 15, min = 0, max = 100),
      numericInput("num2", "Coordinate size",
                   value = 15, min = 0, max = 100),
      numericInput("num3", "Text one",
                   value = 20, min = 0, max = 100),
      sliderInput("size1", "Point text size",
                  min = 0, max = 25,
                  value = 8),
      sliderInput("size2", "Points size",
                  min = 0, max = 25,
                  value = 8),
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
      plotOutput("detectfig",width = "100%",height = '700px')
    )
  )
)



server <- function(input, output, session) {
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
               header = T,
               row.names = 1
    )
  })
  

  
  output$detectfig <- renderPlot({
    library(ggplot2)
    library(ade4)   # 用于计算PcoA
    library(vegan) 
    
    df = filedata()
    dfGroup = metadata()
    
    print(is.null(df) | is.null(dfGroup))
    
    if(is.null(df) | is.null(dfGroup)){
      warning("Please upload files!")
    } 
    else{
      if (input$Type == 'y') {
        df=t(df)
        df.dist = vegdist(df,method='euclidean')    #基于euclidean距离
        pcoa =  dudi.pco(df.dist,
                         scannf = F,   # 一种逻辑值，指示是否应该显示特征值条形图
                         nf=2)         # 保留几个维度的坐标信息
        
        # 整理绘图所需的数据
        data = pcoa$li
        data$name = rownames(data)
        data$group = dfGroup$Group
        
        # 绘图, p1 为全局变量
        p1 <<- ggplot(data,aes(x = A1,
                               y = A2,
                               color = group,
                               group = group,
                               fill = group
        ))+
          geom_point(size = input$size2)+
          theme_classic()+
          geom_vline(xintercept = 0, color = 'gray', size = 0.4) +   # 在0处添加垂直线条
          geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
          stat_ellipse(aes(x=A1,    # 添加置信区间圈
                           y=A2,
          ),
          geom = "polygon",
          level = 0.95,
          alpha=0.4)+
          geom_text(                # 添加文本标签
            aes(label=name),   
            vjust=1.5,            
            size=input$size1,
            color = "black"
          )+
          labs(  # 更改x与y轴坐标为pcoa$eig/sum(pcoa$eig)
            x = paste0("PCoA1 (",as.character(round(pcoa$eig[1] / sum(pcoa$eig) * 100,2)),"%)"),
            y = paste0("PCoA2 (",as.character(round(pcoa$eig[2] / sum(pcoa$eig) * 100,2)),"%)")
          ) + 
          theme(legend.title=element_text(size=input$num1),
                axis.text=element_text(size=input$num2),
                axis.title.x = element_text(size=input$num3),
                axis.title.y = element_text(size=input$num3),
                legend.text = element_text(size=input$num2))
        
        print(p1)
        
      } else if (input$Type == 'n') {
        df=t(df)
        df.dist = vegdist(df,method='euclidean')    #基于euclidean距离
        pcoa =  dudi.pco(df.dist,
                         scannf = F,   # 一种逻辑值，指示是否应该显示特征值条形图
                         nf=2)         # 保留几个维度的坐标信息
        
        # 整理绘图所需的数据
        data = pcoa$li
        data$name = rownames(data)
        data$group = dfGroup$Group
        
        # 绘图, p2 为全局变量
        p2 <<- ggplot(data,aes(x = A1,
                               y = A2,
                               color = group,
                               group = group,
                               fill = group
        ))+
          geom_point(size = input$size2)+
          theme_classic()+
          geom_vline(xintercept = 0, color = 'gray', size = 0.4) +   # 在0处添加垂直线条
          geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
          geom_text(                # 添加文本标签
            aes(label=name),   
            vjust=1.5,            
            size=input$size1,
            color = "black"
          )+
          labs(  # 更改x与y轴坐标为pcoa$eig/sum(pcoa$eig)
            x = paste0("PCoA1 (",as.character(round(pcoa$eig[1] / sum(pcoa$eig) * 100,2)),"%)"),
            y = paste0("PCoA2 (",as.character(round(pcoa$eig[2] / sum(pcoa$eig) * 100,2)),"%)")
          ) + 
          theme(legend.title=element_text(size=input$num1),
                axis.text=element_text(size=input$num2),
                axis.title.x = element_text(size=input$num3),
                axis.title.y = element_text(size=input$num3),
                legend.text = element_text(size=input$num2))
        
        print(p2)
      }
    }
  
    ## 下载图片写法
    output$plotDown <- downloadHandler(
      filename = function(){
        paste0(input$file1, '.',input$extPlot)
      },
      content = function(file){
        if(input$Type=='y'){
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
        } else if(input$Type=='n'){
          if(input$extPlot == 'pdf'){
            pdf(file)
          }else if(input$extPlot == 'png'){
            png(file)
          }else{
            jpeg(file)
          }
          # 打印全局变量 p2 (ggplot对象)
          print(p2)
          dev.off()
        }
      }
    )
  })
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50327)

