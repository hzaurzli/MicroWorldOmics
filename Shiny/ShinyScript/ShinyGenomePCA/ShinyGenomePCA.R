library(shiny)

ui <- fluidPage(
  titlePanel("Genome PCA"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Step 1: Choose alignment File",
                accept = c(
                  "aln",
                  "fasta/aln,text/plain",
                  ".aln")
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
      plotOutput("detectfig",width = "100%",height = '500px')
    )
  )
)



server <- function(input, output, session) {
  alignment <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    fasta2genlight(infile$datapath,
                   chunkSize = 10,
                   parallel = F)
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
    library(adegenet)

    
    flu = alignment()
    dfGroup = metadata()
    
    print(is.null(flu) | is.null(dfGroup))
    
    if(is.null(flu) | is.null(dfGroup)){
      warning("Please upload files!")
    } 
    else{
      if (input$Type == 'y') {
        df.pca<-glPca(flu,nf=3)  
        df.pca.scores<-as.data.frame(df.pca$scores)  
        df.pca.scores$id = row.names(df.pca.scores)
        dfGroup$id = row.names(df.pca.scores)
        pca = na.omit(merge(df.pca.scores,dfGroup,
                            by = 'id',all = T))
        pca = pca[!duplicated(pca$id),]
        
        if(ncol(pca) == 5){
          p1 <<- ggplot(data=pca,
                        aes(x = PC1, y = PC2,
                            color = pca[,4],
                            shape = pca[,5]))+
            geom_point(size=input$size1)+
            labs(
              x = paste0("PC1"),
              y = paste0("PC2")
            ) + 
            stat_ellipse(aes(fill = type_1),
                         type="norm",geom="polygon",alpha=0.2,color=NA) + 
            theme_classic() +
            theme(legend.title=element_text(size=input$num1),
                  axis.text=element_text(size=input$num2),
                  axis.title.x = element_text(size=input$num3),
                  axis.title.y = element_text(size=input$num3),
                  legend.text = element_text(size=input$num2)) +
            labs(color = colnames(pca)[4], 
                 shape = colnames(pca)[5])
            
          
          print(p1)
        }
          
        else if(ncol(dfGroup) == 4){
          p1 <<- ggplot(data=pca,
                        aes(x = PC1, y = PC2,
                            color = pca[,4]))+
            geom_point(size=input$size1)+
            labs(
              x = paste0("PC1"),
              y = paste0("PC2")
            ) + 
            theme_classic() +
            theme(legend.title=element_text(size=input$num1),
                  axis.text=element_text(size=input$num2),
                  axis.title.x = element_text(size=input$num3),
                  axis.title.y = element_text(size=input$num3),
                  legend.text = element_text(size=input$num2)) +
            labs(color = colnames(pca)[4], 
                 shape = colnames(pca)[5])
          
          print(p1)
        }
        
        else{
          warning("Too many categories, at most two categories are supported!")
        }
        
      } else if (input$Type == 'n') {
          df.pca<-glPca(flu,nf=3)  
          df.pca.scores<-as.data.frame(df.pca$scores)  
          df.pca.scores$id = row.names(df.pca.scores)
          pca = na.omit(merge(df.pca.scores,ssuit_info_new,
                              by = 'id',all = T))
          pca = pca[!duplicated(pca$id),]
          
          if(ncol(dfGroup) == 5){
            p2 <<- ggplot(data=pca,
                          aes(x = PC1, y = PC2,
                              color = pca[,4],
                              shape = pca[,5]))+
              geom_point(size=input$size1)+
              labs(
                x = paste0("PC1"),
                y = paste0("PC2")
              ) + 
              stat_ellipse(aes(fill = type_1),
                           type="norm",geom="polygon",alpha=0.2,color=NA) + 
              theme_classic() +
              theme(legend.title=element_text(size=input$num1),
                    axis.text=element_text(size=input$num2),
                    axis.title.x = element_text(size=input$num3),
                    axis.title.y = element_text(size=input$num3),
                    legend.text = element_text(size=input$num2)) +
              labs(color = colnames(pca)[4], 
                   shape = colnames(pca)[5])
            
            print(p2)
          }
          
          else if(ncol(dfGroup) == 4){
            p2 <<- ggplot(data=pca,
                          aes(x = PC1, y = PC2,
                              color = pca[,4]))+
              geom_point(size=input$size1)+
              labs(
                x = paste0("PC1"),
                y = paste0("PC2")
              ) + 
              theme_classic() +
              theme(legend.title=element_text(size=input$num1),
                    axis.text=element_text(size=input$num2),
                    axis.title.x = element_text(size=input$num3),
                    axis.title.y = element_text(size=input$num3),
                    legend.text = element_text(size=input$num2)) +
              labs(color = colnames(pca)[4], 
                   shape = colnames(pca)[5])
            
            print(p2)
          }
          else{
            warning("Too many categories, at most two categories are supported!")
          }
 
      }
    }
    
    ## 下载图片写法
    output$plotDown <- downloadHandler(
      filename = function(){
        paste0(input$dataset, '.',input$extPlot)
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
runApp(app, port = 50328)


