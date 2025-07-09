library(shiny)
library(igraph)
library(WGCNA)
library(RColorBrewer)
library(flashClust)



ui <- tagList(
  fluidPage(
    titlePanel("Differential network analysis based on DiffCoEx"),
    sidebarLayout(
      sidebarPanel(
        # uiOutput 做上传文件的 ui, 对应后面的 output$file1
        uiOutput('file1'),
        uiOutput('file2'),
        numericInput('num','Step 3: User defined parameter for soft thresholding',
                     value = 6,min = 1,max = 10000),
        actionButton('reset', 'RESET'),
        actionButton('start', 'START'),
        hr(),
        h4(strong('Download all results (large data be patient)')),
        downloadButton("downloadData", "Download dissimilarity table"),
        br(),
        br(),
        downloadButton("downloadData1", "Download correlation matrix"),
        br(),
        br(),
        downloadButton("downloadData2", "Download ASV/OTU modules"),
        br(),
        br(),
        downloadButton("plotDown", "Download figure"),
        br(),
        br(),
        radioButtons('extPlot', 'Plot output format',
                     choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
        hr(),
        h5('Developer:'),
        h6('Small runze (shiny app)'),
        br(),
        h5('Github: '),
        h6('https://github.com/hzaurzli (Small runze)'),
      ),
      mainPanel(
        h4("Dissimilarity coefficient T values and ASV/OTU dendrogram and module"),
        shinycssloaders::withSpinner(
          dataTableOutput("table")
        ),
        br(),
        h5(strong("Larger T values represent that the weights (regulatory relationships) of ASV/OTU pairs are less different between the two treatments!!!")),
        hr(),
        h4("The correlation matrix divides modules according to the degree of difference"),
        shinycssloaders::withSpinner(
          dataTableOutput("table_module"),
        ),
        hr(),
        plotOutput(outputId = "detectfig", width = "90%")
      )
    )
  )
)


server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  otuC1_dat <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T,row.names = 1)
  })
  
  
  otuC2_dat <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T,row.names = 1)
  })
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 1: Choose otu abundance matrices (Treatment 1)",
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
      fileInput("file2", "Step 2: Choose otu abundance matrices (Treatment 2)",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  output$table <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  output$table_module <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  observeEvent(input$start, {
    
    otuC1_dat = as.matrix(otuC1_dat())
    otuC2_dat = as.matrix(otuC2_dat())
    
    if(is.null(otuC1_dat) | is.null(otuC2_dat) ){
      warning("Please upload files!")
    }
    else{
      # 定义 β 值
      beta1=input$num
      
      name1 = colnames(otuC1_dat)
      name2 = colnames(otuC2_dat)
      
      name_i1 = name1[order(name1)]
      name_i2 = name2[order(name2)]
      
      if (!identical(name_i1, name_i2)) {
        message('Error: different column names for documents 1 and 2')
      }
      else{
        datC1 = otuC1_dat
        datC2 = otuC2_dat
        
        output$table <- renderDataTable({
          AdjMatC1<-sign(cor(datC1,method="spearman"))*(cor(datC1,method="spearman"))^2
          AdjMatC2<-sign(cor(datC2,method="spearman"))*(cor(datC2,method="spearman"))^2
          
          corC1 = cor(datC1,method="spearman")
          corC2 = cor(datC2,method="spearman")
          
          corC1_dat = corC1
          corC2_dat = corC2
          
          corC1_dat[upper.tri(corC1_dat,diag = T)] <- 0
          corC2_dat[upper.tri(corC2_dat,diag = F)] <- 0
          
          cor_dat = corC1_dat + corC2_dat
          cor_final <<- cor_dat
          
          diag(AdjMatC1)<-0
          diag(AdjMatC2)<-0
          collectGarbage()
          # 计算相异矩阵 T
          dissTOMC1C2=TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta1/2))
          collectGarbage()
          
          geneTreeC1C2 <<- flashClust(as.dist(dissTOMC1C2), method = "average");
          
          dynamicModsHybridC1C2 = cutreeDynamic(dendro = geneTreeC1C2, distM = dissTOMC1C2,method="hybrid",cutHeight=.996,deepSplit = T, pamRespectsDendro = FALSE,minClusterSize = 20);
          
          #Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
          dynamicColorsHybridC1C2 = labels2colors(dynamicModsHybridC1C2)
          
          #the next step merges clusters which are close (see WGCNA package documentation)
          mergedColorC1C2<-mergeCloseModules(rbind(datC1,datC2),dynamicColorsHybridC1C2,cutHeight=.2)$color
          colorh1C1C2 <<- mergedColorC1C2
          
          
          item = names(table(colorh1C1C2))
          item_order = c()
          ann = c()
          num = 1
          for (i in item) {
            item_order = c(item_order,which(colorh1C1C2 == i))
            ann = c(ann,rep(paste('M',num,sep = ''),
                            length(which(colorh1C1C2 == i))))
            num = num + 1
          }
          
          datC1 = datC1[,item_order]
          datC2 = datC2[,item_order]
          
          ann_row = data.frame(ann)
          row.names(ann_row) = colnames(datC1)
          ann_row_final <<- ann_row
          
          corC1 = cor(datC1,method="spearman")
          corC2 = cor(datC2,method="spearman")
          
          corC1_dat = corC1
          corC2_dat = corC2
          
          corC1_dat[upper.tri(corC1_dat,diag = T)] <- 0
          corC2_dat[lower.tri(corC2_dat,diag = F)] <- 0
          
          cor_dat = corC1_dat + corC2_dat
          cor_final <<- cor_dat
          
          
          myAdjacencyMatrix = as.matrix(dissTOMC1C2)
          colnames(myAdjacencyMatrix) = name1
          g  <- graph.adjacency(myAdjacencyMatrix,weighted=TRUE)
          df <- get.data.frame(g)
          colnames(df) = c('from','to','T value')
          df_dat <<- df
        },options = list(pageLength = 10))
      }
      
      
      output$table_module <- renderDataTable({
          ann_row_final
      },options = list(pageLength = 10))
    
      
      output$detectfig <- renderPlot({
        pheatmap::pheatmap(cor_final,cluster_cols = F,
                           cluster_rows = F,show_rownames = F,
                           show_colnames = F,border_color = NA,
                           annotation_row = ann_row_final,
                           annotation_legend = T)
        
      })
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
  
  
  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(cor_final,file,row.names = T,quote = F)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(ann_row_final,file,row.names = T,quote = F)
    }
  )
  
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
      pheatmap::pheatmap(cor_final,cluster_cols = F,
                         cluster_rows = F,show_rownames = F,
                         show_colnames = F,border_color = NA,
                         annotation_row = ann_row_final,
                         annotation_legend = T)
      dev.off()
    }
  )
  
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50758)
