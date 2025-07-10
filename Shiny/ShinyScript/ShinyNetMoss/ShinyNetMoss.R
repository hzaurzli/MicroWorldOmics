library(shinycssloaders)
library(shinythemes)
library(ggplot2)
library(shiny)
library(DT)
library(NetMoss2)
library(rsparcc)
library(shinyFiles)


button_color_css <- "
#DivCompClear, #FinderClear, #EnterTimes{
/* Change the background color of the update button
to blue. */
background: DodgerBlue;

/* Change the text size to 15 pixels. */
font-size: 15px;
}"


ui <- fluidPage(
  
  #Navbar structure for UI
  navbarPage("NetMoss2 for microbiomarker finding", 
             theme = shinytheme("lumen"),
             # tabPanel 是导航栏换页
             tabPanel("For multiple files", fluid = TRUE,
                      tags$style(button_color_css),
                      # Sidebar layout with a input and output definitions
                      sidebarLayout(
                        sidebarPanel(
                          h3(strong('NetMoss2 multiple files')),
                          hr(),
                          fluidRow(
                            column(12,
                                   uiOutput('meta'),
                                   uiOutput('text1')
                            ),
                            column(6,
                                   h5('Step2: Case input directory'),
                                   shinyDirButton('folder1', 
                                                  'Case directory', 
                                                  FALSE),
                                   br(),
                                   br(),
                                   uiOutput('text2'),
                                   br(),
                                   sliderInput("NetMoss_Score", "NetMoss score",
                                               value = 0.3, min = 0, max = 1),
                                   actionButton('reset', 'RESET'),
                            ),
                            column(6,
                                   h5('Step3: Ctrl input directory'),
                                   shinyDirButton('folder2', 
                                                  'Control directory', 
                                                  FALSE),
                                   br(),
                                   br(),
                                   uiOutput('text3'),
                                   br(),
                                   sliderInput("p_adj", "Markers padj cutoff",
                                               value = 0.05, min = 0, max = 1),
                                   actionButton('start', 'START'),
                            ),
                            column(9,
                                   hr(),
                                   h4('Download all results'),
                                   downloadButton("downloadTable", "Download markers table"),
                                   br(),
                                   br(),
                                   downloadButton("downloadGraph", "Download network Graph"),
                                   br(),
                                   br(),
                                   radioButtons('extPlot', 'Plot output format',
                                                choices = c("PNG"='png', 
                                                            'PDF'='pdf',
                                                            'JPEG'='jpeg'), inline = T),
                                   hr(),
                                   h4('Markers verification (after finding markers)'),
                                   sliderInput("train_num", "Verification training number",
                                               value = 20, min = 5, max = 200),
                                   actionButton('start_v', 'START'),
                                   br(),
                                   br(),
                                   downloadButton("downloadGraph_v", "Download roc Graph"),
                                   br(),
                                   br(),
                                   radioButtons('extPlot_v', 'Plot output format',
                                                choices = c("PNG"='png', 
                                                            'PDF'='pdf',
                                                            'JPEG'='jpeg'), inline = T),
                                   
                            )
                          ),
                        ),
                        mainPanel(
                          titlePanel("Markers table and network graph"),
                          withSpinner(dataTableOutput(outputId = "table")),
                          plotOutput("detectfig",width = "100%"),
                          uiOutput('title'),
                          withSpinner(plotOutput("detectfig_v",width = "100%"))
                        )
                      )
             ),
             
             tabPanel("For single file", fluid = TRUE,
                      tags$style(button_color_css),
                      # Sidebar layout with a input and output definitions
                      sidebarLayout(
                        sidebarPanel(
                          h3(strong('NetMoss2 single file')),
                          hr(),
                          fluidRow(
                            uiOutput('meta1'),
                            uiOutput('file1'),
                            uiOutput('file2'),
                            
                            sliderInput("NetMoss_Score1", "NetMoss score",
                                        value = 0.3, min = 0, max = 1),
                            sliderInput("p_adj1", "Markers padj cutoff",
                                        value = 0.05, min = 0, max = 1),
                            actionButton('reset1', 'RESET'),
                            actionButton('start1', 'START'),
                            hr(),
                            h4('Download all results'),
                            downloadButton("downloadTable1", "Download markers table"),
                            br(),
                            br(),
                            downloadButton("downloadGraph1", "Download network Graph"),
                            br(),
                            br(),
                            radioButtons('extPlot1', 'Plot output format',
                                         choices = c("PNG"='png', 
                                                     'PDF'='pdf',
                                                     'JPEG'='jpeg'), inline = T),
                            
                            hr(),
                            h4('Markers verification (after finding markers)'),
                            sliderInput("train_num1", "Verification training number",
                                        value = 20, min = 5, max = 200),
                            actionButton('start_v1', 'START'),
                            br(),
                            br(),
                            downloadButton("downloadGraph_v1", "Download roc Graph"),
                            br(),
                            br(),
                            radioButtons('extPlot_v1', 'Plot output format',
                                         choices = c("PNG"='png', 
                                                     'PDF'='pdf',
                                                     'JPEG'='jpeg'), inline = T),
                            
                          ),
                        ),
                        mainPanel(
                          titlePanel("Markers table and network graph"),
                          withSpinner(dataTableOutput(outputId = "table1")),
                          plotOutput("detectfig1",width = "100%"),
                          uiOutput('title1'),
                          withSpinner(plotOutput("detectfig_v1",width = "100%"))
                        )
                      )
             ),
             
  ),
)


server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  path_get <- path.expand('~')
  print(path_get)
  
  # 选择文件夹
  shinyDirChoose(input, 'folder1', roots=c(path = path_get))
  shinyDirChoose(input, 'folder2', roots=c(path = path_get))
  
  
  metadata_r <- reactive({
    infile <- input$meta
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  
  observeEvent(input$reset, {
    values$meta <- NULL
    output$meta <- renderUI({
      fileInput("meta", "Step 1: Choose metadata",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$text1 <- NULL
    output$text1 <- renderUI({
      textInput("text1", "Data root path (put your data in this path)",value = path_get)
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$text2 <- NULL
    output$text2 <- renderUI({
      textInput("text2", "Directory path",value = "")
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$text3 <- NULL
    output$text3 <- renderUI({
      textInput("text3", "Directory path",value = "")
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$folder1, {
    if(length(input$folder1) == 1){
      print('ok')
    }
    else{
      path_vec = c()
      path_list = input$folder1$path
      for (i in 1:length(path_list)) {
        path_vec = c(path_vec,path_list[[i]][1])
      }
      path_path = path_vec[1]
      for (i in 2:length(path_vec)) {
        path_path = paste(path_path, path_vec[i],sep = '/')
      }
      path_get <- path.expand('~')
      path_all <- paste(path_get,path_path,sep = '')
      
      output$text2 <- renderUI({
        textInput("text2", "Directory path",value = path_all)
      })
    }
  }, ignoreNULL = F)
  
  
  observeEvent(input$folder2, {
    if(length(input$folder2) == 1){
      print('ok')
    }
    else{
      path_vec = c()
      path_list = input$folder2$path
      for (i in 1:length(path_list)) {
        path_vec = c(path_vec,path_list[[i]][1])
      }
      path_path = path_vec[1]
      for (i in 2:length(path_vec)) {
        path_path = paste(path_path, path_vec[i],sep = '/')
      }
      path_get <- path.expand('~')
      path_all <- paste(path_get,path_path,sep = '')
      
      output$text3 <- renderUI({
        textInput("text3", "Directory path",value = path_all)
      })
    }
  }, ignoreNULL = F)
  
  
  output$table <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  output$detectfig <- renderPlot({
    NULL
  })
  
  output$detectfig_v <- renderPlot({
    NULL
  })
  
  
  observeEvent(input$reset, {
    values$detectfig <- NULL
    output$detectfig <- renderPlot({
      NULL
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$detectfig_v <- NULL
    output$detectfig_v <- renderPlot({
      NULL
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$start, {
    
    metadata = metadata_r()
    
    if(is.null(metadata)){
      warning("Please upload files!")
    }
    else{
      if(!(nzchar(input$text1) | nzchar(input$text2))){
        warning("Please upload files!")
      }
      else{
        case_dir = input$text2
        control_dir = input$text3
        
        case_dir_final <<- case_dir
        control_dir_final <<- control_dir
        
        gg = strsplit(case_dir,'/')[[1]]
        gg = gg[-length(gg)]
        
        net_case_dir = paste(paste(gg, collapse = "/"),'net_case_dir',sep = '/')
        net_control_dir = paste(paste(gg, collapse = "/"),'net_control_dir',sep = '/')
        
        output$table <- renderDataTable({
          
          if(file.exists(net_case_dir)){
            unlink(net_case_dir)
            print('delete')
          }
          
          if(file.exists(net_control_dir)){
            unlink(net_control_dir)
          }
          
          set.seed(123)
          #construct networks  ####if files exist, skip
          library(rsparcc)
          netBuild(case_dir = case_dir,
                   control_dir = control_dir,
                   method = "sparcc")
          
          #calculate NetMoss score
          nodes_result = NetMoss(case_dir = case_dir,    
                                 control_dir = control_dir,    
                                 net_case_dir = net_case_dir,   
                                 net_control_dir = net_control_dir) 
          
          nodes_result_plot <<- nodes_result
          result = nodes_result[[1]]
          
          netmoss_score = input$NetMoss_Score
          padj = input$p_adj
          marker = data.frame(result[which(result$p.adj < padj),])
          marker = data.frame(marker[which(marker$NetMoss_Score > netmoss_score),])
          rownames(marker) = marker$taxon_names
          
          marker_df = data.frame(marker)
          row.names(marker_df) = NULL
          marker_final <<- marker
          
          metadata[,2] = as.factor(metadata[,2])
          metadata_final <<- metadata
          
          print(dim(marker))
          print(dim(metadata_final))
          
          final <<- data.frame(marker_df)
        }, options = list(pageLength = 6))
        
        output$detectfig <- renderPlot({
          #plot networks
          pp = getwd()
          netPlot(nodes_result_plot)
          unlink(paste(pp,'NetMoss_score.pdf',sep = '/'))
        })
      }
    }
    
  })
  
  
  observeEvent(input$start_v, {
    output$detectfig_v <- renderPlot({
      
      trainnum = input$train_num
      
      myROC = netROC(case_dir = case_dir_final,
                     control_dir = control_dir_final,
                     marker = marker_final,
                     metadata = metadata_final,
                     plot.roc = FALSE, 
                     train.num = trainnum) 
      
      
      tp.fp = myROC[[1]]
      
      output$title <- renderUI({
        h5('ROC curve')
      })
      
      p2 = ggplot() +
        geom_path(aes(FPR, TPR), data = tp.fp) +
        coord_equal() +
        annotate(
          "text",
          x = .75,
          y = .15,
          label = paste("AUC =", round(WeightedAUC(tp.fp), 2))
        ) +
        labs(x = "False positive fraction", y = "True positive fraction") +
        theme_bw()
      
      pr <<- p2
      print(p2)
    })
  })
  
  output$downloadTable <- downloadHandler(
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
      
      pp = getwd()
      netPlot(nodes_result_plot)
      unlink(paste(pp,'NetMoss_score.pdf',sep = '/'))
      dev.off()
    }
  )
  
  
  output$downloadGraph_v <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot_v)
    },
    content = function(file){
      if(input$extPlot_v == 'pdf'){
        pdf(file)
      }else if(input$extPlot_v == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      
      print(pr)
      dev.off()
    }
  )
  
  #####################################
  #####################################
  values <- reactiveValues(
    file = NULL
  )
  
  
  
  metadata_r1 <- reactive({
    infile <- input$meta1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  
  observeEvent(input$reset1, {
    values$meta1 <- NULL
    output$meta1 <- renderUI({
      fileInput("meta1", "Step 1: Choose metadata",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  casedata <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T
    )
  })
  
  ctrldata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T
    )
  })
  
  
  observeEvent(input$reset1, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 2: Choose case matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset1, {
    values$file <- NULL
    output$file2 <- renderUI({
      fileInput("file2", "Step 3: Choose control matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table1 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  output$detectfig1 <- renderPlot({
    NULL
  })
  
  output$detectfig_v1 <- renderPlot({
    NULL
  })
  
  
  observeEvent(input$reset1, {
    values$detectfig1 <- NULL
    output$detectfig1 <- renderPlot({
      NULL
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset1, {
    values$detectfig_v1 <- NULL
    output$detectfig_v1 <- renderPlot({
      NULL
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset1, {
    values$detectfig1 <- NULL
    output$detectfig1 <- renderPlot({
      NULL
    })
  }, ignoreNULL = F)
  
  
  
  observeEvent(input$start1, {
    case_dir = casedata()
    control_dir = ctrldata()
    
    case_dir_final <<- case_dir
    control_dir_final <<- control_dir
    
    metadata = metadata_r1()
    
    if(is.null(case_dir) | is.null(control_dir) | is.null(metadata)){
      warning("Please upload files!")
    }
    else{
      
      output$table1 <- renderDataTable({
        infile <- input$file1
        
        hh = strsplit(infile$datapath,'/')[[1]]
        hh = hh[-length(hh)]
        
        #construct networks  ####if files exist, skip
        library(rsparcc)
        ww = netBuild(case_dir = case_dir,
                      control_dir = control_dir,
                      method = "sparcc")
        
        path_s = strsplit(ww[1],' ')[[1]][5]
        
        file.copy(paste(path_s,'d_net.txt',sep = '/'),
                  paste(paste(hh, collapse = "/"),'d_net.txt',sep = '/'))
        
        file.copy(paste(path_s,'h_net.txt',sep = '/'),
                  paste(paste(hh, collapse = "/"),'h_net.txt',sep = '/'))
        
        
        Sys.sleep(5)
        unlink(paste(path_s,'d_net.txt',sep = '/'))
        unlink(paste(path_s,'h_net.txt',sep = '/'))
        net_case_dir = read.table(paste(paste(hh, collapse = "/"),'d_net.txt',sep = '/'),
                                  header = T,sep = '\t', row.names = 1)
        
        net_control_dir = read.table(paste(paste(hh, collapse = "/"),'h_net.txt',sep = '/'),
                                     header = T,sep = '\t', row.names = 1)
        #calculate NetMoss score
        nodes_result = NetMoss(case_dir = case_dir,    
                               control_dir = control_dir,    
                               net_case_dir = net_case_dir,   
                               net_control_dir = net_control_dir)
        
        nodes_result_plot1 <<- nodes_result
        result = nodes_result[[1]]
        
        
        netmoss_score = input$NetMoss_Score1
        padj = input$p_adj1
        marker = data.frame(result[which(result$p.adj < padj),])
        marker = data.frame(marker[which(marker$NetMoss_Score > netmoss_score),])
        rownames(marker) = marker$taxon_names
        
        
        marker_df = data.frame(marker)
        row.names(marker_df) = NULL
        marker_final <<- marker
        
        metadata[,2] = as.factor(metadata[,2])
        metadata_final <<- metadata
        
        print(dim(marker))
        print(dim(metadata_final))
        
        final1 <<- data.frame(marker_df)
      }, options = list(pageLength = 9))
      
      output$detectfig1 <- renderPlot({
        #plot networks
        pp = getwd()
        netPlot(nodes_result_plot1)
        unlink(paste(pp,'NetMoss_score.pdf',sep = '/'))
      })
      
    }
  })
  
  
  observeEvent(input$start_v1, {
    output$detectfig_v1 <- renderPlot({
      
      trainnum = input$train_num1
      
      myROC = netROC(case_dir = case_dir_final,
                     control_dir = control_dir_final,
                     marker = marker_final,
                     metadata = metadata_final,
                     plot.roc = FALSE, 
                     train.num = trainnum) 
      
      
      tp.fp = myROC[[1]]
      
      output$title1 <- renderUI({
        h5('ROC curve')
      })
      
      p2 = ggplot() +
        geom_path(aes(FPR, TPR), data = tp.fp) +
        coord_equal() +
        annotate(
          "text",
          x = .75,
          y = .15,
          label = paste("AUC =", round(WeightedAUC(tp.fp), 2))
        ) +
        labs(x = "False positive fraction", y = "True positive fraction") +
        theme_bw()
      
      pr <<- p2
      print(p2)
    })
  })
  
  output$downloadTable1 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final1,file,row.names = T,quote = F)
    }
  )
  
  
  output$downloadGraph1 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot1)
    },
    content = function(file){
      if(input$extPlot1 == 'pdf'){
        pdf(file)
      }else if(input$extPlot1 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      
      pp = getwd()
      netPlot(nodes_result_plot1)
      unlink(paste(pp,'NetMoss_score.pdf',sep = '/'))
      dev.off()
    }
  )
  
  
  output$downloadGraph_v1 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot_v1)
    },
    content = function(file){
      if(input$extPlot_v1 == 'pdf'){
        pdf(file)
      }else if(input$extPlot_v1 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      
      print(pr)
      dev.off()
    }
  )
  
}  

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 51235)
