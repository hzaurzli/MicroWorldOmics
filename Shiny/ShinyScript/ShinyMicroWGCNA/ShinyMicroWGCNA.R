library(dplyr)
library(ggplot2)
library(shiny)
library(DT)
library(ggrepel)
library(tidyr)
library(shinycssloaders)
library(shinythemes)
library(WGCNA)
options(stringsAsFactors = FALSE)


button_color_css <- "
#DivCompClear, #FinderClear, #EnterTimes{
/* Change the background color of the update button
to blue. */
background: DodgerBlue;

/* Change the text size to 15 pixels. */
font-size: 15px;
}"


# Define UI
ui <- fluidPage(
  
  #Navbar structure for UI
  navbarPage("Microbiome ASV/OTU WGCNA", 
     theme = shinytheme("lumen"),
     ##################
     navbarMenu("Step by step",
        tabPanel("Step1: soft threshold filtering", fluid = TRUE,
             tags$style(button_color_css),
             # Sidebar layout with a input and output definitions
             sidebarLayout(
               sidebarPanel(
                 h3(strong('Soft threshold filtering')),
                 hr(),
                 fluidRow(uiOutput('file1'),
                          uiOutput('file2'),
                          sliderInput("range", "Range:",
                                      min = 0, max = 50,
                                      value = c(1,20)),
                          actionButton('reset1', 'RESET'),
                          actionButton('start1', 'START'),
                          hr(),
                          h4('Download all results'),
                          downloadButton("downloadGraph1_1", "Download scale independence graph"),
                          br(),
                          br(),
                          downloadButton("downloadGraph1_2", "Download mean connectivity graph"),
                          br(),
                          br(),
                          radioButtons('extPlot1', 'Plot output format',
                                       choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                 ),
               ),
               
               mainPanel(
                 titlePanel("Scale independence and Mean connectivity"),
                 withSpinner(plotOutput(outputId = "detectfig1")),
                 hr()
               )
            )
        ),
        
        tabPanel("Step2: One-step network construction", fluid = TRUE,
           tags$style(button_color_css),
           # Sidebar layout with a input and output definitions
           sidebarLayout(
             sidebarPanel(
               h3(strong('One-step network construction')),
               hr(),
               fluidRow(uiOutput('file3'),
                        uiOutput('file4'),
                        h4('Parameters for one-step network'),
                        column(6,
                          selectInput('cortype', 'CorType', choices = c('pearson','bicor'),selected = 'pearson'),
                          selectInput('tomtype', 'TOMType', choices = c('unsigned','signed','signed Nowick','unsigned 2','signed 2','signed Nowick 2','none'),selected = 'unsigned'),
                          numericInput('maxblocksize','MaxBlockSize',
                                       value = 6000,min = 0,max = 50000),
                          numericInput('minmodulesize','MinModuleSize',
                                       value = 30,min = 0,max = 50000),
                          actionButton('reset2', 'RESET'),
                          actionButton('start2', 'START'),
                         
                        ),
                        column(6,
                          numericInput('power','Power',
                                      value = 12,min = 0,max = 50000), 
                          selectInput('networktype', 'NetworkType', choices = c('unsigned','signed','signed hybrid'), selected = 'unsigned'),
                          numericInput('reassignthreshold','Reassign Threshold',
                                       value = 0,min = 0,max = 50000),
                          numericInput('mergecutheight','Merge cut height',
                                       value = 0.25,min = 0,max = 10), 
                        ),
                        column(9,
                          hr(),
                          h4('Download all results'),
                          downloadButton("downloadGraph2", "Download graph"),
                          br(),
                          br(),
                          radioButtons('extPlot2', 'Plot output format',
                                       choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                        )
               ),
             ),
             
             mainPanel(
               titlePanel("One-step network construction"),
               withSpinner(plotOutput(outputId = "detectfig2")),
               hr()
             )
           )
        ),
        
        tabPanel("Step3: Modules and traits relationship", fluid = TRUE,
           tags$style(button_color_css),
           # Sidebar layout with a input and output definitions
           sidebarLayout(
             sidebarPanel(
               h3(strong('Modules and traits relationship')),
               hr(),
               fluidRow(uiOutput('file5'),
                        uiOutput('file6'),
                        
                        h4('Parameters for one-step network'),
                        column(6,
                               selectInput('cortype1', 'CorType', choices = c('pearson','bicor'),selected = 'pearson'),
                               selectInput('tomtype1', 'TOMType', choices = c('unsigned','signed','signed Nowick','unsigned 2','signed 2','signed Nowick 2','none'),selected = 'unsigned'),
                               numericInput('maxblocksize1','MaxBlockSize',
                                            value = 6000,min = 0,max = 50000),
                               numericInput('minmodulesize1','MinModuleSize',
                                            value = 30,min = 0,max = 50000),
                               actionButton('reset3', 'RESET'),
                               actionButton('start3', 'START'),
                               
                        ),
                        column(6,
                               numericInput('power1','Power',
                                            value = 12,min = 0,max = 50000), 
                               selectInput('networktype1', 'NetworkType', choices = c('unsigned','signed','signed hybrid'), selected = 'unsigned'),
                               numericInput('reassignthreshold1','Reassign Threshold',
                                            value = 0,min = 0,max = 50000),
                               numericInput('mergecutheight1','Merge cut height',
                                            value = 0.25,min = 0,max = 10), 
                        ),
                        column(9,
                              hr(),
                              h4('Download all results'),
                              downloadButton("downloadTable", "Download relationship table"),
                              br(),
                              br(),
                              downloadButton("downloadGraph3", "Download relationship graph"),
                              br(),
                              br(),
                              radioButtons('extPlot3', 'Plot output format',
                                           choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                        )
                    ),
             ),
             
             mainPanel(
               titlePanel("Modules and traits relationship matirx and graph"),
               withSpinner(dataTableOutput(outputId = "table")),
               hr(),
               plotOutput(outputId = "detectfig3")
             )
           )
        ),
        
        tabPanel("Step4: Export modules", fluid = TRUE,
           tags$style(button_color_css),
           # Sidebar layout with a input and output definitions
           sidebarLayout(
             sidebarPanel(
               h3(strong('Export modules')),
               hr(),
               fluidRow(
                 uiOutput('file7'),
                 uiOutput('file8'),
                 
                 h4('Parameters for one-step network'),
                 column(6,
                        selectInput('cortype2', 'CorType', choices = c('pearson','bicor'),selected = 'pearson'),
                        selectInput('tomtype2', 'TOMType', choices = c('unsigned','signed','signed Nowick','unsigned 2','signed 2','signed Nowick 2','none'),selected = 'unsigned'),
                        numericInput('maxblocksize2','MaxBlockSize',
                                     value = 6000,min = 0,max = 50000),
                        numericInput('minmodulesize2','MinModuleSize',
                                     value = 30,min = 0,max = 50000),
                        actionButton('reset4', 'RESET'),
                        actionButton('start4', 'START'),
                        
                 ),
                 column(6,
                        numericInput('power2','Power',
                                     value = 12,min = 0,max = 50000), 
                        selectInput('networktype2', 'NetworkType', choices = c('unsigned','signed','signed hybrid'), selected = 'unsigned'),
                        numericInput('reassignthreshold2','Reassign Threshold',
                                     value = 0,min = 0,max = 50000),
                        numericInput('mergecutheight2','Merge cut height',
                                     value = 0.25,min = 0,max = 10), 
                 ),
                 column(9,
                        hr(),
                        h4('Download all results'),
                        uiOutput('file9'),
                        uiOutput('file10'),
                        downloadButton("downloadTable1", "Download module edges table"),
                        br(),
                        br(),
                        downloadButton("downloadTable2", "Download module nodes table")
                 )
               ),
             ),
             
             mainPanel(
               titlePanel("Modules weight matrix"),
               dataTableOutput(outputId = "table1"),
               hr(),
               withSpinner(plotOutput(outputId = "detectfig4"))
             )
           )
        ),
        
        tabPanel("Step5: Network heatmap plot (All ASV/OTU)", fluid = TRUE,
           tags$style(button_color_css),
           # Sidebar layout with a input and output definitions
           sidebarLayout(
             sidebarPanel(
               h3(strong('Network heatmap plot (All ASV/OTU)')),
               hr(),
               fluidRow(
                  uiOutput('file11'),
                  uiOutput('file12'),
                  
                  h4('Parameters for one-step network'),
                  column(6,
                         selectInput('cortype3', 'CorType', choices = c('pearson','bicor'),selected = 'pearson'),
                         selectInput('tomtype3', 'TOMType', choices = c('unsigned','signed','signed Nowick','unsigned 2','signed 2','signed Nowick 2','none'),selected = 'unsigned'),
                         numericInput('maxblocksize3','MaxBlockSize',
                                      value = 6000,min = 0,max = 50000),
                         numericInput('minmodulesize3','MinModuleSize',
                                      value = 30,min = 0,max = 50000),
                         actionButton('reset5', 'RESET'),
                         actionButton('start5', 'START'),
                         
                  ),
                  column(6,
                         numericInput('power3','Power',
                                      value = 12,min = 0,max = 50000), 
                         selectInput('networktype3', 'NetworkType', choices = c('unsigned','signed','signed hybrid'), selected = 'unsigned'),
                         numericInput('reassignthreshold3','Reassign Threshold',
                                      value = 0,min = 0,max = 50000),
                         numericInput('mergecutheight3','Merge cut height',
                                      value = 0.25,min = 0,max = 10), 
                  ),
                  column(9,
                         hr(),
                         h4('Download all results'),
                         downloadButton("downloadGraph5", "Download graph"),
                         br(),
                         br(),
                         radioButtons('extPlot5', 'Plot output format',
                                      choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                  )
               ),
             ),
             
             mainPanel(
               titlePanel("Network heatmap graph"),
               withSpinner(plotOutput(outputId = "detectfig5"))
             )
           )
        ),
  
     ),
  ),
)


server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  otudata <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  
  sampledata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  observeEvent(input$reset1, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "First: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  output$detectfig1 <- renderPlot({
    return(NULL)
  })
  
  observeEvent(input$reset1, {
    values$file <- NULL
    output$file2 <- renderUI({
      fileInput("file2", "Second: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  

  observeEvent(input$start1, {
    range = input$range
    range_seq = seq(range[1],range[2])
    expro = otudata()
    samples = sampledata()
    
    if(is.null(expro) | is.null(samples)){
      warning("Please upload files!")
    }
    else{
      output$detectfig1 <- renderPlot({
        
        datExpr = as.data.frame(t(expro))
        nGenes = ncol(datExpr)
        nSamples = nrow(datExpr)
        
        ##软阈值筛选##
        powers <<- range_seq
        sft <<- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        par(mfrow = c(1,2))
        cex1 = 0.9
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
             main = paste("Scale independence"));
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels=powers,cex=cex1,col="red");
        abline(h=0.90,col="red")
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
             xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
             main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
      })
    }
    
  })
  
  
  output$downloadGraph1_1 <- downloadHandler(
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
      
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
           main = paste("Scale independence"));
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           labels=powers,cex=cex1,col="red");
      abline(h=0.90,col="red")
      dev.off()
    }
  )
  
  output$downloadGraph1_2 <- downloadHandler(
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
      
      plot(sft$fitIndices[,1], sft$fitIndices[,5],
           xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
           main = paste("Mean connectivity"))
      text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
      
      dev.off()
    }
  )
  
  
  ###############################
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  otudata2 <- reactive({
    infile <- input$file3
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  
  sampledata2 <- reactive({
    infile <- input$file4
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  observeEvent(input$reset2, {
    values$file <- NULL
    output$file3 <- renderUI({
      fileInput("file3", "First: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$detectfig2 <- renderPlot({
    return(NULL)
  })
  
  
  observeEvent(input$reset2, {
    values$file <- NULL
    output$file4 <- renderUI({
      fileInput("file4", "Second: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$start2, {
    
    expro = otudata2()
    samples = sampledata2()
    
    if(is.null(expro) | is.null(samples)){
      warning("Please upload files!")
    }
    else{
      output$detectfig2 <- renderPlot({
        
        datExpr = as.data.frame(t(expro))
        nGenes = ncol(datExpr)
        nSamples = nrow(datExpr)
        
        power_dat = input$power
        cortype_dat = input$cortype
        maxblocksize_dat = input$maxblocksize
        tomtype_dat = input$tomtype
        networktype_dat = input$networktype
        minmodulesize_dat = input$minmodulesize
        reassignthreshold_dat = input$reassignthreshold
        mergecutheight_dat = input$mergecutheight
        
        print(power_dat)
        print(cortype_dat)
        print(maxblocksize_dat)
        print(tomtype_dat)
        print(networktype_dat)
        print(minmodulesize_dat)
        print(reassignthreshold_dat)
        print(mergecutheight_dat)
        
        net = blockwiseModules(datExpr, 
                               power = power_dat,
                               corType = cortype_dat,
                               maxBlockSize = maxblocksize_dat,
                               TOMType = tomtype_dat,
                               networkType = networktype_dat,
                               minModuleSize = minmodulesize_dat,
                               reassignThreshold = reassignthreshold_dat, 
                               mergeCutHeight = mergecutheight_dat,
                               numericLabels = TRUE, 
                               pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = paste(getwd(),"asv-otu-TOM",sep = '/'),
                               verbose = 3)
        
        unlink(paste(getwd(),"asv-otu-TOM-block.1.RData",sep = '/'), recursive=TRUE)
        
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        table(moduleColors)
        MEs = net$MEs;
        geneTree = net$dendrograms[[1]];
        
        mergedColors = labels2colors(net$colors)
        # Plot the dendrogram and the module colors underneath
        plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
        
        })
    }
  })
  

  output$downloadGraph2 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot2)
    },
    content = function(file){
      if(input$extPlot2 == 'pdf'){
        pdf(file)
      }else if(input$extPlot2 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      
      
      dev.off()
    }
  )
  
  
  ###############################
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  otudata3 <- reactive({
    infile <- input$file5
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  
  sampledata3 <- reactive({
    infile <- input$file6
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  observeEvent(input$reset3, {
    values$file <- NULL
    output$file5 <- renderUI({
      fileInput("file5", "First: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$detectfig3 <- renderPlot({
    return(NULL)
  })
  
  
  output$table <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  observeEvent(input$reset3, {
    values$file <- NULL
    output$file6 <- renderUI({
      fileInput("file6", "Second: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$start3, {
    
    expro = otudata3()
    samples = sampledata3()
    
    if(is.null(expro) | is.null(samples)){
      warning("Please upload files!")
    }
    else{
      output$detectfig3 <- renderPlot({
        
        datExpr = as.data.frame(t(expro))
        nGenes = ncol(datExpr)
        nSamples = nrow(datExpr)
        
        power_dat = input$power1
        cortype_dat = input$cortype1
        maxblocksize_dat = input$maxblocksize1
        tomtype_dat = input$tomtype1
        networktype_dat = input$networktype1
        minmodulesize_dat = input$minmodulesize1
        reassignthreshold_dat = input$reassignthreshold1
        mergecutheight_dat = input$mergecutheight1
        
        print(power_dat)
        print(cortype_dat)
        print(maxblocksize_dat)
        print(tomtype_dat)
        print(networktype_dat)
        print(minmodulesize_dat)
        print(reassignthreshold_dat)
        print(mergecutheight_dat)
        
        net = blockwiseModules(datExpr, 
                               power = power_dat, 
                               corType = cortype_dat,
                               maxBlockSize = maxblocksize_dat,
                               TOMType = tomtype_dat,
                               networkType = networktype_dat,
                               minModuleSize = minmodulesize_dat,
                               reassignThreshold = reassignthreshold_dat, 
                               mergeCutHeight = mergecutheight_dat,
                               numericLabels = TRUE, 
                               pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = paste(getwd(),"asv-otu-TOM",sep = '/'),
                               verbose = 3)
        
        unlink(paste(getwd(),"asv-otu-TOM-block.1.RData",sep = '/'), recursive=TRUE)
        
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        table(moduleColors)
        MEs = net$MEs;
        geneTree = net$dendrograms[[1]];
        
        
        mergedColors = labels2colors(net$colors)
   
        moduleLabelsAutomatic = net$colors
        moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
        moduleColorsWW = moduleColorsAutomatic
        MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
        MEsWW = orderMEs(MEs0)
        modTraitCor = cor(MEsWW, samples, use = "p")
        colnames(MEsWW)
        modlues=MEsWW
        
        row.names(modTraitCor) = paste("ME",seq(1,nrow(modTraitCor)),sep = '_')
        modTraitP = corPvalueStudent(modTraitCor, nSamples)
        textMatrix <<- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
        dim(textMatrix) = dim(modTraitCor)
        
        dat_cor <<- modTraitCor
        dat_p <<- modTraitP
        modlues_dat <<- modlues
        samples_dat <<- samples
        
        labeledHeatmap(Matrix = dat_cor, 
                       xLabels = colnames(samples_dat), 
                       yLabels = row.names(dat_cor), 
                       cex.lab = 0.5,  
                       yColorWidth=0.01, 
                       xColorWidth = 0.01,
                       ySymbols = colnames(modlues_dat), 
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50), 
                       textMatrix = textMatrix, 
                       setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
                       main = paste("Module-trait relationships"))
        
      })
      
      output$table <- renderDataTable({
        a = c(dat_cor[,1],dat_cor[,2])
        b = c(dat_p[,1],dat_p[,2])
        c = c(rep(colnames(dat_cor)[1],nrow(dat_cor)), 
              rep(colnames(dat_cor)[2],nrow(dat_cor)))
        d = c(row.names(dat_cor), row.names(dat_cor))
        
        for (i in 3:ncol(dat_cor)) {
          a = c(a,dat_cor[,i])
          b = c(b,dat_p[,i])
          c = c(c,rep(colnames(dat_cor)[i],nrow(dat_cor)))
          d = c(d, row.names(dat_cor))
        }
        
        dat_table = data.frame(Cor = a, P_val = b, Group = c, Module = d)
        final <<- dat_table
        
      })
    }
  })
  
  
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final,file,row.names = T,quote = F)
    }
  )
  
  
  output$downloadGraph3 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot3)
    },
    content = function(file){
      if(input$extPlot3 == 'pdf'){
        pdf(file)
      }else if(input$extPlot3 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      
      labeledHeatmap(Matrix = dat_cor, 
                     xLabels = colnames(samples_dat), 
                     yLabels = row.names(dat_cor), 
                     cex.lab = 0.5,  
                     yColorWidth=0.01, 
                     xColorWidth = 0.01,
                     ySymbols = colnames(modlues_dat), 
                     colorLabels = FALSE,
                     colors = blueWhiteRed(50), 
                     textMatrix = textMatrix, 
                     setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
                     main = paste("Module-trait relationships"))
      
      dev.off()
    }
  )
  
  
  
  ###############################
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  otudata4 <- reactive({
    infile <- input$file7
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  
  sampledata4 <- reactive({
    infile <- input$file8
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  observeEvent(input$reset4, {
    values$file <- NULL
    output$file7 <- renderUI({
      fileInput("file7", "First: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table1 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  observeEvent(input$reset4, {
    values$file <- NULL
    output$file8 <- renderUI({
      fileInput("file8", "Second: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$file9 <- renderUI({
    selectInput(inputId ="site_module", 
                label ="Select module", 
                choices = c('No module to select'))
  })
  
  output$file10 <- renderUI({
    sliderInput("weight_cutoff", "Weight cutoff:",
                min = 0, max = 1,
                value = 0.02) 
  })
  
  
  output$detectfig4 <- renderPlot({
    return(NULL)
  })
  
  
  
  observeEvent(input$start4, {
    
    expro = otudata4()
    samples = sampledata4()
    
    if(is.null(expro) | is.null(samples)){
      warning("Please upload files!")
    }
    else{
      output$detectfig4 <- renderPlot({
        
        datExpr = as.data.frame(t(expro))
        nGenes = ncol(datExpr)
        nSamples = nrow(datExpr)
        datExpr_dat <<- datExpr
        
        power_dat = input$power2
        cortype_dat = input$cortype2
        maxblocksize_dat = input$maxblocksize2
        tomtype_dat = input$tomtype2
        networktype_dat = input$networktype2
        minmodulesize_dat = input$minmodulesize2
        reassignthreshold_dat = input$reassignthreshold2
        mergecutheight_dat = input$mergecutheight2
        
        print(power_dat)
        print(cortype_dat)
        print(maxblocksize_dat)
        print(tomtype_dat)
        print(networktype_dat)
        print(minmodulesize_dat)
        print(reassignthreshold_dat)
        print(mergecutheight_dat)
        
        net = blockwiseModules(datExpr, 
                               power = power_dat, 
                               corType = cortype_dat,
                               maxBlockSize = maxblocksize_dat,
                               TOMType = tomtype_dat,
                               networkType = networktype_dat,
                               minModuleSize = minmodulesize_dat,
                               reassignThreshold = reassignthreshold_dat, 
                               mergeCutHeight = mergecutheight_dat,
                               numericLabels = TRUE, 
                               pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = paste(getwd(),"asv-otu-TOM",sep = '/'),
                               verbose = 3)
        
        unlink(paste(getwd(),"asv-otu-TOM-block.1.RData",sep = '/'), recursive=TRUE)
        
        net_dat <<- net
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        table(moduleColors)
        
        moduleLabelsAutomatic = net$colors
        moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
        moduleColorsWW = moduleColorsAutomatic
        MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
        MEsWW = orderMEs(MEs0)
        MEsWW_dat <<- MEsWW
        modTraitCor = cor(MEsWW, samples, use = "p")
        colnames(MEsWW)
        modlues=MEsWW
        
        TOM = TOMsimilarityFromExpr(datExpr, power = power_dat)
        TOM_dat <<- TOM
        moduleLabels = net_dat$colors
        moduleColors = labels2colors(net_dat$colors)
        table(moduleColors)
        MEs = net$MEs;
        geneTree = net$dendrograms[[1]];
        
        mergedColors = labels2colors(net$colors)
        
        moduleLabelsAutomatic = net$colors
        moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
        moduleColorsWW = moduleColorsAutomatic
        MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
        MEsWW = orderMEs(MEs0)
        modTraitCor = cor(MEsWW, samples, use = "p")
        colnames(MEsWW)
        modlues=MEsWW
        
        row.names(modTraitCor) = paste("ME",seq(1,nrow(modTraitCor)),sep = '_')
        modTraitP = corPvalueStudent(modTraitCor, nSamples)
        textMatrix <<- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
        dim(textMatrix) = dim(modTraitCor)
        
        dat_cor <<- modTraitCor
        dat_p <<- modTraitP
        modlues_dat <<- modlues
        samples_dat <<- samples
        
        index_dat = data.frame(id = paste('ME',seq(1:length(colnames(MEsWW_dat))),sep = '_'),
                               color = colnames(MEsWW_dat))
        index_dat[,2] = gsub("ME","",index_dat[,2])
        
        index_select_dat <<- index_dat
        
     
        labeledHeatmap(Matrix = dat_cor, 
                       xLabels = colnames(samples_dat), 
                       yLabels = row.names(dat_cor), 
                       cex.lab = 0.5,  
                       yColorWidth=0.01, 
                       xColorWidth = 0.01,
                       ySymbols = colnames(modlues_dat), 
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50), 
                       textMatrix = textMatrix, 
                       setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
                       main = paste("Module-trait relationships"))
      })
      
      
      # 顺序执行程序,先运行output$detectfig4,在运行output$file9
      output$file9 <- renderUI({
        selectInput(inputId ="site_module", label ="Select module", 
                    choices = index_select_dat[,1],
                    selected = index_select_dat[1,1])
      })
    }
  })
  
  
  observeEvent(input$site_module, {
    print(input$site_module)
    if (input$site_module != 'No module to select') {
      output$table1 <- renderDataTable({
        
        modules_index = as.character(input$site_module)
        modules_index_dat <<- modules_index
        
        moduleColors_dat <<- labels2colors(net_dat$colors)
        id_index = which(index_select_dat[,1]==modules_index)
        id_index_dat <<- id_index
        modules = index_select_dat[id_index,2]
        
        # Select module probes选择模块探测
        probes = names(datExpr_dat)
        inModule = is.finite(match(moduleColors_dat, modules))
        
        modProbes = probes[inModule]
        
        modTOM = TOM_dat[inModule, inModule]
        dimnames(modTOM) = list(modProbes, modProbes)
        modTOM_dat <<- modTOM
        modProbes_dat <<- modProbes
        moduleColors_dat1 <<- moduleColors_dat[inModule]
        # Export the network into edge and node list files Cytoscape can read
        weight_cutoff_dat <<- input$weight_cutoff
        cyt = exportNetworkToCytoscape(modTOM_dat,weighted = TRUE,
                                       threshold = weight_cutoff_dat,
                                       nodeNames = modProbes_dat,
                                       #altNodeNames = modGenes,
                                       nodeAttr = moduleColors_dat1)
        
        
        unlink(paste(getwd(),paste("AS-green-FPKM-One-step-CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),sep = '/'), recursive=TRUE)
        unlink(paste(getwd(),paste("AS-green-FPKM-One-step-CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),sep = '/'), recursive=TRUE)
        
        final <<- cyt$edgeData
      })
    }
  })
  
  
  output$downloadTable1 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      cyt = exportNetworkToCytoscape(modTOM_dat,weighted = TRUE,
                                     threshold = weight_cutoff_dat,
                                     nodeNames = modProbes_dat,
                                     #altNodeNames = modGenes,
                                     nodeAttr = moduleColors_dat1)
      
      write.csv(cyt$edgeData,file,row.names = T,quote = F)
      
      
    }
  )
  
  output$downloadTable2 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      cyt = exportNetworkToCytoscape(modTOM_dat,weighted = TRUE,
                                     threshold = weight_cutoff_dat,
                                     nodeNames = modProbes_dat,
                                     #altNodeNames = modGenes,
                                     nodeAttr = moduleColors_dat1)
      
      write.csv(cyt$nodeData,file,row.names = T,quote = F)
    }
  )
  
  ###############################
  values <- reactiveValues(
    file = NULL
  )
  
  
  otudata5 <- reactive({
    infile <- input$file11
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  
  sampledata5 <- reactive({
    infile <- input$file12
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  observeEvent(input$reset5, {
    values$file <- NULL
    output$file11 <- renderUI({
      fileInput("file11", "First: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$detectfig5 <- renderPlot({
    return(NULL)
  })
  
  
  observeEvent(input$reset5, {
    values$file <- NULL
    output$file12 <- renderUI({
      fileInput("file12", "Second: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$start5, {
    
    expro = otudata5()
    samples = sampledata5()
    
    if(is.null(expro) | is.null(samples)){
      warning("Please upload files!")
    }
    else{
      output$detectfig5 <- renderPlot({
        
        datExpr = as.data.frame(t(expro))
        nGenes = ncol(datExpr)
        nSamples = nrow(datExpr)
        
        power_dat = input$power3
        cortype_dat = input$cortype3
        maxblocksize_dat = input$maxblocksize3
        tomtype_dat = input$tomtype3
        networktype_dat = input$networktype3
        minmodulesize_dat = input$minmodulesize3
        reassignthreshold_dat = input$reassignthreshold3
        mergecutheight_dat = input$mergecutheight3
        
        print(power_dat)
        print(cortype_dat)
        print(maxblocksize_dat)
        print(tomtype_dat)
        print(networktype_dat)
        print(minmodulesize_dat)
        print(reassignthreshold_dat)
        print(mergecutheight_dat)
        
        net = blockwiseModules(datExpr, 
                               power = power_dat,
                               corType = cortype_dat,
                               maxBlockSize = maxblocksize_dat,
                               TOMType = tomtype_dat,
                               networkType = networktype_dat,
                               minModuleSize = minmodulesize_dat,
                               reassignThreshold = reassignthreshold_dat, 
                               mergeCutHeight = mergecutheight_dat,
                               numericLabels = TRUE, 
                               pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = paste(getwd(),"asv-otu-TOM",sep = '/'),
                               verbose = 3)
        
        unlink(paste(getwd(),"asv-otu-TOM-block.1.RData",sep = '/'), recursive=TRUE)
        
        moduleLabels = net$colors
        moduleColors <<- labels2colors(net$colors)
        table(moduleColors)
        MEs = net$MEs
        geneTree <<- net$dendrograms[[1]];
        
        mergedColors = labels2colors(net$colors)
        
        ## 可视化基因网络
        dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = power_dat);
        # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
        plotTOM <<- dissTOM^7;
        diag(plotTOM) = NA
        TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot (All ASV/OTU)")
        
      })
    }
  })
  
  
  output$downloadGraph5 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot5)
    },
    content = function(file){
      if(input$extPlot5 == 'pdf'){
        pdf(file)
      }else if(input$extPlot5 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      TOMplot(plotTOM, geneTree, moduleColors, 
              main = "Network heatmap plot (All ASV/OTU)")
      
      dev.off()
    }
  )
  
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50700)


