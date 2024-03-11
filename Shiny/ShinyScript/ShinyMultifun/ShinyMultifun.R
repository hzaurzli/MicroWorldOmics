library(dplyr)
library(ggplot2)
library(shiny)
library(DT)
library(ggrepel)
library(tidyr)
library(shinycssloaders)
library(shinythemes)
library(ggClusterNet)
library(phyloseq)
library(tidyverse)
library(ape)
library(sna)
library(tidyfst)


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
  navbarPage("Microbiome comprehensive analysis", 
             theme = shinytheme("lumen"),
   # tabPanel 是导航栏换页
   tabPanel("Microbial diversity", fluid = TRUE,
            tags$style(button_color_css),
            # Sidebar layout with a input and output definitions
            sidebarLayout(
              sidebarPanel(
                h3(strong('Microbial diversity')),
                hr(),
                fluidRow(uiOutput('file1'),
                         uiOutput('file2'),
                         uiOutput('file3'),
                                         
                         actionButton('reset', 'RESET'),
                         actionButton('start', 'START'),
                         hr(),
                         h4('Download all results'),
                         downloadButton("downloadTable", "Download table"),
                         br(),
                         br(),
                         radioButtons('extPlot', 'Plot output format',
                                      choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                         downloadButton("downloadGraph", "Download graph"),
                         ),
                ),
              mainPanel(
                titlePanel("Microbial diversity table and graph"),
                withSpinner(dataTableOutput(outputId = "table")),
                hr(),
                plotOutput(outputId = "detectfig")
                )
            )
        ),
       ##################
     navbarMenu("High-dimensional analysis",
        tabPanel("DCA analysis", fluid = TRUE,
           tags$style(button_color_css),
           # Sidebar layout with a input and output definitions
           sidebarLayout(
             sidebarPanel(
               h3(strong('DCA analysis')),
               hr(),
               fluidRow(uiOutput('file4'),
                        uiOutput('file5'),
                        uiOutput('file6'),
                    
                        actionButton('reset1', 'RESET'),
                        actionButton('start1', 'START'),
                        hr(),
                        h4('Download all results'),
                        downloadButton("downloadTable1", "Download table"),
                        br(),
                        br(),
                        radioButtons('extPlot1', 'Plot output format',
                                     choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                        downloadButton("downloadGraph1", "Download graph"),
                        
                     ),
             ),
             mainPanel(
               titlePanel("DCA analysis table and graph"),
               withSpinner(dataTableOutput(outputId = "table1")),
               hr(),
               plotOutput(outputId = "detectfig1")
             )
           )
        ),
        
        tabPanel("CCA analysis", fluid = TRUE, 
           tags$style(button_color_css),
           # Sidebar layout with a input and output definitions
           sidebarLayout(
             sidebarPanel(
               h3(strong('CCA analysis')),
               hr(),
               fluidRow(uiOutput('file7'),
                        uiOutput('file8'),
                        uiOutput('file9'),
                        
                        actionButton('reset2', 'RESET'),
                        actionButton('start2', 'START'),
                        hr(),
                        h4('Download all results'),
                        downloadButton("downloadTable2", "Download table"),
                        br(),
                        br(),
                        radioButtons('extPlot2', 'Plot output format',
                                     choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                        downloadButton("downloadGraph2", "Download graph"),
               ),
             ),
             mainPanel(
               titlePanel("CCA analysis table and graph"),
               withSpinner(dataTableOutput(outputId = "table2")),
               hr(),
               plotOutput(outputId = "detectfig2")
             )
           )
        ),
        tabPanel("RDA analysis", fluid = TRUE, 
           tags$style(button_color_css),
           # Sidebar layout with a input and output definitions
           sidebarLayout(
             sidebarPanel(
               h3(strong('RDA analysis')),
               hr(),
               fluidRow(uiOutput('file10'),
                        uiOutput('file11'),
                        uiOutput('file12'),
                        
                        
                        actionButton('reset3', 'RESET'),
                        actionButton('start3', 'START'),
                        hr(),
                        h4('Download all results'),
                        downloadButton("downloadTable3", "Download table"),
                        br(),
                        br(),
                        radioButtons('extPlot3', 'Plot output format',
                                     choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                        downloadButton("downloadGraph3", "Download graph"),
               ),
             ),
             mainPanel(
               titlePanel("RDA analysis table and graph"),
               withSpinner(dataTableOutput(outputId = "table3")),
               hr(),
               plotOutput(outputId = "detectfig3")
             )
           )
        ),
        tabPanel("MDS analysis", fluid = TRUE, 
                 tags$style(button_color_css),
                 # Sidebar layout with a input and output definitions
                 sidebarLayout(
                   sidebarPanel(
                     h3(strong('MDS analysis')),
                     hr(),
                     fluidRow(uiOutput('file13'),
                              uiOutput('file14'),
                              uiOutput('file15'),
                              
                              
                              actionButton('reset4', 'RESET'),
                              actionButton('start4', 'START'),
                              hr(),
                              h4('Download all results'),
                              downloadButton("downloadTable4", "Download table"),
                              br(),
                              br(),
                              radioButtons('extPlot4', 'Plot output format',
                                           choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                              downloadButton("downloadGraph4", "Download graph"),),
                   ),
                   mainPanel(
                     titlePanel("MDS analysis table and graph"),
                     withSpinner(dataTableOutput(outputId = "table4")),
                     hr(),
                     plotOutput(outputId = "detectfig4")
                   )
                 )
        ),
        tabPanel("PCoA analysis", fluid = TRUE, 
                 tags$style(button_color_css),
                 # Sidebar layout with a input and output definitions
                 sidebarLayout(
                   sidebarPanel(
                     h3(strong('PCoA analysis')),
                     hr(),
                     fluidRow(uiOutput('file16'),
                              uiOutput('file17'),
                              uiOutput('file18'),
                              
                              
                              actionButton('reset5', 'RESET'),
                              actionButton('start5', 'START'),
                              hr(),
                              h4('Download all results'),
                              downloadButton("downloadTable5", "Download table"),
                              br(),
                              br(),
                              radioButtons('extPlot5', 'Plot output format',
                                           choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                              downloadButton("downloadGraph5", "Download graph"),),
                   ),
                   mainPanel(
                     titlePanel("PCoA analysis table and graph"),
                     withSpinner(dataTableOutput(outputId = "table5")),
                     hr(),
                     plotOutput(outputId = "detectfig5")
                   )
                 )
        ),
        tabPanel("PCA analysis", fluid = TRUE, 
                 tags$style(button_color_css),
                 # Sidebar layout with a input and output definitions
                 sidebarLayout(
                   sidebarPanel(
                     h3(strong('PCA analysis')),
                     hr(),
                     fluidRow(uiOutput('file19'),
                              uiOutput('file20'),
                              uiOutput('file21'),
                              
                              
                              actionButton('reset6', 'RESET'),
                              actionButton('start6', 'START'),
                              hr(),
                              h4('Download all results'),
                              downloadButton("downloadTable6", "Download table"),
                              br(),
                              br(),
                              radioButtons('extPlot6', 'Plot output format',
                                           choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                              downloadButton("downloadGraph6", "Download graph"),),
                   ),
                   mainPanel(
                     titlePanel("PCA analysis table and graph"),
                     withSpinner(dataTableOutput(outputId = "table6")),
                     hr(),
                     plotOutput(outputId = "detectfig6")
                   )
                 )
        ),
        tabPanel("LDA analysis", fluid = TRUE, 
                 tags$style(button_color_css),
                 # Sidebar layout with a input and output definitions
                 sidebarLayout(
                   sidebarPanel(
                     h3(strong('LDA analysis')),
                     hr(),
                     fluidRow(uiOutput('file22'),
                              uiOutput('file23'),
                              uiOutput('file24'),
                              
                              
                              actionButton('reset7', 'RESET'),
                              actionButton('start7', 'START'),
                              hr(),
                              h4('Download all results'),
                              downloadButton("downloadTable7", "Download table"),
                              br(),
                              br(),
                              radioButtons('extPlot7', 'Plot output format',
                                           choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                              downloadButton("downloadGraph7", "Download graph"),),
                   ),
                   mainPanel(
                     titlePanel("LDA analysis table and graph"),
                     withSpinner(dataTableOutput(outputId = "table7")),
                     hr(),
                     plotOutput(outputId = "detectfig7")
                   )
                 )
        ),
        
        tabPanel("NMDS analysis", fluid = TRUE, 
                 tags$style(button_color_css),
                 # Sidebar layout with a input and output definitions
                 sidebarLayout(
                   sidebarPanel(
                     h3(strong('NMDS analysis')),
                     hr(),
                     fluidRow(uiOutput('file25'),
                              uiOutput('file26'),
                              uiOutput('file27'),
                              
                              
                              actionButton('reset8', 'RESET'),
                              actionButton('start8', 'START'),
                              hr(),
                              h4('Download all results'),
                              downloadButton("downloadTable8", "Download table"),
                              br(),
                              br(),
                              radioButtons('extPlot8', 'Plot output format',
                                           choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                              downloadButton("downloadGraph8", "Download graph"),),
                   ),
                   mainPanel(
                     titlePanel("NMDS analysis table and graph"),
                     withSpinner(dataTableOutput(outputId = "table8")),
                     hr(),
                     plotOutput(outputId = "detectfig8")
                   )
                 )
              ),
        ),
   
      ##################
      navbarMenu("Relative abundance analysis",
           tabPanel("Relative abundance bar plot", fluid = TRUE,
                    tags$style(button_color_css),
                    # Sidebar layout with a input and output definitions
                    sidebarLayout(
                      sidebarPanel(
                        h3(strong('Relative abundance bar plot')),
                        hr(),
                        fluidRow(
                          uiOutput('file28'),
                          uiOutput('file29'),
                          uiOutput('file30'),
                          
                          selectInput("class_levels", "Levels:", 
                                      choices = c("Kingdom","Phylum",
                                                  "Class","Order","Family",
                                                  "Genus","Species"),
                                      selected = 'Genus' 
                          ),
                          
                          actionButton('reset9', 'RESET'),
                          actionButton('start9', 'START'),
                          hr(),
                          h4('Download all results'),
                          downloadButton("downloadTable9", "Download table"),
                          br(),
                          br(),
                          radioButtons('extPlot9', 'Plot output format',
                                       choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                          downloadButton("downloadGraph9", "Download graph")
                        
                          
                        ),
                      ),
                      mainPanel(
                        titlePanel("Relative abundance table and graph"),
                        withSpinner(dataTableOutput(outputId = "table9")),
                        hr(),
                        plotOutput(outputId = "detectfig9")
                      )
                    )
                ),
           
           tabPanel("Relative abundance ternary diagrams", fluid = TRUE,
                    tags$style(button_color_css),
                    # Sidebar layout with a input and output definitions
                    sidebarLayout(
                      sidebarPanel(
                        h3(strong('Relative abundance ternary diagrams')),
                        hr(),
                        fluidRow(
                          uiOutput('file31'),
                          uiOutput('file32'),
                          uiOutput('file33'),
                          
                          selectInput("class_levels1", "Levels:", 
                                      choices = c("Kingdom","Phylum",
                                                  "Class","Order","Family",
                                                  "Genus","Species"),
                                      selected = 'Phylum' 
                          ),
                          
                          actionButton('reset10', 'RESET'),
                          actionButton('start10', 'START'),
                          hr(),
                          h4('Download all results'),
                          radioButtons('extPlot10', 'Plot output format',
                                       choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                          downloadButton("downloadGraph10", "Download graph")
                          
                          
                        ),
                      ),
                      mainPanel(
                        titlePanel("Relative abundance ternary diagrams"),
                        hr(),
                        textOutput(outputId = "error_text"),
                        plotOutput(outputId = "detectfig10"),
                      )
                    )
           ),
           
           tabPanel("Relative abundance circlize", fluid = TRUE,
                tags$style(button_color_css),
                # Sidebar layout with a input and output definitions
                sidebarLayout(
                  sidebarPanel(
                    h3(strong('Relative abundance circlize')),
                    hr(),
                    fluidRow(
                      uiOutput('file34'),
                      uiOutput('file35'),
                      uiOutput('file36'),
                      
                      selectInput("class_levels2", "Levels:", 
                                  choices = c("Kingdom","Phylum",
                                              "Class","Order","Family",
                                              "Genus","Species"),
                                  selected = 'Phylum' 
                      ),
                      sliderInput("size", label = "Label text size:",
                                  min = 0, max = 1.5, value = 0.8),
                      br(),
                      actionButton('reset11', 'RESET'),
                      actionButton('start11', 'START'),
                      hr(),
                      h4('Download all results'),
                      radioButtons('extPlot11', 'Plot output format',
                                   choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                      downloadButton("downloadGraph11", "Download graph")
                      
                      
                    ),
                  ),
                  mainPanel(
                    titlePanel("Relative abundance circlize"),
                    hr(),
                    plotOutput(outputId = "detectfig11"),
                  )
                )
       ),
    ),

     ##################
     navbarMenu("Microbial differential analysis",
            tabPanel("Wilcox.test", fluid = TRUE,
                     tags$style(button_color_css),
                     # Sidebar layout with a input and output definitions
                     sidebarLayout(
                       sidebarPanel(
                         h3(strong('Wilcox.test')),
                         hr(),
                         fluidRow(
                           uiOutput('file37'),
                           uiOutput('file38'),
                           uiOutput('file39'),
                           
                           tags$h5("Step4: Choose comparative group"),
                           uiOutput('file40'),
                           uiOutput('file41'),
                           br(),
                           sliderInput("Pval", "P value cut-off",
                                       min = 0, max = 0.1,
                                       value = 0.05),
                           br(),
                           actionButton('reset12', 'RESET'),
                           actionButton('start12', 'START'),
                           hr(),
                           h4('Download all results'),
                           downloadButton("downloadTable12", "Download table"),
                    
                         ),
                       ),
                       mainPanel(
                         titlePanel("Wilcox.test result table"),
                         withSpinner(dataTableOutput(outputId = "table12")),
                       )
                     )
                  ),
              
                    
                    tabPanel("DESeq2", fluid = TRUE,
                             tags$style(button_color_css),
                             # Sidebar layout with a input and output definitions
                             sidebarLayout(
                               sidebarPanel(
                                 h3(strong('DESeq2')),
                                 hr(),
                                 fluidRow(
                                   uiOutput('file42'),
                                   uiOutput('file43'),
                                   uiOutput('file44'),
                                   
                                   tags$h5("Step4: Choose comparative group"),
                                   uiOutput('file45'),
                                   uiOutput('file46'),
                                   br(),
                                   sliderInput("Pval1", "P value cut-off",
                                               min = 0, max = 0.1,
                                               value = 0.05),
                                   br(),
                                   actionButton('reset13', 'RESET'),
                                   actionButton('start13', 'START'),
                                   hr(),
                                   h4('Download all results'),
                                   downloadButton("downloadTable13", 
                                                  "Download table"),
                                   
                                 ),
                               ),
                               mainPanel(
                                 titlePanel("DESeq2 result table"),
                                 withSpinner(dataTableOutput(outputId = "table13")),
                               )
                             )
                    ),
                    
                    tabPanel("edgeR", fluid = TRUE,
                             tags$style(button_color_css),
                             # Sidebar layout with a input and output definitions
                             sidebarLayout(
                               sidebarPanel(
                                 h3(strong('edgeR')),
                                 hr(),
                                 fluidRow(
                                   uiOutput('file47'),
                                   uiOutput('file48'),
                                   uiOutput('file49'),
                                   
                                   tags$h5("Step4: Choose comparative group"),
                                   uiOutput('file50'),
                                   uiOutput('file51'),
                                   br(),
                                   sliderInput("Pval2", "P value cut-off",
                                               min = 0, max = 0.1,
                                               value = 0.05),
                                   br(),
                                   actionButton('reset14', 'RESET'),
                                   actionButton('start14', 'START'),
                                   hr(),
                                   h4('Download all results'),
                                   downloadButton("downloadTable14", 
                                                  "Download table"),
                                   
                                 ),
                               ),
                               mainPanel(
                                 titlePanel("edgeR result table"),
                                 withSpinner(dataTableOutput(outputId = "table14")),
                               )
                             )
                    ),
   
          ),
           
          # tabPanel 是导航栏换页
          tabPanel("Microbial community network", fluid = TRUE,
            tags$style(button_color_css),
            # Sidebar layout with a input and output definitions
            sidebarLayout(
              sidebarPanel(
                h3(strong('Microbial community network')),
                hr(),
                fluidRow(uiOutput('file52'),
                         uiOutput('file53'),
                         uiOutput('file54'),
                         
                         selectInput("class_levels4", "Levels:", 
                                     choices = c("Kingdom","Phylum",
                                                 "Class","Order","Family",
                                                 "Genus","Species"),
                                     selected = 'Phylum' 
                         ),
                         
                         actionButton('reset15', 'RESET'),
                         actionButton('start15', 'START'),
                         hr(),
                         h4('Download all results'),
                         downloadButton("downloadTable15", "Download table"),
                         br(),
                         br(),
                         radioButtons('extPlot15', 'Plot output format',
                                      choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                         downloadButton("downloadGraph15", "Download graph"),
                ),
              ),
              mainPanel(
                titlePanel("Microbial community network table and graph"),
                withSpinner(dataTableOutput(outputId = "table15")),
                hr(),
                plotOutput(outputId = "detectfig15",width = "80%")
              )
            )
          ),
   
    ),
)

# Define server
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
  
  taxdata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata <- reactive({
    infile <- input$file3
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
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
      fileInput("file2", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$file <- NULL
    output$file3 <- renderUI({
      fileInput("file3", "Step 3: Choose metadata file",
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
    
    index = c("Shannon","Inv_Simpson",
              "Pielou_evenness","Simpson_evenness" ,
              "Richness" ,"Chao1","ACE" )
    
    
    otu = as.matrix(otudata())
    sample = sampledata()
    tax = as.matrix(taxdata())
    
    output$table <- renderDataTable({
      OTU = otu_table(otu, taxa_are_rows = TRUE)
      TAX = tax_table(tax)
      physeq_p = phyloseq(OTU, TAX)
      
      random_tree = rtree(ntaxa(physeq_p), 
                          rooted=TRUE, 
                          tip.label=taxa_names(physeq_p))
      sampledata = sample_data(sample)
      ps = phyloseq(OTU, TAX, sampledata, random_tree)
      
      
      samplesize = min(phyloseq::sample_sums(ps))
      ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
      mapping = phyloseq::sample_data(ps11)
      ps11 = phyloseq::filter_taxa(ps11, function(x) sum(x ) >0 , TRUE); ps11
      mapping$Group = as.factor(mapping$Group)
      count = as.data.frame(t(ggClusterNet::vegan_otu(ps11)))
      alpha=vegan::diversity(count, "shannon")
      x = t(count)
      
      print('11111')
      
      Shannon = vegan::diversity(x)
      Shannon
      Inv_Simpson <- vegan::diversity(x, index = "invsimpson")
      Inv_Simpson
      S <- vegan::specnumber(x);S
      S2 = rowSums(x>0)
      Pielou_evenness <- Shannon/log(S)
      Simpson_evenness <- Inv_Simpson/S
      est <- vegan::estimateR(x)
      Richness <- est[1, ]
      Chao1 <- est[2, ]
      ACE <- est[4, ]
      report = cbind(Shannon, Inv_Simpson, Pielou_evenness, Simpson_evenness,
                     Richness, Chao1,ACE)
      head(report)
      index = merge(mapping,report , by="row.names",all=F)
      sel = c(match("Inv_Simpson",colnames(index)),
              match("Pielou_evenness",colnames(index)),
              match("Simpson_evenness",colnames(index)),
              match("Richness",colnames(index)),
              match("Chao1",colnames(index)),
              match("ACE",colnames(index)),
              match("Shannon",colnames(index)))
      
      
      n = length(sel) + 3
      data = cbind(data.frame(ID = 1:length(index$Group),group = index$Group),index[sel])
      
      print('wwwww')
      
      result = EasyStat::MuiKwWlx2(data = data,num = c(3:(n -1)))
      result1 <<- EasyStat::FacetMuiPlotresultBox(data = data,num = c(3:(n -1)),
                                                  result = result,
                                                  sig_show ="abc",
                                                  ncol = 3 )
      
      
      
      table = result1[[1]]$data
      colnames(table)[3] = c('value')
      final <<- table[,-c(1,4,5)]
    })
    
    output$detectfig <- renderPlot({
      p1_1 <<- result1[[1]] +
        ggplot2::guides(fill = guide_legend(title = NULL))
      p1_1
    })
  })
  
  ## 下载图片写法
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
      # 打印全局变量 p1 (ggplot对象)
      print(p1_1)
      dev.off()
    }
  )
  
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final,file,row.names = T,quote = F)
    }
  )
  
  ##################### DCA
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata1 <- reactive({
    infile <- input$file4
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata1 <- reactive({
    infile <- input$file5
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata1 <- reactive({
    infile <- input$file6
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset1, {
    values$file <- NULL
    output$file4 <- renderUI({
      fileInput("file4", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset1, {
    values$file <- NULL
    output$file5 <- renderUI({
      fileInput("file5", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset1, {
    values$file <- NULL
    output$file6 <- renderUI({
      fileInput("file6", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table1 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start1, {

    otu = as.matrix(otudata1())
    sample = sampledata1()
    tax = as.matrix(taxdata1())
    
    
    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    
    ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    
    #-DCA
    method = "DCA"
    ordi = phyloseq::ordinate(ps1_rela, method="DCA", distance="bray")
    points = ordi$rproj[,1:2]
    colnames(points) = c("x", "y")
    eig = ordi$evals^2
    table = points
    colnames(table) = c("DCA_1","DCA_2")
    
    g = sample_data(ps)$Group %>% unique() %>% length()
    n = sample_data(ps)$Group%>% length()
    o = n/g
    ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    unif = phyloseq::distance(ps, method="bray")
    # adonis
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$Group = as.factor(map$Group)
    colbar = length(levels(map$Group))
    
    points = cbind(points, map[match(rownames(points), rownames(map)), ])
    points$ID = row.names(points)

    output$table1 <- renderDataTable({
      final1 <<- table
    })
    
    output$detectfig1 <- renderPlot({
      p = ggplot(points, aes(x=x, y=y, fill = Group)) +
        geom_point(alpha=.7, size=5, pch = 21) +
        labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
             y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
             title=title1) +
        stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      p1 <<- p + ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
      p1
    })
  })
  
  
  ## 下载图片写法
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
      # 打印全局变量 p1 (ggplot对象)
      print(p1)
      dev.off()
    }
  )
  
  output$downloadTable1 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final1,file,row.names = T,quote = F)
    }
  )
  
  ##################### CCA
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata2 <- reactive({
    infile <- input$file7
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata2 <- reactive({
    infile <- input$file8
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata2 <- reactive({
    infile <- input$file9
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset2, {
    values$file <- NULL
    output$file7 <- renderUI({
      fileInput("file7", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset2, {
    values$file <- NULL
    output$file8 <- renderUI({
      fileInput("file8", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset2, {
    values$file <- NULL
    output$file9 <- renderUI({
      fileInput("file9", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table2 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start2, {
    
    otu = as.matrix(otudata2())
    sample = sampledata2()
    tax = as.matrix(taxdata2())
    
    
    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    
    ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    
    #-CCA
    method = "CCA"
    ordi = phyloseq::ordinate(ps1_rela, method="CCA", distance="bray")
    points = ordi$CA$u[,1:2]
    colnames(points) = c("x", "y")
    eig = ordi$CA$eig^2
    table = points
    colnames(table) = c("CCA_1","CCA_2")
    
    g = sample_data(ps)$Group %>% unique() %>% length()
    n = sample_data(ps)$Group%>% length()
    o = n/g
    ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    unif = phyloseq::distance(ps, method="bray")
    # adonis
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$Group = as.factor(map$Group)
    colbar = length(levels(map$Group))
    
    points = cbind(points, map[match(rownames(points), rownames(map)), ])
    points$ID = row.names(points)
    
    output$table2 <- renderDataTable({
      final2 <<- table
    })
    
    output$detectfig2 <- renderPlot({
      p = ggplot(points, aes(x=x, y=y, fill = Group)) +
        geom_point(alpha=.7, size=5, pch = 21) +
        labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
             y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
             title=title1) +
        stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      p2 <<- p + ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
      p2
    })
  })
  
  
  ## 下载图片写法
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
      # 打印全局变量 p1 (ggplot对象)
      print(p2)
      dev.off()
    }
  )
  
  output$downloadTable2 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final2,file,row.names = T,quote = F)
    }
  )
  
  
  ##################### RDA
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata3 <- reactive({
    infile <- input$file10
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata3 <- reactive({
    infile <- input$file11
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata3 <- reactive({
    infile <- input$file12
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset3, {
    values$file <- NULL
    output$file10 <- renderUI({
      fileInput("file10", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset3, {
    values$file <- NULL
    output$file11 <- renderUI({
      fileInput("file11", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset3, {
    values$file <- NULL
    output$file12 <- renderUI({
      fileInput("file12", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table3 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start3, {
    
    otu = as.matrix(otudata3())
    sample = sampledata3()
    tax = as.matrix(taxdata3())
    
    
    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    
    ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    
    #-RDA
    method = "RDA"
    ordi = phyloseq::ordinate(ps1_rela, method="RDA", distance="bray")
    points = ordi$CA$u[,1:2]
    colnames(points) = c("x", "y")
    eig = ordi$CA$eig
    table = points
    colnames(table) = c("RDA_1","RDA_2")
    
    g = sample_data(ps)$Group %>% unique() %>% length()
    n = sample_data(ps)$Group%>% length()
    o = n/g
    ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    unif = phyloseq::distance(ps, method="bray")
    # adonis
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$Group = as.factor(map$Group)
    colbar = length(levels(map$Group))
    
    points = cbind(points, map[match(rownames(points), rownames(map)), ])
    points$ID = row.names(points)
   
    output$table3 <- renderDataTable({
      final3 <<- table
    })
    
    output$detectfig3 <- renderPlot({
      p = ggplot(points, aes(x=x, y=y, fill = Group)) +
        geom_point(alpha=.7, size=5, pch = 21) +
        labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
             y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
             title=title1) +
        stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      p3 <<- p + ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
      p3
    })
  })
  
  
  ## 下载图片写法
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
      # 打印全局变量 p1 (ggplot对象)
      print(p3)
      dev.off()
    }
  )
  
  output$downloadTable3 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final3,file,row.names = T,quote = F)
    }
  )
  
  
  
  ##################### MDS
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata4 <- reactive({
    infile <- input$file13
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata4 <- reactive({
    infile <- input$file14
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata4 <- reactive({
    infile <- input$file15
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset4, {
    values$file <- NULL
    output$file11 <- renderUI({
      fileInput("file11", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset4, {
    values$file <- NULL
    output$file14 <- renderUI({
      fileInput("file14", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset4, {
    values$file <- NULL
    output$file15 <- renderUI({
      fileInput("file15", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table4 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start4, {
    
    otu = as.matrix(otudata4())
    sample = sampledata4()
    tax = as.matrix(taxdata4())
    
  
    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    
    ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    
    #-MDS
    method = "MDS"
    ordi = phyloseq::ordinate(ps1_rela, method="MDS", distance="bray")
    points = ordi$vectors[,1:2]
    colnames(points) = c("x", "y")
    eig = ordi$values[,1]
    table = points
    colnames(table) = c("MDS_1","MDS_2")
    
    g = sample_data(ps)$Group %>% unique() %>% length()
    n = sample_data(ps)$Group%>% length()
    o = n/g
    ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    unif = phyloseq::distance(ps, method="bray")
    # adonis
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$Group = as.factor(map$Group)
    colbar = length(levels(map$Group))
    
    points = cbind(points, map[match(rownames(points), rownames(map)), ])
    points$ID = row.names(points)
    
    output$table4 <- renderDataTable({
      final4 <<- table
    })
    
    output$detectfig4 <- renderPlot({
      p = ggplot(points, aes(x=x, y=y, fill = Group)) +
        geom_point(alpha=.7, size=5, pch = 21) +
        labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
             y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
             title=title1) +
        stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      p4 <<- p + ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
      p4
    })
  })
  
  
  ## 下载图片写法
  output$downloadGraph4 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot4)
    },
    content = function(file){
      if(input$extPlot4 == 'pdf'){
        pdf(file)
      }else if(input$extPlot4 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      # 打印全局变量 p1 (ggplot对象)
      print(p4)
      dev.off()
    }
  )
  
  output$downloadTable4 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final4,file,row.names = T,quote = F)
    }
  )
  
  
  ##################### PCoA
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata5 <- reactive({
    infile <- input$file16
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata5 <- reactive({
    infile <- input$file17
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata5 <- reactive({
    infile <- input$file18
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset5, {
    values$file <- NULL
    output$file16 <- renderUI({
      fileInput("file16", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset5, {
    values$file <- NULL
    output$file17 <- renderUI({
      fileInput("file17", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset5, {
    values$file <- NULL
    output$file18 <- renderUI({
      fileInput("file18", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table5 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start5, {
    
    otu = as.matrix(otudata5())
    sample = sampledata5()
    tax = as.matrix(taxdata5())
    
    
    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    
    ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    
    #-PCoA
    method = "PCoA"
    unif = phyloseq::distance(ps1_rela , method="bray", type="samples")
    #这里请记住pcoa函数
    pcoa = stats::cmdscale(unif, k=2, eig=T)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y")
    eig = pcoa$eig
    table = points
    colnames(table) = c("PCoA_1","PCoA_2")
    
    g = sample_data(ps)$Group %>% unique() %>% length()
    n = sample_data(ps)$Group%>% length()
    o = n/g
    ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    unif = phyloseq::distance(ps, method="bray")
    # adonis
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$Group = as.factor(map$Group)
    colbar = length(levels(map$Group))
    
    points = cbind(points, map[match(rownames(points), rownames(map)), ])
    points$ID = row.names(points)

    output$table5 <- renderDataTable({   
      final5 <<- table
    })
    
    output$detectfig5 <- renderPlot({
      p = ggplot(points, aes(x=x, y=y, fill = Group)) +
        geom_point(alpha=.7, size=5, pch = 21) +
        labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
             y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
             title=title1) +
        stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      p5 <<- p + ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
      p5
    })
  })
  
  
  ## 下载图片写法
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
      # 打印全局变量 p1 (ggplot对象)
      print(p5)
      dev.off()
    }
  )
  
  output$downloadTable5 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final5,file,row.names = T,quote = F)
    }
  )
  
  
  ##################### PCA
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata6 <- reactive({
    infile <- input$file19
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata6 <- reactive({
    infile <- input$file20
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata6 <- reactive({
    infile <- input$file21
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset6, {
    values$file <- NULL
    output$file19 <- renderUI({
      fileInput("file19", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset6, {
    values$file <- NULL
    output$file20 <- renderUI({
      fileInput("file20", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset6, {
    values$file <- NULL
    output$file21 <- renderUI({
      fileInput("file21", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table6 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start6, {
    
    otu = as.matrix(otudata6())
    sample = sampledata6()
    tax = as.matrix(taxdata6())
    

    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    
    ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    
    #-PCA
    method = "PCA"
    otu_table = as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela)))
    otu.pca = stats::prcomp(t(otu_table), scale.default = TRUE)
    points = otu.pca$x[,1:2]
    colnames(points) = c("x", "y")
    eig=otu.pca$sdev
    eig=eig*eig
    table = points
    colnames(table) = c("PCA_1","PCA_2")
    
    g = sample_data(ps)$Group %>% unique() %>% length()
    n = sample_data(ps)$Group%>% length()
    o = n/g
    ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    unif = phyloseq::distance(ps, method="bray")
    # adonis
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$Group = as.factor(map$Group)
    colbar = length(levels(map$Group))
    
    points = cbind(points, map[match(rownames(points), rownames(map)), ])
    points$ID = row.names(points)

    output$table6 <- renderDataTable({
      final6 <<- table
    })
    
    output$detectfig6 <- renderPlot({
      p = ggplot(points, aes(x=x, y=y, fill = Group)) +
        geom_point(alpha=.7, size=5, pch = 21) +
        labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
             y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
             title=title1) +
        stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      p6 <<- p + ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
      p6
    })
  })
  
  
  ## 下载图片写法
  output$downloadGraph6 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot6)
    },
    content = function(file){
      if(input$extPlot6 == 'pdf'){
        pdf(file)
      }else if(input$extPlot6 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      # 打印全局变量 p1 (ggplot对象)
      print(p6)
      dev.off()
    }
  )
  
  output$downloadTable6 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final6,file,row.names = T,quote = F)
    }
  )
  
  ##################### LDA
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata7 <- reactive({
    infile <- input$file22
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata7 <- reactive({
    infile <- input$file23
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata7 <- reactive({
    infile <- input$file24
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset7, {
    values$file <- NULL
    output$file22 <- renderUI({
      fileInput("file22", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset7, {
    values$file <- NULL
    output$file23 <- renderUI({
      fileInput("file23", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset7, {
    values$file <- NULL
    output$file24 <- renderUI({
      fileInput("file24", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table7 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start7, {
    
    otu = as.matrix(otudata7())
    sample = sampledata7()
    tax = as.matrix(taxdata7())
    
 
    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    
    ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    
    #-LDA
    method = "LDA"
    otu_table = as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela)))
    data = t(otu_table)
    data = as.data.frame(data)
    data = scale(data, center = TRUE, scale = TRUE)
    dim(data)
    data1 = data[,1:10]
    map = as.data.frame(sample_data(ps1_rela))
    model = MASS::lda(data, map$Group)
    ord_in = model
    axes = c(1:2)
    points = data.frame(predict(ord_in)$x[, axes])
    colnames(points) = c("x", "y")
    eig= ord_in$svd^2
    
    table = points
    colnames(table) = c("LDA_1","LDA_2")
    
    g = sample_data(ps)$Group %>% unique() %>% length()
    n = sample_data(ps)$Group%>% length()
    o = n/g
    ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    unif = phyloseq::distance(ps, method="bray")
    # adonis
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$Group = as.factor(map$Group)
    colbar = length(levels(map$Group))
    
    points = cbind(points, map[match(rownames(points), rownames(map)), ])
    points$ID = row.names(points)

    output$table7 <- renderDataTable({
      final7 <<- table
    })
    
    output$detectfig7 <- renderPlot({
      p = ggplot(points, aes(x=x, y=y, fill = Group)) +
        geom_point(alpha=.7, size=5, pch = 21) +
        labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
             y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
             title=title1) +
        stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      p7 <<- p + ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
      p7
    })
  })
  
  
  ## 下载图片写法
  output$downloadGraph7 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot7)
    },
    content = function(file){
      if(input$extPlot7 == 'pdf'){
        pdf(file)
      }else if(input$extPlot7 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      # 打印全局变量 p1 (ggplot对象)
      print(p7)
      dev.off()
    }
  )
  
  output$downloadTable7 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final7,file,row.names = T,quote = F)
    }
  )
  
  
  ##################### NMDS
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata8 <- reactive({
    infile <- input$file25
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata8 <- reactive({
    infile <- input$file26
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata8 <- reactive({
    infile <- input$file27
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset8, {
    values$file <- NULL
    output$file25 <- renderUI({
      fileInput("file25", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset8, {
    values$file <- NULL
    output$file26 <- renderUI({
      fileInput("file26", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset8, {
    values$file <- NULL
    output$file27 <- renderUI({
      fileInput("file27", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table8 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start8, {
    
    otu = as.matrix(otudata8())
    sample = sampledata8()
    tax = as.matrix(taxdata8())
    

    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    
    ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    
    #-NMDS
    method = "NMDS"
    ordi = phyloseq::ordinate(ps1_rela, method="NMDS", distance="bray")
    points = ordi$points[,1:2]
    colnames(points) = c("x", "y")
    
    table = points
    colnames(table) = c("NMDS_1","NMDS_2")
    
    g = sample_data(ps)$Group %>% unique() %>% length()
    n = sample_data(ps)$Group%>% length()
    o = n/g
    ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    unif = phyloseq::distance(ps, method="bray")
    # adonis
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$Group = as.factor(map$Group)
    colbar = length(levels(map$Group))
    
    points = cbind(points, map[match(rownames(points), rownames(map)), ])
    points$ID = row.names(points)
  
    output$table8 <- renderDataTable({
      final8 <<- table
    })
    
    output$detectfig8 <- renderPlot({
      p = ggplot(points, aes(x=x, y=y, fill = Group)) +
        geom_point(alpha=.7, size=5, pch = 21) +
        labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
             y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
             title=title1) +
        stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      p8 <<- p + ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
      p8
    })
  })
  
  
  ## 下载图片写法
  output$downloadGraph8 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot8)
    },
    content = function(file){
      if(input$extPlot8 == 'pdf'){
        pdf(file)
      }else if(input$extPlot8 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      # 打印全局变量 p1 (ggplot对象)
      print(p8)
      dev.off()
    }
  )
  
  output$downloadTable8 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final8,file,row.names = T,quote = F)
    }
  )
  
  
  ##################### rb
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata9 <- reactive({
    infile <- input$file28
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata9 <- reactive({
    infile <- input$file29
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata9 <- reactive({
    infile <- input$file30
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset9, {
    values$file <- NULL
    output$file28 <- renderUI({
      fileInput("file28", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset9, {
    values$file <- NULL
    output$file29 <- renderUI({
      fileInput("file29", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset9, {
    values$file <- NULL
    output$file30 <- renderUI({
      fileInput("file30", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table9 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  
  observeEvent(input$start9, {
    
    otu = as.matrix(otudata9())
    sample = sampledata9()
    tax = as.matrix(taxdata9())

    output$table9 <- renderDataTable({
      colnames(tax) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  
      OTU = otu_table(otu, taxa_are_rows = TRUE)
      TAX = tax_table(tax)
      physeq_p = phyloseq(OTU, TAX)
      
      random_tree = rtree(ntaxa(physeq_p), 
                          rooted=TRUE, 
                          tip.label=taxa_names(physeq_p))
      sampledata = sample_data(sample)
      ps = phyloseq(OTU, TAX, sampledata, random_tree)
      
      j= as.character(input$class_levels)
      psdata <- ggClusterNet::tax_glom_wt(ps = ps,ranks = j)
      psdata = psdata%>% phyloseq::transform_sample_counts(function(x) {x/sum(x)})
      otu = phyloseq::otu_table(psdata)
      tax = phyloseq::tax_table(psdata)
      
      for (i in 1:dim(tax)[1]) {
        if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:10])) {
          tax[i,j] =tax[i,j]
        } else {
          tax[i,j]= "others"
        }
      }
      phyloseq::tax_table(psdata)= tax
      Taxonomies <- psdata %>%phyloseq::psmelt()
      Taxonomies$Abundance = Taxonomies$Abundance * 100
      colnames(Taxonomies) <- gsub(j,"aa",colnames(Taxonomies))
      data = c()
      i = 2
      for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
        a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
        b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
        c <- Taxonomies %>%
          dplyr::filter(Group == a)
        c$Abundance <- c$Abundance/b
        data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
        if (i == 1) {
          table = data
        }
        if (i != 1) {
          table = rbind(table,data)}}
      Taxonomies = table
      by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
      zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance), sd(Abundance))
      iris_groups<- dplyr::group_by(Taxonomies, aa)
      cc<- dplyr::summarise(iris_groups, sum(Abundance))
      head(cc)
      colnames(cc)= c("aa","allsum")
      cc<- dplyr::arrange(cc, desc(allsum))
      head(zhnagxu2)
      colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
      zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
      zhnagxu3 = zhnagxu2
      Taxonomies_x = plyr::ddply(zhnagxu3,"group", summarize,label_sd = cumsum(Abundance),label_y = cumsum(Abundance) - 0.5*Abundance)
      head( Taxonomies_x )
      Taxonomies_x = cbind(as.data.frame(zhnagxu3),as.data.frame(Taxonomies_x)[,-1])
      Taxonomies_x$label = Taxonomies_x$aa
      Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))

      Taxonomies_x <<- Taxonomies_x[,c(1,2,3,4)]
      colnames(Taxonomies_x)[1] = 'Item'
      final9 <<- Taxonomies_x
    })
    
    observeEvent(input$class_levels, {
      output$detectfig9 <- renderPlot({
        p9 <<- ggplot(Taxonomies_x , aes(x =  group, y = Abundance, fill = aa, order = aa)) +
          geom_bar(stat = "identity",width = 0.5,color = "black") +
          theme(axis.title.x = element_blank()) +
          theme(legend.text=element_text(size=6)) +
          scale_y_continuous(name = "Relative abundance (%)") +
          guides(fill = guide_legend(title = j)) +
          labs(x="",y="Relative abundance (%)",title= "") + scale_fill_hue()
        p9
      })
    })
  })
  
  
  ## 下载图片写法
  output$downloadGraph9 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot9)
    },
    content = function(file){
      if(input$extPlot9 == 'pdf'){
        pdf(file)
      }else if(input$extPlot9 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      # 打印全局变量 p1 (ggplot对象)
      print(p9)
      dev.off()
    }
  )
  
  output$downloadTable9 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final9,file,row.names = T,quote = F)
    }
  )
  
  
  ##################### td
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata10 <- reactive({
    infile <- input$file31
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata10 <- reactive({
    infile <- input$file32
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata10 <- reactive({
    infile <- input$file33
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset10, {
    values$file <- NULL
    output$file31 <- renderUI({
      fileInput("file31", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset10, {
    values$file <- NULL
    output$file32 <- renderUI({
      fileInput("file32", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset10, {
    values$file <- NULL
    output$file33 <- renderUI({
      fileInput("file33", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  
  observeEvent(input$start10, {
    
    otu = as.matrix(otudata10())
    sample = sampledata10()
    tax = as.matrix(taxdata10())
    
    if(length(unique(sample[,1])) != 3){
      text = 'Error: Number of Group types not equal to 3'
      
      output$error_text <- renderPrint({
        text
      })
    }
    else{
      colnames(tax) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
      
      OTU = otu_table(otu, taxa_are_rows = TRUE)
      TAX = tax_table(tax)
      physeq_p = phyloseq(OTU, TAX)
      
      random_tree = rtree(ntaxa(physeq_p), 
                          rooted=TRUE, 
                          tip.label=taxa_names(physeq_p))
      sampledata = sample_data(sample)
      ps = phyloseq(OTU, TAX, sampledata, random_tree)
      
      ps_rela = phyloseq::transform_sample_counts(ps,function(x) x / sum(x) )
      otu = ggClusterNet::vegan_otu(ps_rela) %>% as.data.frame()
      iris.split <- split(otu,as.factor(as.factor(phyloseq::sample_data(ps)$Group)))
      iris.apply <- lapply(iris.split,function(x)colSums(x[]))
      iris.combine <- do.call(rbind,iris.apply)
      ven2 = t(iris.combine) %>% as.data.frame()
      A <- combn(colnames(ven2),3)
      ven2$mean = rowMeans(ven2)
      tax = ggClusterNet::vegan_tax(ps)
      otutax = cbind(ven2,tax)
      
      j <<- which(colnames(otutax) == input$class_levels1)
      otutax[,j][otutax[,j] == ""] = "Unknown"
      otutax <<- otutax 
      i= 1
      x <<- A[1,i]
      y <<- A[2,i]
      z <<- A[3,i]
      
      output$detectfig10 <- renderPlot({
        p10 <<- ggtern::ggtern(data=otutax,
                               aes_string(x = x, 
                                          y = y, 
                                          z = z,
                                          color = colnames(otutax)[j],
                                          size ="mean")) +
          geom_point() + 
          theme_void()
        
        print(p10)
      })
      
    }
    
  })
  
  
  ## 下载图片写法
  output$downloadGraph10 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot10)
    },
    content = function(file){
      if(input$extPlot10 == 'pdf'){
        pdf(file)
      }else if(input$extPlot10 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      # 打印全局变量 p1 (ggplot对象)
      print(p10)
      dev.off()
    }
  )
  
  
  ##################### cl
  
  values <- reactiveValues(
    file = NULL
  )
  
  otudata11 <- reactive({
    infile <- input$file34
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata11 <- reactive({
    infile <- input$file35
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata11 <- reactive({
    infile <- input$file36
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset11, {
    values$file <- NULL
    output$file34 <- renderUI({
      fileInput("file34", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset11, {
    values$file <- NULL
    output$file35 <- renderUI({
      fileInput("file35", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset11, {
    values$file <- NULL
    output$file36 <- renderUI({
      fileInput("file36", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  
  observeEvent(input$start11, {
    
    otu = as.matrix(otudata11())
    sample = sampledata11()
    tax = as.matrix(taxdata11())
    
    colnames(tax) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    
    OTU = otu_table(otu, taxa_are_rows = TRUE)
    TAX = tax_table(tax)
    physeq_p = phyloseq(OTU, TAX)
    
    random_tree = rtree(ntaxa(physeq_p), 
                        rooted=TRUE, 
                        tip.label=taxa_names(physeq_p))
    sampledata = sample_data(sample)
    ps = phyloseq(OTU, TAX, sampledata, random_tree)
    
    ps_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
    ps_P <- ps_rela %>%
      ggClusterNet::tax_glom_wt(rank = as.character(input$class_levels2))
    ps_P
    otu_P = as.data.frame((ggClusterNet::vegan_otu(ps_P)))
    tax_P = as.data.frame(ggClusterNet::vegan_tax(ps_P))
    sub_design <<- as.data.frame(phyloseq::sample_data(ps_P))
    count2 =   otu_P
    iris.split <- split(count2,as.factor(sub_design$Group))
    iris.apply <- lapply(iris.split,function(x)colSums(x[]))
    iris.combine <- do.call(rbind,iris.apply)
    ven2 = t(iris.combine)
    lev = input$class_levels2
    
    Taxonomies <- ps %>%
      ggClusterNet::tax_glom_wt(rank = as.character(lev)) %>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x)} )%>%
      phyloseq::psmelt() %>%
      dplyr::arrange( Phylum)
    iris_groups<- dplyr::group_by(Taxonomies, Phylum)
    ps0_sum <- dplyr::summarise(iris_groups, mean(Abundance), sd(Abundance))
    ps0_sum[is.na(ps0_sum)] <- 0
    
    colnames(ps0_sum) = c("ID","mean","sd")
    ps0_sum <- dplyr::arrange(ps0_sum,desc(mean))
    ps0_sum$mean <- ps0_sum$mean *100
    ps0_sum <- as.data.frame(ps0_sum)
    
    top_P = ps0_sum$ID[1:10]
    otu_P = as.data.frame(t(otu_P))
    otu_tax = merge(ven2,tax_P,by = "row.names",all = F)
    
    otu_tax[,lev] = as.character(otu_tax[,lev])
    otu_tax[,lev][is.na(otu_tax[,lev])] = "others"
    for (i in 1:nrow(otu_tax)) {
      if(otu_tax[,lev] [i] %in% top_P){otu_tax[,lev] [i] = otu_tax[,lev] [i]}
      else if(!otu_tax[,lev] [i] %in% top_P){otu_tax[,lev] [i] = "others"}
    }
    otu_tax[,lev] = as.factor(otu_tax[,lev])
    head(otu_tax)
    otu_mean = otu_tax[as.character(unique(sub_design$Group))]
    head(otu_mean)
    row.names(otu_mean) = row.names(otu_tax)
    iris.split <- split(otu_mean,as.factor(otu_tax[,lev]))
    iris.apply <- lapply(iris.split,function(x)colSums(x[]))
    iris.combine <- do.call(rbind,iris.apply)
    mer_otu_mean <<- t(iris.combine)
    
    mi_sam <<- RColorBrewer::brewer.pal(9,"Set1")
    mi_tax <<- colorRampPalette(RColorBrewer::brewer.pal(9,"Set3"))(length(row.names(mer_otu_mean)))
    
    
    output$detectfig11 <- renderPlot({
      grid.col = NULL
      grid.col[as.character(unique(sub_design$Group))] = mi_sam
      
      grid.col[row.names(mer_otu_mean)] = mi_tax
      
      circlize::circos.par(gap.degree = c(rep(2, nrow(mer_otu_mean)-1), 10, rep(2, ncol(mer_otu_mean)-1), 10),
                           start.degree = 180)
      circlize::chordDiagram(mer_otu_mean,
                             directional = F,
                             diffHeight = 0.06,
                             grid.col = grid.col,
                             reduce = 0,
                             transparency = 0.5,
                             annotationTrack =c("grid", "axis"),
                             preAllocateTracks = 2
      )
      cex_size <<- input$size
      circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
        circlize::circos.text(circlize::CELL_META$xcenter, 
                              circlize::CELL_META$ylim[1], 
                              circlize::CELL_META$sector.index,
                              cex = cex_size,
                              facing = "clockwise", 
                              niceFacing = TRUE, 
                              adj = c(0, 0.5))}, 
        bg.border = NA)# here set bg.border to NA is important
      circlize::circos.clear()
    })
      
  })
  
  
  ## 下载图片写法
  output$downloadGraph11 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot11)
    },
    content = function(file){
      if(input$extPlot11 == 'pdf'){
        pdf(file)
      }else if(input$extPlot11 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      # 画出 circlize, 提供下载对象
      grid.col = NULL
      grid.col[as.character(unique(sub_design$Group))] = mi_sam
      
      grid.col[row.names(mer_otu_mean)] = mi_tax
      
      
      circlize::circos.par(gap.degree = c(rep(2, nrow(mer_otu_mean)-1), 10, 
                                          rep(2, ncol(mer_otu_mean)-1), 10),
                           start.degree = 180)
      circlize::chordDiagram(mer_otu_mean,
                             directional = F,
                             diffHeight = 0.06,
                             grid.col = grid.col,
                             reduce = 0,
                             transparency = 0.5,
                             annotationTrack =c("grid", "axis"),
                             preAllocateTracks = 2
      )
      
      circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
        circlize::circos.text(circlize::CELL_META$xcenter, 
                              circlize::CELL_META$ylim[1], 
                              circlize::CELL_META$sector.index,
                              cex = cex_size,
                              facing = "clockwise", 
                              niceFacing = TRUE, 
                              adj = c(0, 0.5))}, 
        bg.border = NA)# here set bg.border to NA is important
      circlize::circos.clear()
      
      dev.off()
    }
  )
  
  
  ##################### wc
  
  values <- reactiveValues(
    file = NULL
  )
  
  choose1 <- reactiveValues(
    list_choose1 = NULL
  )
  
  otudata12 <- reactive({
    infile <- input$file37
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata12 <- reactive({
    infile <- input$file38
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata12 <- reactive({
    infile <- input$file39
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset12, {
    values$file <- NULL
    output$file37 <- renderUI({
      fileInput("file37", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset12, {
    values$file <- NULL
    output$file38 <- renderUI({
      fileInput("file38", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset12, {
    values$file <- NULL
    output$file39 <- renderUI({
      fileInput("file39", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table12 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  # 填充复选框内容, 即先读取metadata的表头,当触发input_metadata()信号后执行下面
  observeEvent(sampledata12(), {

    sample = sampledata12()
    id.g = sample[,1] %>% unique() %>% as.character() %>% combn(2)
    
    for (i in 1:ncol(id.g)) {
      choose1$list_choose1[i] = paste(id.g[1,i], 'vs', id.g[2,i], sep = '_')
    }
    print(choose1$list_choose1)
  })
  
  
  output$res1 <- renderPrint({
    input$check1
  })
  
  
  observeEvent(!is.null(input$check1), {
    group_val <<- input$check1
  }, ignoreNULL = F)
  
  
  # 这里的input$reset记作点击的次数,若无则input$reset=0
  observeEvent(input$reset12, {
    values$file <- NULL
    # 点击过执行下面这句
    if (input$reset12[1] > 0) {
      choose1$list_choose1 = c("None")
    }
    output$file40 <- renderUI({
      if (is.null(choose1$list_choose1)) {
        selectInput("check1", "Choose", c("None"))
      }
      else{
        selectInput("check1", "Choose", choose1$list_choose1)
      }
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset12, {
    values$file <- NULL
    output$file41 <- renderUI({
      verbatimTextOutput(outputId = "res1")
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$start12, {
    
    otu = as.matrix(otudata12())
    sample = sampledata12()
    tax = as.matrix(taxdata12())
    
    output$table12 <- renderDataTable({
      cutoff = as.numeric(input$Pval)
      
      OTU = otu_table(otu, taxa_are_rows = TRUE)
      TAX = tax_table(tax)
      physeq_p = phyloseq(OTU, TAX)
      
      random_tree = rtree(ntaxa(physeq_p), 
                          rooted=TRUE, 
                          tip.label=taxa_names(physeq_p))
      sampledata = sample_data(sample)
      ps = phyloseq(OTU, TAX, sampledata, random_tree)
      
      
      map= sample_data(ps)
      
      id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
      check_li = c()
      li = c(strsplit(group_val, "_")[[1]][1],strsplit(group_val, "_")[[1]][3])
      for (id in 1:ncol(id.g)) {
        check_li[id] = identical(li,id.g[,id])
      }
      
      i = which(check_li == TRUE)
      ASV_table = ps %>%
        scale_micro(method = "sampling") %>%
        subset_samples.wt("Group", id.g[,i]) %>%
        vegan_otu() %>% 
        t() %>%
        as.data.frame()
     
      groupings <- sample[which(sample[,1] %in% id.g[,i]),]
      groupings$ID = row.names(groupings)
      
      pvals <- apply(ASV_table, 1, function(x) wilcox.test(x ~ groupings$Group, exact=F)$p.value)
      dat <<- pvals %>% as.data.frame()
      
      colnames(dat) = "P_val"
      tab.d12 = dat %>%
        rownames_to_column(var = "id") %>%
        dplyr::select(id,P_val) %>%
        dplyr::filter(P_val < cutoff) %>%
        dplyr::rename(
          OTU = id
        )  %>%
        dplyr::mutate(group = " wilcox.test.rare")

    
      final12 <<- data.frame(tab.d12)
    })
    
  })
  
  
  output$downloadTable12 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final12,file,row.names = T,quote = F)
    }
  )
  
  ##################### deseq2
  
  values <- reactiveValues(
    file = NULL
  )
  
  choose2 <- reactiveValues(
    list_choose2 = NULL
  )
  
  otudata13 <- reactive({
    infile <- input$file42
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata13 <- reactive({
    infile <- input$file43
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata13 <- reactive({
    infile <- input$file44
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset13, {
    values$file <- NULL
    output$file42 <- renderUI({
      fileInput("file42", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset13, {
    values$file <- NULL
    output$file43 <- renderUI({
      fileInput("file43", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset13, {
    values$file <- NULL
    output$file44 <- renderUI({
      fileInput("file44", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table13 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  # 填充复选框内容, 即先读取metadata的表头,当触发input_metadata()信号后执行下面
  observeEvent(sampledata13(), {
    
    sample = sampledata13()
    id.g = sample[,1] %>% unique() %>% as.character() %>% combn(2)
    
    for (i in 1:ncol(id.g)) {
      choose2$list_choose2[i] = paste(id.g[1,i], 'vs', id.g[2,i], sep = '_')
    }
    print(choose2$list_choose2)
  })
  
  
  observeEvent(!is.null(input$check2), {
    random_val <<- input$check2
  }, ignoreNULL = F)
  
  
  output$res2 <- renderPrint({
    input$check2
  })
  
  
  observeEvent(!is.null(input$check2), {
    group_val <<- input$check2
  }, ignoreNULL = F)
  
  
  # 这里的input$reset记作点击的次数,若无则input$reset=0
  observeEvent(input$reset13, {
    values$file <- NULL
    # 点击过执行下面这句
    if (input$reset13[1] > 0) {
      choose2$list_choose2 = c("None")
    }
    output$file45 <- renderUI({
      if (is.null(choose2$list_choose2)) {
        selectInput("check2", "Choose", c("None"))
      }
      else{
        selectInput("check2", "Choose", choose2$list_choose2)
      }
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset13, {
    values$file <- NULL
    output$file46 <- renderUI({
      verbatimTextOutput(outputId = "res2")
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$start13, {
    
    otu = as.matrix(otudata13())
    sample = sampledata13()
    tax = as.matrix(taxdata13())
    
    output$table13 <- renderDataTable({
      cutoff = as.numeric(input$Pval1)
      
      OTU = otu_table(otu, taxa_are_rows = TRUE)
      TAX = tax_table(tax)
      physeq_p = phyloseq(OTU, TAX)
      
      random_tree = rtree(ntaxa(physeq_p), 
                          rooted=TRUE, 
                          tip.label=taxa_names(physeq_p))
      sampledata = sample_data(sample)
      ps = phyloseq(OTU, TAX, sampledata, random_tree)
      
      
      map= sample_data(ps)
      
      id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
      check_li = c()
      li = c(strsplit(group_val, "_")[[1]][1],strsplit(group_val, "_")[[1]][3])
      for (id in 1:ncol(id.g)) {
        check_li[id] = identical(li,id.g[,id])
      }
      
      i = which(check_li == TRUE)
      
      ASV_table = ps %>%
        filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
        subset_samples.wt("Group", id.g[,i]) %>%
        vegan_otu() %>% t() %>%
        as.data.frame()
      
      groupings <- sample[which(sample[,1] %in% id.g[,i]),]
      groupings$ID = row.names(groupings)
      
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = ASV_table,
                                            colData=groupings,
                                            design = ~ Group)
      dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")
      res <- DESeq2::results(dds_res, tidy=TRUE, format="DataFrame")
      rownames(res) <- res$row
      res <- res[,-1]
      
      tab.d5 = res %>%
        rownames_to_column(var = "id") %>%
        dplyr::select(id,log2FoldChange,padj) %>%
        dplyr::filter(padj < cutoff) %>%
        dplyr::rename(
          OTU = id,
          log2FoldChange = log2FoldChange,
          p = padj
        )  %>%
        dplyr::mutate(group = "DESeq2")
      
      colnames(tab.d5)[3] = 'P_val'
      
      final13 <<- data.frame(tab.d5)
    })
    
  })
  
  
  output$downloadTable13 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final13,file,row.names = T,quote = F)
    }
  )
  
  
  ##################### edgeR
  
  values <- reactiveValues(
    file = NULL
  )
  
  choose3 <- reactiveValues(
    list_choose3 = NULL
  )
  
  otudata14 <- reactive({
    infile <- input$file47
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata14 <- reactive({
    infile <- input$file48
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata14 <- reactive({
    infile <- input$file49
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset14, {
    values$file <- NULL
    output$file47 <- renderUI({
      fileInput("file47", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset14, {
    values$file <- NULL
    output$file48 <- renderUI({
      fileInput("file48", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset14, {
    values$file <- NULL
    output$file49 <- renderUI({
      fileInput("file49", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table14 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  # 填充复选框内容, 即先读取metadata的表头,当触发input_metadata()信号后执行下面
  observeEvent(sampledata14(), {
    
    sample = sampledata14()
    id.g = sample[,1] %>% unique() %>% as.character() %>% combn(2)
    
    for (i in 1:ncol(id.g)) {
      choose3$list_choose3[i] = paste(id.g[1,i], 'vs', id.g[2,i], sep = '_')
    }
    print(choose3$list_choose3)
  })
  
  
  observeEvent(!is.null(input$check3), {
    random_val <<- input$check3
  }, ignoreNULL = F)
  
  
  output$res3 <- renderPrint({
    input$check3
  })
  
  
  observeEvent(!is.null(input$check3), {
    group_val <<- input$check3
  }, ignoreNULL = F)
  
  
  # 这里的input$reset记作点击的次数,若无则input$reset=0
  observeEvent(input$reset14, {
    values$file <- NULL
    # 点击过执行下面这句
    if (input$reset14[1] > 0) {
      choose3$list_choose3 = c("None")
    }
    output$file50 <- renderUI({
      if (is.null(choose3$list_choose3)) {
        selectInput("check3", "Choose", c("None"))
      }
      else{
        selectInput("check3", "Choose", choose3$list_choose3)
      }
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset14, {
    values$file <- NULL
    output$file51 <- renderUI({
      verbatimTextOutput(outputId = "res3")
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$start14, {
    
    otu = as.matrix(otudata14())
    sample = sampledata14()
    tax = as.matrix(taxdata14())
    
    output$table14 <- renderDataTable({
      cutoff = as.numeric(input$Pval2)
      
      OTU = otu_table(otu, taxa_are_rows = TRUE)
      TAX = tax_table(tax)
      physeq_p = phyloseq(OTU, TAX)
      
      random_tree = rtree(ntaxa(physeq_p), 
                          rooted=TRUE, 
                          tip.label=taxa_names(physeq_p))
      sampledata = sample_data(sample)
      ps = phyloseq(OTU, TAX, sampledata, random_tree)
      
      
      map= sample_data(ps)
      
      id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
      check_li = c()
      li = c(strsplit(group_val, "_")[[1]][1],strsplit(group_val, "_")[[1]][3])
      for (id in 1:ncol(id.g)) {
        check_li[id] = identical(li,id.g[,id])
      }
      
      i = which(check_li == TRUE)
      
      phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
        require("edgeR")
        require("phyloseq")
        # Enforce orientation.
        if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
        x = as(otu_table(physeq), "matrix")
        # Add one to protect against overflow, log(0) issues.
        x = x + 1
        # Check `group` argument
        if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
          # Assume that group was a sample variable name (must be categorical)
          group = get_variable(physeq, group)
        }
        # Define gene annotations (`genes`) as tax_table
        taxonomy = tax_table(physeq, errorIfNULL=FALSE)
        if( !is.null(taxonomy) ){
          taxonomy = data.frame(as(taxonomy, "matrix"))
        }
        # Now turn into a DGEList
        y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
        # Calculate the normalization factors
        z = calcNormFactors(y, method=method)
        # Check for division by zero inside `calcNormFactors`
        if( !all(is.finite(z$samples$norm.factors)) ){
          stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
        }
        # Estimate dispersions
        return(estimateTagwiseDisp(estimateCommonDisp(z)))
      }
      
      
      phylo <- ps %>%
        subset_samples.wt("Group", id.g[,i])
      
      test <- phyloseq_to_edgeR(physeq = phylo, group="Group")
      
      et = exactTest(test)
      
      tt = topTags(et, n=nrow(test$table), adjust.method="fdr", sort.by="PValue")
      res <- tt@.Data[[1]]
      head(res)
      
      tab.d4 = res %>%
        rownames_to_column(var = "id") %>%
        dplyr::select(id,logFC,FDR) %>%
        dplyr::filter(FDR < cutoff) %>%
        dplyr::rename(
          OTU = id,
          p = FDR
        )  %>%
        dplyr::mutate(group = "edgeR")
      
      colnames(tab.d4)[3] = 'P_val'
      
      final14 <<- data.frame(tab.d4)
    })
    
  })
  
  
  output$downloadTable14 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final14,file,row.names = T,quote = F)
    }
  )
  
  
  ################### network
  values <- reactiveValues(
    file = NULL
  )
  
  otudata15 <- reactive({
    infile <- input$file52
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T,
             row.names = 1
    )
  })
  
  taxdata15 <- reactive({
    infile <- input$file53
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  sampledata15 <- reactive({
    infile <- input$file54
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = T, 
             row.names = 1
    )
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset15, {
    values$file <- NULL
    output$file52 <- renderUI({
      fileInput("file52", "Step 1: Choose abundance matrix",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset15, {
    values$file <- NULL
    output$file53 <- renderUI({
      fileInput("file53", "Step 2: Choose taxonomy file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset15, {
    values$file <- NULL
    output$file54 <- renderUI({
      fileInput("file54", "Step 3: Choose metadata file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table15 <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  observeEvent(input$start15, {
    
    otu = as.matrix(otudata15())
    sample = sampledata15()
    tax = as.matrix(taxdata15())
    
    output$table15 <- renderDataTable({
      colnames(tax) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
      
      OTU = otu_table(otu, taxa_are_rows = TRUE)
      TAX = tax_table(tax)
      physeq_p = phyloseq(OTU, TAX)
      
      random_tree = rtree(ntaxa(physeq_p), 
                          rooted=TRUE, 
                          tip.label=taxa_names(physeq_p))
      sampledata = sample_data(sample)
      ps = phyloseq(OTU, TAX, sampledata, random_tree)
      
      result = corMicro(ps = ps,
                        N = 150,
                        method.scale = "TMM",
                        r.threshold=0.8,
                        p.threshold=0.05,
                        method = "spearman"
                        
                        
      )
      
      cor = result[[1]]
      ps_net = result[[3]]
      otu_table = ps_net %>%
        vegan_otu() %>%
        t() %>%
        as.data.frame()
      tax_table = ps_net %>%
        vegan_tax() %>%
        as.data.frame()
      netClu = data.frame(ID = row.names(tax_table),group = rep(1,length(row.names(tax_table)))[1:length(row.names(tax_table))] )
      netClu$group = as.factor(netClu$group)
      
      result2 = PolygonClusterG (cor = cor,nodeGroup =netClu )
      
      node = result2[[1]]
      nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
      edge = edgeBuild(cor = cor,node = node)
      edge_dat <<- edge
      
      class_type = input$class_levels4
      id = which(colnames(nodes) == class_type)
      nodes_dat = nodes[,c(1,2,id,ncol(nodes))]
      colnames(nodes_dat)[3] = c('taxonomy')
      nodes_dat <<- nodes_dat
      
      final15 <<- edge[,c(3,4,5,8)]
    })
    
    observeEvent(input$class_levels4, {
      output$detectfig15 <- renderPlot({
        p15 <<- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                        data = edge_dat, size = 0.5) +
          geom_point(aes(X1, X2, fill = taxonomy, size = mean),pch = 21, data = nodes_dat) +
          scale_colour_brewer(palette = "Set1") +
          scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
          theme(panel.background = element_blank()) +
          theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
          theme(legend.background = element_rect(colour = NA)) +
          theme(panel.background = element_rect(fill = "white",  colour = NA)) +
          theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
        p15
      })
    })
  })
  
  ## 下载图片写法
  output$downloadGraph15 <- downloadHandler(
    filename = function(){
      paste0(input$dataset, '.',input$extPlot15)
    },
    content = function(file){
      if(input$extPlot15 == 'pdf'){
        pdf(file)
      }else if(input$extPlot15 == 'png'){
        png(file)
      }else{
        jpeg(file)
      }
      # 打印全局变量 p1 (ggplot对象)
      print(p15)
      dev.off()
    }
  )
  
  output$downloadTable15 <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(final15,file,row.names = T,quote = F)
    }
  )
  
  
}

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50700)