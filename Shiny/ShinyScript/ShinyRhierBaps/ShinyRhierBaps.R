library(shiny)
library(shinycssloaders)
library(shinythemes)
library(rhierbaps)
library(ggtree)
library(phytools)
library(ape)


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
  navbarPage("BAPs clustering method", 
   theme = shinytheme("lumen"),
   # tabPanel 是导航栏换页
   tabPanel("BAPs cluster calculation", fluid = TRUE,
        tags$style(button_color_css),
        # Sidebar layout with a input and output definitions
        sidebarLayout(
          sidebarPanel(
            h3(strong('BAPs cluster calculation')),
            hr(),
            fluidRow(uiOutput('file1'),
                     
                     numericInput("max_depth", "Max depth",
                                  value = 2, min = 2, max = 1000),
                     
                     numericInput("n_pops", "Max number of populations",
                                  value = 20, min = 1, max = 1000),
                     radioButtons("inSelect", "Run until the algorithm converges to a local optimum",
                                  c("Yes", "No"), 
                                  selected = "Yes",
                                  inline = TRUE),
                     actionButton('reset', 'RESET'),
                     actionButton('start', 'START'),
                     hr(),
                     h4('Download all results'),
                     downloadButton("downloadTable", "Download table")
            ),
          ),
          mainPanel(
            titlePanel("BAPs cluster table"),
            withSpinner(dataTableOutput(outputId = "table")),
          )
        )
      ),
   
   tabPanel("BAPs cluster graph", fluid = TRUE,
            tags$style(button_color_css),
            # Sidebar layout with a input and output definitions
            sidebarLayout(
              sidebarPanel(
                h3(strong('BAPs cluster graph')),
                hr(),
                fluidRow(uiOutput('file2'),
                         uiOutput('file3'),
                         numericInput("max_depth1", "Max depth",
                                      value = 2, min = 2, max = 1000),
                         
                         numericInput("n_pops1", "Max number of populations",
                                      value = 20, min = 1, max = 1000),
                         radioButtons("inSelect1", "Run until the algorithm converges to a local optimum",
                                      c("Yes", "No"), 
                                      selected = "Yes",
                                      inline = TRUE),
                         radioButtons("inSelect2", "Ignore branch length",
                                      c("Yes", "No"), 
                                      selected = "No",
                                      inline = TRUE),
                         actionButton('reset1', 'RESET'),
                         actionButton('start1', 'START'),
                         br(),
                         br(),
                         uiOutput('file4'),
                         hr(),
                         h4('Download all results'),
                         radioButtons('extPlot', 'Plot output format',
                                      choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                         downloadButton("downloadGraph", "Download graph"),
                ),
              ),
              mainPanel(
                titlePanel("BAPs cluster phylogenetic tree"),
                withSpinner(plotOutput("detectfig")),
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
  
  
  snp_matrix <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    load_fasta(infile$datapath)
  })
  
  
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "First: Choose isolates fasta",
                accept = c(
                  "aln",
                  "fasta/aln,text/plain",
                  ".fa")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  observeEvent(input$start, {
    
    snp_matrix <- snp_matrix()
    
    if(is.null(snp_matrix)){
      warning("Please upload files!")
    } 
    else {
      
      max.depth = input$max_depth
      n.pops = input$n_pops
      
      output$table <- renderDataTable({
        if(input$inSelect == 'Yes') {
          hb.results <- hierBAPS(snp_matrix, max.depth = max.depth, 
                                 n.pops = n.pops, n.extra.rounds = Inf, 
                                 quiet = TRUE)
          
          final <<- hb.results$partition.df
        }
        else{
          hb.results <- hierBAPS(snp.matrix, 
                                 max.depth = max.depth, 
                                 n.pops = n.pops, 
                                 quiet = TRUE)
          
          final <<- hb.results$partition.df
        }
      },options = list(pageLength = 10))
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
  
  
  ###########################
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  snp_matrix1 <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    load_fasta(infile$datapath)
  })
  
  
  
  read_newick <- reactive({
    infile <- input$file3
    if (is.null(infile)){
      return(NULL)      
    }
    phytools::read.newick(infile$datapath)
  })
  

  observeEvent(input$reset1, {
    values$file <- NULL
    output$file2 <- renderUI({
      fileInput("file2", "First: Choose isolates fasta",
                accept = c(
                  "aln",
                  "fasta/aln,text/plain",
                  ".fa")
      )
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset1, {
    values$file <- NULL
    output$file3 <- renderUI({
      fileInput("file3", "First: Choose tree (newick)",
                accept = c(
                  "newick","nwk",
                  ".treefile,text/plain",
                  ".tree")
      )
    })
  }, ignoreNULL = F)
  
  
  output$detectfig <- renderPlot({
    NULL
  })
  

  output$file4 <- renderUI({
    # 自动添加色块,有多少个分类色块就添加多少个色块
    selectInput(inputId ="select_levels", 
                label ="Select Levels", 
                choices = c("NULL"))
  })
  
  
  observeEvent(input$start1, {
    
    snp_matrix <- snp_matrix1()
    tree_file <- read_newick()
    
    if(is.null(snp_matrix) | is.null(tree_file)){
      warning("Please upload files!")
    } 
    else {
      
      max.depth = input$max_depth1
      n.pops = input$n_pops1
      
      output$detectfig <- renderPlot({
        if(input$inSelect1 == 'Yes') {
          hb.results <- hierBAPS(snp_matrix, max.depth = max.depth, 
                                 n.pops = n.pops, n.extra.rounds = Inf, 
                                 quiet = TRUE)
          
          output$file4 <- renderUI({
            # 自动添加色块,有多少个分类色块就添加多少个色块
            selectInput(inputId ="select_levels", 
                        label ="Select Levels", 
                        choices = colnames(hb.results$partition.df)[-1])
          })
          
          
          if (input$inSelect2 == 'No') {
            colnam = input$select_levels
            df = hb.results$partition.df
            colnam = 'level 2'
            id = which(colnames(df) == colnam)
            colnames(df)[id] = 'BAPs cluster'
            
            gg <- ggtree(tree_file, layout = "circular")
            gg <- gg %<+% df
            gg <- gg + geom_tippoint(aes(color = factor(`BAPs cluster`)))
            g <<- gg
            gg
          }
          else {
            colnam = input$select_levels
            df = hb.results$partition.df
            colnam = 'level 2'
            id = which(colnames(df) == colnam)
            colnames(df)[id] = 'BAPs cluster'
            
            gg <- ggtree(tree_file, layout = "circular", 
                         branch.length = "none")
            gg <- gg %<+% df
            gg <- gg + geom_tippoint(aes(color = factor(`BAPs cluster`)))
            gg <- gg + theme(legend.position = "right")
            g <<- gg
            gg
          }
          
        }
        else{
          hb.results <- hierBAPS(snp.matrix, 
                                 max.depth = max.depth, 
                                 n.pops = n.pops, 
                                 quiet = TRUE)
          
          output$file4 <- renderUI({
            # 自动添加色块,有多少个分类色块就添加多少个色块
            selectInput(inputId ="select_levels", 
                        label ="Select Levels", 
                        choices = colnames(hb.results$partition.df)[-1])
          })
          
          
          if (input$inSelect2 == 'No') {
            colnam = input$select_levels
            df = hb.results$partition.df
            colnam = 'level 2'
            id = which(colnames(df) == colnam)
            colnames(df)[id] = 'BAPs cluster'

            gg <- ggtree(tree_file, layout = "circular")
            gg <- gg %<+% df
            gg <- gg + geom_tippoint(aes(color = factor(`BAPs cluster`)))
            g <<- gg
            gg
          }
          else {
            colnam = input$select_levels
            df = hb.results$partition.df
            colnam = 'level 2'
            id = which(colnames(df) == colnam)
            colnames(df)[id] = 'BAPs cluster'
            
            gg <- ggtree(tree_file, layout = "circular", branch.length = "none")
            gg <- gg %<+% df
            gg <- gg + geom_tippoint(aes(color = factor(`BAPs cluster`)))
            gg <- gg + theme(legend.position = "right")
            g <<- gg
            gg
          }
        }
      })
    }
  })
  
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
      print(g)
      dev.off()
    }
  )
  
}  

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50909)