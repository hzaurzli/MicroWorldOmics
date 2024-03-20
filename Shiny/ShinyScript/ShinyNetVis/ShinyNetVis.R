library(shiny)
library(shinycssloaders)
library(colourpicker)
library(ggClusterNet)
library(tidyverse)
library(phyloseq)
library(ape)
library(sna)
library(igraph)
library(ggrepel)


dropdownButton <- function(label = "", status = c("default", "primary", "success", "info", "warning", "danger"), ..., width = NULL) {
  
  status <- match.arg(status)
  # dropdown button content
  html_ul <- list(
    class = "dropdown-menu",
    style = if (!is.null(width)) 
      paste0("width: ", validateCssUnit(width), ";"),
    lapply(X = list(...), FUN = tags$li, style = "margin-left: 10px; margin-right: 10px;")
  )
  # dropdown button apparence
  html_button <- list(
    class = paste0("btn btn-", status," dropdown-toggle"),
    type = "button", 
    `data-toggle` = "dropdown"
  )
  html_button <- c(html_button, list(label))
  html_button <- c(html_button, list(tags$span(class = "caret")))
  # final result
  tags$div(
    class = "dropdown",
    do.call(tags$button, html_button),
    do.call(tags$ul, html_ul),
    tags$script(
      "$('.dropdown-menu').click(function(e) {
      e.stopPropagation();
});")
  )
}

ui <- fluidPage(
  titlePanel("Network visualization"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(uiOutput('file1'),
               uiOutput('file2'),
               uiOutput('file3'),
               
               
               column(6,
                      div(style = "font-size: 11px", 
                          numericInput('N_val', 'Filter by abundance',  
                                       value = 400,min = 0,max = 500000)
                      ),
                      div(style = "font-size: 11px", 
                          selectInput('method_val', 'Correlation calculation',  
                                      choices = c("spearman", "kendall"), 
                                      selected = "spearman")
                      ),
                      div(style = "font-size: 11px",
                          numericInput('p_val','Significance cutoff',
                                       value = 0.05,min = 0,max = 1)
                      ),
                      div(style = "font-size: 11px",
                          numericInput('R_val','Pvalue calculate times',
                                       value = 10,min = 0,max = 500000)
                      ),
                      br(),
                      actionButton('reset', 'RESET')
               ),
               column(6,
                      div(style = "font-size: 11px",
                          numericInput('r_val', 'Correlation cutoff', 
                                       value = 0.6, min = 0, max = 1)
                      ),
                      div(style = "font-size: 11px",
                          uiOutput('file4'),    
                      ),
                      div(style = "font-size: 11px",
                          numericInput('maxnode_val','Max node size',
                                       value = 2,min = 0,max = 500000)
                      ),
                      div(style = "font-size: 11px",
                          numericInput('step_val','Random sampling times',
                                       value = 100,min = 0,max = 500000)
                      ),
                      br(),
                      actionButton('start', 'START')
               ),
               column(9,
                      hr(),
                      h5(strong('Draw network figure')),
                      div(style = "font-size: 11px",
                          colourInput("col", 'Boundary color', "grey40"),
                      ),
                      div(style = "font-size: 11px",
                          strong("Choose different levels colors (points colors)"),
                          br(),
                          br(),
                          dropdownButton(
                            label = "Choose colors", status = "default", width = 50,
                            uiOutput('file5')
                          )
                      ),
                      br(),
                      actionButton('draw', 'Draw'),
                      hr(),
                      h5(strong('Download all results')),
                      downloadButton("downloadData", "Download Table"),
                      br(),
                      br(),
                      downloadButton("downloadGraph", "Download Graph"),
                      br(),
                      br(),
                      radioButtons('extPlot', 'Plot output format',
                                   choices = c("PNG"='png', 
                                               'PDF'='pdf',
                                               'JPEG'='jpeg'), inline = T),
                      hr()
               ),
      ),
      
    ),
    mainPanel(align="center",
              h4("Network weight matrix and graph"),
              shinycssloaders::withSpinner(
                dataTableOutput("table")
              ),
              hr(),
              plotOutput("detectfig",width = "75%")
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  choose <- reactiveValues(
    list_choose = NULL
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
  
  
  observeEvent(taxdata(), {
    taxdata = taxdata()
    col_metadata = taxdata
    choose$list_choose <<- colnames(col_metadata)
    print(choose$list_choose)
  })
  
  
  observeEvent(input$reset, {
    values$file <- NULL
    # 点击过执行下面这句
    if (input$reset[1] > 0) {
      choose$list_choose = c("None")
    }
    output$file4 <- renderUI({
      if (is.null(choose$list_choose)) {
        selectInput('fill_val', 'Fill coulor of node',  
                    choices = c("None"))
      }
      else{
        selectInput('fill_val', 'Fill coulor of node',  
                    choices = choose$list_choose,
                    selected = choose$list_choose[2])
      }
    })
  }, ignoreNULL = F)
  
  
  
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
  
  
  output$detectfig <- renderPlot({
    return(NULL)
  })
  
  
  observeEvent(input$start, {
    
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
      
      
      N_df = input$N_val
      r_df = input$r_val
      p_df = input$p_val
      maxnode_df = input$maxnode_val
      method_df = input$method_val
      fill_df = input$fill_val
      step_df = input$step_val
      R_df = input$R_val	
      
      tab.r = network.pip(
        ps = ps,
        N = N_df,
        r.threshold = r_df,
        p.threshold = p_df,
        maxnode = maxnode_df,
        method = method_df,
        fill = fill_df,
        step = step_df,
        R = R_df,
        
        big = TRUE,
        select_layout = FALSE,
        layout_net = "model_maptree2",
        label = TRUE,
        lab = "elements",
        group = "Group",
        size = "igraph.degree",
        zipi = TRUE,
        ram.net = TRUE,
        clu_method = "cluster_fast_greedy",
        ncpus = 1
      )
      
      
      dat = tab.r[[2]]
      node = dat$net.cor.matrix$node
      edge = dat$net.cor.matrix$edge
      colnames(node)[which(colnames(node) == fill_df)] = "Class_levels"
      
      node_dat <<- node
      edge_dat <<- edge
      color_dat <<- unique(node$Class_levels)
      
      x <- colors()
      y <- sample(x,length(color_dat),replace = F)
      y_dat <<- y
      
      output$file5 <- renderUI({
        # 自动添加色块,有多少个分类色块就添加多少个色块
        lapply(1:length(color_dat), function(i) {
          colourpicker::colourInput(paste0("col_", i), paste0(color_dat[i]), y[i])
        })
      })
      
      
      final_tmp = data.frame()
      Type = c()
      for (j in 1:length(dat$net.cor.matrix$cortab)) {
        myAdjacencyMatrix = as.matrix(dat$net.cor.matrix$cortab[[j]])
        g  <- graph.adjacency(myAdjacencyMatrix,weighted=TRUE)
        final_df <- get.data.frame(g)
        final_tmp <- rbind(final_tmp,final_df)
        type_tmp <- rep(names(dat$net.cor.matrix$cortab)[j],
                        nrow(final_df))
        Type <- c(Type,type_tmp)
        print(nrow(final_df))
      }
      
      final <<- data.frame(final_tmp,Type)
      
    },options = list(pageLength = 10))
  })
  
  
  observeEvent(input$draw, {
    output$detectfig <- renderPlot({
      
      col_df = input$col
      
      fill_color = c()
      for (n in 1:length(color_dat)) {
        # 运行input$col_1,input$col_2这样的命令,取对应的色块
        cols <- paste0("input$col_", n)
        print(cols)
        cols <- eval(parse(text = cols))
        print(cols)
        if(is.null(cols)){
          fill_color <- append(fill_color, y_dat[n])
        } 
        else {
          fill_color <- append(fill_color, cols)
        }
        
      }
      print(fill_color)
      
      p <- ggplot() + 
        geom_segment(aes(x = X1, y = Y1, 
                         xend = X2, yend = Y2,
                         color = cor),
                     data = edge_dat, 
                     size = 0.03,alpha = 0.1) +
        geom_point(aes(X1, X2,
                       fill = Class_levels,
                       size = igraph.degree),
                   pch = 21, data = node_dat,color = col_df) +
        facet_wrap(.~label,scales="free_y",nrow = 2) +
        
        scale_colour_manual(values = c("#6D98B5","#D48852")) +
        scale_fill_manual(values = fill_color) +
        scale_size(range = c(0.8, 5)) +
        scale_x_continuous(breaks = NULL) +
        scale_y_continuous(breaks = NULL) +
        theme(panel.background = element_blank(),
              plot.title = element_text(hjust = 0.5)
        ) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank()
        ) +
        theme(legend.background = element_rect(colour = NA)) +
        theme(panel.background = element_rect(fill = "white",  
                                              colour = NA)) +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major = element_blank())
      p1 <<- p
      p
    })
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(df,file,row.names = T,quote = F)
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
      # 打印全局变量 p1 (ggplot对象)
      print(p1)
      dev.off()
    }
  )
  
}

# Run the application
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50905)
