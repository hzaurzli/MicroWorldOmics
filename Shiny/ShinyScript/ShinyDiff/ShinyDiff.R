library(shiny)
library(Maaslin2)

## 这个函数作 renderUI 和 uiOutput 无法在菜单栏翻页的时候进行渲染,只能在单页中使用
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


ui <- tagList(
      fluidPage(
      tags$h2("Maaslin2 differential abundance"),
      sidebarLayout(
        sidebarPanel(
          uiOutput('file1'),
          uiOutput('file2'),
          tags$h5(strong("Step3: Choose fix effects")),
          uiOutput('file3'),
          br(),
          uiOutput('file4'),
          tags$h5(strong("Step4: Choose random effects")),
          uiOutput('file5'),
          br(),
          uiOutput('file6'),
          actionButton('reset', 'RESET'),
          hr(),
          radioButtons("inSelect", "Standardize",
                       c("TRUE", "FALSE"),
                       selected = "FALSE",
                       inline = TRUE),
          actionButton('start', 'START',
                       style = "color: white; 
                           background-color: blue; 
                           position: relative;
                           height: 35px;
                           width: 70px;
                           text-align:center;
                           text-indent: -2px;
                           border-radius: 6px;
                           border-width: 2px"),
          downloadButton("downloadData", "Download result")
        ),
        mainPanel(
          h4("Differential abundance table"),
          br(),
          shinycssloaders::withSpinner(dataTableOutput("table"))
      )
    )
  )
)


server <- function(input, output, session) {
  values <- reactiveValues(
    file = NULL
  )
  
  choose <- reactiveValues(
    list_choose = NULL
  )
  
  input_data <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    as.character(infile$datapath)
  })
  
  
  input_metadata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    as.character(infile$datapath)
  })
  
  # 填充复选框内容, 即先读取metadata的表头,当触发input_metadata()信号后执行下面
  observeEvent(input_metadata(), {
    input_metadata = input_metadata()
    col_metadata = read.table(input_metadata, header=T, sep="\t")
    choose$list_choose <<- colnames(col_metadata)
    print(choose$list_choose)
  })
  
  
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 1: Choose OTU abundance matrices",
                accept=c("text/tsv",
                         "text/plain",
                         ".tsv")
      )
    })
  }, ignoreNULL = F)
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file2 <- renderUI({
      fileInput("file2", "Step 2: Choose metadata table",
                accept=c("text/tsv",
                         "text/plain",
                         ".tsv")
      )
    })
  }, ignoreNULL = F)

  
  # 这里的input$reset记作点击的次数,若无则input$reset=0
  observeEvent(input$reset, {
    values$file <- NULL
    # 点击过执行下面这句
    if (input$reset[1] > 0) {
      choose$list_choose = c("None")
    }
    output$file3 <- renderUI({
      if (is.null(choose$list_choose)) {
        dropdownButton(
          label = "Fix effects", status = "default", width = 80,
          checkboxGroupInput(inputId = "check1", 
                             label = "Choose", 
                             choices = c("None"))
        )
      }
      else{
        dropdownButton(
          label = "Fix effects", status = "default", width = 80,
          checkboxGroupInput(inputId = "check1", 
                             label = "Choose", 
                             choices = choose$list_choose)
        )
      }
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$file <- NULL
    output$file4 <- renderUI({
      verbatimTextOutput(outputId = "res1")
    })
  }, ignoreNULL = F)
  
  
  # 这里的input$reset记作点击的次数,若无则input$reset=0
  observeEvent(input$reset, {
    values$file <- NULL
    # 点击过执行下面这句
    if (input$reset[1] > 0) {
      choose$list_choose = c("None")
    }
    output$file5 <- renderUI({
      if (is.null(choose$list_choose)) {
        dropdownButton(
          label = "Random effects", status = "default", width = 80,
          checkboxGroupInput(inputId = "check2", 
                             label = "Choose", 
                             choices = c("None"))
        )
      }
      else{
        dropdownButton(
          label = "Random effects", status = "default", width = 80,
          checkboxGroupInput(inputId = "check2", 
                             label = "Choose", 
                             choices = choose$list_choose)
        )
      }
    })
  }, ignoreNULL = F)
  
  
  observeEvent(input$reset, {
    values$file <- NULL
    output$file6 <- renderUI({
      verbatimTextOutput(outputId = "res2")
    })
  }, ignoreNULL = F)
  
  # 获得多选框所选中的值
  observeEvent(!is.null(input$check1), {
    fix_val <<- input$check1
  }, ignoreNULL = F)
  
  # 获得多选框所选中的值
  observeEvent(!is.null(input$check2), {
    random_val <<- input$check2
  }, ignoreNULL = F)
  
  
  output$res1 <- renderPrint({
    input$check1
  })
  
  output$res2 <- renderPrint({
    input$check2
  })
  
  # 初始状态下加载空表,否则一直loading
  output$table <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  # 观察是否有点击start的信号,有则进行下面步骤
  observeEvent(input$start, {
    input_data = input_data()
    input_metadata = input_metadata()
    if(is.null(input_data) | is.null(input_metadata) | is.null(fix_val) | is.null(random_val)){
      warning("Please upload files!")
    } 
    else{
      output$table <- renderDataTable({
        fit_data <<- Maaslin2(input_data, input_metadata, 
                              paste(getwd(),'demo_output',sep = '/'),
                              fixed_effects = fix_val,
                              random_effects = random_val,
                              standardize = input$inSelect)
        
        unlink(paste(getwd(),'demo_output',sep = '/'), recursive=TRUE)
        final <<- fit_data$results
      }, options = list(pageLength = 7))
    }
  })
  
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
runApp(app, port = 50368)