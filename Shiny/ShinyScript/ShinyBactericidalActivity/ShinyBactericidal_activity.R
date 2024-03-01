library(shiny)

ui <- tagList(
  fluidPage(
    titlePanel("BAE (Bactericidal activity efficiency (Slope!))"),
    sidebarLayout(
      sidebarPanel(
        # uiOutput 做上传文件的 ui, 对应后面的 output$file1
        uiOutput('file1'),
        # 对应后面的 output$file2
        uiOutput('file2'),
        
        actionButton('reset', 'RESET'),
        hr(),
        downloadButton("downloadData", "Download"),
        hr(),
        h5('Developer:'),
        h6('Small runze (shiny app)'),
        br(),
        h5('Github: '),
        h6('https://github.com/hzaurzli (Small runze)'),
        br(),
        h5('Cition: A standardized approach for accurate quantification of murein hydrolase activity in high-throughput assays'),
        h6('')
      ),
      mainPanel(
        h4("BAE table"),
        br(),
        br(),
        shinycssloaders::withSpinner(
          dataTableOutput("table")
        )
      )
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  bae <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T)
  })
  
  concentration <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T)
  })
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 1: Choose bacterial quantity table (OD600, > 5 time points)",
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
      fileInput("file2", "Step 2: Choose drugs concentration table",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table <- renderDataTable({
    
    dat <- bae()
    concentration_table <- concentration()
    
    if(is.null(dat) | is.null(concentration_table)){
      warning("Please upload files!")
    } 
    else{
      result = data.frame()
      final = data.frame()
      
      find_peaks <- function (x, m = 3){
        shape <- diff(sign(diff(x, na.pad = FALSE)))
        pks <- sapply(which(shape < 0), FUN = function(i){
          z <- i - m + 1
          z <- ifelse(z > 0, z, 1)
          w <- i + m + 1
          w <- ifelse(w < length(x), w, length(x))
          if(all(x[c(z : i, (i + 2) : w)] < x[i + 1])) return(i + 1) else return(numeric(0))
        })
        pks <- unlist(pks)
        pks
      }
      
      
      for (j in 2:ncol(dat)) {
        for (i in 5:nrow(dat)) {
          dat_test = dat[1:i,]
          res = lm(dat_test[,j] ~ dat_test[,1])
          res_sum = summary(res)
          coef = res_sum$coefficients[2,1]
          p_val = res_sum$coefficients[2,4]
          R2 = res_sum$adj.r.squared
          result[i-4,1] = colnames(dat)[j]
          result[i-4,2] = coef
          result[i-4,3] = p_val
          result[i-4,4] = R2
          colnames(result) = c('group','coef','p_val','r2')
        }
        
        if(is.null(find_peaks(result[,4]))){
          id = which.max(result[,4])
          group = result[id,1]
          rate = result[id,2]
          p_val = result[id,3]
          R2 = result[id,4]
        }
        else{
          id = find_peaks(result[,4],m = 1)
          vec = c()
          num = 1
          for (m in id) {
            vec[num] = result[m,4]
            num = num + 1
          }
          max_val = max(vec)
          id = which(result[,4] == max_val)
          group = result[id,1]
          rate = result[id,2]
          p_val = result[id,3]
          R2 = result[id,4]
        }
        final[j-1,1] = group
        final[j-1,2] = rate
        final[j-1,3] = p_val
        final[j-1,4] = R2
        colnames(final) = c('Group','BAE_val','P_val','R2')
      }
      
      concentration_table = read.csv('D:/Documents/Desktop/concentration.csv',header = T)
      final[,2] = final[,2] / concentration_table[,2]
      finalTable <<- final
    }
  }, options = list(pageLength = 10))
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(linkList,file,row.names = T,quote = F)
    }
  )
  
  
}

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50632)

