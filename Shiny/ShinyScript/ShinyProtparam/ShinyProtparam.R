library(tidyverse)
library(dplyr)
library(Biostrings)
library(Peptides)  
library(shiny)
library(shinydashboard)
library(waiter)


header <- dashboardHeader(
  title = "Protein properties"
)

sidebar <- dashboardSidebar(
  fileInput("file1", "Choose FASTA File",
            multiple = TRUE,
            accept = c("text/fa",
                       "text/comma-separated-values,text/plain",
                       ".fa")),
  useWaiter(), # include dependencies
  actionButton("show", "Update"),
  # Horizontal line ----
  tags$hr(),
  selectInput("dataset", "Choose a dataset:",
              choices = c("properties")),
  downloadButton("downloadData", "Download")
)




body <- dashboardBody(
  title = "Examples of DataTables",
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        'input.dataset === "properties"',
        checkboxGroupInput("show_vars", "Columns in properties to show:",
                           c("ID","length","molecular_weight",
                             "instability","hydrophobicity",
                             "aliphatic","pI","charge"), 
                           selected = c("ID","length","molecular_weight",
                                        "instability","hydrophobicity","aliphatic",
                                        "pI","charge")
        )
      )
    ),
    mainPanel(
      tabsetPanel(
        id = 'dataset',
        tabPanel("properties", DT::dataTableOutput("mytable1"))
      )
    )
  )
)


ui <- dashboardPage(header, sidebar,body)

take_fa <- function(dat){
  inFile <- dat
  fa <- readAAStringSet(inFile$datapath)
  
  table = data.frame(fa) %>%
    rownames_to_column("name") %>%
    mutate("length" = Peptides::lengthpep(seq = fa)) %>% 
    mutate("molecular_weight" = mw(seq = fa)) %>%
    mutate("instability" = instaIndex(seq = fa)) %>%
    mutate("hydrophobicity" = hydrophobicity(seq = fa)) %>%
    mutate("aliphatic" = aIndex(seq = fa)) %>%     
    mutate("pI" = pI(seq = fa)) %>% 
    mutate("charge" = charge(seq = fa)) %>%
    as_tibble()
  
  table = as.data.frame(table[,c(1,3,4,5,6,7,8,9)])
  colnames(table) = c("ID","length","molecular_weight",
                      "instability","hydrophobicity","aliphatic",
                      "pI","charge")
  

  return(table)
}


server <- function(input, output) {
  
  observeEvent(input$show, {
    
    waiter_show( # show the waiter
      html = spin_fading_circles() # use a spinner
    )
    
    dataTable = take_fa(input$file1)
    
    
    output$mytable1 <- DT::renderDataTable({
      DT::datatable(dataTable[, input$show_vars, drop = FALSE],
                    extensions = 'Buttons',
                    options = list(                                                     
                      fixedColumns = TRUE,
                      autoWidth = TRUE,
                      ordering = TRUE,
                      dom = 'Bliftsp',
                      buttons = c('copy', 'csv', 'excel')
                    ),)
    })
    
    #output$dataset <- dataTable
    
    waiter_hide() # hide the waiter
    
    # put downloadevent into loading event
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$dataset, ".csv", sep = "")
      },
      content = function(file) {
        
        write.csv(dataTable, file, row.names = FALSE)
      }
    )
  })
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50326)