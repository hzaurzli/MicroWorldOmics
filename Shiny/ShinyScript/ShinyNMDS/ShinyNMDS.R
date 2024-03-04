#### NMDS plot ####

library(shiny)
library(shinythemes)
library(tidyverse)
library(dplyr)
library(vegan)
library(tools)
library(shinycssloaders)
library(ggrepel)
library(RColorBrewer)

#### ui ####
ui <- fluidPage(
  theme = shinytheme("sandstone"),
  titlePanel("NMDS app"),
  sidebarLayout(
    sidebarPanel(
      fileInput('otu', 
                'Step1: Choose OTU table CSV',
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      tags$hr(style = "border-color: black;"),
      fileInput('samples', 
                'Step2: Choose sample list',
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      tags$hr(style = "border-color: black;"),
      sliderInput("percent_treshold", 
                  "Step3: Filter OTUs by percentage per sample", 
                  min = 0, 
                  max = 100, c(3), 
                  post = "%", 
                  step = 0.5),
      numericInput("no_samples", 
                   "Number of samples with percentage >= upper value", 
                   value = 1, 
                   min = 1, 
                   step = 1),
      helpText("Max number of samples:"),
      textOutput("sample_range"),
      tags$hr(style = "border-color: black;"),
      uiOutput("grouping_factor"),
      radioButtons("factor_select",
                   "Colours by", 
                   c("Factor" = "Factor", "Numeric" = "Values"), 
                   inline = T,
                   selected = "Factor"),
      uiOutput("label_factor"),
      downloadButton("downloadMultivar", 
                     "Table ready for NMDS"),
      tags$hr(style = "border-color: black;"),
      radioButtons("dissimilarity",
                   "Dissimilarity matrix",
                   c("Hellinger distance" = "hell", "Bray-Curtis" = "bray"),
                   inline = T,
                   selected = "hell"
      ),
      uiOutput("fitted"),
      tags$br(),
      downloadButton("downloadMultivarFinal", 
                     "Table with NMDS results for plotting"),
      tags$br(),
      tags$br(),
      downloadButton("downloadPlotFinal", 
                     "Download final plot as PDF"),
      tags$hr(style = "border-color: black;"),
      tags$br()
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("NMDS",
                 h5(textOutput("caption1")),
                 tableOutput("contents1") %>% withSpinner(type = getOption("spinner.type", default = 4)),
                 h5(textOutput("caption2")),
                 tableOutput("contents2") %>% withSpinner(type = getOption("spinner.type", default = 4)),
                 plotOutput("contents3") %>% withSpinner(type = getOption("spinner.type", default = 4))
        ),
        tabPanel("About",
                 h4("Plots for fast insight into community data"),
                 p("Tested on real data from the paper", a("Tláskal et al., 2017.", href = "https://academic.oup.com/femsec/article-abstract/93/12/fix157/4604780", target = "_blank"), "App is producing same results as metaMDS and envfit functions from the vegan package alone."),
                 p("Please note that apps hosted for free on shinyapps.io are limited to 1GB of memory. Therefore loading of larger OTU tables may take a while. If server disconnects after upload try to decrease size of excel file by e.g. deleting of singleton OTUs."),
                 p("packages:", 
                   p(a("tidyverse", href = "https://www.tidyverse.org/", target="_blank")), 
                   p(a("vegan", href = "https://cran.r-project.org/web/packages/vegan/index.html", target="_blank")), 
                   p(a("ggrepel", href = "https://github.com/slowkow/ggrepel", target="_blank")), 
                   p(a("shinycssloaders", href = "https://github.com/andrewsali/shinycssloaders", target="_blank")),
                   p(a("openxlsx", href = "https://github.com/awalker89/openxlsx", target="_blank")))
        )
      )
    )
  )
)


#### server ####  
server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)
  
  # OTUs
  dataset_otu <- reactive({
    
    infile = input$otu  
    
    if (is.null(infile))
      return(NULL)
    
    read.csv(infile$datapath)
    # readxl::read_excel(infile$datapath, sheet = input$sheet_otu, col_names = input$header_otu)
    
  })
  
  # samples
  dataset_samples <- reactive({
    infile = input$samples  
    
    if (is.null(infile))
      return(NULL)
    
    read.csv(infile$datapath, header = T)
    # readxl::read_excel(infile$datapath, sheet = input$sheet_samples, col_names = input$header_samples)
  })
  
  # total number of samples
  samples_count <- reactive({
    number = nrow(dataset_samples())
  })
  
  output$sample_range <- renderText({ 
    samples_count() 
  })
  
  # text and table
  output$caption1 <- renderText({
    "first 5 rows of OTU table are displayed"
  })
  
  output$contents1 <- renderTable({
    head(dataset_otu(), 5)
  })
  
  # text and table
  output$caption2 <- renderText({
    "first 5 rows of sample table are displayed"
  })
  output$contents2 <- renderTable({
    head(dataset_samples(), 5)
  })
  
  # ggplot grouping factor
  ggplot_factor <- reactive({
    factor <- dataset_samples()[,input$grouping_factor_input, drop = FALSE] 
    colnames(factor)<- "ggplot_factor" 
    factor 
  })
  
  # filter percentage
  filtered_titles <- reactive({
    otus_percent <- dataset_otu()
    tbl_df(otus_percent) %>% 
      gather(sample, per, (2:ncol(otus_percent))) %>% 
      group_by_at(c(1,2)) %>%  
      filter(per >= input$percent_treshold) %>%
      ungroup() %>% 
      group_by_at(c(1)) %>% 
      dplyr::summarise(treshold_count = n()) %>%
      filter(treshold_count >= input$no_samples) %>% 
      select(c(1))
  })
  
  # filtr multivar OTUs
  otus_multivar <- reactive({
    filtered_titles_list <- filtered_titles()
    otus_percent <- dataset_otu()
    tbl_df(otus_percent) %>% 
      gather(sample, per, (2:ncol(otus_percent))) %>%
      right_join(filtered_titles_list) %>% 
      spread(sample, per) 
  })
  
  # vegan matrix 
  otus_multivar_for_plot <- reactive({
    dataset_samples <- dataset_samples()
    otus_multivar <- otus_multivar()
    dataset_samples
    otus_multivar <- gather(otus_multivar, sample, perc, 2:ncol(otus_multivar)) %>% 
      spread(1, perc) %>% # this function orders NMDS input table according to sample names, must be consistent with dataset_samples() order
      arrange((match(sample, dataset_samples$sample_name))) %>% # ensures same order of NMDS input as dataset_samples() table but needs "sample_name" as column name 
      tibble::column_to_rownames(var = "sample")
  })
  
  # vegan matrix download
  output$downloadMultivar <- downloadHandler(
    filename = function() {
      paste("vegan_ready", input$otu, sep = "_")
    },
    content = function(file) {
      write.xlsx(otus_multivar_for_plot(), file, colNames = TRUE, rowNames = TRUE)
    })
  
  # NMDS without envfit
  mdsord <- reactive({
    if(input$dissimilarity == "hell") {
      set.seed(31)
      mdsord <- metaMDS(comm = decostand(otus_multivar_for_plot(), "hellinger"), distance = "euclidean", trace = FALSE, k = 2, trymax = 200, autotransform = FALSE)
    } else {
      set.seed(31)
      mdsord <- metaMDS(comm = otus_multivar_for_plot(), distance = "bray", trace = FALSE, k = 2, trymax = 200, autotransform = FALSE)
    }
    NMDS_data <- dataset_samples()
    ggplot_factor <- as.data.frame(ggplot_factor())
    NMDS_x <- mdsord$points[ ,1]  
    NMDS_y <- mdsord$points[ ,2]
    nmds_stress <- round(mdsord$stress, digits = 3)
    NMDS_data_final <- cbind(NMDS_data, NMDS_x, NMDS_y, ggplot_factor, data.frame(nmds_stress))
  })
  
  # NMDS envfit included, important are same parametres and set.seed
  mdsord_fitted <- reactive({
    if(input$dissimilarity == "hell") {
      set.seed(31)
      mdsord <- metaMDS(comm = decostand(otus_multivar_for_plot(), "hellinger"), distance = "euclidean", trace = FALSE, k = 2, trymax = 200, autotransform = FALSE)
    } else {
      set.seed(31)
      mdsord <- metaMDS(comm = otus_multivar_for_plot(), distance = "bray", trace = FALSE, k = 2, trymax = 200, autotransform = FALSE)
    }
    if(is.null(input$fitted_factors)){
    } else {
      set.seed(31)
      fitted_plot <- envfit(mdsord, fitted_df(), permutations = 999, arrow.mul = 1) 
      envfit_scores <- as.data.frame(scores(fitted_plot, display = "vectors"))
      envfit_scores <- cbind.data.frame(envfit_scores, env.variables = rownames(envfit_scores), stringsAsFactors = FALSE)
    }
  })
  
  # NMDS final points and variables table download
  output$downloadMultivarFinal <- downloadHandler(
    filename = function() {
      paste("final_positions", input$otu, sep = "_")
    },
    content = function(file) {
      l <- list("nmds_points" = mdsord(), "variables_score" = mdsord_fitted())
      write.xlsx(l, file)
    })
  
  output$grouping_factor <- renderUI({
    selectInput("grouping_factor_input", "Grouping factor",
                colnames(dataset_samples()),
                selected = NULL)
  })
  
  
  label_df <- reactive({
    label_variables <- dataset_samples()[,input$label_factor_input, drop = FALSE] 
    label_variables <- unlist(label_variables) # data.frame to atomic vector which is needed for geom_text 
  })
  
  
  fitted_df <- reactive({
    variables <- dataset_samples()[,input$fitted_factors, drop = FALSE] 
    variables 
  })
  
  
  # ggplot NMDS, points are from NMDS without envfit, arrows for env variables are from NMDS with envfit
  mdsord_final <- reactive({
    mdsord <- mdsord()
    mdsord_fitted <- mdsord_fitted()
    stress <- unique(mdsord$nmds_stress) # for stress value of NMDS
    # coloured by factor or value
    if (input$factor_select == "Factor") {
      mdsord$ggplot_factor <- as.factor(mdsord$ggplot_factor)
    } else {
      mdsord$ggplot_factor <- as.numeric(mdsord$ggplot_factor)
    }
    # if env variables are available
    if(is.null(input$fitted_factors)){
      ggplot(data = mdsord, aes(y = NMDS_y, x = NMDS_x)) +
        geom_point(aes(colour = ggplot_factor), show.legend = TRUE, size = 4.5) +
        {if (is.factor(mdsord$ggplot_factor)== TRUE) {
          scale_colour_brewer(palette = "Paired", type = "div") # color in the case of discrete values    
        } else {
          scale_color_viridis_c() # color in the case of continuous values
        }} +
        annotate("text", x = (0+max(mdsord$NMDS_x)), y = (0+min(mdsord$NMDS_y)), label = paste("stress\n", stress), size = 3.5) +
        theme_bw() + 
        ggtitle("NMDS plot")
    } else {
      ggplot(data = mdsord, aes(y = NMDS_y, x = NMDS_x)) +
        geom_point(aes(colour = ggplot_factor), show.legend = TRUE, size = 4.5) +
        {if (is.factor(mdsord$ggplot_factor)== TRUE) {
          scale_colour_brewer(palette = "Paired", type = "div") # color in the case of discrete values    
        } else {
          scale_color_viridis_c() # color in the case of continuous values
        }} +
        annotate("text", x = (0+max(mdsord$NMDS_x)), y = (0+min(mdsord$NMDS_y)), label = paste("stress\n", stress), size = 3.5) +
        theme_bw() +
        ggtitle("NMDS plot") +
        geom_segment(data = mdsord_fitted(),
                     aes(x = 0, xend = 1.2*NMDS1, y = 0, yend = 1.2*NMDS2),
                     arrow = arrow(length = unit(0.25, "cm")), colour = "#556b2f", size = 0.7) +
        geom_text(data = mdsord_fitted(),
                  aes(x = 1.2*NMDS1, y = 1.2*NMDS2, label = env.variables),
                  size = 6,
                  hjust = -0.3)
    }
  })
  
  output$contents3 <- renderPlot({
    mdsord_final()
  })
  
  # download NMDS
  output$downloadPlotFinal <- downloadHandler(
    filename = function() { 
      paste(tools::file_path_sans_ext(input$otu), '.pdf', sep='') 
    },
    content = function(file) {
      ggsave(file, plot = mdsord_final(), device = "pdf", dpi = 300, height = 210, width = 297, units = "mm")
    }
  )
  
}

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50635)

