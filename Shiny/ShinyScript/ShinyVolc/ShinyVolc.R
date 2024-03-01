library(shiny)
library(shinyBS)
library(dplyr)
library(ggplot2)
library(metricsgraphics)
library(RColorBrewer)


ui <- fluidPage(
  titlePanel("Volcano Plot"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Step 1: Choose data File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      hr(),
      radioButtons('sep', 'Separator',
                   c(Tab='\t',
                     Comma=','
                   ),
                   selected=','),
      checkboxInput("gene_names", "Show gene names", value = FALSE),
      hr(),
      h4("Axes"),
      sliderInput("lfcr", "Log2(Fold-Change) Range:", 
                  -10, 10, value = c(-2.5, 2.5), 
                  step=0.1, animate=FALSE),
      sliderInput("lo", "-Log10(P-Value):", 
                  0, 15, value = 4, step=0.05),
      hr(),
      h4("Cut-offs Selection"),
      sliderInput("hl", "P-Value Threshold:",
                  1, 6, value = 1.30, step=0.1),
      verbatimTextOutput('conversion'),
      sliderInput("vl", "log2(FC) Threshold:", 
                  0,2, value = 0.8, step=0.1),
      hr(),
      downloadButton('downloadData', 'Download Selected DE genes list'),
      hr(),
      downloadButton('downloadPlot', 'Download Volcano Plot (PDF)')  
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("ggPlot", plotOutput("ggplot", width = "720px", height = "720px")),
                  tabPanel("Cut-off Selected", dataTableOutput("tableOut"))
      )
    )
  )
)


server <- function(input, output, session) {    
    
  data <- reactive({
	 infile <- input$file1
	 if (is.null(infile)){
	   return(NULL)    
	 }
	 read.csv(infile$datapath,
			  quote='"',
			  stringsAsFactors=FALSE
	 )
	})
	  
  output$plot <- renderPlot({ 
    dat <- data();
    mask <-  with(dat, -log10(as.numeric(dat$P.Value))>input$hl & abs(dat$logFC)>input$vl)
    cols <- ifelse(mask, "red", "black")
    plot(as.numeric(dat$logFC), -log10(as.numeric(dat$P.Value)),
         xlim=input$lfcr, ylim=range(0,input$lo),
         xlab="log2(Fold-change)", ylab="-log10(P.Value)",
         cex = 0.35, pch = 16, col = cols)
    abline(h=input$hl, col="red")
    abline(v=-input$vl, col="blue")
    abline(v=input$vl, col="blue")
    tmp <- dat[-log10(as.numeric(dat$P.Value))>input$hl & abs(dat$logFC)>input$vl,]
    if(input$gene_names) try(text(tmp$logFC, -log10(tmp$P.Value), tmp$ID))
  })
  
  #     output$ggplot <- renderPlot({ 
  ggplotInput <- reactive({ 
    dat <- data();
    dat2 <- data.frame(x=as.numeric(dat$logFC), y=-log10(as.numeric(dat$P.Value)), ID=dat$ID)
    p <- ggplot(dat2, aes(x, y, label= ID)) + geom_point() +
      geom_vline(xintercept = input$vl, color = "blue") + #add vertical line
      geom_vline(xintercept = -input$vl, color = "blue") + #add vertical line
      geom_hline(yintercept = input$hl, color = "red") +  #add vertical line
      labs(x="log2(Fold-change)", y="-log10(P.Value)") + 
      scale_x_continuous("log2(Fold-change)", limits = input$lfcr) +
      scale_y_continuous("-log10(P.Value)", limits = range(0,input$lo)) + theme_bw()
    
    tmp <- dat[-log10(as.numeric(dat$P.Value))>input$hl & abs(dat$logFC)>input$vl,]
    
    q <- p + annotate("text", x=tmp$logFC, y=-log10(tmp$P.Value), 
                      label=tmp$ID, size=-log10(as.numeric(tmp$P.Value)), 
                      vjust=-0.1, hjust=-0.1)
    
    #         if(input$gene_names) print(q) else print(p)
    if(input$gene_names) q else p
  })
  
  
  output$ggplot <- renderPlot({
    print(ggplotInput())
  })
  
  output$conversion <- renderPrint(10^-(input$hl))
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste(gsub(".csv","", input$file1), '_selected.csv', sep='') 
    },
    content = function(file) {
      dat <-  data.frame(data());
      write.csv(dat[-log10(as.numeric(dat$P.Value))>input$hl & abs(dat$logFC)>input$vl
                    ,c("ID","logFC","P.Value")], 
                file, 
                row.names=FALSE,
                quote=FALSE)
    }
  )
  
  
  output$downloadPlot<- downloadHandler(
    filename <- function() {
      paste(gsub(".csv","", input$file1), Sys.Date(),'.pdf',sep='')
    },
    content = function(file) {
      pdf(file)
      print(ggplotInput())
      dev.off()
    }
  )
  
  output$tableOut <- renderDataTable({
    dat <-  data.frame(data())
    dat[-log10(as.numeric(dat$P.Value))>input$hl & abs(dat$logFC)>input$vl,c("ID","logFC","P.Value")]
  }) 
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50631)