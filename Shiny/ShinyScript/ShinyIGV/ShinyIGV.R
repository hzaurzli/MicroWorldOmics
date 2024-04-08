library(shiny)
#BiocManager::install("igvShiny")
library(igvShiny)
library(htmlwidgets)
library(colourpicker)
library(GenomicAlignments)
library(dplyr)
# library(biomaRt)

# Set options
options(shiny.maxRequestSize=100*1024^3) #100G
options(encoding = "UTF-8")
options(stringsAsFactors = FALSE)
rm(list=ls())

#----------------------------------------------------------------------------------------------------
# 内置参考基因组
stock.genomes <- sort(get_css_genomes())
#stock.genomes
# [1] "ASM294v2"        "ASM985889v3"     "bosTau8"         "bosTau9"         "canFam3"        
# [6] "canFam4"         "canFam5"         "ce11"            "chm13v1.1"       "danRer10"       
# [11] "danRer11"        "dm3"             "dm6"             "dmel_r5.9"       "galGal6"        
# [16] "GCA_000182895.1" "GCA_003086295.2" "GCA_011100615.1" "GCF_001433935.1" "GCF_016699485.2"
# [21] "gorGor4"         "gorGor6"         "hg18"            "hg19"            "hg38"           
# [26] "hg38_1kg"        "hs1"             "macFas5"         "mm10"            "mm39"           
# [31] "mm9"             "NC_016856.1"     "panPan2"         "panTro4"         "panTro5"        
# [36] "panTro6"         "rn6"             "rn7"             "sacCer3"         "susScr11"       
# [41] "tair10" 
#----------------------------------------------------------------------------------------------------
ui = shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("genomeChooser", "Choose stock igv genome:", stock.genomes, selected = "hg38"),
      hr(),
      textInput("roi", label="Input gene name or location:", placeholder="gene or chrN:start-end"),
      actionButton("searchButton", "Search"),
      hr(),
      h5("Load local custom data:"),
      fileInput("bam","Choose bam file", accept = ".bam"),
      actionButton("addBamLocalFileButton", "Add BAM Track"),
      hr(),
      fileInput("bedgraph","Choose bedgraph file", accept = ".bedgraph"),
      colourInput("bedgraph_col", "Select colour", "blue"),
      actionButton("addBedGraphTrackButton", "Add BedGraph Track"),
      downloadButton(
        outputId = "downloadBedgraph",
        label = "DemoBedGraphFile",
        class = NULL
      ),
      hr(),
      fileInput("bed","Choose bed file", accept = ".bed"),
      colourInput("bed_col", "Select colour", "red"),
      actionButton("addBedTrackButton", "Add Bed Track"),
      downloadButton(
        outputId = "downloadBed",
        label = "DemoBedFile",
        class = NULL
      ),
      hr(),
      fileInput("gwas","Choose gwas data file", accept = c(".txt","tsv")),
      actionButton("addGwasTrackButton", "Add GWAS Track"),
      downloadButton(
        outputId = "downloadGwas",
        label = "DemoGWASFile",
        class = NULL
      ),
      hr(),
      h5("Remove user added tracks:"),
      actionButton("removeUserTracksButton", "Remove User Tracks"),
      hr(),
      h5("Get or clean current gene location:"),
      actionButton("getChromLocButton", "Get Region"),
      actionButton("clearChromLocButton", "Clear Region"),
      div(style = "background-color: white; width: 140px; height:20px; padding-left: 5px;
                   margin-top: 30px; border: 1px solid blue; font-size: 10px;",
          htmlOutput("chromLocDisplay")),
      hr(),
      width = 3
    ),
    mainPanel(
      igvShinyOutput('igvShiny_0'),
      width = 9
    )
  ) # sidebarLayout
))
#----------------------------------------------------------------------------------------------------
server = function(input, output, session) {
  
  observeEvent(input$genomeChooser, ignoreInit=TRUE, {
    newGenome <- input$genomeChooser
    #printf("new genome: %s", newGenome)
    genomeSpec <- parseAndValidateGenomeSpec(genomeName=newGenome, initialLocus="all")
    output$igvShiny_0 <- renderIgvShiny(
      igvShiny(genomeSpec)
    )
  })
  
  observeEvent(input$searchButton, {
    searchString = isolate(input$roi)
    #printf("--- search: %s", searchString)
    if(nchar(searchString) > 0)
      showGenomicRegion(session, id = "igvShiny_0", searchString)
  })
  
  observeEvent(input$addBedGraphTrackButton, {
    #printf("---- addBedgraphTrack")
    #printf("current working directory: %s", getwd())
    bedgraphFile <- input$bedgraph
    trackName <- basename(bedgraphFile$datapath)
    tbl.bedgraph <- read.table(bedgraphFile$datapath, header = T, sep = "\t", check.names = F)
    showGenomicRegion(session, id="igvShiny_0", "chr1:1,234,756-1,323,470")
    loadBedGraphTrack(session, id="igvShiny_0", trackName=trackName, tbl=tbl.bedgraph,
                      color=input$bedgraph_col, autoscale=TRUE)
  })
  
  observeEvent(input$addBedTrackButton, {
    #printf("---- addBedTrack")
    #printf("current working directory: %s", getwd())
    bedFile <- input$bed
    trackName <- basename(bedFile$datapath)
    tbl.bed <- read.table(bedFile$datapath, header = T, sep = "\t", check.names = F)
    showGenomicRegion(session, id="igvShiny_0", "chr1:1,234,756-1,323,470")
    loadBedTrack(session, id="igvShiny_0", trackName=trackName, tbl=tbl.bed, color=input$bed_col);
  })
  
  observeEvent(input$addGwasTrackButton, {
    #printf("---- addGWASTrack")
    #printf("current working directory: %s", getwd())
    gwasFile <- input$gwas
    trackName <- basename(gwasFile$datapath)
    tbl.gwas <- read.table(gwasFile$datapath, header = T, sep = "\t", check.names = F)
    showGenomicRegion(session, id="igvShiny_0", "chr19:45,248,108-45,564,645")
    loadGwasTrack(session, id="igvShiny_0", trackName=trackName, tbl=tbl.gwas, deleteTracksOfSameName=FALSE)
  })
  
  observeEvent(input$addBamLocalFileButton, {
    #printf("---- addBamLocalFileButton")
    bamFile <- input$bam$datapath
    trackName <- basename(bamFile)
    x <- readGAlignments(bamFile, param = Rsamtools::ScanBamParam(what="seq"))
    showGenomicRegion(session, id="igvShiny_0", "chr21:10,399,581-10,405,329")
    loadBamTrackFromLocalData(session, id = "igvShiny_0", trackName = trackName,
                              data = x, displayMode = "squished")
  })
  
  observeEvent(input$removeUserTracksButton, {
    #printf("---- removeUserTracks")
    removeUserAddedTracks(session, id = "igvShiny_0")
  })
  
  observeEvent(input$trackClick, {
    #printf("--- trackclick event")
    x <- input$trackClick
    print(x)
  })
  
  observeEvent(input[["igv-trackClick"]], {
    #printf("--- igv-trackClick event")
    x <- input[["igv-trackClick"]]
    attribute.name.positions <- grep("name", names(x))
    attribute.value.positions <- grep("value", names(x))
    attribute.names <- as.character(x)[attribute.name.positions]
    attribute.values <- as.character(x)[attribute.value.positions]
    tbl <- data.frame(name=attribute.names,
                      value=attribute.values,
                      stringsAsFactors=FALSE)
    dialogContent <- renderTable(tbl)
    html <- HTML(dialogContent())
    showModal(modalDialog(html, easyClose=TRUE))
  })
  
  observeEvent(input$getChromLocButton, {
    # printf("--- getChromLoc event")
    # sends message to igv.js in browser; currentGenomicRegion.<id> event sent back
    # see below for how that can be captured and displayed
    getGenomicRegion(session, id = "igvShiny_0")
  })
  
  observeEvent(input$clearChromLocButton, {
    output$chromLocDisplay <- renderText({" "})
  })
  
  observeEvent(input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]], {
    newLoc <- input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]]
    #printf("new chromLocString: %s", newLoc)
    output$chromLocDisplay <- renderText({newLoc})
  })
  
  ##ü TODO add fasta and bam files to inst
  genomes <- c("hg38", "hg19", "mm10", "tair10", "rhos")
  loci <- c("chr5:88,466,402-89,135,305",  "chr1:7,426,231-7,453,241", "MEF2C", "Mef2c",
            "1:7,432,931-7,440,395", "NC_007494.2:370,757-378,078",
            "chr1:6,575,383-8,304,088")
  
  output$igvShiny_0 <- renderIgvShiny({
    cat("--- starting renderIgvShiny\n");
    genomeOptions <- parseAndValidateGenomeSpec(genomeName="hg38", initialLocus="chr5:134,099,595-134,235,369")
    x <- igvShiny(genomeOptions,
                  displayMode="SQUISHED",
                  tracks=list()
    )
    cat("--- ending renderIgvShiny\n");
    return(x)
  })
  
  output$downloadBedgraph <- downloadHandler(
    filename = function() {
      paste("Example_bedgraph_file",".bedgraph", sep="")
    },
    content = function(file) {
      tbl.bedgraph <- read.table("./demo_data/demo_test.bedgraph",
                                 header = T,
                                 sep = "\t", check.names = F,
                                 stringsAsFactors = F)
      write.table(tbl.bedgraph,
                  file = file,
                  sep = "\t",
                  quote = F,
                  row.names = F)
    }
  )
  
  output$downloadBed <- downloadHandler(
    filename = function() {
      paste("Example_bed_file",".bed", sep="")
    },
    content = function(file) {
      tbl.bed <- read.table("./demo_data/demo_test.bed",
                            header = T,
                            sep = "\t", check.names = F,
                            stringsAsFactors = F)
      write.table(tbl.bed,
                  file = file,
                  sep = "\t",
                  quote = F,
                  row.names = F)
    }
  )
  
  output$downloadGwas <- downloadHandler(
    filename = function() {
      paste("Example_gwas_file",".txt", sep="")
    },
    content = function(file) {
      tbl.gwas <- read.table("./demo_data/demo_test_gwas.txt",
                             header = T,
                             sep = "\t", check.names = F,
                             stringsAsFactors = F)
      write.table(tbl.gwas,
                  file = file,
                  sep = "\t",
                  quote = F,
                  row.names = F)
    }
  )
  
} # server
#------------------------------------------------------------------------------------------------------------------------
# Run the application
app = shinyApp(ui = ui, server = server)
runApp(app, port = 55368)