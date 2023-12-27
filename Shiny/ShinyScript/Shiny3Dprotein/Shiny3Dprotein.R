library(shiny)
library(r3dmol)


ui <- fluidPage(
  titlePanel("Shiny 3D protein"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Step 1: Choose 3D structure (pdb) File",
                accept = c(
                  "pdb",
                  "PDB/pdb,text/plain",
                  ".pdb")
      ),
      hr(),
      ## 适用于少量文本
      textInput("Beta", "Beta sheet color", value = '#636efa'),
      ## 适用于少量文本
      textInput("alpha", "Alpha helix color", value = '#ff7f0e'),
      ## 适用于少量文本
      textInput("cartoon", "Cartoon color", value = '#00cc96'),

      hr(),
      h5('Developer:'),
      h6('Small runze (shiny app)'),
      br(),
      h5('Github: '),
      h6('https://github.com/hzaurzli (Small runze)')
    ),
    mainPanel(
      r3dmolOutput("detectfig",width = "100%",height = '500px')
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  # 上传文件后自动读取文件
  pdbread <-  reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(bio3d::read.pdb(system.file("extdata", 
                                         "1LNIA_decoy2_180.pdb",
                                         package="TMscoreAlign")))      
    }
    bio3d::read.pdb(infile$datapath)
  })

  output$detectfig <- renderR3dmol({

    pdb = pdbread()
    
    outfile <- tempfile(fileext = '.pdb')
    
    tryCatch(
      bio3d::write.pdb(pdb = pdb, xyz = pdb$xyz, file = outfile), 
      error=function(e){
        bio3d::write.pdb(pdb = pdb, xyz = pdb$xyz, file = outfile)
      },finnally = {
        bio3d::write.pdb(pdb = pdb, xyz = pdb$xyz, file = outfile)
      }
    )
    
    if(is.null(pdb)){
      warning("Please upload files!")
    } 
    else{
      r3dmol(
        viewer_spec = m_viewer_spec(
          cartoonQuality = 10, # 图形质量
          lowerZoomLimit = 50, # 缩放下限
          upperZoomLimit = 350 # 缩放上限
        )
      ) %>%
        # 添加模型
        m_add_model(data = outfile, format = "pdb") %>%  
        # 模型缩放到整体
        m_zoom_to() %>%
        # 设置 Cartoon 样式，并且颜色为绿色
        m_set_style(style = m_style_cartoon(color = input$cartoon)) %>% 
        # 设置 beta-折叠为蓝紫色
        m_set_style(sel = m_sel(ss = 's'),                 
                    style = m_style_cartoon(color = input$beta, arrows = TRUE)) %>% 
        # 设置 alpha-螺旋为橙色
        m_set_style(sel = m_sel(ss = 'h'),
                    style = m_style_cartoon(color = input$alpha)) %>%
        # 初始角度按Y轴旋转90度
        m_rotate(angle = 90, axis = 'y') %>%
        # 旋转动画
        m_spin()
    
    }
  })
}

# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50330)