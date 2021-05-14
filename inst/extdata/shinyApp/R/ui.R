library(shiny); library(shinydashboard); library(shinydashboardPlus); library(yaml); library(plotly); library(visNetwork); library(DT); library(shinyWidgets); library(shinyBS); library(shinyjs)
source('R/function.R')
lis.cfg <- yaml.load_file('config/config.yaml')
tit <- sub('^(title|width):', '', lis.cfg$title)


js <- "function openFullscreen(elem) {
  if (elem.requestFullscreen) {
    elem.requestFullscreen();
  } else if (elem.mozRequestFullScreen) { /* Firefox */
    elem.mozRequestFullScreen();
  } else if (elem.webkitRequestFullscreen) { /* Chrome, Safari and Opera */
    elem.webkitRequestFullscreen();
  } else if (elem.msRequestFullscreen) { /* IE/Edge */
    elem.msRequestFullscreen();
  }
}"


data_ui <- function(id, deg = FALSE) {
  ns <- NS(id)
  # if (deg == FALSE) tabPanel("Primary Visualization", value='primary',
  if (deg==FALSE) box(width = 12, title = "Data (replicates aggregated)", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      navbarPage('Parameters:',
      tabPanel("Basic",
      fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%', '1%', '10%', '1%', '10%'), '',
      radioButtons(inputId=ns('scale'), label='Scale by', choices=c('No', 'Row', 'Column'), selected='No', inline=TRUE), '', 
      radioButtons(inputId=ns('log'), label='Log/exp-transform', choices=c("No", "log2", "exp2"), selected='No', inline=TRUE),
      bsTooltip(id=ns('log'), title="No: original values in uploaded data are used.", placement = "bottom", trigger = "hover")
      ))), # tabPanel
      tabPanel("Filter",
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '30%', '1%', '24%', '1%', '24%'), '',
      numericInput(inputId=ns("A"), label="Threshold (A) to exceed", value=0), '',
      numericInput(inputId=ns("P"), label="Proportion (P) of samples with values >= A", value=0), '',
      numericInput(inputId=ns("CV1"), label="Min coefficient of variation (CV1)", value=-10^4), '', 
      numericInput(inputId=ns("CV2"), label="Max coefficient of variation (CV2)", value=10^4)
      )), actionButton(inputId=ns('fil.but'), label="Submit"), verbatimTextOutput(ns("fil.par")) 
      ),
      tabPanel("Re-order columns", column(12, uiOutput(ns('col.order'))))
      ), # navbarPage 
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dt")), ""))
     # )
  ) else if (deg == TRUE) {
    box(width = 12, title = "Data (with replicates)", closable = FALSE, solidHeader = TRUE, 
      collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '30%', '1%', '24%', '1%', '24%'), '',
      numericInput(inputId=ns("A"), label="Threshold (A) to exceed", value=0), '',
      numericInput(inputId=ns("P"), label="Proportion (P) of samples with values >= A", value=0), '',
      numericInput(inputId=ns("CV1"), label="Min coefficient of variation (CV1)", value=-10^4), '', 
      numericInput(inputId=ns("CV2"), label="Max coefficient of variation (CV2)", value=10^4)
      )), actionButton(inputId=ns('fil.but'), label="Submit"), verbatimTextOutput(ns("fil.par")),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dt")),  ""))
    )
  }
}

upload_ui <- function(id) {
   ns <- NS(id)
   tabPanel("Data sets", value='dataSets',
      box(width = 12, title = 'Data & aSVGs', closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      
      h4(strong("Step1: choose custom or default data sets")),
      fluidRow(splitLayout(cellWidths=c('1%', '30%'), '', 
      selectInput(ns("fileIn"), NULL, c('none', 'customData'), 'none')
      )), br(),
      fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '10%'), '', h4(strong("Step 2: upload custom data")), '', actionButton(ns("cusHelp"), "Help", icon = icon('question-circle')))),
      fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '18%', '1%', '25%', '1%', '25%'), '',
      fileInput(ns("geneInpath"), "2A: upload formatted data matrix", accept=c(".txt", ".csv"), multiple=FALSE), '',
      radioButtons(inputId=ns('dimName'), label='2B: is column or row gene?', choices=c("None", "Row", "Column"), selected='None', inline=TRUE), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "columns in the data matrix corresponds with "rows" in the targets file respectively.'),
      fileInput(ns("target"), "2C (optional): upload targets file for columns", accept=c(".txt", ".csv"), multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "rows" in the data matrix corresponds with "rows" in the row metadata file respectively.'),
      fileInput(ns("met"), "2D (optional): upload metadata file for rows", accept=c(".txt", ".csv"), multiple=FALSE))

      )),
      h4(strong("Step 3: upload custom aSVG(s)")),
      fluidRow(splitLayout(cellWidths=c('1%', '27%', '1%', '28%'), '',
      tags$div(class='tp', span(class='tpt', 'The data is matched with a single aSVG file.'),
      fileInput(ns("svgInpath1"), "3A: upload one aSVG file", accept=".svg", multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'The data is matched with multiple aSVG files (e.g. developmental stages).'),
      fileInput(ns("svgInpath2"), "3B (optional): upload multiple aSVG files", accept=".svg", multiple=TRUE))
      )),
      bsTooltip(id='svgInpath', title = "The data is matched with a single aSVG file.", placement = "bottom", trigger = "hover"),
      div(style = "font-size: 10px; padding: 0px 0px; margin:0%",
      fluidRow(splitLayout(cellWidths=c('1%', '27%', '1%', '28%', '1%', '28%'), '',
      downloadButton(ns("dld.sgl"), "Example1: data matched with a single aSVG"), '',
      downloadButton(ns("dld.mul"), "Example2: data matched with multiple aSVGs"), '',
      downloadButton(ns("dld.st"), "Example3:  spatiotemporal data-aSVG")
      ))), br(), 

      h4(strong("Additional files")), 
      fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '35%'), '',
      tags$div(class='tp', span(class='tpt', 'Upload a config file in ".yaml" format.'),
      fileInput(ns("config"), "Upload a config file (optional)", accept=".yaml", multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'The batched data sets will be listed under "Step 1".'),
      fileInput(ns("tar"), "Upload batched data, aSVGs in two separate tar files (optional)", accept=c(".tar"), multiple=TRUE))
      )),
      div(style = "font-size: 10px; padding: 0px 0px; margin:0%",
      fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '35%'), '',
      downloadButton(ns("dld.cfg"), "Example config file"), '', downloadButton(ns("dld.bat"), "Example data/aSVGs in batch")
      )))
      )
      
      )# tabPanel
  
}

shm_ui <- function(id, data.ui) {
  ns <- NS(id)
  tabPanel("Spatial Heatmap", value = 'shm2', icon = icon('image'),
    br(),
    # list(
    # width = ifelse(input$lgdTog %% 2 == 0, 9, 12), 
      # boxPad(color = NULL, title = NULL, solidHeader = FALSE, 
    # Append matrix heatmap, network with SHMs.    
    do.call(tabsetPanel, append(list(type = "pills", id='shmMhNet', selected = "shm1",
    # tabsetPanel(type = "pills", id=NULL, selected="shm1",
  
      tabPanel(title="Image", value='shm1',  
      navbarPage('Parameters:',
      tabPanel("Basic",
      fluidRow(splitLayout(cellWidths=c('0.5%', '9%', '0.5%', '6%', '0.5%', '10%', '0.5%', '10%', '0.5%', '10%', '0.5%', '10%', '0.5%', '13%', '0.5%', '13%'), '',  
      div(style = "margin-top: 23px;",
      actionButton(ns("fs"), "Full screen", onclick = "openFullscreen(document.getElementById('barSHM'))")), '',
      div(title = 'Number of columns for the subplots.', 
      numericInput(inputId=ns('col.n'), label='Columns', value=2, min=1, max=Inf, step=1, width=150)), '',
      div(title='Data column: by the column order in data matrix.',
      selectInput(inputId = ns("genCon"), label = "Display by", choices = c("Gene"="gene", "Condition"="con", "Data column"="none"), selected = '', multiple = FALSE, width = NULL)), '',
      numericInput(inputId=ns('scale.shm'), label='Scale plots', value='', min=0.1, max=Inf, step=0.1, width=150), '',
      numericInput(inputId=ns('scrollH'), label='Scrolling height', value=450, min=1, max=Inf, step=5, width=140), '',
      numericInput(inputId=ns('title.size'), label='Title size', value='', min=0.1, max=Inf, step=0.5, width=150), '',
      div(style = "margin-top: 23px;",
      dropdownButton(inputId=ns('dropdown'), label='Color key', circle=FALSE, icon= icon("sliders"), status='primary', inline=FALSE, width=250,
      fluidRow(splitLayout(cellWidths=c('1%', '60%', '35%'), '', textInput(ns("color"), "Color scheme", '', placeholder=paste0('Eg: ', ''), width=200),
      actionButton(ns("col.but"), "Confirm", icon=NULL, style = "margin-top: 24px;"))), 
      radioButtons(inputId=ns('cs.v'), label='Color key based on', choices=c("Selected rows", "All rows"), selected='', inline=TRUE)
      )), '', 
      actionButton(ns("lgdTog"), "Toggle legend plot", icon=icon('toggle-on'), class="btn-xs", title = "Click to toggle legend plot", style='margin-top: 24px')
      )), # fluidRow
      # bsPopover(id=ns('genCon'), title="Data column: by the column order in data matrix.", placement = "top", trigger = "hover"),
      textOutput(ns('h.w.c')), textOutput(ns('msg.col')),
      fluidRow(splitLayout(cellWidths=c('0.5%', '99.5%'), '',   checkboxGroupInput(inputId=ns("tis"), label="Select features to be transparent", choices='', selected='', inline=TRUE)))
      ), # tabPanel
      tabPanel("Value legend",
      column(9, offset=0, style='padding-left:50px; padding-right:80px; padding-top:0px; padding-bottom:5px',
      fluidRow(splitLayout(cellWidths=c('1%', '32%', '1%', '13%', '1%', '17%', '1%', '16%', '1%', '28%'), '', 
      actionButton(ns("val.lgd"), "Add/Remove", icon=icon("refresh"), style = "margin-top: 24px;"), '',  
      numericInput(inputId=ns('val.lgd.row'), label='Rows', value='', min=1, max=Inf, step=1, width=150), '',
      numericInput(inputId=ns('val.lgd.key'), label='Key size', value='', min=0.0001, max=1, step=0.01, width=150), '',
      numericInput(inputId=ns('val.lgd.text'), label='Text size', value='', min=0.0001, max=Inf, step=1, width=140), '',
      radioButtons(inputId=ns('val.lgd.feat'), label='Include features', choices=c('No', 'Yes'), selected='', inline=TRUE)
      ))) # column
      ), # tabPanel
      tabPanel("Shape outline",
      splitLayout(cellWidths=c('1%', '15%', '1%', '13%'), '', 
      selectInput(ns('line.color'), label='Line color', choices=c('grey70', 'black', 'red', 'green', 'blue'), selected=''), '', 
      numericInput(inputId=ns('line.size'), label='Line size', value='', min=0.05, max=Inf, step=0.05, width=150) 
      )), # tabPanel
     tabPanel("Download",

#     tags$div(title="Download the spatial heatmaps and legend plot.",
      # h1(strong("Download paramters:"), style = "font-size:20px;"),
      fluidRow(splitLayout(cellWidths=c('0.7%', '25%', '1%', '15%', '1%', '16%', '1%', '15%', '1%', '15%'), '',
      radioButtons(inputId=ns('ext'), label='File type', choices=c('NA', "png", "jpg", "pdf"), selected='', inline=TRUE), '', 
      numericInput(inputId=ns('res'), label='Resolution (dpi)', value='', min=10, max=Inf, step=10, width=150), '',
      radioButtons(inputId=ns('lgd.incld'), label='Include legend plot', choices=c('Yes', 'No'), selected='', inline=TRUE), '', 
      numericInput(inputId=ns('lgd.size'), label='Legend plot size', value='', min=-1, max=Inf, step=0.1, width=140), '',
      downloadButton(ns("dld.shm"), "Download", style = "margin-top: 24px;")
      )), # fluidRow
      bsTooltip(id=ns('ext'), title="Select a file type to download.", placement = "bottom", trigger = "hover")

      #fluidRow(splitLayout(cellWidths=c('18%', '1%', '25%', '1%', '18%') 
      # tags$div(title="Alegend plot.",
      # ), '', 
      #)) # fluidRow 
      ), # tabPanel 
     
     # navbarMenu("More", # Create dropdown menu.
     tabPanel("Relative size",
       tags$div(title="Only applicable in multiple aSVGs.",
       numericInput(inputId=ns('relaSize'), label='Relative sizes', value='', min=0.01, max=Inf, step=0.1, width=140))
      ), # tabPanel

      tabPanel(title="Re-match features", value='rematch',
        column(12, fluidRow(splitLayout(cellWidths=c('0.2%', "40%", '20%', "30%"), '',
          uiOutput(ns('svg'), style = 'margin-left:-5px'), '', 
          actionButton(ns("match"), "Confirm re-matching", icon=icon("refresh"), style="color: #fff; background-color:#3498DB;border-color: #2e6da4;margin-top: 24px;")
        ))), verbatimTextOutput(ns('msg.match')),
        column(12, uiOutput(ns('ft.match')))
      )
      #) # navbarMenu
      ), # navbarPage 
 
    verbatimTextOutput(ns('msg.shm')), uiOutput(ns('shm.ui')), data.ui
    ), # tabPanel 

      tabPanel(title='Animation', value='shm2', 
      fluidRow(splitLayout(cellWidths=c('1.5%', '15%', '1%', '16%', '1%', '12%', '1%', '12%'), '', 
      radioButtons(inputId=ns("ggly.but"), label="Show animation", choices=c("Yes", "No"), selected='', inline=TRUE), '',
      uiOutput(ns('tran.t')), '', uiOutput(ns('anm.scale')), '', uiOutput(ns('dld.anm.but'))
      )), textOutput(ns('tran')), uiOutput(ns('sld.fm')),
      # The input ids should be unique, so no legend plot parameters are added here.
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "61%", "30%"), "", plotOutput(ns("bar2")), htmlOutput(ns("ggly")), plotOutput(ns("lgd2"))))
      ),

      tabPanel(title='Video', value='shm3',
      navbarPage('Parameters:',
      tabPanel("Basic",
      fluidRow(splitLayout(cellWidths=c('1%', '8%', '1%', '10%', '1%', '13%', '1%', '10%', '1%', '8%', '2%', '18%'), '',
      numericInput(inputId=ns('vdo.key.row'), label='Key rows', value=2, min=1, max=Inf, step=1, width=270), '',
      numericInput(inputId=ns('vdo.key.size'), label='Key size', value=0.04, min=0.01, max=Inf, step=0.1, width=270), '',
      radioButtons(inputId=ns("vdo.val.lgd"), label="Key value", choices=c("Yes", "No"), selected='No', inline=TRUE), '', 
      radioButtons(inputId=ns("vdo.label"), label="Feature label", choices=c("Yes", "No"), selected='No', inline=TRUE), '',
      numericInput(inputId=ns('vdo.lab.size'), label='Label size', value=2, min=0, max=Inf, step=0.5, width=150), '',
      uiOutput(ns('video.dim'))
      )) # fluidRow
      ),
      tabPanel("More",
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '13%'), '', 
      numericInput(inputId=ns('vdo.itvl'), label='Transition time (s)', value=1, min=0.1, max=Inf, step=1, width=270), '',
      numericInput(inputId=ns('vdo.res'), label='Resolution (dpi)', value=400, min=1, max=1000, step=5, width=270)
      )) # fluidRow
      )
      ), # navbarPage
      textOutput(ns('tran.vdo')), htmlOutput(ns('ffm')), 
      radioButtons(inputId=ns("vdo.but"), label="Show/update video", choices=c("Yes", "No"), selected='No', inline=TRUE),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", uiOutput(ns('video')), ""))
      ) # tabPanel

      ), network_ui(ns('net')) )) # append, do.call

    #  ) # list

  )
}

network_ui <- function(id) {
  ns <- NS(id)
  list(
  tabPanel("Matrix Heatmap", id = 'mhm', icon = icon('border-all'),
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '10%', '1%', '20%', '1%', '7%', '1%', '20%', '1%', '10%'), '', 
      radioButtons(inputId=ns('measure'), label="Measure:", choices=c('correlation', 'distance'), selected='correlation', inline=TRUE), '', 
      radioButtons(inputId=ns("cor.abs"), label="Cor.absolute:", choices=c('No', 'Yes'), selected='No', inline=TRUE), '', 
      radioButtons(inputId=ns("thr"), label="Select by:", choices=c('proportion'='p', 'number'='n', 'value'='v'), selected='p', inline=TRUE), '',
      numericInput(inputId=ns('mhm.v'), label='Cutoff: ', value=0.2, min=-Inf, max=Inf, step=NA, width=NULL), '',
      radioButtons(inputId=ns("mat.scale"), label="Scale by:", choices=c("No", "Column", "Row"), selected='No', inline=TRUE), '',
      # radioButtons(inputId="mhm.but", label="Show plot:", choices=c("Yes", "No"), selected='No', inline=TRUE)
      actionButton(ns("mhm.but"), "Update", icon=icon("refresh"), style="color: #fff; background-color:purple;border-color: #2e6da4")
      )),
      plotlyOutput(ns("HMly"))
  ),
  tabPanel("Interactive Network", id = 'net', icon = icon('project-diagram'), 
      fluidRow( # If column widths are not integers, columns are vertically aligned.
      column(2, offset=0, style='padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      textInput(ns("color.net"), "Color scheme:", 'yellow,orange,red', placeholder='Eg: yellow,orange,red', width=150), actionButton(ns("col.but.net"), "Go", icon=icon("refresh"))
      ),
      column(1, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      selectInput(inputId=ns("net.type"), label="Network type:", choices=c('signed', 'unsigned', 'signed hybrid', 'distance'), selected='signed', width=100)
      ),
      column(1, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      numericInput(ns("min.size"), "Minmum module size:", value=15, min=15, max=5000, width=150)
      ),
      column(2, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      tags$span(style="color:brown;font-weight:NULL", selectInput(ns("gen.sel"), "Select a target gene:", c("None"), selected='None', width=150))
      ),
      column(2, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      selectInput(ns("ds"),"Module splitting sensitivity level:", 3:2, selected=3, width=170)
      ),
      column(1, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      selectInput(ns("adj.in"), "Adjacency threshold (the smaller, the more edges):", sort(seq(0, 1, 0.002), decreasing=TRUE), selected=1, width=200)
      ),
      column(2, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      tags$span(style="color:brown", numericInput(ns("max.edg"), "Maximun edges (too many edges may crash the app):", value=50, min=1, max=500, width=200)),
      htmlOutput(ns("edge"))
      ),
      column(1, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      actionButton(ns("cpt.nw"), "Update", icon=icon("refresh"), style="color: #fff; background-color:purple;border-color: #2e6da4")
      )
      ),
      fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput(ns("bar.net")), visNetworkOutput(ns("vis")), ""))
      ) # tabPanel("Interactive Network"
  )
}


deg_ui <- function(id) { 
  ns <- NS(id)
  tabPanel('Spatial Enrichment', value='deg',
    data_ui(ns("datDEG"), deg = TRUE),
    box(title='Spatial Enrichment', status="primary", solidHeader=TRUE, collapsible=TRUE, width=12,
      navbarPage(NULL, 
        tabPanel(title="Parameters",
          dataTableOutput(ns("dt.vs3")), br(),
          fluidRow(splitLayout(cellWidths=c('1%', '30%', '1%', '30%', '1%', '12%', '1%', '12%'), '',
            uiOutput(ns('ssg.sam')), '', uiOutput(ns('ssg.con')), '', 
            selectInput(ns("sam.con"), "Compare by", c('feature', 'factor', 'feature__factor'), 'feature'), '',
            uiOutput(ns('ssg.tar'))
           )),
          fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%', '1%', '20%', '1%', '10%', '1%', '10%'), '',
            selectInput(ns('ssg.meth'), label='Select methods', choices=c('edgeR', 'limma', 'DESeq2', 'ROKU', 'distinct'), selected=c('edgeR', 'limma'), multiple=TRUE), '',
            selectInput(ns("edg.lim.nor"), "Normalize-edgeR/limma", c("CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "none"), 'CNF-TMM'), '',
            selectInput(ns("rok.dis.nor"), "Normalize-ROKU/distinct", c("CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "ESF", "VST", "rlog", "none"), 'CNF-TMM'), '',
            numericInput(ns('ssg.fc'), 'Log2-fold change', 1, min=0, max=Inf, step=1), '',
            numericInput(ns('ssg.fdr'), 'FDR', 0.05, min=0, max=1, step=0.01)
          )),
          actionButton(ns("ssg.update"), "Update", icon=icon("refresh"), style="color: #fff; background-color:purple;border-color: #2e6da4") 
        ), # tabPanel
        tabPanel(title="Results in plots", value='w.meth',
          dataTableOutput(ns("dt.vs1")), br(),
          strong('Over-Expressed'), br(), column(12, column(8, plotOutput(ns("upset1"))), column(4, plotOutput(ns("ovl1")))),
          strong('Under-Expressed'), br(), column(12, column(8, plotOutput(ns("upset2"))), column(4, plotOutput(ns("ovl2")))),
        ),
        tabPanel(title="Link with spatial heatmap", value='dt.deg', 
          dataTableOutput(ns("dt.vs2")), br(),
          downloadButton(ns("dld.ssg.tab"), "Download"), dataTableOutput(ns("dt.deg"))
        ) 
      ) # navbarPage
    ) # box(title='Spatial Enrich'
  )
}


ui <- function(request) {
  dashboardPage(
 
  # includeCSS("style.css"),
  dashboardHeader(title = NULL, titleWidth = 0),
  dashboardSidebar(collapsed = TRUE, disable = TRUE, width = 0, sidebarMenu() ),
  controlbar = dashboardControlbar(id = "right.bar", collapsed = FALSE, overlay = FALSE, width = 51),
  dashboardBody(
   # tags$head(HTML('<title>spatialHeatmap</title>')),
   useShinyjs(), 
   includeCSS("www/style.css"),
   tags$head(tags$link(rel="stylesheet", type="text/css", href="style.css")),
   
   tags$head(tags$script(src = "javascript.js"), tags$script(HTML(js))),
   includeScript(path = "www/javascript.js"),
   tags$script(src="javascript.js"),
   tags$head(HTML("<script type='text/javascript' src='javascript.js'></script>")),
   HTML('<head>
         <link rel="stylesheet" type="text/css" href="style.css">
         <script type="text/javascript" src="www/javascript.js"></script>
         </head>
   '),

   # Place title on the right of dashboard header.
   tags$script(HTML('
     $(document).ready(function() {
       $("header").find("nav").append(\'<span class="myClass">spatialHeatmap</span>\');
     })
   ')),
    fluidRow(
 
      fluidRow(
      column(6, textInput(inputId='search', label='Search by gene IDs (e.g. ENSG00000000971,ENSG00000001617):', value='', placeholder='Muliple IDs must only be separated by space or comma.', width='100%')),
      column(1, actionButton('search.but', 'Submit'), align = "center", style = "margin-bottom: 2px;", style = "margin-top: 20px;")
      ),
      # tabsetPanel(type = "pills", id=NULL, selected="primary", data_ui('dat') ),
      # tabsetPanel(type = "pills", id = 'shm.sup', selected = "shm2", shm_ui('shmAll'))
      do.call(tabsetPanel, append(list(type = "pills", id = 'shm.sup', selected="shm2", upload_ui('upl'), shm_ui('shmAll', data_ui('dat')), deg_ui('deg')),
        list(tabPanel("About", value='about',
          box(width = 12, title = "", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
        'Coming soon!'
          )
        ))
    ))

    )
 ) # dashboardBody
)

}
