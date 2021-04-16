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


data_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Primary Visualization", value='primary',
  box(width = 12, title = "Data (replicates aggregated)", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      navbarPage('Parameters:',
      tabPanel("Basic",
      fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%', '1%', '10%', '1%', '10%'), '',
      radioButtons(inputId=ns('scale'), label='Scale by', choices=c('No', 'Row', 'Column'), selected='No', inline=TRUE), '', 
      #log_exp_ui('data')
      radioButtons(inputId=ns('log'), label='Log/exp-transform', choices=c("No", "log2", "exp2"), selected='No', inline=TRUE)
      # dropdownButton(inputId='drdnFil', label='Filter', circle=FALSE, icon=NULL, status='primary')
      ))), # tabPanel
      tabPanel("Filter",
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '30%', '1%', '24%', '1%', '24%'), '',
      numericInput(inputId=ns("A"), label="Value (A) to exceed", value=0), '',
      numericInput(inputId=ns("P"), label="Proportion (P) of samples with values >= A", value=0), '',
      numericInput(inputId=ns("CV1"), label="Min coefficient of variation (CV1)", value=-10^4), '', 
      numericInput(inputId=ns("CV2"), label="Max coefficient of variation (CV2)", value=10^4)
      )), actionButton(inputId=ns('fil.but'), label="Submit"), verbatimTextOutput(ns("fil.par")) 
      ),
      tabPanel("Re-order columns",
      column(12, uiOutput(ns('col.order'))) 
      )
      ), # navbarPage 
      fluidRow(column(1, offset=0, style='padding-left:10px; padding-right:10px; padding-top:0px; padding-bottom:5px',
      bsTooltip(id='drdnFil', title="Filter the data by a proportion of a threthold value across all samples (pOverA) and coefficient of variance (CV).", placement = "bottom", trigger = "hover"),
      bsTooltip(id='P', title="Filter the data by a proportion of a threthold value across all samples (pOverA) and coefficient of variance (CV).", placement = "bottom", trigger = "hover")
      )
      #column(1, offset=0, style='padding-left:10px; padding-right:10px; padding-top:0px; padding-bottom:5px',
      #dropdownButton(inputId='drdn.scale', label='Transform', circle=FALSE, icon=NULL, status='primary', width=400
      #))
      # column(10, offset=0, style='padding-left:40px; padding-right:10px; padding-top:0px; padding-bottom:5px', uiOutput('col.order'))
      ),
      fluidRow(
      column(6, textInput(inputId=ns('search'), label='Search:', value='', placeholder='Muliple IDs must only be separated by space or comma.', width='100%')),
      column(1, actionButton(ns('search.but'), 'Submit'), align = "center", style = "margin-bottom: 2px;", style = "margin-top: 20px;")
      ),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dt")), ""))
      )
  )
}

upload_ui <- function(id) {
   ns <- NS(id)
   tabPanel("Upload", value='upload',
      box(width = 12, title = 'Upload data & aSVGs', closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      
      h4(strong("Step1: choose custom or default data sets")),
      fluidRow(splitLayout(cellWidths=c('1%', '30%'), '', 
      selectInput(ns("fileIn"), NULL, c('none', 'customData', 'customComputedData'), 'none')
      )), br(),
      h4(strong("Step 2: upload custom data")),
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

shm_ui <- function(id) {
  ns <- NS(id)
  tabPanel("Spatial Heatmap", value = 'shm2', icon = icon('image'),
    p(tags$head(tags$style(HTML('#lgdB{background-color:#d2d6de;}'))), actionButton(ns("lgdB"), "Toggle Legend Plot", icon=icon('toggle-on'), class="btn-xs", title = "Click to toggle legend plot"), hover = TRUE),
     uiOutput(ns('shm.ui')), uiOutput(ns('lgd.ui'))
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
  tabPanel('Spatial Enrich', value='deg',
    box(width = 12, title = "Data (with replicates)", closable = FALSE, solidHeader = TRUE, 
      collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dt.rep")), ""))
    ),
    box(title='Spatial Enrich', status="primary", solidHeader=TRUE, collapsible=TRUE, width=12,
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

shinyUI(dashboardPage(
 
  # includeCSS("style.css"),
  dashboardHeader(title = NULL, titleWidth = 0),
  dashboardSidebar(collapsed = TRUE, disable = TRUE, width = 0, sidebarMenu() ),
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
      tabsetPanel(type = "pills", id=NULL, selected="primary",
        upload_ui('upl'), data_ui('dat'), deg_ui('deg'),
        tabPanel("About", value='about',
          box(width = 12, title = "", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
        'Coming soon!'
          )
        )
      ),
      do.call(tabsetPanel, append(list(type = "pills", id = 'shm.sup', selected = "shm2", shm_ui('shmAll')), network_ui('net')))
    )
 ) # dashboardBody
)) # shinyUI(dashboardPage(
