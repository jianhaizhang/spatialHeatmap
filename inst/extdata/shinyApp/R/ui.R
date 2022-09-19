# The Shiny modules (e.g. search_ui) are temporarily placed in this file only for debugging purpose, and will be moved to independent files in the R folder after the App development is completed.

library(shiny); library(shinydashboard); library(shinydashboardPlus); library(yaml); library(plotly); library(visNetwork); library(DT); library(shinyWidgets); library(shinyBS); library(shinyjs)
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

# A module can only have a "ui" element without ther "server".
search_ui <- function(id, lab=label) {
  # Applying "ns" on an id: this id is prefixed with the module id in HTML output.
  ns <- NS(id) 
  fluidRow(splitLayout(cellWidths=c('1%', '13%', '1%', '85%'), '',
  radioButtons(inputId=ns('sch.mode'), label='Search mode', choices=c('Single', 'Multiple'), selected='Multiple', inline=TRUE), '',
  uiOutput(ns('sch.box'))
  ))
}


data_ui <- function(id, dim.ui=NULL, tailor.ui=NULL, deg=FALSE) {
  ns <- NS(id)
  # if (deg == FALSE) tabPanel("Primary Visualization", value='primary',
  if (deg==FALSE) box(width = 12, title = "Data (replicates aggregated)", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
    tabsetPanel(type="pills", id=ns('dtab.shm'), selected='dTabAll',
      tabPanel('Complete', value='dTabAll',
      navbarPage('Parameters:',
      tabPanel("Basic",
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '15%', '1%', '15%', '1%', '15%', '1%', '15%'), '',
      actionButton(ns("dat.all.but"), "Confirm selection", style='margin-top:24px'), '',
      selectInput(ns("normDat"), "Normalize", c('None'='none', "CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "ESF", "VST", "rlog"), selected='none'), '',
      selectInput(inputId=ns('log'), label='Log/exp-transform', choices=c("No", 'Log2'="log2", 'Exp2'="exp2"), selected='No'), '', 
      selectInput(inputId=ns('scaleDat'), label='Scale by', choices=c('No'='No', 'Row'='Row', 'Selected'='Selected', 'All'='All'), selected='Row'), '',
      numericInput(ns('page'), label='Page height', value=300, min=50, max=Inf, step=50, width=150)
      )),
      bsTooltip(id=ns('normDat'), title="CNF: calcNormFactors in edgeR. <br/> ESF: estimateSizeFactors in DESeq2. <br/> VST: varianceStabilizingTransformation in DESeq2. <br/> rlog: regularized log in DESeq2.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('scaleDat'), title="Row: scale each row independently. <br/> Selected: scale across all selected genes as a whole. <br/> All: scale across all genes as a whole.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('log'), title="No: original values in uploaded data are used. <br/> Log2: transform data to log2-scale. <br/> Exp2: transform data to the power of 2.", placement = "top", trigger = "hover")
      ), # tabPanel
      tabPanel("Filter",
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '30%', '1%', '24%', '1%', '24%'), '',
      numericInput(inputId=ns("A"), label="Threshold (A) to exceed", value=0), '',
      numericInput(inputId=ns("P"), label="Proportion (P) of samples with values >= A", value=0), '',
      numericInput(inputId=ns("CV1"), label="Min coefficient of variation (CV1)", value=-10^4), '', 
      numericInput(inputId=ns("CV2"), label="Max coefficient of variation (CV2)", value=10^4)
      )), actionButton(inputId=ns('fil.but'), label="Submit"), verbatimTextOutput(ns("fil.par")) 
      ),
      tabPanel("Threshold",
      fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%', '1%', '10%'), '', 
      textInput(ns("sig.max"), "Signal threshold (max)", '', placeholder=('Default to max signal.'), width=200), '',
      textInput(ns("sig.min"), "Signal threshold (min)", '', placeholder=('Default to min signal.'), width=200), '',
      actionButton(ns("sig.but"), "Confirm", icon=NULL, style = "margin-top: 24px;")
      )), br(), span(uiOutput(ns('msg.sig.thr')), style='color:red')  
      ),
      tabPanel("Re-order columns", column(12, uiOutput(ns('col.order'))))
      ), # navbarPage 
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dtAll")), ""))
     ), # tabPanel('Complete'
      tabPanel('Selected', value='dTabSel',
        actionButton(ns("tran.scale.but.sel"), "Transform/scale data", style='margin-top:10px'),
        fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dtSel")), ""))
      ),
      tabPanel('Selected Profile', value='dTabSelProf',
        actionButton(ns("tran.scale.but.prof"), "Transform/scale data", style='margin-top:10px'),
        fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", plotOutput(ns("selProf")), ""))
      ), 
      tabPanel('Single-cell metadata', value='dTabScell',
        column(12, id='colTailorUI', tailor.ui),
        column(12, id='colDimUI', dim.ui) 
      )
   ) # tabsetPanel(

  ) else if (deg == TRUE) {
    box(width = 12, title = "Data (with replicates)", closable = FALSE, solidHeader = TRUE, 
      collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '30%', '1%', '24%', '1%', '24%'), '',
      numericInput(inputId=ns("A"), label="Threshold (A) to exceed", value=0), '',
      numericInput(inputId=ns("P"), label="Proportion (P) of samples with values >= A", value=0), '',
      numericInput(inputId=ns("CV1"), label="Min coefficient of variation (CV1)", value=-10^4), '', 
      numericInput(inputId=ns("CV2"), label="Max coefficient of variation (CV2)", value=10^4)
      )), actionButton(inputId=ns('fil.but'), label="Submit"), verbatimTextOutput(ns("fil.par")),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dtRep")),  ""))
    )
  }
}

upload_ui <- function(id) {
   ns <- NS(id)
   tabPanel(title="Datasets", value='landing',
   navbarPage('', 
     tabPanel(title="Gallery", value='gallery',
      fluidRow(
        column(4, id='brainHum', style='text-align:center', uiOutput(ns('brain.hum'))),
        column(4, id='mouse', style='text-align:center', uiOutput(ns('mouse'))),
        column(4, id='chicken', style='text-align:center', uiOutput(ns('chicken')))
      ),
      br(), p(em('Arabidopsis thaliana'), style='font-size:18px'),
      fluidRow(
        column(4, id='organArab', style='text-align:center', uiOutput(ns('organ.arab'))),
        column(4, id='shootArab', style='text-align:center', uiOutput(ns('shoot.arab'))),
        column(4, id='rootArab', style='text-align:center', uiOutput(ns('root.arab')))
      ),
        column(4, id='stageArab', style='text-align:center', uiOutput(ns('stage.arab'))),
        column(4, style='width:100%'), column(4, style='width:100%'),
      fluidRow(
        column(4, id='clpRice', style='text-align:center', uiOutput(ns('clp.rice'))),
        column(4, ), column(4, )
      ),
     ), # tabPanel(title="Gallery",
     tabPanel(title="Data & aSVGs", value='datSVG',
      h4(strong("Step 1: choose custom or default data sets")),
      fluidRow(splitLayout(cellWidths=c('1%', '30%'), '', 
      selectInput(ns("fileIn"), label=NULL, choices=c('customBulkData'), selected='')
      )),
      uiOutput(ns('bulk.sce')), uiOutput(ns('svg.upl')),
      #fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '10%'), '', h4(strong("Step 2: upload custom data")), '', actionButton(ns("cusHelp"), "Help", icon = icon('question-circle')))),
      #fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '18%', '1%', '25%', '1%', '25%'), '',
      #fileInput(ns("geneInpath"), "2A: upload formatted data matrix", accept=c(".txt", ".csv"), multiple=FALSE), '',
      #radioButtons(inputId=ns('dimName'), label='2B: is column or row gene?', choices=c("None", "Row", "Column"), selected='None', inline=TRUE), '',
      #tags$div(class='tp', span(class='tpt', 'Ensure "columns in the data matrix corresponds with "rows" in the targets file respectively.'),
      #fileInput(ns("target"), "2C (optional): upload targets file for columns", accept=c(".txt", ".csv"), multiple=FALSE)), '',
      #tags$div(class='tp', span(class='tpt', 'Ensure "rows" in the data matrix corresponds with "rows" in the row metadata file respectively.'),
      #fileInput(ns("met"), "2D (optional): upload metadata file for rows", accept=c(".txt", ".csv"), multiple=FALSE))
      #)),
      #h4(strong('Single-Cell Data')), 
      #fileInput(ns("sglCell"), "", accept=c(".rds"), multiple=FALSE),
      #h5(strong("Step 3: upload custom aSVG(s)")),
      #fluidRow(splitLayout(cellWidths=c('1%', '27%', '1%', '28%'), '',
      #tags$div(class='tp', span(class='tpt', 'The data is matched with a single aSVG file.'),
      #fileInput(ns("svgInpath1"), "3A: upload one aSVG file", accept=".svg", multiple=FALSE)), '',
      #tags$div(class='tp', span(class='tpt', 'The data is matched with multiple aSVG files (e.g. developmental stages).'),
      #fileInput(ns("svgInpath2"), "3B (optional): upload multiple aSVG files", accept=".svg", multiple=TRUE))
      #)),
      bsTooltip(id='svgInpath', title = "The data is matched with a single aSVG file.", placement = "bottom", trigger = "hover"),
      div(style = "font-size: 10px; padding: 0px 0px; margin:0%",
      fluidRow(splitLayout(cellWidths=c('1%', '21%', '1%', '21%', '1%', '24%', '1%', '24%'), '',
      downloadButton(ns("dld.sgl"), "Example1: data & a single aSVG"), '',
      downloadButton(ns("dld.mul"), "Example2: data & multiple aSVGs"), '',
      downloadButton(ns("dld.st"), "Example3: spatiotemporal data & aSVG"), '',
      downloadButton(ns("dld.covis"), "Example4: co-visualization data & aSVG")
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
   ) # tabPanel(title="Data & aSVGs",
   ) # tabsetPanel(selected="gallery",
   ) # tabPanel(title="Datasets", 
  
}

# Match spatial features between data and aSVG.
match_ui <- function(id) { 
  ns <- NS(id)
  list(
  column(12, fluidRow(splitLayout(cellWidths=c('1%', '40%', '5%', '10%', '3%', '10%'), '',
    uiOutput(ns('svgs'), style = 'margin-left:5px'), '',
    uiOutput(ns('match.but'), style = 'margin-top:23px'), '',
    uiOutput(ns('match.reset'), style = 'margin-top:23px') 
    ))), verbatimTextOutput(ns('msg.match')),
  column(12, uiOutput(ns('ft.match')))
  )
}

shm_ui <- function(id, data.ui, search.ui) {
  ns <- NS(id)
  tabPanel("Spatial Heatmap", value='shmPanelAll', icon=NULL,
    div(style='margin-top:10px'),
    # list(
    # width = ifelse(input$lgdTog %% 2 == 0, 9, 12), 
      # boxPad(color = NULL, title = NULL, solidHeader = FALSE, 
    # Append matrix heatmap, network with SHMs.    
    do.call(tabsetPanel, append(list(type="pills", id=ns('shmMhNet'), selected="shm1",
    # tabsetPanel(type = "pills", id=NULL, selected="shm1", 
      tabPanel(title="Image", value='shm1',
      column(12, search.ui, style='z-index:5'),  
      navbarPage('Parameters:', id=ns('shmPar'),
      tabPanel("Basic", value='basic', 
      fluidRow(splitLayout(cellWidths=c('0.5%', '8%', '0.5%', '8%', '0.5%', '8%', '0.5%', '9%', '0.5%', '11%', '0.5%', '8%', '0.5%', '8%', '0.5%', '8%'), '',  
      actionButton(ns("fs"), "Full screen", onclick = "openFullscreen(document.getElementById('barSHM'))"), '',
      div(title = 'Number of columns for the subplots.',
      dropdownButton(inputId=ns('colDrop'), label='Columns', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width='100%',
        sliderInput(ns("col.n"), "", min=1, max=50, step=1, value=2, width='100%')
      )
      ), '',
      div(title='Data column: by the column order in data matrix.',
      dropdownButton(inputId=ns('disDrop'), label='Display by', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250,
        radioButtons(inputId=ns('genCon'), label='', choices=c("Gene"="gene", "Condition"="con", "Data column"="none"), selected='', inline=FALSE, width='100%')
      )), '',
      dropdownButton(inputId=ns('scaleDrop'), label='Scale plots', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250,
        sliderInput(ns("scale.shm"), "", min=0.1, max=10, step=0.1, value=1, width='100%')
      ), '',
      dropdownButton(inputId=ns('scroDrop'), label='Scrolling height', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250,
        sliderInput(ns("scrollH"), "", min=50, max=10000, step=50, value=450, width='100%')
      ), '',
      dropdownButton(inputId=ns('titDrop'), label='Title size', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250,
        sliderInput(ns("title.size"), "", min=0, max=100, step=0.5, value=12, width='100%')
      ), '',
      dropdownButton(inputId=ns('dropdown'), label='Color key', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250,
      fluidRow(splitLayout(cellWidths=c('1%', '60%', '35%'), '', textInput(ns("color"), "Color scheme", '', placeholder=paste0('Eg: ', ''), width=200),
      actionButton(ns("col.but"), "Confirm", icon=NULL, style = "margin-top: 24px;"))), 
      radioButtons(inputId=ns('cs.v'), label='Color key based on', choices=c("Selected rows", "All rows"), selected='', inline=TRUE)
      ), '', 
      dropdownButton(inputId=ns('togDrop'), label='Horizontal layout', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250,
        sliderInput(ns("togSld"), "Adjust horizontal layout", min=0, max=1, step=0.05, value=0.67, width='100%')
      )
      )), # fluidRow
      # bsPopover(id=ns('genCon'), title="Data column: by the column order in data matrix.", placement = "top", trigger = "hover"),
      textOutput(ns('h.w.c')), textOutput(ns('msg.col')), div(style='margin-top:10px') 
      ), # tabPanel
      tabPanel("Transparency",
        fluidRow(splitLayout(cellWidths=c('0.5%', '99.5%'), '',  checkboxGroupInput(inputId=ns("tis"), label="Select features to be transparent", choices='', selected='', inline=TRUE)))  
      ),
      tabPanel("Value legend",
      fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '10%', '1%', '10%', '1%', '10%', '1%', '10%'), '', 
      numericInput(inputId=ns('val.lgd.row'), label='Rows', value='', min=1, max=Inf, step=1, width=150), '',
      numericInput(inputId=ns('val.lgd.key'), label='Key size', value='', min=0.0001, max=1, step=0.01, width=150), '',
      numericInput(inputId=ns('val.lgd.text'), label='Text size', value='', min=0.0001, max=Inf, step=1, width=140), '',
      radioButtons(inputId=ns('val.lgd.feat'), label='Include features', choices=c('No', 'Yes'), selected='', inline=TRUE), '',
      actionButton(ns("val.lgd"), "Add/Remove", icon=icon("sync"), style = "margin-top: 24px;") 
      ))
      ), # tabPanel
      tabPanel("Shape outline",
      splitLayout(cellWidths=c('1%', '15%', '1%', '13%'), '', 
      selectInput(ns('line.color'), label='Line color', choices=c('grey70', 'black', 'red', 'green', 'blue'), selected=''), '', 
      numericInput(inputId=ns('line.size'), label='Line size', value='', min=0.05, max=Inf, step=0.05, width=150) 
      )), # tabPanel
     tabPanel("Download",

#     tags$div(title="Download the spatial heatmaps and legend plot.",
      # h1(strong("Download paramters:"), style = "font-size:20px;"),
      fluidRow(splitLayout(cellWidths=c('0.7%', '15%', '1%', '10%', '1%', '12%', '1%', '12%', '1%', '8%', '1%', '8%'), '',
      radioButtons(inputId=ns('ext'), label='File type', choices=c("jpg", "png", "pdf"), selected='jpg', inline=TRUE), '', 
      numericInput(inputId=ns('res'), label='Resolution (dpi)', value='', min=10, max=Inf, step=10, width=150), '',
      radioButtons(inputId=ns('lgd.incld'), label='Include legend plot', choices=c('Yes', 'No'), selected='', inline=TRUE), '', 
      numericInput(inputId=ns('lgd.size'), label='Legend plot size', value='', min=-1, max=Inf, step=0.1, width=140), '',
      actionButton(ns("dld.but"), "Confirm", icon=icon("sync"), style = "margin-top: 24px;"), '',
      # downloadButton(ns("dld.shm"), "Download", style = "margin-top: 24px;")
      uiOutput(ns('dldBut'))
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
        # column(12, fluidRow(splitLayout(cellWidths=c('1%', "40%", '10%', "30%"), '',
        #  uiOutput(ns('svg'), style = 'margin-left:-5px'), '', 
        #  actionButton(ns("match"), "Confirm re-matching", icon=icon("sync"))
        #))), verbatimTextOutput(ns('msg.match')),
        #column(12, uiOutput(ns('ft.match')))
        match_ui(ns('rematch'))
      ),
      tabPanel("Raster image",
       tags$div(title="Overlay templates of raster images with spatial heatmap images.",
         splitLayout(cellWidths=c('1%', '10%', '1%', '10%', '1%', '10%'), '', 
           selectInput(ns('raster'), label='Superimposing', choices=c('Yes', 'No'), selected='Yes'), '',
           selectInput(ns('coal'), label='Black-white', choices=c('Yes', 'No'), selected='No'), '',
           div(style='margin-top:25px',
           dropdownButton(inputId=ns('dpwAlpOver'), label='Alpha', circle=FALSE, icon=icon("", verify_fa = FALSE), status='primary', inline=FALSE, width=300,
           fluidRow(splitLayout(cellWidths=c('1%', '60%', '35%'), '',
           sliderInput(ns('alpOver'), "", min=0, max=1, step=0.05, value=1, width='100%'),
           actionButton(ns("alpOverBut"), "Confirm", icon=icon("sync"), style='margin-top:31px')))
           ))
         )
       )
      ), # tabPanel
      tabPanel("Single cell", value='scellTab', 
        tags$div(title="Single cell",
        fluidRow(splitLayout(cellWidths=c('1%', '7%', '1%', '15%', '1%', '10%'), '',  
        selectInput(ns('profile'), label='Profile', choices=c('No', 'Yes'), selected='No'), '', 
        selectInput(ns('dims'), label='Dimentionality reduction', choices=c('TSNE', 'PCA', 'UMAP'), selected='UMAP'), '', 
        uiOutput(ns('tarCell'))
        ))
       )
      ) # tabPanel
      #) # navbarMenu
      ), # navbarPage 
 
    verbatimTextOutput(ns('msgSHM')), uiOutput(ns('shm.ui')), data.ui
    ), # tabPanel 

      tabPanel(title='Interactive', value='interTab',
      navbarPage('', id=ns('interNav'),
      tabPanel('Plot', value='interPlot',
        fluidRow(splitLayout(cellWidths=c("1%", "13%", '5%', "80%"), '',
        actionButton(ns("ggly.but"), "Click to show/update", icon=icon("sync"), style="color:#fff; background-color:#499fe9;border-color:#2e6da4"), '',
        uiOutput(ns('sld.fm'))
        )),
        # The input ids should be unique, so no legend plot parameters are added here.
        fluidRow(splitLayout(cellWidths=c("1%", "7%", "61%", "30%"), "", plotOutput(ns("bar2")), htmlOutput(ns("ggly")), plotOutput(ns("lgd2"))))
      ),
      tabPanel('Parameters',
        fluidRow(splitLayout(cellWidths=c('1.5%', '16%', '1%', '12%', '1%', '12%'), '',
          numericInput(ns('t'), label='Transition time (s)', value=2, min=0.1, max=Inf, step=NA, width=270), '',
          numericInput(ns('scale.ly'), label='Scale plot', value=1, min=0.1, max=Inf, step=0.1, width=170), '',
          downloadButton(ns("dld.anm"), "Download", style = "margin-top: 24px;") 
         )), textOutput(ns('tran'))
      )) # navbarPage
      ),
      tabPanel(title='Video', value='vdoTab',
      navbarPage('', id=ns('vdoNav'),
      tabPanel('Video', value='video',
      actionButton(ns("vdo.but"), "Click to show/update", icon=icon("sync"), style="color:#fff; background-color:#499fe9;border-color:#2e6da4"),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", uiOutput(ns('video')), ""))
      ),
      tabPanel("Parameters",
      fluidRow(splitLayout(cellWidths=c('1%', '8%', '1%', '10%', '1%', '13%', '1%', '10%', '1%', '8%', '2%', '14%', '1%', '8%'), '',
      numericInput(inputId=ns('vdo.key.row'), label='Key rows', value=2, min=1, max=Inf, step=1, width=270), '',
      numericInput(inputId=ns('vdo.key.size'), label='Key size', value=0.04, min=0.01, max=Inf, step=0.1, width=270), '',
      radioButtons(inputId=ns("vdo.val.lgd"), label="Key value", choices=c("Yes", "No"), selected='No', inline=TRUE), '', 
      radioButtons(inputId=ns("vdo.label"), label="Feature label", choices=c("Yes", "No"), selected='No', inline=TRUE), '',
      numericInput(inputId=ns('vdo.lab.size'), label='Label size', value=2, min=0, max=Inf, step=0.5, width=150), '',
      selectInput(ns("vdo.dim"), label="Fixed dimension", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected='640x480', width=110), '',
      numericInput(inputId=ns('vdo.bar.width'), label='Color bar width', value=0.1, min=0.01, max=0.9, step=0.03, width=270)
      )), # fluidRow

      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '13%'), '', 
      numericInput(inputId=ns('vdo.itvl'), label='Transition time (s)', value=1, min=0.1, max=Inf, step=1, width=270), '',
      numericInput(inputId=ns('vdo.res'), label='Resolution (dpi)', value=400, min=1, max=1000, step=5, width=270)
      )), # fluidRow
      textOutput(ns('tran.vdo')), htmlOutput(ns('ffm')) 
      )
      ) # navbarPage
      ) # tabPanel
      ), network_ui(ns('net')) )) # append, do.call

    #  ) # list

  )
}

network_ui <- function(id) {
  ns <- NS(id)
  list(
  tabPanel(id='cluster', "Cluster", value='clus', icon=NULL,
      fluidRow(splitLayout(style='margin-top:3px;margin-bottom:3px', cellWidths=c('1%', '15%', '1%', '15%'), '', 
      dropdownButton(inputId=ns('dpbMea'), label='Similarity/Dissimilarity', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      radioButtons(inputId=ns('measure'), label="Measure", choices=c('Correlation', 'Distance'), selected='Correlation', inline=TRUE, width='100%'), 
      div(title='Only applicable when "Correlation" is selected.',
      radioButtons(inputId=ns("cor.abs"), label="Absolute correlation", choices=c('No', 'Yes'), selected='No', inline=TRUE, width='100%')
      )
      ), '',
      dropdownButton(inputId=ns('dpbThr'), label='Subset', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=400,
        div(title='Subset the most similar neighbors at a cutoff by top proportion, specific number, or similarity/dissimilairty value.',
        radioButtons(inputId=ns("thr"), label="By", choices=c('Proportion'='p', 'Number'='n', 'Similarity/Dissimilarity'='v'), selected='p', inline=TRUE, width='100%'),
        numericInput(inputId=ns('mhm.v'), label='Cutoff', value=0.2, min=-Inf, max=Inf, step=0.1, width='100%')
        )
      )
      )),
      bsTooltip(id=ns('dpbMea'), title='Measure to subset the nearest neighbors for the target gene (s) in spatial heatmaps.', placement = "top", trigger = "hover"),
      navbarPage('', id=ns('clusNav'),
      tabPanel(strong('Matrix Heatmap (MHM)'), value='mhmPlot', 
      actionButton(ns("mhm.but"), "Click to show/update", icon=icon("sync"), style="color: #fff; background-color:#499fe9;border-color:#2e6da4"), plotlyOutput(ns("HMly"))
      ),
      tabPanel('Parameter (MHM)', value='mhmPar',
      fluidRow(splitLayout(cellWidths=c('1%', '20%'), '', 
      radioButtons(inputId=ns("mat.scale"), label="Scale by", choices=c("No", "Column", "Row"), selected='No', inline=TRUE)
      # radioButtons(inputId="mhm.but", label="Show plot:", choices=c("Yes", "No"), selected='No', inline=TRUE)
      ))
      ), # tabPanel('Plot', value='mhmPar'
    tabPanel(strong("Interactive Network (NET)"), value='netPlot', icon=NULL, 
      actionButton(ns("cpt.nw"), "Click to show/update", icon=icon("sync"), style="color: #fff; background-color:#499fe9;border-color:#2e6da4"),
      fluidRow(splitLayout(cellWidths=c("1%", "5%", "92%", "2%"), "", plotOutput(ns("bar.net")), visNetworkOutput(ns("vis")), ""))
    ),
    tabPanel("Parameter (NET)", value='netPar', icon=NULL, 
      #fluidRow( # If column widths are not integers, columns are vertically aligned.
      fluidRow(splitLayout(style='margin-top:3px;margin-bottom:3px', cellWidths=c('0.5%', '9%', '0.5%', '10%', '0.5%', '14%', '0.5%', '15%', '0.5%', '11%'), '',
      dropdownButton(inputId=ns('dpwNetTar'), label='Target gene', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250,
        selectInput(ns("gen.sel"), "", c("None"), selected='None', width='100%')
      ), '',
      dropdownButton(inputId=ns('dpwNetType'), label='Network type', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      selectInput(inputId=ns("net.type"), label="", choices=c('signed', 'unsigned', 'signed hybrid', 'distance'), selected='signed', width='100%')
      ), '',
      dropdownButton(inputId=ns('dpwModSize'), label='Module size/splitting', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
        numericInput(ns("min.size"), "Min module size", value=15, min=15, max=5000, width='100%'),
        selectInput(ns("ds"), HTML("Module splitting sensitivity level <br/> (Larger value results in more <br/> modules with smaller sizes)"), 3:2, selected=3, width='100%')
      ), '',
      dropdownButton(inputId=ns('dpwNetAdj'), label='Adjacency/edges', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=300,
        selectInput(ns("adj.in"), HTML("Adjacency threshold <br/> (the smaller, the more edges)"), sort(seq(0, 1, 0.002), decreasing=TRUE), selected=1, width='100%'),
        numericInput(ns("max.edg"), HTML("Maximun edges <br/> (too many edges may crash the app)"), value=50, min=1, max=500, width='100%'), htmlOutput(ns("edge"))
      ), '',
      dropdownButton(inputId=ns('dpwColNet'), label='Color key', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=300,
      fluidRow(splitLayout(cellWidths=c('1%', '60%', '35%'), '',
      textInput(ns("color.net"), "Color scheme", 'yellow,orange,red', placeholder='Eg: yellow,orange,red', width='100%'),
      actionButton(ns("col.but.net"), "Confirm", icon=icon("sync"), style='margin-top:24px'))) # fluidRow(splitLayout
      )
      )) # fluidRow(splitLayout
     )
    ) # navbarPage('', id=ns('mhmNav')
  )
  ) # list
}


deg_ui <- function(id) { 
  ns <- NS(id)
  tabPanel('Spatial Enrichment', value='deg',
    data_ui(ns("datDEG"), deg=TRUE),
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
            selectInput(ns('ssg.meth'), label='Select methods', choices=c('edgeR', 'limma', 'DESeq2', 'distinct'), selected=c('edgeR', 'limma'), multiple=TRUE), '',
            selectInput(ns("edg.lim.nor"), "Normalize-edgeR/limma", c("CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "none"), 'CNF-TMM'), '',
            selectInput(ns("rok.dis.nor"), "Normalize-distinct", c("CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "ESF", "VST", "rlog", "none"), 'CNF-TMM'), '',
            numericInput(ns('ssg.fc'), 'Log2-fold change', 1, min=0, max=Inf, step=1), '',
            numericInput(ns('ssg.fdr'), 'FDR', 0.05, min=0, max=1, step=0.01)
          )),
          actionButton(ns("ssg.update"), "Update", icon=icon("sync"), style="color: #fff; background-color:purple;border-color: #2e6da4") 
        ), # tabPanel
        tabPanel(title="Results in plots", value='w.meth',
          dataTableOutput(ns("dt.vs1")), br(),
          strong('Over-Expressed'), br(), column(12, column(8, plotOutput(ns("upset1"))), column(4, plotOutput(ns("ovl1")))),
          strong('Under-Expressed'), br(), column(12, column(8, plotOutput(ns("upset2"))), column(4, plotOutput(ns("ovl2")))),
        ),
        tabPanel(title="Link with spatial heatmap", value='dt.deg', 
          dataTableOutput(ns("dt.vs2")), br(),
          column(12, search_ui(ns('deg')), style='z-index:5'),  
          column(12, dataTableOutput(ns("dt.deg"))),
          downloadButton(ns("dld.ssg.tab"), "Download")
        ) 
      ) # navbarPage
    ) # box(title='Spatial Enrich'
  )
}

# Dimension reduction in single-cell data.
dim_ui <- function(id) { 
  ns <- NS(id); uiOutput(ns('dim.ui'))
  #list(
  #  column(12,
  #  fluidRow(splitLayout(cellWidths=c('1%', '48%', '1%', '48%'), '',
  #    plotOutput(ns('tsne')), '', plotOutput(ns('pca'))
  #  )), dataTableOutput(ns('scell.cdat'))
  #)
  #)
}

# Dimension reduction in single-cell data.
tailor_match_ui <- function(id) { 
  ns <- NS(id)
  col1 <- column(6, id='dimCellBut', 
  fluidRow(splitLayout(cellWidths=c('1%', '50%', '1%', '20%', '1%', '5%'), '',
    selectInput(ns('dimCell'), label='Cells before co-clusterings', choices=c("UMAP", "PCA", "TSNE")), '', 
    actionButton(ns("coclusPlotBut"), label='Covisualizing results', style='margin-bottom:-60px'), ''
  ))
  );
  col2 <- column(6, id='selBlkButCan',
  fluidRow(splitLayout(cellWidths=c('1%', '40%', '1%', '25%', '1%', '12%'), '', uiOutput(ns('selBlk')), '',
  actionButton(ns("selBlkBut"), label='Final confirmation', style='margin-bottom:-60px'), '',
  actionButton(ns("selBlkCancel"), label='Reset', style='margin-bottom:-60px'), '',
  actionButton(ns("tailorHelp"), "Help", style='margin-bottom:-60px', icon = icon('question-circle'))
  ))
  ) 
  list( # Inside the list: "," is the separator. Outside the list ";" is the separator.
  col1, col2,
  # if (hide==TRUE) hidden(col1), if (hide==TRUE) hidden(col2), 
  uiOutput(ns('dim.ui'))
  )
}

covis_man_ui <- function(id) { 
  ns <- NS(id)
    tabsetPanel(type = "pills", id=ns('tabSetCell'), selected="datCell",
      tabPanel(title="Parameters", value='parMan', br(),
        actionButton(ns("parManBut"), "Update"), 
        h5(strong('Quality Control')),
        fluidRow(splitLayout(cellWidths=c('1%', '10%','1%', '33%'), '',
          numericInput(ns('cntThr'), label='Min counts', value=0, min=0, max=Inf, step=50, width=150), '',
          numericInput(ns('nmads'), label='Min median absolute deviations (MADs)', value=3, min=1, max=Inf, step=1, width=150)
        )),
        bsTooltip(ns('cntThr'), title = "perCellQCMetrics", placement = "bottom", trigger = "hover"),
        bsTooltip(ns('nmads'), title = "isOutlier", placement = "bottom", trigger = "hover"),  
        
        h5(strong('Normalization')),
        fluidRow(splitLayout(cellWidths=c('1%', '8%','1%', '8%'), '',
          numericInput(ns('minSize'), label='Min size', value=100, min=5, max=Inf, step=50, width=150), '',
          numericInput(ns('maxSize'), label='Max size', value=3000, min=10, max=Inf, step=100, width=150)
        )),
        bsTooltip(ns('minSize'), title = "quickCluster", placement = "bottom", trigger = "hover"),
        bsTooltip(ns('maxSize'), title = "computeSumFactors", placement = "bottom", trigger = "hover"),

        h5(strong('Variance Modelling')),
        fluidRow(splitLayout(cellWidths=c('1%', '23%','1%', '23%'), '',
          numericInput(ns('hvgN'), label='Top variable genes by number', value=3000, min=5, max=Inf, step=50, width=150), '',
          numericInput(ns('hvgP'), label='Top variable genes by proportion', value=0.1, min=0.0001, max=1, step=0.1, width=150)
        )),
        bsTooltip(ns('hvgN'), title = "getTopHVGs", placement = "bottom", trigger = "hover"),
        bsTooltip(ns('hvgP'), title = "getTopHVGs", placement = "bottom", trigger = "hover"),

        h5(strong('Dimension reduction')),
        fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%'), '',
        numericInput(ns('minRank'), label='Min PCs', value=5, min=2, max=Inf, step=1, width=150), '',
        numericInput(ns('maxRank'), label='Max PCs', value=50, min=3, max=Inf, step=5, width=150)
        )),
        bsTooltip(ns('minRank'), title = "denoisePCA", placement = "bottom", trigger = "hover"),
        bsTooltip(ns('maxRank'), title = "denoisePCA", placement = "bottom", trigger = "hover"),
        fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%'), '',
        numericInput(ns('ncomT'), label='t-SNE dimensions to obtain', value=2, min=2, max=Inf, step=1, width=150), '',
        numericInput(ns('ntopT'), label='Number of top HVGs', value=500, min=50, max=Inf, step=100, width=150)
        )),
        bsTooltip(ns('ncomT'), title = "runTSNE", placement = "bottom", trigger = "hover"),
        bsTooltip(ns('ntopT'), title = "runTSNE", placement = "bottom", trigger = "hover"),
        fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%', '1%', '20%'), '',
        numericInput(ns('ncomU'), label='UMAP dimensions to obtain', value=2, min=2, max=Inf, step=1, width=150), '',
        numericInput(ns('ntopU'), label='Number of top HVGs', value=500, min=50, max=Inf, step=100, width=150), '',
        numericInput(ns('pcs'), label='Number of PCs', value=50, min=5, max=Inf, step=5, width=150)
        )), span(textOutput(ns('msg.umap')), style='color:red'),
        bsTooltip(ns('ncomU'), title = "runUMAP", placement = "bottom", trigger = "hover"),
        bsTooltip(ns('ntopU'), title = "runUMAP", placement = "bottom", trigger = "hover"),
        fluidRow(splitLayout(cellWidths=c('1%', '13%', '1%', '5%'), '',
        h5(strong('Build neighbor graphs'))
        )),
        fluidRow(splitLayout(cellWidths=c('1%', '20%'), '',
        selectInput(ns('nn.graph'), label='Method', choices=c('buildSNNGraph', 'buildKNNGraph'), selected='buildSNNGraph')
        )), uiOutput(ns('graph.par')),
        fluidRow(splitLayout(cellWidths=c('1%', '13%', '1%', '5%'), '',
          h5(strong('Clustering')),  
        )),
        fluidRow(splitLayout(cellWidths=c('1%', '20%'), '',
        selectInput(ns('scell.cluster'), label='Method', choices=c('cluster_walktrap', 'cluster_fast_greedy', 'cluster_leading_eigen'), selected='cluster_walktrap')
        )), uiOutput(ns('clus.par'))
 
    ), # tabPanel(title="Parameters"
      tabPanel(title="Data Table", value='datCell', br(),
      box(id='bulk', title = "Bulk Data", width = 12, closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      dataTableOutput(ns("datCovisBlk"))
      ),
      box(id='cell', title = "Cell Data", width = 12, closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      dataTableOutput(ns("datCell"))
      )

      ), # navbarPage tabPanel 
      tabPanel("Dimensionality Reduction", value='dimred',
      navbarPage('', id=ns('dimredNav'), 
      tabPanel('Plot', dim_ui(ns('dim'))),
     )), #tabPanel("Dimensionality Reduction", navbarPage
      tabPanel("Matching cells and bulk", value='manualMatch',
      #navbarPage('Parameters:',
      #tabPanel(title="Match spatial features", value='matchCell',
        match_ui(ns('rematchCell'))
      #) # navbarMenu
      #) # navbarPage 
      ) #tabPanel
  ) 
}
scell_ui <- function(id) { 
  ns <- NS(id)
  tabPanel("Co-visualization", value='scell', icon=NULL,
    fluidRow(splitLayout(cellWidths=c('1%', '15%','1%', '15%'), '',
      # uiOutput(ns('methCovis')), '', uiOutput(ns('direc'))
      selectInput(ns('methCovis'), label='Methods', choices=c('Annotation/manual'='man', 'Automatic'='auto'), selected='man'), '',
      selectInput(ns('direc'), label='Mapping direction', choices=cho <- c('Cell2Bulk'='toBulk', 'Bulk2Cell'='toCell'), selected='toBulk')
    )),
    covis_man_ui(ns('covisMan')),

    tabsetPanel(type = "pills", id=ns('tabSetAuto'), selected="autoMatch",
      tabPanel("Automatic", value='autoMatch',  
      navbarPage('', id='autoMatNav', selected='tailor',  
      tabPanel('Result',
        # h5(strong('Subset assignments')),  
        fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '10%'), '',
        numericInput(ns('asgThr'), label='Assignment threshold', value=0, min=0, max=1, step=0.1, width=150), '',
        actionButton(ns("coclusAsg"), "Confirm", style='margin-top:24px')
        )),
        dataTableOutput(ns("resAsg")) 
      ),


      tabPanel('Parameters',
        h3(strong('Initial filtering bulk and cells')),  
        fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '10%', '1%', '18%', '1%', '18%', '1%', '10%', '1%', '10%'), '',
        numericInput(ns('initFilBlkP'), label='P', value=0.05, min=0, max=1, step=0.1, width=150), '',
        numericInput(ns('initFilBlkA'), label='A', value=5, min=0, max=1000, step=5, width=150), '',
        numericInput(ns('initFilBlkCV1'), label='Min coefficient of variation (CV1)', value=0.05, min=-1000, max=1000, step=0.1, width=200), '',
        numericInput(ns('initFilBlkCV2'), label='Max coefficient of variation (CV2)', value=100, min=-1000, max=1000, step=0.1, width=200), '',
        numericInput(ns('initFilPGen'), label='P by gene', value=0.05, min=0, max=1, step=0.1, width=150), '',
        numericInput(ns('initFilPCell'), label='P by cell', value=0.01, min=0, max=1, step=0.1, width=150)
        )), 
        h3(strong('Normalizing bulk and cells')),  
        fluidRow(splitLayout(cellWidths=c('1%', '15%'), '',
        selectInput(ns('normCoclus'), label='Method', choices=c('computeSumFactors'='fct', 'CPM'='cpm'), selected='fct', width=200)
        )), 
        h3(strong('Secondary filtering bulk and cells')),  
        fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '10%', '1%', '18%', '1%', '18%', '1%', '10%', '1%', '10%'), '',
        numericInput(ns('filBlkP'), label='P', value=0.1, min=0, max=1, step=0.1, width=150), '',
        numericInput(ns('filBlkA'), label='A', value=1, min=0, max=1000, step=5, width=150), '',
        numericInput(ns('filBlkCV1'), label='Min coefficient of variation (CV1)', value=0.1, min=-1000, max=1000, step=0.1, width=200), '',
        numericInput(ns('filBlkCV2'), label='Max coefficient of variation (CV2)', value=200, min=-1000, max=1000, step=0.1, width=200), '',
        numericInput(ns('filPGen'), label='P by gene', value=0.01, min=0, max=1, step=0.1, width=150), '',
        numericInput(ns('filPCell'), label='P by cell', value=0.1, min=0, max=1, step=0.1, width=150)
        )), 
        h3(strong('Clustering')),  
        fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '13%'), '',
        selectInput(ns('clusDim'), label='Reducing dimensionality', choices=c('PCA', 'UMAP'), selected='PCA'), '',
        selectInput(ns('clusMeth'), label='Building graph', choices=c('buildSNNGraph'='snn', 'buildKNNGraph'='knn'), selected='knn')
        # selectInput(ns('clusSim'), label='Similarity', choices=c('Spearman', 'Pearson'), selected='Spearman')
        )),
        h5(strong('Refining cell clusters')),
        fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '10%'), '',
        numericInput(ns('clusSimT'), label='Similarity threshold', value=0.2, min=0, max=1, step=0.1, width=150), '',
        numericInput(ns('clusSimP'), label='Similarity proportion', value=0.8, min=0, max=1, step=0.1, width=150)
        )),
        h5(strong('Co-clustering bulk and cells')),
        fluidRow(splitLayout(cellWidths=c('1%', '10%'), '',
        numericInput(ns('clusDimN'), label='Top dims', value=13, min=5, max=50, step=1, width=150)
        )), 
        actionButton(ns("coclusPar"), "Confirm parameters", style='margin-top:1px')
      ), #tabPanel 
     tabPanel("Tailoring matching", value='tailor',
         tailor_match_ui(ns('tailor'))
     ) #tabPanel("Tailoring matching", navbarPage     
     )) # tabPanel("Auto-matching", navbarPage(

   ) # tabsetPanel(
  ) # tabPanel("Single Cell",
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
       $("header").find("nav").append(\'<span class="mainTitle">spatialHeatmap</span>\');
     })
   ')),
    fluidRow( 
      do.call(tabsetPanel, append(list(type = "pills", id = 'shm.sup', selected="landing", upload_ui('upl'), shm_ui('shmAll', data_ui('dat', dim_ui('datDim'), tailor_match_ui('datDimTailor')), search_ui('sear')), deg_ui('deg'), scell_ui('scell')),
        list(tabPanel("About", value='about',
          if (0) box(width = 12, title = "", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
          ),
        includeHTML("instruction/about.html")
        ))
    )) # do.call append

    )
 ) # dashboardBody
)

}
