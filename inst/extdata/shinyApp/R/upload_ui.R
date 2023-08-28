# Module for uploading files.
upload_ui <- function(id) {
   ns <- NS(id)
   tabPanel(title="Datasets", value='dataset',
   navbarPage('', selected='datSVG', 
     tabPanel(title="Data & aSVGs", value='datSVG',
      fluidRow(splitLayout(cellWidths=c('10px', '430px', '1px', '70px'), '', 
      h4(strong("Step1: upload custom or select default datasets")),
      actionButton(ns("dathelp"), 'Help', style=paste0('margin-top:10px;', hp))
      )),
      fluidRow(splitLayout(cellWidths=c('1px', '430px'), '', 
      selectInput(ns("fileIn"), label=NULL, choices=('Choosing a dataset'='none'), selected='none')
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
      bsTooltip(id='svgInpath', title = "The data is matched with a single aSVG file.", placement = "bottom", trigger = "hover")

      #h4(strong("Additional files")), 
      #fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '35%'), '',
      #tags$div(class='tp', span(class='tpt', 'Upload a config file in ".yaml" format.'),
      #fileInput(ns("config"), "Upload a config file (optional)", accept=".yaml", multiple=FALSE)), '',
      #tags$div(class='tp', span(class='tpt', 'The batched data sets will be listed under "Step 1".'),
      #fileInput(ns("tar"), "Upload batched data, aSVGs in two separate tar files (optional)", accept=c(".tar"), multiple=TRUE))
      #)),
      #div(style = "font-size: 10px; padding: 0px 0px; margin:0%",
      #fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '35%'), '',
      #downloadButton(ns("dld.cfg"), "Example config file"), '', downloadButton(ns("dld.bat"), "Example data/aSVGs in batch")
      # )))
   ), # tabPanel(title="Data & aSVGs",
   tabPanel(title="Example datasets", value='exDat',
     h4(strong("Download pre-configured assay/image datasets.")),
     fluidRow(splitLayout(cellWidths=c('10px', '250px', '1px', '255px', '1px', '285px', '1px', '295px'), '',
     downloadButton(ns("dld.sgl"), "Example1: data & a single aSVG"), '',
     downloadButton(ns("dld.mul"), "Example2: data & multiple aSVGs"), '',
     downloadButton(ns("dld.st"), "Example3: multi-variable data & aSVG"), '',
     downloadButton(ns("dld.covis"), "Example4: co-visualization data & aSVG")
     )) # br(), 
   ), tabPanel(title=span('Help', style=hp), value='help', htmlOutput(ns('help'))) 
  ) # tabsetPanel(selected="gallery",
   ) # tabPanel(title="Datasets", 
}
