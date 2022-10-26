# Module for uploading files.
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
      downloadButton(ns("dld.st"), "Example3: multi-dimensional data & aSVG"), '',
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
