# Module for processing data.
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
      )
      #tabPanel('Single-cell metadata', value='dTabScell',
      #  column(12, id='colTailorUI', tailor.ui),
      #  column(12, id='colDimUI', dim.ui) 
      #)
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
