# Module for processing data.
data_ui <- function(id, dim.ui=NULL, tailor.ui=NULL) {
  ns <- NS(id)
  # if (deg == FALSE) tabPanel("Primary Visualization", value='primary',
  box(width = 12, title = "Data Table", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      navbarPage('', id=ns('settNav'), selected='dat', 
      tabPanel("Settings", value='sett', 
      actionButton(ns("run"), "Run", icon=icon("sync"), style = run.col),
      div(id=ns('submsg'),
      fluidRow(splitLayout(cellWidths=c('12px', '70px', '1px', '70px', '1px', '95px', '1px', '89px'), '',
        numericInput(ns('r1'), label='Row start', value=1, min=1, max=Inf, step=1), '',
        numericInput(ns('r2'), label='Row end', value=100, min=2, max=Inf, step=1), '',
        numericInput(ns('c1'), label='Column start', value=1, min=1, max=Inf, step=1), '',
        numericInput(ns('c2'), label='Column end', value=20, min=2, max=Inf, step=1)
      ))),
      fluidRow(splitLayout(cellWidths=c('5px', '120px', '1px', '150px', '15px', '424px', '30px', '202px', '1px', '100px'), '',
      selectInput(ns("normDat"), "1. Normalize", c('None', "CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "ESF", "VST", "rlog"), selected='CNF-TMM'), '',
      selectInput(ns('log'), label='2. Log/exp-transform', choices=c("No", 'Log2'="log2", 'Exp2'="exp2"), selected='No'), '',
      div(id=ns('filter'),
      fluidRow(splitLayout(cellWidths=c('1px', '100px', '1px', '120px', '1px', '100px', '1px', '100px'), '',
      numericInput(ns("A"), label="3.a Cutoff (A)", value=0), '',
      numericInput(ns("P"), label="3.b Proportion (P)", value=0), '',
      numericInput(ns("CV1"), label="3.c CV1", value=-10^4), '',
      numericInput(ns("CV2"), label="3.d CV2", value=10^4)
      ))), '', 
      div(id=ns('thr'),
      fluidRow(splitLayout(cellWidths=c('1px', '100px', '1px', '100px'), '', 
      numericInput(ns("sig.min"), "4.a Min value", value=-10^4), '',
      numericInput(ns("sig.max"), "4.b Max value", value=10^4)
      ))), '', 
      selectInput(ns('scl'), label='5. Scale by', choices=c('No'='No', 'Row'='Row', 'Selected'='Selected', 'All'='All'), selected='No')
      )), div(style='height:200px'),
      bsTooltip(id=ns('submsg'), title="Subsetting the data matrix for display only, not for downstream analysis.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('normDat'), title="Output: log2 scale. <br/> CNF: calcNormFactors (edgeR). <br/> ESF: estimateSizeFactors (DESeq2). <br/> VST: varianceStabilizingTransformation (DESeq2). <br/> rlog: regularized log (DESeq2).", placement = "right", trigger = "hover"),
      bsTooltip(id=ns('scl'), title="Row: scale each row independently. <br/> Selected: scale across all selected rows as a whole. <br/> All: scale across all rows as a whole.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('log'), title="No: skipping this step. <br/> Log2: transform data to log2-scale. <br/> Exp2: transform data to the power of 2.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('filter'), title="Rows passing the following filtering will remain: <br/> 1. Expression values >= A across >= P of all samples <br/> 2. Coefficient of variation (CV) is between CV1 and CV2.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('thr'), title='Values > "Max value": set to "Max value". <br/> Values < "Min value": set to "Min value".', placement = "top", trigger = "hover")
      ), # tabPanel

      tabPanel("Complete", value='dat',
      fluidRow(splitLayout(cellWidths=c('5px', '130px', '1px', '115px', '1px', '80px', '1px', '80px', '1px', '170px'), '',
      actionButton(ns("selRow"), "Confirm selection", style=run.top), '',
      actionButton(ns("deSel"), "Deselect rows", style='margin-top:24px'), '',
      numericInput(ns('page'), label='Page height', value=300, min=50, max=Inf, step=50, width=150), '',
      selectInput(ns('spk'), label='Sparklines', choices=c('No', 'Yes'), selected='No'), '',
      selectInput(ns('datIn'), label='Input data', choices=c('Complete'='all'), selected='all') 
      )),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dtAll")), ""))
      ), # tabPanel
      tabPanel('Selected', value='dTabSel',
        fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dtSel")), "")),
        fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", plotOutput(ns("selProf")), ""))
      )
      # navbarPage 
      #tabPanel('Single-cell metadata', value='dTabScell',
      #  column(12, id='colTailorUI', tailor.ui),
      #  column(12, id='colDimUI', dim.ui) 
      #)
   ) # tabsetPanel(
  ) 
}
