# Module for co-visualization through annotation/manual methods.
covis_man_ui <- function(id) { 
  ns <- NS(id)
    tabsetPanel(type = "pills", id=ns('tabSetCell'), selected="datCell",
      tabPanel(title="Settings", value='parMan', br(),
      actionButton(ns("parManBut"), "Run", style=run.col),
      h5(strong('Jointly normalizing bulk and cell data')),
      fluidRow(splitLayout(cellWidths=c('12px', '250px'), '',
        selectInput(ns('norm'), label='Method', choices=c('computeSumFactors'='fct', 'CPM'='cpm'), selected='fct')
      )),
      h5(strong('Filtering')),
      div(id=ns('filBlk'),
      fluidRow(splitLayout(cellWidths=c('12px', '70px', '1px', '70px', '1px', '72px', '1px', '72px'), '',
        numericInput(ns('filBlkP'), label='P', value=0.1, min=0, max=1, step=0.1), '',
        numericInput(ns('filBlkA'), label='A', value=1, min=0, max=1000, step=5), '',
        numericInput(ns('filBlkCV1'), label='CV1 (min)', value=0.1, min=-1000, max=1000, step=0.1), '',
        numericInput(ns('filBlkCV2'), label='CV2 (max)', value=200, min=-1000, max=1000, step=0.1)
      ))),
      div(id=ns('filCell'),
      fluidRow(splitLayout(cellWidths=c('12px', '72px', '1px', '98px', '1px', '98px'), '',
        numericInput(ns('cutoff'), label='Cutoff', value=1, min=-Inf, max=Inf, step=0.5), '',
        numericInput(ns('filPGen'), label='P in gene (P1)', value=0.01, min=0, max=1, step=0.1), '',
        numericInput(ns('filPCell'), label='P in cell (P2)', value=0.1, min=0, max=1, step=0.1)
      ))),
      bsTooltip(ns('filBlk'), title = "Filtering bulk data: rows passing the following filtering will remain. <br/> 1. Expression  values >= A across >= P of all samples <br/> 2. Coefficient of variation (CV) is between CV1 and CV2.", placement = "bottom",      trigger = "hover"),
      bsTooltip(ns('filCell'), title = "Filtering cell data: <br/> 1. genes with expression values >= Cutoff across >= P1 cells    will remain <br/> 2. cells with expression values >= Cutoff across >= P2 genes will remain.", placement = "bottom", trigger =      "hover"),
      h5(strong('Dimension reduction')),
      fluidRow(splitLayout(cellWidths=c('12px', '112px', '1px', '112px'), '',
        numericInput(ns('minRank'), label='Min dimensions', value=5, min=2, max=Inf, step=1), '',
        numericInput(ns('maxRank'), label='Max dimensions', value=50, min=3, max=Inf, step=5)
      )) 
    ), # tabPanel(title="Parameters"
      tabPanel(title="Data Table", value='datCell', br(),
      actionButton(ns("subdat"), "Run", icon=icon("sync"), style = run.col),
      div(id=ns('submsg'), 
      fluidRow(splitLayout(cellWidths=c('12px', '70px', '1px', '70px', '1px', '95px', '1px', '89px'), '',
        numericInput(ns('r1'), label='Row start', value=1, min=1, max=Inf, step=1), '',
        numericInput(ns('r2'), label='Row end', value=500, min=2, max=Inf, step=1), '',
        numericInput(ns('c1'), label='Column start', value=1, min=1, max=Inf, step=1), '',
        numericInput(ns('c2'), label='Column end', value=20, min=2, max=Inf, step=1)
      ))),
      bsTooltip(id=ns('submsg'), title="Subsetting the data matrix for display only, not for downstream analysis.", placement =    "top", trigger = "hover"),
      box(id='datall', title = "", width = 12, closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      dataTableOutput(ns("datall"))
      )
      ), # navbarPage tabPanel 
      tabPanel("Dimension Reduction", value='dimred',
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
