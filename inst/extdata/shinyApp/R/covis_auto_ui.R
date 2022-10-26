# Module for co-visualization through automatic method.
covis_auto_ui <- function(id) { 
  ns <- NS(id)
  tabsetPanel(type = "pills", id=ns('tabSetCellAuto'), selected="datCell",
    tabPanel(title="Parameters", value='parAuto', br(),
      actionButton(ns("parAutoBut"), "Update", style='margin-top:1px'),
      h5(strong('Normalizing bulk and cell data')),  
      fluidRow(splitLayout(cellWidths=c('1%', '15%'), '',
        selectInput(ns('normCoclus'), label='Method', choices=c('computeSumFactors'='fct', 'CPM'='cpm'), selected='fct', width=200)
      )), 
      h5(strong('Filtering')),  
      fluidRow(splitLayout(cellWidths=c('1%', '6%', '1%', '6%', '1%', '28%', '1%', '28%', '1%', '6%', '1%', '6%'), '',
        numericInput(ns('filBlkP'), label='P', value=0.1, min=0, max=1, step=0.1, width=150), '',
        numericInput(ns('filBlkA'), label='A', value=1, min=0, max=1000, step=5, width=150), '',
        numericInput(ns('filBlkCV1'), label='Min coefficient of variation (CV1)', value=0.1, min=-1000, max=1000, step=0.1, width=200), '',
        numericInput(ns('filBlkCV2'), label='Max coefficient of variation (CV2)', value=200, min=-1000, max=1000, step=0.1, width=200), '',
        numericInput(ns('filPGen'), label='P in gene', value=0.01, min=0, max=1, step=0.1, width=150), '',
        numericInput(ns('filPCell'), label='P in cell', value=0.1, min=0, max=1, step=0.1, width=150)
      )), 
      bsTooltip(ns('filBlkP'), title = "Filtering bulk data", placement = "bottom", trigger = "hover"),
      bsTooltip(ns('filBlkA'), title = "Filtering bulk data", placement = "bottom", trigger = "hover"),
      bsTooltip(ns('filBlkCV1'), title = "Filtering bulk data", placement = "bottom", trigger = "hover"),
      bsTooltip(ns('filBlkCV2'), title = "Filtering bulk data", placement = "bottom", trigger = "hover"),
      bsTooltip(ns('filPGen'), title = "Filtering cell data", placement = "bottom", trigger = "hover"),
      bsTooltip(ns('filPCell'), title = "Filtering cell data", placement = "bottom", trigger = "hover"),
      h5(strong('Dimension reduction')),
      fluidRow(splitLayout(cellWidths=c('1%', '11%', '1%', '11%'), '',
        numericInput(ns('minRank'), label='Min dimensions', value=5, min=2, max=Inf, step=1, width=150), '',
        numericInput(ns('maxRank'), label='Max dimensions', value=50, min=3, max=Inf, step=5, width=150)
      )),
      h5(strong('Co-clustering')),
      fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '14%', '1%', '13%', '1%', '10%'), '',
        selectInput(ns('dimSel'), label='Dimensions for building graphs', choices=c('PCA', 'UMAP'), selected='PCA'), '',
        selectInput(ns('graphMeth'), label='Building graphs', choices=c('buildSNNGraph'='snn', 'buildKNNGraph'='knn'), selected='knn'), '',
        selectInput(ns('clusMeth'), label='Detecting clusters', choices=c('cluster_walktrap'='wt', 'cluster_fast_greedy'='fg', 'cluster_leading_eigen'='le'), selected='wt'), '',
        numericInput(ns('asgThr'), label='Similarity cutoff', value=0, min=0, max=1, step=0.1, width=150)
      )),
      bsTooltip(ns('asgThr'), title = "Assignments with similarities above this cutoff are retained.", placement = "bottom", trigger = "hover")
    ), # tabPanel(title="Parameters"
      
    tabPanel(title="Data Table", value='datCell', br(),
      box(id='bulk', title = "Bulk Data", width = 12, closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      dataTableOutput(ns("datCovisBlk"))
      ),
      box(id='cell', title = "Cell Data", width = 12, closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      dataTableOutput(ns("datCovisCell"))
      )
      ), # navbarPage tabPanel 
      tabPanel("Results", value='result', dim_ui(ns('dim'))),
     tabPanel("Tailoring assignments", value='tailor',
       tailor_match_ui(ns('tailor'))
     ) #tabPanel
  ) 
}
