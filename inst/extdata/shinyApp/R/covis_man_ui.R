# Module for co-visualization through annotation/manual methods.
covis_man_ui <- function(id) { 
  ns <- NS(id)
    tabsetPanel(type = "pills", id=ns('tabSetCell'), selected="dimred",
      tabPanel(title="Settings", value='parMan', br(),
      actionButton(ns("parManBut"), "Run", style=run.col),
      h5(strong('Jointly normalizing bulk and cell data')),
      fluidRow(splitLayout(cellWidths=c('12px', '250px'), '',
        selectInput(ns('norm'), label='Method', choices=c('computeSumFactors (scran)'='fct', 'CPM (Counts/million)'='cpm'), selected='fct')
      )),
      h5(strong('Filtering')),
      fluidRow(id=ns('filBlk'), style='width:350px', splitLayout(cellWidths=c('12px', '97px', '1px', '70px', '1px', '72px', '1px', '72px'), '',
        numericInput(ns('filBlkP'), label='Proportion (P)', value=0.1, min=0, max=1, step=0.1), '',
        numericInput(ns('filBlkA'), label='Cutoff (A)', value=1, min=0, max=1000, step=5), '',
        numericInput(ns('filBlkCV1'), label='CV1 (min)', value=0.1, min=-1000, max=1000, step=0.1), '',
        numericInput(ns('filBlkCV2'), label='CV2 (max)', value=200, min=-1000, max=1000, step=0.1)
      )),
      fluidRow(id=ns('filCell'), style='width:300px', splitLayout(cellWidths=c('12px', '72px', '1px', '98px', '1px', '98px'), '',
        numericInput(ns('cutoff'), label='Cutoff (A)', value=1, min=-Inf, max=Inf, step=0.5), '',
        numericInput(ns('filPGen'), label='P in gene (P1)', value=0.01, min=0, max=1, step=0.1), '',
        numericInput(ns('filPCell'), label='P in cell (P2)', value=0.1, min=0, max=1, step=0.1)
      )),
      bsTooltip(ns('filBlk'), title = "Filtering bulk data: <br/> 1. Genes with expression values > A across at least P of all samples, <br/> 2. Coefficient of variance (CV) is between CV1 and CV2.", placement = "top", trigger = "hover"),
      bsTooltip(ns('filCell'), title = "Filtering cell data: <br/> 1. Genes with expression values > A across at least P1 of cells, <br/> 2. cells with expression values > A across at least P2 of genes.", placement = "top", trigger ="hover"),
      h5(strong('Dimension reduction')),
      fluidRow(id=ns('dims'), style='width:270px', splitLayout(cellWidths=c('12px', '112px', '1px', '112px'), '',
        numericInput(ns('minRank'), label='Min dimensions', value=5, min=2, max=Inf, step=1), '',
        numericInput(ns('maxRank'), label='Max dimensions', value=50, min=3, max=Inf, step=5)
      )), 
      bsTooltip(ns('dims'), title = "Mix and max principal components (PCs) to retain in PCA, which will be used as input for UMAP and TSNE.", placement = "top", trigger = "hover")
    ), # tabPanel(title="Parameters"
    tabPanel(title="Experiment and Image Data", value='datCell',
      dat_all_ui(ns('dat'))
    ),
      tabPanel(span(id=ns('scTab'), "Single-Cell Data"), value='dimred',
      bsTooltip(ns('scTab'), title = "Exploring single-cell data before tissue-cell matching and co-visualization.", placement = "top", trigger = "hover"),
      tabPanel('Plot', dim_ui(ns('dim'), coclus=FALSE))
     ), #tabPanel("Dimensionality Reduction", navbarPage
      tabPanel(span(id=ns('matTab'), "Matching Cells and Bulk"), value='manualMatch',
      bsTooltip(ns('matTab'), title = "Matching cells and bulk tissues for co-visualization.", placement = "top", trigger = "hover"),
      #navbarPage('Parameters:',
      #tabPanel(title="Match spatial features", value='matchCell',
        match_ui(ns('rematchCell'))
      #) # navbarMenu
      #) # navbarPage 
      ), #tabPanel
      tabPanel(span('Help', style=hp.txt), value='help', htmlOutput(ns('help')))
  ) 
}
