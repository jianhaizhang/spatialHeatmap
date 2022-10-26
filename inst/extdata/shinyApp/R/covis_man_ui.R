# Module for co-visualization through annotation/manual methods.
covis_man_ui <- function(id) { 
  ns <- NS(id)
    tabsetPanel(type = "pills", id=ns('tabSetCell'), selected="datCell",
      tabPanel(title="Parameters", value='parMan', br(),
        actionButton(ns("parManBut"), "Update"), 
        h4(strong('Bulk tissues')), br(),
        fluidRow(splitLayout(cellWidths=c('1%', '10%'), '',
        selectInput(ns("normBlk"), "Normalization", c('None'='none', "CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "ESF", "VST", "rlog"), selected='VST')
        )), 
        h4(strong('Single Cells')), br(),
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
        )), 
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
