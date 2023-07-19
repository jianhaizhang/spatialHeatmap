# Module for co-visualization through automatic method.
covis_auto_ui <- function(id) { 
  ns <- NS(id)
  tabsetPanel(type = "pills", id=ns('tabSetCellAuto'), selected="datCell",
    tabPanel(title="Settings", value='parAuto', br(),
      actionButton(ns("parAutoBut"), "Run", style=run.col),
      h5(strong('Jointly normalizing bulk and cell data')),  
      fluidRow(splitLayout(cellWidths=c('12px', '250px'), '',
        selectInput(ns('normCoclus'), label='Method', choices=c('computeSumFactors (scran)'='fct', 'CPM (counts/million)'='cpm'), selected='fct')
      )), 
      h5(strong('Filtering')),  
      fluidRow(id=ns('filBlk'), style='width:350px', splitLayout(cellWidths=c('12px', '97px', '1px', '70px', '1px', '72px', '1px', '72px'), '',
        numericInput(ns('filBlkP'), label='Proportion (P)', value=0.1, min=0, max=1, step=0.1), '',
        numericInput(ns('filBlkA'), label='Cutoff (A)', value=1, min=0, max=1000, step=5), '',
        numericInput(ns('filBlkCV1'), label='CV1 (min)', value=0.1, min=-1000, max=1000, step=0.1), '',
        numericInput(ns('filBlkCV2'), label='CV2 (max)', value=200, min=-1000, max=1000, step=0.1)
      )),
      fluidRow(id=ns('filCell'), style='width:350px', splitLayout(cellWidths=c('12px', '72px', '1px', '98px', '1px', '98px'), '',
        numericInput(ns('cutoff'), label='Cutoff (A)', value=1, min=-Inf, max=Inf, step=0.5), '',
        numericInput(ns('filPGen'), label='P in gene (P1)', value=0.01, min=0, max=1, step=0.1), '',
        numericInput(ns('filPCell'), label='P in cell (P2)', value=0.1, min=0, max=1, step=0.1)
      )),
      bsTooltip(ns('filBlk'), title = "Filtering bulk data: <br/> 1. Genes with expression values > A across at least P of all samples, <br/> 2. Coefficient of variation (CV) is between CV1 and CV2.", placement = "top", trigger = "hover"),
      bsTooltip(ns('filCell'), title = "Filtering cell data: <br/> 1. genes with expression values > A across at least P1 of cells will remain, <br/> 2. cells with expression values > A across at least P2 of genes", placement = "top", trigger = "hover"),
      h5(strong('Joint dimension reduction for bulk and single-cell data')),
      fluidRow(id=ns('dims'), style='width:250px', splitLayout(cellWidths=c('12px', '60px', '1px', '60px', '1px', '80px'), '',
        numericInput(ns('minRank'), label='Min PCs', value=14, min=2, max=Inf, step=1), '',
        numericInput(ns('maxRank'), label='Max PCs', value=50, min=3, max=Inf, step=5), '',
        selectInput(ns('dimSel'), label='Method', choices=c('PCA', 'UMAP'), selected='PCA')
      )),
      bsTooltip(ns('dims'), title = "Mix and max principal components (PCs) to retain in PCA, which will be used as input for UMAP.", placement = "top", trigger = "hover"),
      h5(strong('Co-clustering on joint dimensions')),
      fluidRow(id=ns('coclusTp'), style='width:450px', splitLayout(cellWidths=c('12px', '200px', '1px', '210px', '1px', '130px'), '',
        selectInput(ns('graphMeth'), label='Building graphs', choices=c('buildSNNGraph (scran)'='snn', 'buildKNNGraph (scran)'='knn'), selected='knn'), '',
        selectInput(ns('clusMeth'), label='Detecting clusters', choices=c('cluster_walktrap (igraph)'='wt', 'cluster_fast_greedy (igraph)'='fg', 'cluster_leading_eigen (igraph)'='le'), selected='wt'), '',
        numericInput(ns('asgThr'), label='Similarity cutoff (A)', value=0, min=0, max=1, step=0.1, width=150)
      )),
      bsTooltip(ns('coclusTp'), title = "Co-clustering is performed on the top joint dimensions: <br/> 1. A graph is built where nodes are cells (or tissues) and edges are connections between nearest neighbors, <br/> 2. The graph is partitioned to obtain clusters, <br/> 3. In each cluster, cells are assigned to tissues with a nearest-neighbor approach based on Spearman correlation coefficient (similarity). <br/> A: To retain robust bulk-cell assignments, set a similarity cutoff.", placement = "top", trigger = "hover")
    ), # tabPanel(title="Parameters" 
    tabPanel(title="Experiment and Image Data", value='datCell',
      dat_all_ui(ns('dat'))
    ),  
     tabPanel(span(id=ns('resTab'), "Results"), value='result', 
     dim_ui(ns('dim')),
     bsTooltip(ns('resTab'), title = "Exploring co-clustering results before co-visualization.", placement = "top", trigger = "hover")
     ),
     tabPanel("Tailoring Assignments", value='tailor',
       tailor_match_ui(ns('tailor'))
     ), #tabPanel
     tabPanel(span('Help', style=hp.txt), value='help', htmlOutput(ns('help')))
  ) 
}
