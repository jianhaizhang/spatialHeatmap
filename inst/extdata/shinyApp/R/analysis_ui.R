
net_ui <- function(id) {
  ns <- NS(id)
  list(
  # Result and settings of the same module should be in the same module ui.
  navbarPage('', id=ns('cluNetNav'),
  tabPanel(strong('Result'), value='result',
    fluidRow(splitLayout(cellWidths=c('1%', '8%'), '', 
      uiOutput(ns('dld'))
    )),
    fluidRow(splitLayout(cellWidths=c("1%", "5%", "92%", "2%"), "", plotOutput(ns("bar.net")), visNetworkOutput(ns("vis")), ""))
  ),
  tabPanel(strong('Settings'), value='set',
    fluidRow(splitLayout(style='margin-top:3px;margin-bottom:3px', cellWidths=c('1px', '145px', '1px', '155px', '1px', '132px', '1px', '100px'), '',
    dropdownButton(inputId=ns('dpwNetType'), label='Adjacency matrix', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      selectInput(inputId=ns("net.type"), label="", choices=c('Signed (S)'='signed', 'Unsigned (U)'='unsigned', 'Signed hybrid (H)'='signed hybrid', 'Distance (D)'='distance'), selected='signed', width='100%')
    ), '',
    dropdownButton(inputId=ns('dpwModSize'), label='Module partitioning', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      numericInput(ns("min.size"), "Min module size", value=15, min=15, max=5000, width='100%'),
      selectInput(ns("ds"), HTML("Module partitioning strigency <br/> (Larger value results in more <br/> modules with smaller sizes)"), 3:2, selected=3, width='100%')
    ), '',
    dropdownButton(inputId=ns('dpwNetAdj'), label='Edges to show', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=300,
      selectInput(ns("adj.in"), HTML("Adjacency threshold <br/> (the smaller, the more edges to show)"), sort(seq(0, 1, 0.002), decreasing=TRUE), selected=1, width='100%'),
      numericInput(ns("max.edg"), HTML("Maximun edges <br/> (too many edges may crash the app)"), value=50, min=1, max=500, width='100%'), htmlOutput(ns("edge"))
    ), '',
    dropdownButton(inputId=ns('dpwColNet'), label='Color key', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=300,
      fluidRow(splitLayout(cellWidths=c('1%', '60%', '35%'), '',
      textInput(ns("color.net"), "Color scheme", 'yellow,orange,red', placeholder='Eg: yellow,orange,red', width='100%'),
      actionButton(ns("col.but.net"), "Confirm", icon=icon("sync"), style='margin-top:24px'))) # fluidRow(splitLayout
    )
    )), # fluidRow(splitLayout
    bsTooltip(ns('dpwNetType'), title="S: negative correlations (cor) beteen genes are maintained negative; U: absolute cor is used; H: negative cor is treated 0; D: use distance measure.", placement = "top", trigger = "hover"),
    bsTooltip(ns('dpwModSize'), title="A module is a group of genes sharing highly similar expression patterns.", placement = "top", trigger = "hover"), 
    bsTooltip(ns('dpwNetAdj'), title="Too many edges (adjacency between genes) may crash the app, so this option limits edges to show in the graph.", placement = "top", trigger = "hover") 
  ) 
  )
  )
}

kmean_ui <- function(id) {
  ns <- NS(id)
  list(
  # Result and settings of the same module should be in the same module ui.
  navbarPage('', id=ns('cluNetNav'),
  tabPanel(strong('Result'), value='result',
    fluidRow(splitLayout(cellWidths=c('1%', '8%'), '', 
      uiOutput(ns('dld'))
    )),
    # Too identical outputs: Shiny does not run. 
    # uiOutput(ns('dld'), style='margin-top:24px'),
    column(id='elb', width=6,
      plotOutput(ns("elbow"), hover = hoverOpts(id = ns("plot_hover"))) %>% spsComps::bsTooltip(title= paste0(as.character(HTML('The optimal k is usually at/near the elbow in the line graph.<br/>')), as.character(a(href='https://www.analyticsvidhya.com/blog/2021/01/in-depth-intuition-of-k-means-clustering-algorithm-in-machine-learning/', target='blank', 'More.'))), placement='top', html=TRUE, click_inside=TRUE)
    ),
    column(id='clus', width=6, plotOutput(ns("cluster"), hover = hoverOpts(id = ns("plot_hover")))),
    # bsTooltip(id=ns('elbow'), title='The optimal k is usually at/near the elbow in the line graph.', placement = "top", trigger = "hover"),
    bsTooltip(id=ns('cluster'), title='The cluster containing the query is indicated by the suffix ".query".', placement = "top", trigger = "hover")
  ),
  tabPanel(strong('Settings'), value='set',
  fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '10%', '1%', '10%', '1%', '10%'), '', 
    numericInput(ns('seed'), label='Seed', value=10, min=1, max=Inf, step=1, width='100%'), '',
    numericInput(ns('maxK'), label='Max k', value=10, min=2, max=Inf, step=1, width='100%'), '',
    numericInput(ns('k'), label='Choosen k', value=2, min=2, max=Inf, step=1, width='100%'), '',
    selectInput(ns("dimred"), label="Dimension reduction", choices=c("TSNE", "PCA", "UMAP"))
  )),
  bsTooltip(id=ns('maxK'), title='Run kmeans-clustering with k (number of clusters) ranging from 2 to "Max k" to find the optimal k (elbow method).', placement = "top", trigger = "hover"),
  bsTooltip(id=ns('k'), title='Chosen k with the elbow method', placement = "top", trigger = "hover"),
  bsTooltip(id=ns('dimred'), title='Visualizing clusters.', placement = "top", trigger = "hover")
  )
  )
  )
}

mhm_ui <- function(id) {
  ns <- NS(id)
  list(
  # Result and settings of the same module should be in the same module ui.
  navbarPage('', id=ns('cluNetNav'),
  tabPanel(strong('Result'), value='result',
    fluidRow(splitLayout(cellWidths=c('1%', '8%', '1%', '8%'), '', 
      numericInput(ns('cutH'), label='Cutting height', value=20, min=0, max=Inf, step=0.5, width='100%'), '',
      uiOutput(ns('dld'), style='margin-top:24px')
    )),
    bsTooltip(id=ns('cutH'), title='Cut row dendrograms to get a cluster containing the query.', placement = "right", trigger = "hover"),
    plotlyOutput(ns("HMly"))
  ),
  tabPanel(strong('Settings'), value='set',
  fluidRow(splitLayout(cellWidths=c('1%', '10%'), '', 
    selectInput(inputId=ns("mat.scale"), label="Scale by", choices=c("No", "Row", "Column"), selected='Row')
  ))
  ) 
  )
  )
}

fun_enrich_ui <- function(id) {
  ns <- NS(id)
  list(
   box(id=ns('funEnr'), width = 12, title = "Functional Analysis", closable = FALSE, collapsible=TRUE, solidHeader = TRUE,
   enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '15%', '1%', '10%', '1%', '5%', '1%', '5%'), '',
        selectInput(ns('species'), 'Species', c('None'='none', 'Mouse'='mus', 'Human'='hum', 'Arabidopsis'='arab', 'Zebrafish'='fish', 'Drosophila'='fly'), selected='none', width='100%'), '',
        selectInput(ns('funID'), 'Source ID', c('None'='none', 'UniProt (e.g. Cav2)'='unpr', 'TAIR (e.g. AT1G01030)'='tair', 'Ensembl (e.g. ENSMUSG00000000058)'='ense', 'Entrez (e.g. 12390)'='entr'), selected='none', width='100%'), '',
        selectInput(ns('stpClus'), 'Input cluster/module', c('Step2', 'Step3'), selected='Step2', width='100%'), '',
        actionButton(ns("FunBut"), "Run", icon=icon("sync"), style=run.col), '', uiOutput(ns('dld'), style='margin-top:24px')
      )),
      bsTooltip(ns('stpClus'), title="The cluster/module for enrichemnt.", placement = "top", trigger = "hover"), 
    column(id='go', width=6,
     box(title = "GO Enrichment", status = "primary", solidHeader = TRUE, closable = FALSE, collapsible=TRUE, width=NULL,
      fluidRow(splitLayout(cellWidths=c('2%', '25%', '1%', '15%', '1%', '15%', '1%', '15%'), '',      
      selectInput(ns("ont"), label="Ontology", choices=c('Biological process'='BP', 'Molecular funcion'='MF', 'Cellular component'='CC', 'All'='ALL'), selected='BP', width='100%'), '',
      selectInput(ns("Padj"), label="Adjusting P", choices=c("BH", "BY", "fdr", "none"), selected='BH', width='100%'), '',
      numericInput(ns("Pcut"), label="P cutoff", value=0.05, min=0, max=1, width='100%'), '',
      numericInput(ns("terms"), label="Terms", value=5, min=1, max=Inf, width='100%')
      )),
      plotOutput(ns('go'))
    )
    ),
     column(id='kegg', width=6,
       box(title = "KEGG Enrichment", status = "primary", solidHeader = TRUE, closable = FALSE, collapsible=TRUE, width=NULL,
       fluidRow(splitLayout(cellWidths=c('2%', '15%', '1%', '15%', '1%', '15%'), '', 
       selectInput(ns("PadjKeg"), label="Adjusting P", choices=c("BH", "BY", "fdr", "none"), selected='BH', width='100%'), '',
       numericInput(ns("PcutKeg"), label="P cutoff", value=0.05, min=0, max=1, width='100%'), '',
       numericInput(ns("termsKeg"), label="Terms", value=5, min=1, max=Inf, width='100%')
       )), plotOutput(ns('keg'))
      )
     )
   ) 

  )
}
# Module for matrix heatmap and network analysis.
analysis_ui <- function(id) {
  ns <- NS(id)
  list(
  tabPanel(title="Data Mining", value=ns('clus'), icon=NULL, div(style='margin-top:10px'), 
   box(width = 12, title = "Subsetting/Clustering/Network Analysis", closable = FALSE, collapsible=TRUE, solidHeader = TRUE, status='primary',
    navbarPage('', id=ns('clusNav'), selected=ns('cluNet'),
      tabPanel(strong('Step1. Subsetting the Whole Matrix'), value=ns('subset'),
      fluidRow(splitLayout(style='margin-top:3px;margin-bottom:3px', cellWidths=c('1%', '15%', '1%', '12%', '1%', '8%'), '', 
      div(id=ns('queTp'), selectInput(ns("query"), label="Select a query biomoleclue", choices=c("None"), width='100%')), '',
      selectInput(ns("measure"), label="Measure", choices=c('Correlation'='cor', 'Absolute correlation'='absCor', 'Distance'='dis')), '',
      div(style='margin-top:24px',
      dropdownButton(inputId=ns('dpbThr'), label='Cutoff', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=450,
        radioButtons(inputId=ns("thr"), label="Cutoff metrics", choices=c('Number (N)'='n', 'Proportion (P)'='p', 'Similarity/Dissimilarity (Sim)'='v'), selected='n', inline=TRUE, width='100%'),
        numericInput(inputId=ns('mhm.v'), label='Cutoff value', value=400, min=-Inf, max=Inf, step=0.1, width='100%')
      ))
      )),
      bsTooltip(ns('queTp'), title='Biomolecules visualized in spatial heatmaps are listed here.', placement = "top", trigger = "hover"),
      bsTooltip(id=ns('measure'), title='Subsetting the whole assay matrix to get biomolecules showing similar expression patterns with the query using the "Measure" (similarity) and "Cutoff"', placement = "top", trigger = "hover"),
      bsTooltip(id=ns('dpbThr'), title='Cutoffs: N-a fixed number of most similar biomolecues; P-a fixed proportion of most similar biomolecules; Sim-biomolecules with similarities > a fixed similarity value.', placement = "top", trigger = "hover")
      ), 
      tabPanel(strong('Step2. Clustering/Nework Analysis'), value=ns('cluNet'), br(),
      fluidRow(splitLayout(cellWidths=c('10px', '190px', '1px', '70px'), '', 
        div(style='margin-top:-24px',
        selectInput(ns("method"), label="Select a method", choices=c("Hierarchical clustering"='hcl', "K-means clustering"='kmean', "Network analysis"='net'))), '',
        actionButton(ns("showBut"), "Run", icon=icon("sync"), style=run.col)
      )), 
      bsTooltip(ns('method'), title='Perform clustering/network analysis (WGCNA) on the subset matrix from Step1.', placement = "top", trigger = "hover"),
      div(id=ns('mhm'), mhm_ui(ns('mhm'))), div(id=ns('kmean'), kmean_ui(ns('kmean'))), div(id=ns('net'), net_ui(ns('net')))
      ), 

    tabPanel(span(strong("Step3. Network Analysis (optional)"), title='Optional analysis on the cluster containing the query in Step2.'), value=ns('netOp'), icon=NULL,
      fluidRow(splitLayout(cellWidths=c('1%', '12%'), '', 
        actionButton(ns("showButNet"), "Run", icon=icon("sync"), style=run.col)
      )), 
      div(id=ns('net3'), net_ui(ns('netS3'))) 
    ),
    tabPanel(title=span('Help', style=hp), value='help', htmlOutput(ns('help')))
    ) # navbarPage('', id=ns('mhmNav')
  ), fun_enrich_ui(ns('goKeg'))
  )
  ) # list
}
