
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
    fluidRow(splitLayout(style='margin-top:3px;margin-bottom:3px', cellWidths=c('0.5%', '10%', '0.5%', '14%', '0.5%', '15%', '0.5%', '11%'), '',
    dropdownButton(inputId=ns('dpwNetType'), label='Network type', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      selectInput(inputId=ns("net.type"), label="", choices=c('signed', 'unsigned', 'signed hybrid', 'distance'), selected='signed', width='100%')
    ), '',
    dropdownButton(inputId=ns('dpwModSize'), label='Module splitting', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      numericInput(ns("min.size"), "Min module size", value=15, min=15, max=5000, width='100%'),
      selectInput(ns("ds"), HTML("Module splitting sensitivity level <br/> (Larger value results in more <br/> modules with smaller sizes)"), 3:2, selected=3, width='100%')
    ), '',
    dropdownButton(inputId=ns('dpwNetAdj'), label='Adjacency/edges', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=300,
      selectInput(ns("adj.in"), HTML("Adjacency threshold <br/> (the smaller, the more edges)"), sort(seq(0, 1, 0.002), decreasing=TRUE), selected=1, width='100%'),
      numericInput(ns("max.edg"), HTML("Maximun edges <br/> (too many edges may crash the app)"), value=50, min=1, max=500, width='100%'), htmlOutput(ns("edge"))
    ), '',
    dropdownButton(inputId=ns('dpwColNet'), label='Color key', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=300,
      fluidRow(splitLayout(cellWidths=c('1%', '60%', '35%'), '',
      textInput(ns("color.net"), "Color scheme", 'yellow,orange,red', placeholder='Eg: yellow,orange,red', width='100%'),
      actionButton(ns("col.but.net"), "Confirm", icon=icon("sync"), style='margin-top:24px'))) # fluidRow(splitLayout
    )
    )) # fluidRow(splitLayout
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
    column(id='elbow', width=6, plotOutput(ns("elbow"), hover = hoverOpts(id = ns("plot_hover")))),
    column(id='clus', width=6, plotOutput(ns("cluster"), hover = hoverOpts(id = ns("plot_hover")))),
    bsTooltip(id=ns('elbow'), title='The optimal k is usually at/near the elbow in the line graph.', placement = "top", trigger = "hover")
  ),
  tabPanel(strong('Settings'), value='set',
  fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '10%', '1%', '10%', '1%', '10%'), '', 
    numericInput(ns('seed'), label='Seed', value=10, min=1, max=Inf, step=1, width='100%'), '',
    numericInput(ns('maxK'), label='Max k', value=10, min=2, max=Inf, step=1, width='100%'), '',
    numericInput(ns('k'), label='Choosen k', value=2, min=2, max=Inf, step=1, width='100%'), '',
    selectInput(ns("dimred"), label="Dimension reduction", choices=c("TSNE", "PCA", "UMAP"))
  )),
  bsTooltip(id=ns('maxK'), title='Run kmeans-clustering with k (number of clusters) ranging from 2 to "Max k" to find the optimal k (elbow method).', placement = "top", trigger = "hover"),
  bsTooltip(id=ns('k'), title='Chosen k in the elbow method', placement = "top", trigger = "hover"),
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
   box(width = 12, title = "Functional Analysis", closable = FALSE, collapsible=TRUE, solidHeader = TRUE,
   enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      fluidRow(splitLayout(cellWidths=c('1%', '10%', '1%', '15%', '1%', '10%', '1%', '5%', '1%', '5%'), '',
        selectInput(ns('species'), 'Species', c('None'='none', 'Mouse'='mus', 'Human'='hum', 'Arabidopsis'='arab', 'Zebrafish'='fish', 'Drosophila'='fly'), selected='none', width='100%'), '',
        selectInput(ns('funID'), 'Source ID', c('None'='none', 'UniProt (e.g. Cav2)'='unpr', 'TAIR (e.g. AT1G01030)'='tair', 'Ensembl (e.g. ENSMUSG00000000058)'='ense', 'Entrez (e.g. 12390)'='entr'), selected='none', width='100%'), '',
        selectInput(ns('stpClus'), 'Input cluster/module', c('Step2', 'Step3'), selected='Step2', width='100%'), '',
        actionButton(ns("FunBut"), "Run", icon=icon("sync"), style='margin-top:24px;color:#fff;background-color:#c685c4;border-color:#ddd'), '', uiOutput(ns('dld'), style='margin-top:24px')
      )),
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
      selectInput(ns("query"), label="Select a query biomoleclue", choices=c("None"), width='100%'), '',
      selectInput(ns("measure"), label="Measure", choices=c('Correlation'='cor', 'Absolute correlation'='absCor', 'Distance'='dis')), '',
      div(style='margin-top:24px',
      dropdownButton(inputId=ns('dpbThr'), label='Cutoff', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=400,
        radioButtons(inputId=ns("thr"), label="Cutoff metrics", choices=c('Number'='n', 'Proportion'='p', 'Similarity/Dissimilarity'='v'), selected='n', inline=TRUE, width='100%'),
        numericInput(inputId=ns('mhm.v'), label='Cutoff value', value=200, min=-Inf, max=Inf, step=0.1, width='100%')
      )
      )
      )),
      bsTooltip(id=ns('measure'), title='Measures for subsetting the whole assay matrix.', placement = "top", trigger = "hover"),
      bsTooltip(id=ns('dpbThr'), title='Cutoff metrics and respective values.', placement = "top", trigger = "hover")
      ), 
      tabPanel(strong('Step2. Clustering/Nework Analysis'), value=ns('cluNet'), br(),
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '5%'), '', 
        div(style='margin-top:-24px',
        selectInput(ns("method"), label="Select a method", choices=c("Hierarchical clustering"='hcl', "K-means clustering"='kmean', "Network analysis"='net'))), '',
        actionButton(ns("showBut"), "Run", icon=icon("sync"), style="color:#fff;background-color:#c685c4;border-color:#ddd")
      )), div(id=ns('mhm'), mhm_ui(ns('mhm'))), div(id=ns('kmean'), kmean_ui(ns('kmean'))), div(id=ns('net'), net_ui(ns('net')))
      ), 

    tabPanel(span(strong("Step3. Network Analysis (optional)"), title='Optional analysis on the cluster containing the query in Step2.'), value=ns('netOp'), icon=NULL,
      fluidRow(splitLayout(cellWidths=c('1%', '12%'), '', 
        actionButton(ns("showButNet"), "Run", icon=icon("sync"), style="color:#fff;background-color:#c685c4;border-color:#ddd")
      )), 
      div(id=ns('net3'), net_ui(ns('netS3'))) 
    )
    ) # navbarPage('', id=ns('mhmNav')
  ), fun_enrich_ui(ns('goKeg'))
  )
  ) # list
}
