# Module for processing data.
data_ui <- function(id, dim.ui=NULL, tailor.ui=NULL) {
  ns <- NS(id)
  # if (deg == FALSE) tabPanel("Primary Visualization", value='primary',
  box(id=ns('dat'), width = 12, title = "Experiment and Image Data", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      navbarPage('', id=ns('settNav'), selected='dat',
      tabPanel("Overview", value='over', 
        fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", DTOutput(ns("over")), "")),
      ), 
      tabPanel('Experiment design', value='expDsg',
        div(id=ns('uplRefD'), class='shiny-split-layout',
          div(id=ns('uplRefile'), fileInput(ns("uplRef"), "Upload references", accept=c(".txt", ".csv"), multiple=FALSE) %>% spsComps::bsTooltip(title= paste0(as.character(HTML(' Upload reference experimental variables in a one-column table. Ensure rownames are the same with the table of experiment design (see the example). ')), as.character(a(href='html/shm_shiny_manual.html#21_Static_image', target='blank', 'More.'))), placement='top', html=TRUE, click_inside=TRUE), style='width:270px'),
          div(id=ns('refSelD'), checkboxInput(ns('refSel'), 'Use uploaded references', value = TRUE), style='width:190px;margin-top:22px'),
          div(downloadButton(ns("dldRef"), "Example"), style='width:100px;margin-top:24px')
        ),
        bsTooltip(id=ns('refSelD'), title="Use uploaded references or not.", placement = "top", trigger = "hover"),
        fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("expDsg")), "")),
      ),
      tabPanel("Assay data", value='dat',
      div(class='shiny-split-layout',
        div(id="selRow", actionButton(ns("selRow"), "Plot", style=run.top), style='width:70px'),
        div(id="deSel", actionButton(ns("deSel"), "Deselect rows", style='margin-top:24px'), style='width:115px'),
        div(id="page", numericInput(ns('page'), label='Table height', value=300, min=50, max=Inf, step=50), style='width:85px'),
        div(id="spk", selectInput(ns('spk'), label='Sparklines', choices=c('No', 'Yes'), selected='No'), style='width:80px'),
        div(id=ns("refD"), style='margin-top:24px', 
          dropdownButton(inputId=ns('dpdRef'), label='Reference', circle=FALSE, icon=NULL, status='primary', inline=FALSE, 
          div(selectInput(ns('ref'), label='Use references', choices=c('No', 'Yes'), selected='No'), style='width:105px'), 
          div(selectInput(ns('refLog'), label='Use log2 scale', choices=c('No', 'Yes'), selected='Yes'), style='width:105px')
          ),
        ) %>% spsComps::bsTooltip(title= paste0(as.character(HTML('Use relative expression levels based on references defined in Experiment design. ')), as.character(a(href='html/shm_shiny_manual.html#21_Static_image', target='blank', 'More.'))), placement='top', html=TRUE, click_inside=TRUE),
        div(id=ns('datInD'), selectInput(ns('datIn'), label='Input data', choices=c('Complete'='all'), selected='all'), style='width:170px')
      ),
      bsTooltip(id=ns('ref'), title="Use relative expression levels or not.", placement = "right", trigger = "hover"),
      bsTooltip(id=ns('refLog'), title="Yes: use log2 ratio, e.g. log2(treatment/control). <br/> No: use raw ratio, e.g. treatment/control.", placement = "right", trigger = "hover"),
      bsTooltip(id=ns('selRow'), title="Click to plot spatial heatmaps.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('datInD'), title="Assay data for plotting spatial heatmaps.", placement = "top", trigger = "hover"),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dtAll")), ""))
      ), # tabPanel
      tabPanel("Settings", value='sett', 
      actionButton(ns("run"), "Run", icon=icon("sync"), style = run.col),
      div(id=ns('submsg'), style='width:400px',
      fluidRow(splitLayout(cellWidths=c('12px', '90px', '1px', '90px', '1px', '95px', '1px', '90px'), '',
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
      bsTooltip(id=ns('normDat'), title="CNF: calcNormFactors (edgeR). <br/> ESF: estimateSizeFactors (DESeq2). <br/> VST: varianceStabilizingTransformation (DESeq2). <br/> rlog: regularized log (DESeq2).", placement = "right", trigger = "hover"),
      bsTooltip(id=ns('scl'), title="Row: scale each row independently. <br/> Selected: scale across all selected rows as a whole. <br/> All: scale across all rows as a whole.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('log'), title="No: skipping this step. <br/> Log2: transform data to log2-scale. <br/> Exp2: transform data to the power of 2.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('filter'), title="Rows passing the following filtering will remain: <br/> 1. Expression values > A across at least P of all samples <br/> 2. Coefficient of variation (CV) is between CV1 and CV2.", placement = "top", trigger = "hover"),
      bsTooltip(id=ns('thr'), title='Values > "Max value": set to "Max value". <br/> Values < "Min value": set to "Min value".', placement = "top", trigger = "hover")
      ) # tabPanel
      #tabPanel('Selected', value='dTabSel',
      #  fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput(ns("dtSel")), "")),
      #  fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", plotOutput(ns("selProf")), ""))
      #),
      # navbarPage 
      #tabPanel('Single-cell metadata', value='dTabScell',
      #  column(12, id='colTailorUI', tailor.ui),
      #  column(12, id='colDimUI', dim.ui) 
      #)
   ) # tabsetPanel(
  ) 
}
