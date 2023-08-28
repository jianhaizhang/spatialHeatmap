# Module for plotting SHMs.
deg_ui <- function(id) { 
  ns <- NS(id)
  tabPanel(title=span('Spatial Enrichment', title="This panel is inactive when Co-visualization is enabled."), value='deg',
   # box(title='', status="primary", solidHeader=TRUE, collapsible=TRUE, width=12,
    tabsetPanel(type="pills", id=ns('degAll'), selected="set",
        tabPanel(title="Settings", value='set',
        # dataTableOutput(ns("dtvs3")), br(),
        column(12, style=rec, div(style='margin-top:2px'),
        actionButton(ns("run"), "Run", icon=icon("sync"), style=run.col),
        bsTooltip(ns('run'), title='After clicking, enrichment will be performed for each spatial feature, and results of a certain feature can be queried in the next tab.', placement = "right", trigger = "hover"),
        p('Pre-processing data', style='font-size:18px'),
        fluidRow(splitLayout(cellWidths=c('21px', '423px', '15px', '120px'), '',
        div(id='filter',
        fluidRow(splitLayout(cellWidths=c('1px', '100px', '1px', '123px', '1px', '100px'), '',
        numericInput(ns("A"), label="1.a Cutoff (A)", value=0), '',
        numericInput(ns("P"), label="1.b Proportion (P)", min=0, max=1, step=0.1, value=0), '',
        numericInput(ns("CV1"), label="1.c CV1", value=-10^4), '',
        numericInput(ns("CV2"), label="1.d CV2", value=10^4)
        ))), '',
        div(id=ns('norMethD'), selectInput(ns("norMeth"), "2. Normalize", c("CNF-TMM"='TMM', "CNF-TMMwsp"='TMMwsp', "CNF-RLE"='RLE', "CNF-upperquartile"='upperquartile', 'None'='none'), 'TMM'))
        )),
        bsTooltip(id='filter', title="Rows passing the following filtering will remain: <br/> 1. Expression values > A across at least P of all samples, <br/> 2. Coefficient of variation (CV) is between CV1 and CV2.", placement = "bottom", trigger = "hover"),  
        bsTooltip(ns("norMethD"), title="CNF: calcNormFactors (edgeR)", placement = "top", trigger = "hover"),  
        p('Spatial enrichment', style='font-size:18px'),
        fluidRow(splitLayout(cellWidths=c('10px', '350px', '1px', '350px', '1px', '145px', '1px', '70px', '1px', '110px'), '',
        uiOutput(ns('sam')), '', uiOutput(ns('con')), '', 
        selectInput(ns("comBy"), "Compare across", c('spatial features'='feature', 'variables'='variable'), 'feature') %>% spsComps::bsTooltip(title= paste0(as.character(HTML('Compare across spatial features: variables under the same spatial feature will be treated as replicates. <br/> Compare across variables: spatial featutes under the same variable will be treated as replicates. <br/>')), as.character(a(href='html/shm_shiny_manual.html#3_Spatial_Enrichment',  target='blank', 'More.'))), placement='top', html=TRUE, click_inside=TRUE), '',
        numericInput(ns('outlier'), 'Outliers', 0, min=0, max=Inf, step=1), '',
        selectInput(ns('meth'), label='Select methods', choices=c('edgeR'='edgeR', 'limma-voom'='limma.voom', 'limma'='limma', 'DESeq2'='DESeq2', 'distinct'='distinct'), selected=c('edgeR'), multiple=FALSE)
        )),
        fluidRow(id=ns('fcfdr'), style='width:250px', splitLayout(cellWidths=c('10px', '125px', '1px', '80px'), '',
        numericInput(ns('ssg.fc'), 'Log2-fold change', 1, min=0, max=Inf, step=1), '',
        numericInput(ns('ssg.fdr'), 'FDR', 0.05, min=0, max=1, step=0.01)
        )),
        bsTooltip(id=ns('sam'), title='Selected spatial features will be considered for enrichment.', placement = "top", trigger = "hover"),
        bsTooltip(id=ns('con'), title='Selected experimental variables will be considered for enrichment.', placement = "top", trigger = "hover"),
        # bsTooltip(id=ns('comBy'), title='Compare across spatial features: variables under the same spatial feature will be treated as replicates. <br/> Compare across variables: spatial featutes under the same variable will be treated as replicates.', placement = "top", trigger = "hover"),
        bsTooltip(id=ns('outlier'), title='Allowed outliers in the references.', placement = "top", trigger = "hover"),
        bsTooltip(id=ns('fcfdr'), title='Criteria for selecting enriched or depleted genes.', placement = "top", trigger = "hover"),
        ), column(12, style='margin-top:2px'),
        column(12, style=rec, div(style='margin-top:2px'),
        p('Input data', style='font-size:18px'),
        dat_all_ui(id=ns('dat'))
        )
        ), # tabPanel
        tabPanel(title="Query the Results", value='ovl',
          div(style=rec, 
           fluidRow(splitLayout(cellWidths=c('1%', '12%'), '',
             uiOutput(ns('query'))
           )),
          bsTooltip(ns('query'), title='Enrichment results will be shown for the query.', placement = "top", trigger = "hover"),
          dataTableOutput(ns("dtvs1")),
          bsTooltip(ns('dtvs1'), title='Summary of the query and references. Strings after the "\\_": replicates.', placement = "top", trigger = "hover")
          ), div(style='margin-top:2px'),
          column(12, style=rec, div(style='margin-top:2px'),
          fluidRow(splitLayout(cellWidths=c('1%', '12%'), '',
            selectInput(ns('enrType'), label='Enriched/Depleted', choices=c('Enriched'='up', 'Depleted'='down'), selected=c('up'))
          )),
          bsTooltip(ns('enrType'), title='Overlap of enriched/depleted biomolecules across spatial features (or variables).', placement = "top", trigger = "hover"),
          column(4, plotOutput(ns("upset"))), column(4, plotOutput(ns("matrix"))), column(4, plotOutput(ns("venn"))))
        ),
        tabPanel(title="Results of the Query", value='dtDeg', 
          dataTableOutput(ns("dtvs2")), div(style='margin-top:2px'),
          bsTooltip(ns('dtvs2'), title='Summary of the query and references. Strings after the "\\_": replicates.', placement = "top", trigger = "hover"),
          actionButton(ns("eSHMBut"), "Enrichment SHMs", icon=icon("sync"), style=run.col),
          bsTooltip(ns('eSHMBut'), title='Click to plot enrichment spatial heatmaps.', placement = "top", trigger = "hover"),
          downloadButton(ns("dld.ssg.tab"), "Download"),
          # column(12, search_ui(ns('deg')), style='z-index:5'),  
          column(12, dataTableOutput(ns("dtDeg"))),
          bsTooltip(ns('dtDeg'), title='Results of enriched/depleted (up/down) biomolecules for the query. "total": total references.', placement = "top", trigger = "hover")
        ), tabPanel(title=span('Help', style=hp.txt), value='help', htmlOutput(ns('help'))) 
      ) # navbarPage
    ) # box(title='Spatial Enrich'
    # data_ui(ns("datDEG"), deg=TRUE),
}
