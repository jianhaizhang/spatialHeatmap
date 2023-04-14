# Module for plotting SHMs.
deg_ui <- function(id) { 
  ns <- NS(id)
  tabPanel(title=span('Spatial Enrichment', title="This panel is inactive when Co-visualization is enabled."), value='deg',
    div(style='margin-top:10px'),
    box(title='Spatial Enrichment', status="primary", solidHeader=TRUE, collapsible=TRUE, width=12,
      navbarPage(NULL, id=ns('degAll'), 
        tabPanel(title="Settings", value='set',
        # dataTableOutput(ns("dt.vs3")), br(),
        column(12, style=rec, div(style='margin-top:2px'),
        actionButton(ns("run"), "Run", icon=icon("sync"), style=run.col),
        p('Pre-processing data', style='font-size:18px'),
        fluidRow(splitLayout(cellWidths=c('21px', '423px', '15px', '120px'), '',
        div(id='filter',
        fluidRow(splitLayout(cellWidths=c('1px', '100px', '1px', '123px', '1px', '100px'), '',
        numericInput(ns("A"), label="1.a Cutoff (A)", value=0), '',
        numericInput(ns("P"), label="1.b Proportion (P)", min=0, max=1, step=0.1, value=0), '',
        numericInput(ns("CV1"), label="1.c CV1", value=-10^4), '',
        numericInput(ns("CV2"), label="1.d CV2", value=10^4)
        ))), '',
        selectInput(ns("norMeth"), "2. Normalize", c("CNF-TMM"='TMM', "CNF-TMMwsp"='TMMwsp', "CNF-RLE"='RLE', "CNF-upperquartile"='upperquartile', 'None'='none'), 'TMM')
        )),
        bsTooltip(id='filter', title="Rows passing the following filtering will remain: <br/> 1. Expression values >= A across >= P of all samples <br/> 2. Coefficient of variation (CV) is between CV1 and CV2.", placement = "bottom", trigger = "hover"),  
        p('Spatial enrichment', style='font-size:18px'),
        fluidRow(splitLayout(cellWidths=c('10px', '350px', '1px', '350px', '1px', '145px', '1px', '70px', '1px', '110px'), '',
        uiOutput(ns('ssg.sam')), '', uiOutput(ns('ssg.con')), '', 
        selectInput(ns("comBy"), "Compare across", c('spatial features'='feature', 'variables'='variable'), 'feature'), '',
        numericInput(ns('outlier'), 'Outliers', 0, min=0, max=Inf, step=1), '',
        selectInput(ns('meth'), label='Select methods', choices=c('edgeR'='edgeR', 'limma-voom'='limma', 'DESeq2'='DESeq2', 'distinct'='distinct'), selected=c('edgeR'), multiple=FALSE)
        )),
        fluidRow(splitLayout(cellWidths=c('10px', '125px', '1px', '80px'), '',
        numericInput(ns('ssg.fc'), 'Log2-fold change', 1, min=0, max=Inf, step=1), '',
        numericInput(ns('ssg.fdr'), 'FDR', 0.05, min=0, max=1, step=0.01)
        )),
        bsTooltip(id=ns('comBy'), title='Compare across spatial features: variables under the same spatial feature will be treated as replicates. <br> Compare across variables: spatial featutes under the same variable will be treated as replicates.', placement = "top", trigger = "hover"),
        ), column(12, style='margin-top:2px'),
        column(12, style=rec, div(style='margin-top:2px'),
        p('Input data', style='font-size:18px'),
        dataTableOutput(ns("datAll"))
        )
        ), # tabPanel
        tabPanel(title="Summary of Results", value='ovl',
          div(style=rec, 
           fluidRow(splitLayout(cellWidths=c('1%', '12%'), '',
             uiOutput(ns('query'))
           )),
          dataTableOutput(ns("dt.vs1"))
          ), div(style='margin-top:2px'),
          column(12, style=rec, div(style='margin-top:2px'),
          fluidRow(splitLayout(cellWidths=c('1%', '12%'), '',
            selectInput(ns('enrType'), label='Enriched/Depleted', choices=c('Enriched'='up', 'Depleted'='down'), selected=c('up'))
          )),
          column(4, plotOutput(ns("upset"))), column(4, plotOutput(ns("matrix"))), column(4, plotOutput(ns("venn"))))
        ),
        tabPanel(title="Results of the Query", value='dt.deg', 
          dataTableOutput(ns("dt.vs2")), div(style='margin-top:2px'),
          actionButton(ns("eSHMBut"), "Enrichment SHMs", icon=icon("sync"), style=run.col),
          # column(12, search_ui(ns('deg')), style='z-index:5'),  
          column(12, dataTableOutput(ns("dt.deg"))),
          downloadButton(ns("dld.ssg.tab"), "Download")
        ) 
      ) # navbarPage
    ) # box(title='Spatial Enrich'
    # data_ui(ns("datDEG"), deg=TRUE),

  )
}
