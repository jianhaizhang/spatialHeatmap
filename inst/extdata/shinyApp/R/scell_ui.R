# Module for co-visualization.
scell_ui <- function(id) { 
  ns <- NS(id)
  tabPanel("Co-visualization", value='scell', icon=NULL,
    fluidRow(splitLayout(cellWidths=c('1%', '15%','1%', '15%', '1%', '5%'), '',
      # uiOutput(ns('methCovis')), '', uiOutput(ns('direc'))
      selectInput(ns('methCovis'), label='Methods', choices=c('Annotation/manual'='man', 'Automatic'='auto'), selected='auto'), '',
      selectInput(ns('direc'), label='Mapping direction', choices=cho <- c('Cell2Bulk'='toBulk', 'Bulk2Cell'='toCell'), selected='toBulk'), '',
      actionButton(ns("covisHelp"), "Help", icon = icon('question-circle'), style='margin-top:24px')
    )),
    column(id=ns('covisMan'), width=12, covis_man_ui(ns('covisMan'))), 
    column(id=ns('covisAuto'), width=12, covis_auto_ui(ns('covisAuto')))

  ) # tabPanel("Single Cell",
}
