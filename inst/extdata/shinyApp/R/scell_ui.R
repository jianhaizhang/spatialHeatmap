# Module for co-visualization.
scell_ui <- function(id) { 
  ns <- NS(id)
  tabPanel(title=span('Co-visualization', title="This panel is designed for co-visualization only."), value='scell', icon=NULL,
    fluidRow(splitLayout(cellWidths=c('12px', '300px','1px', '125px'), '',
      # uiOutput(ns('methCovis')), '', uiOutput(ns('direc'))
      selectInput(ns('methCovis'), label='Cell group labels', choices=c('Annotation labels/Manual assignments'='man', 'Automated co-clustering'='auto'), selected='auto'), '',
      selectInput(ns('direc'), label='Mapping direction', choices=cho <- c('Cell-to-bulk'='toBulk', 'Bulk-to-cell'='toCell'), selected='toBulk')
    )),
    column(id=ns('covisMan'), width=12, covis_man_ui(ns('covisMan'))), 
    column(id=ns('covisAuto'), width=12, covis_auto_ui(ns('covisAuto')))
  ) # tabPanel("Single Cell",
}
