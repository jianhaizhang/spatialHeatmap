# Module for co-visualization.
scell_ui <- function(id) { 
  ns <- NS(id)
  tabPanel(title=span(id=ns('covisTab'), "Co-visualization"), value='scell', icon=NULL,
     bsTooltip(ns('covisTab'), title = "This tab is inactivated if single-cell data are not provided.", placement="bottom", trigger="hover"),
    fluidRow(splitLayout(cellWidths=c('12px', '300px','1px', '125px'), '',
      # uiOutput(ns('methCovis')), '', uiOutput(ns('direc'))
      selectInput(ns('methCovis'), label='Cell group labels', choices=c('Annotation (or other) labels'='man', 'Co-clustering labels'='auto'), selected='auto'), '',
      selectInput(ns('direc'), label='Mapping direction', choices=cho <- c('Cell-to-bulk'='toBulk', 'Bulk-to-cell'='toCell'), selected='toBulk')
    )),
    column(id=ns('covisMan'), width=12, covis_man_ui(ns('covisMan'))), 
    column(id=ns('covisAuto'), width=12, covis_auto_ui(ns('covisAuto')))
  ) # tabPanel("Single Cell",
}
