# Match spatial features between data and aSVG.
match_ui <- function(id) { 
  ns <- NS(id)
  list(
  column(12, fluidRow(splitLayout(cellWidths=c('10px', '400px', '1px', '90px', '1px', '80px', '1px', '70px'), '',
    uiOutput(ns('svgs'), style = 'margin-left:0px'), '',
    uiOutput(ns('match.but'), style = 'margin-top:24px'), '',
    uiOutput(ns('match.reset'), style = 'margin-top:24px'), '',
    actionButton(ns("matHelp"), "Help", icon = icon('question-circle'), style='margin-top:24px')
    ))), # verbatimTextOutput(ns('msg.match')),
  column(12, uiOutput(ns('ft.match')))
  )
}
