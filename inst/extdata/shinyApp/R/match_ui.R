# Match spatial features between data and aSVG.
match_ui <- function(id) { 
  ns <- NS(id)
  list(
  column(12, fluidRow(splitLayout(cellWidths=c('1%', '40%', '5%', '10%', '1%', '10%', '1%', '5%'), '',
    uiOutput(ns('svgs'), style = 'margin-left:5px'), '',
    uiOutput(ns('match.but'), style = 'margin-top:23px'), '',
    uiOutput(ns('match.reset'), style = 'margin-top:23px'), '',
    actionButton(ns("matHelp"), "Help", icon = icon('question-circle'), style='margin-top:24px')
    ))), verbatimTextOutput(ns('msg.match')),
  column(12, uiOutput(ns('ft.match')))
  )
}
