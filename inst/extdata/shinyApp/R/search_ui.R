# Module for searching gene ids.
# A module can only have a "ui" element without the "server".
search_ui <- function(id, lab=label) {
  # Applying "ns" on an id: this id is prefixed with the module id in HTML output.
  ns <- NS(id) 
  fluidRow(splitLayout(cellWidths=c('10px', '100px', '15px', '600px'), '',
  radioButtons(inputId=ns('sch.mode'), label='Auto-completion', choices=c('On'='Single', 'Off'='Multiple'), selected='Multiple', inline=TRUE), '',
  uiOutput(ns('sch.box'))
  ))
}
