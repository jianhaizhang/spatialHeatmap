# Module for tailoring co-clustering results.
tailor_match_ui <- function(id) { 
  ns <- NS(id)
  col1 <- column(6, id='dimCellBut', 
  fluidRow(splitLayout(cellWidths=c('1px', '150px', '1px', '115px', '1px'), '',
    selectInput(ns('dimCell'), label='Dimension reduction', choices=c("UMAP", "PCA", "TSNE")), '', 
    actionButton(ns("coclusPlotBut"), label='Co-visualizing', style='margin-top:24px;color:#fff;background-color:#c685c4;border-color:#ddd'), ''
  ))
  );
  col2 <- column(6, id='selBlkButCan', style=rec,
  #div(style="overflow-x: auto; z-index:0",
  fluidRow(splitLayout(cellWidths=c('1px', '205px', '1px', '135px', '1px', '65px', '1px', '105px', '1px', '70px'), '', 
  uiOutput(ns('selBlk')), '',
  actionButton(ns("selBlkBut"), label='Final confirmation', style='margin-bottom:-60px'), '',
  actionButton(ns("selBlkCancel"), label='Reset', style='margin-bottom:-60px'), '',
  downloadButton(ns("dld"), "Download", style = "margin-top: 24px;"), '',
  actionButton(ns("tailorHelp"), "Help", style='margin-bottom:-60px', icon = icon('question-circle'))
  ))
  #)
  ) 
  list( # Inside the list: "," is the separator. Outside the list ";" is the separator.
  col1, col2,
  # if (hide==TRUE) hidden(col1), if (hide==TRUE) hidden(col2), 
  uiOutput(ns('dim.ui'))
  )
}
