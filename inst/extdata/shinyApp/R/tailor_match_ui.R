# Module for tailoring co-clustering results.
tailor_match_ui <- function(id) { 
  ns <- NS(id)
  col1 <- column(6, id='dimCellBut', 
  fluidRow(splitLayout(cellWidths=c('1%', '27%', '1%', '20%', '1%', '5%'), '',
    selectInput(ns('dimCell'), label='Dimension reduction', choices=c("UMAP", "PCA", "TSNE")), '', 
    actionButton(ns("coclusPlotBut"), label='Co-visualizing', style='margin-bottom:-60px'), ''
  ))
  );
  col2 <- column(6, id='selBlkButCan',
  #div(style="overflow-x: auto; z-index:0",
  fluidRow(splitLayout(cellWidths=c('3%', '40%', '1%', '24%', '1%', '10%', '1%', '15%'), '', uiOutput(ns('selBlk')), '',
  actionButton(ns("selBlkBut"), label='Final confirmation', style='margin-bottom:-60px'), '',
  actionButton(ns("selBlkCancel"), label='Reset', style='margin-bottom:-60px'), '',
  actionButton(ns("tailorHelp"), "Help", style='margin-bottom:-60px', icon = icon('question-circle')), '',
  downloadButton(ns("dld"), "Download", style = "margin-top: 24px;")
  ))
  #)
  ) 
  list( # Inside the list: "," is the separator. Outside the list ";" is the separator.
  col1, col2,
  # if (hide==TRUE) hidden(col1), if (hide==TRUE) hidden(col2), 
  uiOutput(ns('dim.ui'))
  )
}
