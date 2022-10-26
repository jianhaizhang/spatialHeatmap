# Module for plots of dimension reduction in co-visualization through annotation/manual methods.
dim_ui <- function(id) { 
  ns <- NS(id); uiOutput(ns('dim.ui'))
  #list(
  #  column(12,
  #  fluidRow(splitLayout(cellWidths=c('1%', '48%', '1%', '48%'), '',
  #    plotOutput(ns('tsne')), '', plotOutput(ns('pca'))
  #  )), dataTableOutput(ns('scell.cdat'))
  #)
  #)
}
