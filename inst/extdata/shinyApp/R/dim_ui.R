# Module for plots of dimension reduction in co-visualization through annotation/manual methods.
dim_ui <- function(id, coclus=NULL) { 
  ns <- NS(id);
  msg.plot <- 'Embedding plot of bulk and single-cell data after co-clustering.' 
  msg.meta <- as.character(HTML(msg.meta.ann))
  if (coclus==TRUE) { 
    msg.plot <- 'Embedding plot of the single-cell data.' 
    msg.meta <- paste0(as.character(HTML(paste0(msg.meta.coclus, "<br/>"))), as.character(a(href='html/shm_shiny_manual.html#coclus', target='blank', 'More.'))) 
  }
  list(
  div(class='shiny-split-layout', 
   div(selectInput(ns('dimMeth'), label='Dimension reduction', choices=c('TSNE', 'UMAP', 'PCA'), selected='TSNE'), style='width:150px'),
   div(id=ns('grpAuto2cellD'), selectInput(ns('grpAuto2cell'), 'Bulk groups', 'sample'), style='width:135px'),
   div(id=ns('grpAuto2blkD'), selectInput(ns('grpAuto2blk'), 'Cell groups', 'assignedBulk'), style='width:135px'),
   div(id=ns('grpAnnD'), selectInput(ns('grpAnn'), 'Cell groups', 'label'), style='width:135px'),
   actionButton(ns('scellRowBut'), 'Confirm row selection', style='margin-top:24px;width:160px'),
   actionButton(ns('scellRowCancelBut'), 'Deselect rows', style='margin-top:24px;width:112px'),
   div(id=ns('coclusD'), selectInput(ns('coclus'), 'Clusters to show', choices='all'), style='width:110px'),
   actionButton(ns('covisBut'), 'Co-visualizing', style=paste0(run.top, 'width:111px')),
   actionButton(ns("covisHelp"), "Help", icon = icon('question-circle'), style='margin-top:24px;width:71px')
   ),
   bsTooltip(ns('grpAuto2cell'), title = 'Replicates in bulk tissues are grouped by tissue labels.', placement = "top", trigger = "hover"),
   bsTooltip(ns('grpAuto2blk'), title = 'Tissue labels assigned to cells as group labels through co-clustering.', placement = "top", trigger = "hover"),
   bsTooltip(ns('grpAnnD'), title = 'Cell group labels obtained from annotation labels, marker genes, etc.', placement = "top", trigger = "hover"),
   bsTooltip(ns('scellRowBut'), title = 'Click to visualize selected cells in the table.', placement = "top", trigger = "hover"),
   bsTooltip(ns('coclusD'), title = 'Show all or a certain cluster (may contain only cells or cells and tissues)', placement="top", trigger = "hover"),
   bsTooltip(ns('covisBut'), title = 'Click for co-visualization plots.', placement="top", trigger = "hover"),
   fluidRow(splitLayout(cellWidths=c('1%', '30%', '1%', '70%'), '', 
     plotOutput(ns('dimPlot')) %>% spsComps::bsTooltip(title=msg.plot, placement='top', html=TRUE, click_inside=TRUE), '', 
     div(dataTableOutput(ns('scellCdat'))) %>% spsComps::bsTooltip(title=msg.meta, placement='left', html=TRUE, click_inside=TRUE)
   )) 
  )
}
