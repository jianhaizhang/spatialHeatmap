
# The Shiny modules (e.g. search_server) are temporarily placed in this file only for debugging purpose, and will be moved to independent files in the R folder after the App development is completed.
options(list(stringsAsFactors=FALSE, shiny.fullstacktrace=TRUE, shiny.maxRequestSize='5G'))

# Every variable in every container should be checked at the beginning. E.g. input$fileIn in reactive({}). These checks will avoid amost all errors/warnings.

desired_bulk_server <- function(input, output, session) {
  observe({
    # validate(need(input$start.but>0, ''))
    withProgress(message="Loading dependencies: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    library(spatialHeatmap); library(SummarizedExperiment); library(shiny); library(shinydashboard); library(shinydashboardPlus); library(ggplot2); library(DT) 
    library(SingleCellExperiment); library(scater)
    incProgress(0.6, detail="please wait ...")
    library(gridExtra); library(grid); library(plotly); library(data.table) 
    # showModal(modal(title = HTML('<center><b>Welcome to spatialHeatmap!</b><center>'), msg = strong('Please wait for the landing page to finish!')))
    incProgress(0.9, detail="please wait...")
   library(shinyWidgets); library(shinyBS); library(shinyjs); library(htmltools)
  })
  })
  observeEvent(input$tailorHelp, {
    showModal(
    div(id='tailorHel', modalDialog(title= HTML('<strong><center>Instructions for Assigning Bulk Tissues</center></strong>'),
      div(style = 'overflow-y:scroll;overflow-x:scroll',
      HTML('
      <ol>
        <li>Mouse over the top of left embedding plot and select "Lasso Select".</li>
        <li>Select target cells on the embedding plot.</li>
        <li>Selected cells are listed on the right table. Select a desired bulk tissue from the dropdown menu.</li>
        <li>If more bulk tissues need to assign, repeat step 2 and 3.</li>
        <li>After all target cells are selected and assigned desired bulk tissues, click "Final confirmation".</li>
      </ol>
      ')
      ))))
    })
  sce <- reactive({
    sgl.cell.ipt <- input$sglCell
    # save(sgl.cell.ipt, file='sgl.cell.ipt')
    if (is.null(sgl.cell.ipt)) return()
    pa <- sgl.cell.ipt$datapath
    withProgress(message="Importing data: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    if (grepl('\\.rds$', pa)) readRDS(pa)
    })
  })

  observeEvent(sce(), ignoreInit=FALSE, {
    if (is.null(sce())) return()
    updateTabsetPanel(session, inputId="asgBulk", selected='asg')
  }) 

  output$coclus <- renderUI({
    sce <- sce(); if (is.null(sce)) return()
    cdat.na <- colnames(colData(sce))
    if (all(c('assignedBulk', 'similarity', 'cluster') %in% cdat.na)) {
      selectInput('coclus', 'Clusters', c('coclusters', sort(unique(sce$cluster))))
    }
  })
  dim.par <- reactiveValues(cocluster.only=TRUE)
  observe({
    sce <- sce(); input$scellRowCancelBut
    grp.sel <- input$coclus
    if (is.null(sce)|is.null(grp.sel)) return()
    cdat.na <- colnames(colData(sce)); row.sel$val <- NULL
    if ('coclusters' %in% grp.sel) {
      dim.par$group.sel <- NULL
      dim.par$cocluster.only <- TRUE
    } else dim.par$group.sel <- grp.sel
  })
  row.sel <- reactiveValues(val=NULL)
  observeEvent(input$scellRowBut, {
    row.sel$val <- input$scellCdat_rows_selected
    dim.par$group.sel <- NULL
  })
  observeEvent(input$scellRowCancelBut, { row.sel$val <- NULL})
  dimred <- reactiveValues()
  output$dimly <- renderPlotly({
    dimCell <- input$dimCell; grp.sel <- dim.par$group.sel
    sce <- sce(); cocluster.only <- dim.par$cocluster.only
    if (is.null(dimCell)|is.null(sce)|is.null(cocluster.only)) return()
    if (!is.null(row.sel$val) & !is.null(grp.sel)) return()
    withProgress(message="Plotting: ", value=0, {
      incProgress(0.3, detail="please wait ...")
      gg <- plot_dim(sce, dim=dimCell, color.by='cluster', group.sel=grp.sel, row.sel=row.sel$val, cocluster.only=cocluster.only)
      dimred$val <- gg
      ggplotly(gg, source='dim', tooltip=c('colour', 'x', 'y')) %>% event_register("plotly_selected")
    })
  })
  df.sel.cell <- reactiveValues(val=NULL, val1=NULL)
  observeEvent(input$selBlkCancel, {
    df.sel.cell$val1 <- df.sel.cell$val <- NULL 
  })
  observeEvent(input$selBulk, {
    cat('Desired bulk for selected cells ... \n')
    sel.blk <- input$selBulk; if (is.null(sel.blk)) return()
    sce <- sce(); if ('none' %in% sel.blk|is.null(sce)) return()
    cat(sel.blk, '\n'); if (!is(dimred$val, 'ggplot')) return()
    event.df <- event_data(event="plotly_selected", source='dim')
    if (!is(event.df, 'data.frame')) return()
    bulk <- subset(sce, , bulkCell=='bulk')
    event.df <- subset(event.df, !key %in% bulk$index)
    event.df <- event.df[, c('x', 'y', 'key')]
    df.val <- df.sel.cell$val
    if (!is.null(df.val)) df.val <- subset(df.val, !key %in% event.df$key)
    event.df$desiredBulk <- sel.blk
    event.df$dimred <- input$dimCell
    df.sel.cell$val <- rbind(event.df, df.val)
    showNotification(paste0('Registered for selected cells: ', sel.blk, '!'), duration=2)
    df.sel.cell$val$key <- as.numeric(df.sel.cell$val$key)
    df.sel.cell$val$dimred <- input$dimCell 
    # df.val <- df.sel.cell$val; save(df.val, file='df.val')
    cat('Done! \n')
  })
  # df.sel.cell$val1 is only used to trigger downstream changes in response to the two buttons. df.sel.cell$val changes every time input$selBulk changes while df.sel.cell$val1 changes only when the two buttons change.
  observeEvent(list(input$selBlkBut), {
    df.sel.cell$val1 <- df.sel.cell$val
  })
  observeEvent(input$scellRowCancelBut, { row.sel$val <- NULL })
  output$scellCdat <- renderDataTable({
    cat('Tailor: colData table ... \n')
    sce <- sce(); if (is.null(sce)) return()
    input$scellRowCancelBut # Deselect rows.
    cdat <- as.data.frame(colData(sce))
    # sf.by from dim_server in scell_server.
    cols <- list(list(targets=seq_len(ncol(cdat)), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    # dom='t' overwrites search box.
    # If "fixedColumns" is used, rows cannot be selected on the fixed part.
    withProgress(message="Rendering metadata table: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    tab <- datatable(cdat, selection=list(mode="multiple", target="row", selected='none'), escape=FALSE, filter="top", extensions=c('Scroller'), plugins = "ellipsis",
      options=list(pageLength=20, lengthMenu=c(10, 20, 50, 100), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=cols), 
      class='cell-border strip hover') %>% formatStyle(0, backgroundColor="white", cursor='pointer')
      cat('Done! \n'); tab
    })
    })

  output$selBlk <- renderUI({
    if (!is(dimred$val, 'ggplot')) return()
    input$selBlkCancel
    # New selection causes 'none'.
    event.df <- event_data(event="plotly_selected", source='dim')
    sce <- sce(); if (is.null(sce)) return()
    blk <- sort(unique(subset(sce, , bulkCell=='bulk')$sample))
    if (is.null(blk)) return()  
    selectInput('selBulk', label='Desired bulk for selected cells', choices=c('none', blk), selected='none')
  })

 output$selCellTab <- renderDataTable({
    sce <- sce(); sel.blk <- input$selBulk  
    if (is.null(sel.blk)|is.null(sce)) return()
    if (!is(dimred$val, 'ggplot')) return()
    event.df <- event_data(event="plotly_selected", source='dim')
    if (!is(event.df, 'data.frame')) return()
    event.df <- event.df[, c('x', 'y', 'key')]
    event.df$desiredBulk <- sel.blk
    event.df$cell <- sce$sample[as.numeric(event.df$key)]
    event.df <- event.df[, c('desiredBulk', 'cell', 'x', 'y', 'key')]
    if (is.null(event.df)) return()
    datatable(event.df, selection='none', escape=FALSE, filter="top", extensions=c('Scroller'), plugins = "ellipsis",
      options=list(pageLength=20, lengthMenu=c(10, 20, 50, 100), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=NULL), 
      class='cell-border strip hover') %>% formatStyle(0, backgroundColor="white", cursor='pointer')
  })

  output$dim.ui <- renderUI({
    row.but <- actionButton('scellRowBut', 'Confirm row selection', style='margin-top:24px')
    row.cancel.but <- actionButton('scellRowCancelBut', 'Deselect rows', style='margin-top:24px')
    lis <- list(
      column(6, 
      fluidRow(splitLayout(cellWidths=c('1%', '97%', '1%'), '',
      # renderUI's output is used in another renderUI on the server side. This section is equivalent to ui, so tsne/umap/pca/scellCdat can be used.
     plotlyOutput('dimly')
     ))
    ),
    column(6, dataTableOutput('selCellTab')),
    column(12, fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '15%', '1%', '10%'), '', row.but, '', row.cancel.but, '', uiOutput('coclus') 
    ))
    ),
    dataTableOutput('scellCdat')
    )
   # The same "row.but" can be used in scell and SHM sections, since these two sections use different ids to distinguish the same module.
   lis 
  })

  output$dld <- downloadHandler(
    filename=function() {
      paste0('selected_cells_with_desired_bulk.', 'txt') 
    },
    content=function(file) { cat("Downloading ... \n")
      df.cell <- df.sel.cell$val
      write.table(df.cell, file, col.names=TRUE, row.names=TRUE, sep='\t'); cat('Done! \n')
    }
  )

 onBookmarked(function(url) { updateQueryString(url) })
 # onBookmarked(updateQueryString)

}





