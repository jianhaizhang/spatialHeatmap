# Module for tailoring co-clustering results.
tailor_match_server <- function(id, sce, sce.upl, section='scell', upl.mod.lis, dat.lis=NULL, session) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns
  observeEvent(input$tailorHelp, {
    showModal(
    div(id = 'tailorIns',
    modalDialog(title = HTML('<strong><center>Tailoring Assignment Results (optional)</center></strong>'),
      div(style = 'overflow-y:scroll;overflow-x:scroll',
      HTML('<img src="image/tailoring.jpg">'),
    ))))
  })
    if (section!='scell') {
      # "hideElement" cannot remove space of the elements while "removeUI" can.
      removeUI(selector="#dimCellBut"); removeUI(selector="#selBlkButCan")
      # If these elements are rendered after calling "hideElement", no elements are hidden.
      # hideElement('dimCell'); hideElement('coclusPlotBut')
      # hideElement('selBlk'); hideElement('selBlkBut')
      # hideElement('selBlkCancel') 
    }
    output$samGrp <- renderUI({
      sce <- sce() 
      covis.type <- sce.upl$covis.type
      if (is.null(sce)|is.null(covis.type)) return()
      cdat.na <- colnames(colData(sce))
      if (all(c('assignedBulk', 'similarity') %in% cdat.na)) {
        if ('toCellAuto' %in% covis.type) {
          selectInput(ns('covisGrp'), 'Bulk groups', 'sample')
        } else if ('toBulkAuto' %in% covis.type) {
          selectInput(ns('covisGrp'), 'Cell groups', 'assignedBulk') 
        }
      }
    })
    output$coclus <- renderUI({
      sce <- sce(); method <- sce.upl$method
      if (is.null(sce)| !'auto' %in% method) return()
      selectInput(ns('coclus'), 'Clusters', c('coclusters', sort(unique(sce$cluster))))
    })
    dim.par <- reactiveValues(cocluster.only=FALSE)
    observe({
      sce <- sce(); method <- sce.upl$method
      input$scellRowCancelBut
      grp <- input$covisGrp; grp.sel <- input$coclus
      if (is.null(sce)|is.null(method)|is.null(grp)|is.null(grp.sel)) return()
      dim.par$grp <- grp; cdat.na <- colnames(colData(sce))
      if (!grp %in% cdat.na) return(); row.sel$val <- NULL
      if ('auto' %in% method) { 
        dim.par$grp <- 'cluster'
        if ('coclusters' %in% grp.sel) {
          dim.par$group.sel <- NULL
          dim.par$cocluster.only <- TRUE
        } else dim.par$group.sel <- grp.sel
      }
    })
    row.sel <- reactiveValues(val=NULL)
    observeEvent(input$scellRowBut, {
      row.sel$val <- input$scellCdat_rows_selected
      dim.par$group.sel <- NULL
  })
    observeEvent(input$scellRowCancelBut, { row.sel$val <- NULL})
    dimred <- reactiveValues()
    output$dimly <- renderPlotly({
      dimCell <- input$dimCell
      grp <- dim.par$grp; grp.sel <- dim.par$group.sel
      sce <- sce(); cocluster.only <- dim.par$cocluster.only
      if (is.null(dimCell)|is.null(sce)|is.null(grp)|is.null(cocluster.only)) return()
      if (!is.null(row.sel$val) & !is.null(grp.sel)) return()
      withProgress(message="Plotting: ", value=0, {
      incProgress(0.3, detail="please wait ...")
      gg <- plot_dim(sce, dim=dimCell, color.by=grp, group.sel=grp.sel, row.sel=row.sel$val, cocluster.only=cocluster.only)
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
    if (sel.blk=='none') return(); cat(sel.blk, '\n')
    if (!is(dimred$val, 'ggplot')) return()
    event.df <- event_data(event="plotly_selected", source='dim')
    if (!is(event.df, 'data.frame')) return()
    bulk <- sce.upl$bulk; if (is.null(bulk)) return()
    event.df <- subset(event.df, !key %in% bulk$index)
    event.df <- event.df[, c('x', 'y', 'key')]
    df.val <- df.sel.cell$val
    if (!is.null(df.val)) df.val <- subset(df.val, !key %in% event.df$key)
    event.df$desiredBulk <- sel.blk
    event.df$dimred <- input$dimCell
    df.sel.cell$val <- rbind(event.df, df.val)
    showNotification(paste0('Registered for selected cells: ', sel.blk, '!'), duration=2)
    df.sel.cell$val$key <- as.numeric(df.sel.cell$val$key)
    # df.val <- df.sel.cell$val; save(df.val, file='df.val')
    cat('Done! \n')
  })
  # df.sel.cell$val1 is only used to trigger downstream changes in response to the two buttons. df.sel.cell$val changes every time input$selBulk changes while df.sel.cell$val1 changes only when the two buttons change.
  observeEvent(list(input$selBlkBut), {
    df.sel.cell$val1 <- df.sel.cell$val
  })
  res.tailor <- reactiveValues()
  observeEvent(list(sce(), df.sel.cell$val1), ignoreNULL=FALSE, {
    cat('Auto-matching: tailoring assignments ... \n')
    sce <- sce(); df.sel.cell <- df.sel.cell$val1
    if (is.null(sce)) return()
    res.tailor$v <- refine_asg(sce.all=sce, df.desired.bulk=df.sel.cell)
    cat('Done! \n'); return(res.tailor)
  })

    output$scellCdat <- renderDataTable({
      cat('Tailor: colData table ... \n')
      sce <- res.tailor$v; if (is.null(sce)) return()
      input$scellRowCancelBut # Deselect rows.
      cdat <- as.data.frame(colData(sce))
      # covisGrp from dim_server in scell_server. 
      if (section=='scell') covisGrp <- input$covisGrp else covisGrp <- dat.lis()$covisGrp
      cols <- list(list(targets=seq_len(ncol(cdat)), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
      # The 1st column is "lable" or "cluster".
      # dom='t' overwrites search box.
      # If "fixedColumns" is used, rows cannot be selected on the fixed part.
      tab <- datatable(cdat, selection=list(mode="multiple", target="row", selected='none'), escape=FALSE, filter="top", extensions=c('Scroller'), plugins = "ellipsis",
      options=list(pageLength=20, lengthMenu=c(10, 20, 50, 100), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=250, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=cols), 
      class='cell-border strip hover') %>% formatStyle(0, backgroundColor="white", cursor='pointer')
      cat('Done! \n'); tab
    })

  output$selBlk <- renderUI({
    input$selBlkCancel
    if (!grepl(na.sgl, upl.mod.lis$ipt$fileIn)) return()
    # New selection causes 'none'.
    if (!is(dimred$val, 'ggplot')) return()
    event.df <- event_data(event="plotly_selected", source='dim')
    sce <- sce(); if (is.null(sce)) return()
    blk <- sort(unique(subset(sce, , bulkCell=='bulk')$sample))
    selectInput(ns('selBulk'), label='Desired bulk for selected cells', choices=c('none', blk), selected='none')
  })

 output$selCellTab <- renderDataTable({
    if (!grepl(na.sgl, upl.mod.lis$ipt$fileIn)) return()
    sce <- sce(); sel.blk <- input$selBulk; 
    if (is.null(sel.blk)|is.null(sce)) return()
    if (!is(dimred$val, 'ggplot')) return()
    event.df <- event_data(event="plotly_selected", source='dim')
    if (!is(event.df, 'data.frame')) return()
    event.df <- event.df[, c('x', 'y', 'key')]
    event.df$desiredBulk <- sel.blk
    event.df$cell <- colData(sce)$sample[as.numeric(event.df$key)]
    event.df <- event.df[, c('desiredBulk', 'cell', 'x', 'y', 'key')]
    if (is.null(event.df)) return()
    datatable(event.df, selection='none', escape=FALSE, filter="top", extensions=c('Scroller'), plugins = "ellipsis", 
    options=list(pageLength=20, lengthMenu=c(10, 20, 50, 100), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=250, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=NULL), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="white", cursor='pointer') 
  })

  output$dim.ui <- renderUI({
   cat('Tailoring: building ui of colData table ... \n')
   if (!grepl(na.sgl, upl.mod.lis$ipt$fileIn)) return()
   row.but <- actionButton(ns('scellRowBut'), 'Confirm row selection', style='margin-top:24px')
   row.cancel.but <- actionButton(ns('scellRowCancelBut'), 'Deselect rows', style='margin-top:24px')
   lis <- list(
     column(6, 
     fluidRow(splitLayout(cellWidths=c('1%', '97%', '1%'), '',
      # renderUI's output is used in another renderUI on the server side. This section is equivalent to ui, so tsne/umap/pca/scellCdat can be used.
      plotlyOutput(ns('dimly')), ''
    ))
    ),
    column(6, style='border-color:#3c8dbc;border-width:1px;border-style:solid;margin-bottom:5px;padding-top:2px;height:420px',
    dataTableOutput(ns('selCellTab'))),
    column(12, fluidRow(splitLayout(cellWidths=c('1%', '12%', '1%', '15%', '1%', '10%', '1%', '10%'), '', uiOutput(ns('samGrp')), '', row.but, '', row.cancel.but, '', uiOutput(ns('coclus'))
    ))
    ),
    dataTableOutput(ns('scellCdat'))
    ); cat('Done! \n')
   # The same "row.but" can be used in scell and SHM sections, since these two sections use different ids to distinguish the same module.
   if (section!='scell') {
     list(column(12, fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '15%'), '', row.but, '', row.cancel.but))),
     dataTableOutput(ns('scellCdat'))
     )
   } else lis 
  })
  output$dld <- downloadHandler(
    filename=function() {
      paste0('selected_cells_with_desired_bulk.', 'txt')
    },
    content=function(file) { cat("Downloading ... \n")
      df.cell <- df.sel.cell$val1
      write.table(df.cell, file, col.names=TRUE, row.names=TRUE, sep='\t'); cat('Done! \n')
    }
  )
  sfBy <- reactiveValues()
  observe({ sfBy$val <- input$covisGrp })
  coclusPlotBut <- reactive({ input$coclusPlotBut })
  selBlkBut <- reactive({ input$selBlkBut })
  selBlkCancel <- reactive({ input$selBlkCancel })
  onBookmark(function(state) { state })
  return(list(res.tailor=res.tailor, covisGrp=sfBy, coclusPlotBut=coclusPlotBut, selBlkBut=selBlkBut, selBlkCancel=selBlkCancel, df.sel.cell=df.sel.cell, row.sel=row.sel))
})}
