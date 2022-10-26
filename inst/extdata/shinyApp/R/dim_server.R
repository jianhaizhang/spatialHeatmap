# Module for plots of dimension reduction in co-visualization through annotation/manual methods.
dim_server <- function(id, sce, sce.upl, section='scell', upl.mod.lis, dat.lis=NULL, session) {
  moduleServer(id, function(input, output, session) {
   ns <- session$ns
   if (section!='scell') {
      hideElement('dimCell'); hideElement('coclusPlotBut')
      hideElement('selBlk'); hideElement('selBlkBut')
      hideElement('selBlkCancel')
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
      } else {
        lab.na <- grep('^label$|^label\\d+', cdat.na, value=TRUE) 
        if (length(lab.na)==0) return()
        selectInput(ns('covisGrp'), 'Cell groups', sort(lab.na))
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
      if (is.null(sce)|is.null(method)|is.null(grp)) return()
      if ('auto' %in% method & is.null(grp.sel)) return()
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
    output$umap <- renderPlot({
      grp <- dim.par$grp; grp.sel <- dim.par$group.sel
      sce <- sce(); cocluster.only <- dim.par$cocluster.only
      if (is.null(sce)|is.null(grp)|is.null(cocluster.only)) return()
      if (!is.null(row.sel$val) & !is.null(grp.sel)) return()
      plot_dim(sce, dim='UMAP', color.by=grp, group.sel=grp.sel, row.sel=row.sel$val, cocluster.only=cocluster.only)
    })
    output$tsne <- renderPlot({
      grp <- dim.par$grp; grp.sel <- dim.par$group.sel
      sce <- sce(); cocluster.only <- dim.par$cocluster.only
      if (is.null(sce)|is.null(grp)|is.null(cocluster.only)) return()
      if (!is.null(row.sel$val) & !is.null(grp.sel)) return()
      plot_dim(sce, dim='TSNE', color.by=grp, group.sel=grp.sel, row.sel=row.sel$val, cocluster.only=cocluster.only)
    })
    output$pca <- renderPlot({
      grp <- dim.par$grp; grp.sel <- dim.par$group.sel
      sce <- sce(); cocluster.only <- dim.par$cocluster.only
      if (is.null(sce)|is.null(grp)|is.null(cocluster.only)) return()
      if (!is.null(row.sel$val) & !is.null(grp.sel)) return()
      plot_dim(sce, dim='PCA', color.by=grp, group.sel=grp.sel, row.sel=row.sel$val, cocluster.only=cocluster.only)
    })

    output$scellCdat <- renderDataTable({
      cat('Single-cell: colData table ... \n')
      sce <- sce(); if (is.null(sce)) return()
      input$scellRowCancelBut # Deselect rows.
      cdat <- as.data.frame(colData(sce))
      # match.lis <- match.mod.lis$val$ft.reorder$ft.rematch
      # covisGrp from dim_server in scell_server. 
      if (section=='scell') covisGrp <- input$covisGrp else covisGrp <- dat.lis()$covisGrp
    cols <- list(list(targets=seq_len(ncol(cdat)), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    sel <- list(mode="multiple", target="row", selected='none')
    if (section!='scell') sel <- 'none'
      # The 1st column is "lable" or "cluster".
      # dom='t' overwrites search box.
      tab <- datatable(cdat, selection=sel, escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
      options=list(pageLength=20, lengthMenu=c(10, 20, 50, 100), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=cols, fixedColumns = list(leftColumns=2)), 
      class='cell-border strip hover') %>% formatStyle(0, backgroundColor="white", cursor='pointer')
      cat('Done! \n'); tab
    })

  output$dim.ui <- renderUI({
   cat('Manual matching: building ui of colData table ... \n')
   if (upl.mod.lis$ipt$fileIn!='customSingleCellData') return()
   row.but <- actionButton(ns('scellRowBut'), 'Confirm row selection', style='margin-top:24px')
   row.cancel.but <- actionButton(ns('scellRowCancelBut'), 'Deselect rows', style='margin-top:24px')
   if ('auto' %in% sce.upl$method) covis.but <- actionButton(ns('covisBut'), 'Co-visualizing', style='margin-top:24px') else covis.but <- NULL
   lis <- list(fluidRow(splitLayout(cellWidths=c('1%', '32%', '1%', '32%', '1%', '32%'), '',
      plotOutput(ns('pca')), '',plotOutput(ns('umap')), '',  plotOutput(ns('tsne'))
    )), br(),
    fluidRow(splitLayout(cellWidths=c('1%', '12%', '1%', '15%', '1%', '10%', '1%', '10%', '1%', '12%'), '', uiOutput(ns('samGrp')), '', row.but, '', row.cancel.but, '', uiOutput(ns('coclus')), '', covis.but
    )), br(),
    dataTableOutput(ns('scellCdat'))
    ); cat('Done! \n')
   if (section!='scell') dataTableOutput(ns('scellCdat')) else lis 
  })
  sfBy <- reactiveValues(); but.covis <- reactiveValues()
  observe({ sfBy$val <- input$covisGrp })
  observe({ but.covis$v <- input$covisBut })
  onBookmark(function(state) { state })
  return(list(covisGrp=sfBy, but.covis=but.covis))
})}
