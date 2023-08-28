# Module for plots of dimension reduction in co-visualization through annotation/manual methods.
dim_server <- function(id, sce, sce.upl, section='scell', upl.mod.lis, dat.lis=NULL, session) {
  moduleServer(id, function(input, output, session) {
   ns <- session$ns
   observeEvent(input$covisHelp, {
     method <- sce.upl$method; if (!check_obj(method)) return()
     if ('auto' %in% method) img <- 'coclus_quick.jpg' else img <- 'ann_quick.jpg' 
     showModal(modal(title = HTML('<center><b>Quick start!</b><center>'), msg = NULL, img=img, img.w="100%"))
    })
   if (section!='scell') {
      hideElement('dimCell'); hideElement('coclusPlotBut')
      hideElement('selBlk'); hideElement('selBlkBut')
      hideElement('selBlkCancel')
    } 
    grp <- reactiveValues()
    observeEvent(list(input$grpAuto2cell, input$grpAuto2blk, input$grpAnn, sce(), sce.upl$covis.type), {
      grpAuto2cell <- input$grpAuto2cell
      grpAuto2blk <- input$grpAuto2blk; grpAnn <- input$grpAnn
      req(any(check_obj(grpAuto2cell), check_obj(grpAuto2blk), check_obj(grpAnn))) 
      sce <- sce(); covis.type <- sce.upl$covis.type
      req(check_obj(list(sce, covis.type)))
      cdat.na <- colnames(colData(sce))
      if (all(c('assignedBulk', 'similarity') %in% cdat.na)) {
        if ('toCellAuto' %in% covis.type) {
          grp$v <- input$grpAuto2cell
        } else if ('toBulkAuto' %in% covis.type) {
          grp$v <- input$grpAuto2blk
        }
      } else grp$v <- input$grpAnn 
    })
    observeEvent(list(sce(), sce.upl$covis.type), {
      sce <- sce(); covis.type <- sce.upl$covis.type
      req(check_obj(list(sce, covis.type)))
      cdat.na <- colnames(colData(sce))
      if (all(c('assignedBulk', 'similarity') %in% cdat.na)) {
        shinyjs::hide(id = "grpAnnD")
        # covisGrp: bulk or cell labels for aggregating. 
        if ('toCellAuto' %in% covis.type) {
          shinyjs::show(id = "grpAuto2cellD")
          shinyjs::hide(id = "grpAuto2blkD")
        } else if ('toBulkAuto' %in% covis.type) {
          shinyjs::show(id = "grpAuto2blkD")
          shinyjs::hide(id = "grpAuto2cellD")
        }
      } else {
        lab.na <- grep('^label$|^label\\d+', cdat.na, value=TRUE) 
        if (length(lab.na)==0) return()
        shinyjs::show(id = "grpAnnD")
        shinyjs::hide(id = "grpAuto2cellD")
        shinyjs::hide(id = "grpAuto2blkD")
        shinyjs::hide(id = "coclusD")
        shinyjs::hide(id = "covisBut")
        updateSelectInput(session, "grpAnn", choices=sort(lab.na))
      }
    })
    output$samGrp <- renderUI({
      sce <- sce() 
      covis.type <- sce.upl$covis.type
      if (is.null(sce)|is.null(covis.type)) return()
      cdat.na <- colnames(colData(sce))
      if (all(c('assignedBulk', 'similarity') %in% cdat.na)) {
        # covisGrp: bulk or cell labels for aggregating. 
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
    observeEvent(list(sce(), sce.upl$method), {
      sce <- sce(); method <- sce.upl$method
      if (is.null(sce)| !'auto' %in% method) return()
      others <- sort(unique(sce$cluster))
      op <- c('all', others)
      updateSelectInput(session, "coclus", choices=op)
    })
    dim.par <- reactiveValues(cocluster.only=FALSE)
    observe({
      sce <- sce(); method <- sce.upl$method
      input$scellRowCancelBut
      grp <- grp$v; grp.sel <- input$coclus
      if (is.null(sce)|is.null(method)|is.null(grp)) return()
      if ('auto' %in% method & is.null(grp.sel)) return()
      dim.par$grp <- grp; cdat.na <- colnames(colData(sce))
      if (!grp %in% cdat.na) return(); row.sel$val <- NULL
      if ('auto' %in% method) { 
        dim.par$grp <- 'cluster'
        if ('all' %in% grp.sel) {
          dim.par$group.sel <- NULL
          dim.par$cocluster.only <- TRUE
        } else dim.par$group.sel <- grp.sel
      } else { 
        dim.par$group.sel <- NULL
        dim.par$cocluster.only <- FALSE
      }
    })
    row.sel <- reactiveValues(val=NULL)
    observeEvent(input$scellRowBut, {
      row.sel$val <- input$scellCdat_rows_selected
      dim.par$group.sel <- NULL
  })
    observeEvent(input$scellRowCancelBut, { row.sel$val <- NULL})
    output$dimPlot <- renderPlot({
      grp <- dim.par$grp; grp.sel <- dim.par$group.sel
      sce <- sce(); cocluster.only <- dim.par$cocluster.only
      dimMeth <- input$dimMeth
      if (!check_obj(list(sce, grp, dimMeth)) | is.null(cocluster.only)) return()
      if (!is.null(row.sel$val) & !is.null(grp.sel)) return()
      req(grp %in% colnames(colData(sce)))
      plot_dim(sce, dim=dimMeth, color.by=grp, group.sel=grp.sel, row.sel=row.sel$val, cocluster.only=cocluster.only, lgd.pos=ifelse(TRUE %in% cocluster.only, 'bottom', 'right'), lgd.l=ifelse(TRUE %in% cocluster.only, -0.07, 0), lgd.r=0.07)
    })

   observeEvent(sce.upl$method, {
     method <- sce.upl$method; req(check_obj(method))
     if ('auto' %in% method) {
       shinyjs::show(id='dtabCo') 
       shinyjs::hide(id='dtabAnn') 
     } else {
       shinyjs::show(id='dtabAnn') 
       shinyjs::hide(id='dtabCo') 
     }
    })
    output$scellCdat <- renderDataTable({
      cat('Single-cell: colData table ... \n')
      sce <- sce(); if (is.null(sce)) return()
      input$scellRowCancelBut # Deselect rows.
      cdat <- as.data.frame(colData(sce))
      # match.lis <- match.mod.lis$val$ft.reorder$ft.rematch
      # covisGrp from dim_server in scell_server. 
      if (section=='scell') covisGrp <- grp$v else covisGrp <- dat.lis()$covisGrp
      if ('auto' %in% sce.upl$method) {
        cna.sel <- c('cluster', 'bulkCell', 'assignedBulk', 'similarity', 'sample', 'index')
        req(all(cna.sel %in% colnames(cdat))); cdat <- cdat[, cna.sel]
      }
      cols <- list(list(targets=seq_len(ncol(cdat)), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
      sel <- list(mode="multiple", target="row", selected='none')
      if (section!='scell') sel <- 'none'
      # The 1st column is "lable" or "cluster".
      # dom='t' overwrites search box.
      tab <- datatable(cdat, selection=sel, escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
      options=list(pageLength=20, lengthMenu=c(10, 20, 50, 100), autoWidth=FALSE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=cols, fixedColumns = list(leftColumns=2)), 
      class='cell-border strip hover') %>% formatStyle(0, backgroundColor="white", cursor='pointer')
      cat('Done! \n'); tab
    })

  sfBy <- reactiveValues(); but.covis <- reactiveValues()
  observe({ sfBy$val <- grp$v })
  observe({ but.covis$v <- input$covisBut })
  onBookmark(function(state) { state })
  return(list(covisGrp=sfBy, but.covis=but.covis))
})}
