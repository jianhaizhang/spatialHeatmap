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
    output$coclus <- renderUI({
      sce <- sce(); method <- sce.upl$method
      if (is.null(sce)| !'auto' %in% method) return()
      others <- sort(unique(sce$cluster))
      op <- setNames(c('coclusters', others), c('all', others))
      selectInput(ns('coclus'), 'Clusters to show', op)
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
    output$dimPlot <- renderPlot({
      grp <- dim.par$grp; grp.sel <- dim.par$group.sel
      sce <- sce(); cocluster.only <- dim.par$cocluster.only
      dimMeth <- input$dimMeth
      if (!check_obj(list(sce, grp, dimMeth)) | is.null(cocluster.only)) return()
      if (!is.null(row.sel$val) & !is.null(grp.sel)) return()  
      plot_dim(sce, dim=dimMeth, color.by=grp, group.sel=grp.sel, row.sel=row.sel$val, cocluster.only=cocluster.only, lgd.pos=ifelse(TRUE %in% cocluster.only, 'bottom', 'right'), lgd.l=ifelse(TRUE %in% cocluster.only, -0.07, 0), lgd.r=0.07)
    })

    output$scellCdat <- renderDataTable({
      cat('Single-cell: colData table ... \n')
      sce <- sce(); if (is.null(sce)) return()
      input$scellRowCancelBut # Deselect rows.
      cdat <- as.data.frame(colData(sce))
      # match.lis <- match.mod.lis$val$ft.reorder$ft.rematch
      # covisGrp from dim_server in scell_server. 
      if (section=='scell') covisGrp <- input$covisGrp else covisGrp <- dat.lis()$covisGrp
      if ('auto' %in% sce.upl$method) cdat <- cdat[, c('cluster', 'bulkCell', 'assignedBulk', 'similarity', 'sample', 'index')]
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

  output$dim.ui <- renderUI({
   cat('Manual matching: building ui of colData table ... \n')
   Sys.sleep(2) # Nested renderUI: if the inner takes longer time than the outer, the outer will not capture the inner output.
   if (!grepl(na.sgl, upl.mod.lis$ipt$fileIn)) return()
   dimMeth <- selectInput(ns('dimMeth'), label='Dimension reduction', choices=c('TSNE', 'UMAP', 'PCA'), selected='TSNE')
   row.but <- actionButton(ns('scellRowBut'), 'Confirm row selection', style='margin-top:24px')
   row.cancel.but <- actionButton(ns('scellRowCancelBut'), 'Deselect rows', style='margin-top:24px')
   covis.but <- NULL; co.lab <- co.but.len <- '0px'
   help.but <- actionButton(ns("covisHelp"), "Help", icon = icon('question-circle'), style='margin-top:24px')
   msg.grp <- 'Cell group labels obtained from annotation labels, marker genes, etc.'
   msg.plot <- 'Embedding plot of the single-cell data.'
   msg.meta <- as.character(HTML(msg.meta.ann))
   if ('auto' %in% sce.upl$method) { 
     covis.but <- actionButton(ns('covisBut'), 'Co-visualizing', style=run.top)
     co.lab <- '110px'; co.but.len <- '111px'
     msg.grp <- 'Tissue labels assigned to cells as group labels through co-clustering.'
     msg.plot <- 'Embedding plot of bulk and single-cell data after co-clustering.'
     msg.meta <- paste0(as.character(HTML(paste0(msg.meta.coclus, "<br/>"))), as.character(a(href='html/shm_shiny_manual.html#coclus', target='blank', 'More.')))
   }
   lis <- list(
     fluidRow(splitLayout(cellWidths=c('10px', '150px', '1px', '135px', '1px', '160px', '1px', '112px', '1px', co.lab, '1px', co.but.len, '1px', '71px'), '', 
     dimMeth, '', uiOutput(ns('samGrp')), '', row.but, '', row.cancel.but, '', uiOutput(ns('coclus')), '', covis.but, '', help.but
   )), # div(style='margin-top:10px'),
     bsTooltip(ns('samGrp'), title = msg.grp, placement = "top", trigger = "hover"),
     bsTooltip(ns('scellRowBut'), title = 'Click this button to visualize selected cells in the table.', placement = "top", trigger = "hover"),
     bsTooltip(ns('coclus'), title = 'Show all or a certain cluster (may contain only cells or cells and tissues)', placement = "top", trigger = "hover"),
     fluidRow(splitLayout(cellWidths=c('1%', '30%', '1%', '70%'), '',
     plotOutput(ns('dimPlot')), '', div(dataTableOutput(ns('scellCdat'))) %>% spsComps::bsTooltip(title=msg.meta, placement='left', html=TRUE, click_inside=TRUE)
     )),
     bsTooltip(ns('dimPlot'), title = msg.plot, placement = "top", trigger = "hover")
     # bsTooltip(ns('scellCdat'), title = msg.meta, placement = "left", trigger = "hover")
    ); cat('Done! \n')
   if (section!='scell') dataTableOutput(ns('scellCdat')) else lis 
  })
  sfBy <- reactiveValues(); but.covis <- reactiveValues()
  observe({ sfBy$val <- input$covisGrp })
  observe({ but.covis$v <- input$covisBut })
  onBookmark(function(state) { state })
  return(list(covisGrp=sfBy, but.covis=but.covis))
})}
