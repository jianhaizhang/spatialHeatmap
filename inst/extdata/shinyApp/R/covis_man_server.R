# Module for co-visualization through annotation/manual methods.
covis_man_server <- function(id, sce.upl, upl.mod.lis, shm.mod.lis, tab, covis.man, session) {
  moduleServer(id, function(input, output, session) {
  cat('Module covis_man_server ... \n')
   ns <- session$ns; 
   observe({
     # sce.man <- reactiveValues(bulk=sce.upl$bulk): if outside "observe", only runs one time when the app is started.
     covis.man$bulk <- sce.upl$bulk
     covis.man$method <- sce.upl$method
     covis.man$covis.type <- sce.upl$covis.type
   })
  observe({ # Bulk data table.
    sce <- sce.upl$bulk; if (is.null(sce)) return()
    withProgress(message="Bulk data table in co-visualization: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    norm.meth <- input$normBlk
    if (is.null(norm.meth)) norm.meth <- 'VST'
    if (grepl('^CNF-', norm.meth)) {
      norm.fun <- 'CNF'
      par.lis <- list(method=sub('.*-', '', norm.meth))
    } else { norm.fun <- norm.meth; par.lis <- NULL }
    covis.man$bulk <- sce <- norm_data(data=sce, norm.fun=norm.fun, parameter.list=par.lis, log2.trans=TRUE)
    output$datCovisBlk <- renderDataTable({
      dat_covis_man(sce, nr=1000, nc=100)
    })
   })
  })

  par.qc <- reactiveValues()
  # Button of 0 cannot trigger observeEvent.
  observeEvent(input$parManBut+1, { 
    par.qc$cnt.thr <- input$cntThr; par.qc$nmads <- input$nmads
  })
  observeEvent(list(sce.upl$cell, par.qc$cnt.thr, par.qc$nmads), {
    cat('Single cell: quality control ... \n')
    cnt.thr <- par.qc$cnt.thr; nmads <- par.qc$nmads
    cell <- sce.upl$cell
    if (is.null(cell)|!is.numeric(cnt.thr)|!is.numeric(nmads)) return()
    withProgress(message="Quality control: ", value=0, {
      incProgress(0.3, detail="in progress ...")
      covis.man$cell.qc <- qc_cell(sce=cell, qc.metric=list(subsets=list(Mt=rowData(cell)$featureType=='mito'), threshold=cnt.thr), qc.filter=list(nmads=nmads))
    }); cat('Done! \n')
  })
  par.norm <- reactiveValues(min.size=100, max.size=3000)
  observeEvent(input$parManBut+1, {
    par.norm$min.size <- input$minSize
    par.norm$max.size <- input$maxSize
  })
  observeEvent(list(covis.man$cell.qc, par.norm$min.size, par.norm$max.size), {
    cat('Single cell: nomalization ... \n')
    cell <- covis.man$cell.qc; if (is.null(cell)) return()
    min.size <- par.norm$min.size; max.size <- par.norm$max.size
    validate(need(round(min.size)==min.size, ''))
    validate(need(round(max.size)==max.size, ''))
    withProgress(message="Normalizing: ", value=0, {
      incProgress(0.3, detail="in progress ...")
      covis.man$cell.nor <- norm_cell(cell, com.sum.fct = list(max.cluster.size = max.size, min.mean = 1), quick.clus = list(min.size = min.size), com = FALSE)
      # covis.man$cell.qc <- NULL; 
      cat('Done! \n') 
    })
  })
  par.var <- reactiveValues(hvg.n=3000, hvg.p=0.1)
  observeEvent(input$parManBut, {
    par.var$hvg.n <- input$hvgN; par.var$hvg.p <- input$hvgP
  })

  observe({
    cat('Top HVGs ... \n')
    cell <- covis.man$cell.nor; if (is.null(cell)) return()
    n <- par.var$hvg.n; p <- par.var$hvg.p
    validate(need(round(n)==n, '')); validate(need(p > 0 & p < 1, ''))
    withProgress(message="Variance modelling: ", value=0, {
    incProgress(0.3, detail="in progress ...")
    if (!'logcounts' %in% assayNames(cell)) {
      covis.man$var <- covis.man$top <- NULL; return()
    }
    covis.man$var <- sce.var <- modelGeneVar(cell)
    incProgress(0.3, detail="in progress ...")
    covis.man$top <- top.hvgs <- getTopHVGs(sce.var, prop=p, n=n)
    cat('Done! \n')
    })
  })

  par.dim <- reactiveValues()
  # Button of 0 cannot trigger observeEvent.
  observeEvent(input$parManBut+1, { 
    par.dim$min.rank <- input$minRank
    par.dim$max.rank <- input$maxRank
    par.dim$ncomT <- input$ncomT; par.dim$ntopT <- input$ntopT
    par.dim$steps <- input$steps
    par.dim$ncomU <- input$ncomU; par.dim$ntopU <- input$ntopU
    par.dim$pcs <- input$pcs
  })
  observe({
    cat('Single cell: reducing dimensions ... \n'); 
    cell <- covis.man$cell.nor; sce.var <- covis.man$var; top <- covis.man$top
    min.rank <- par.dim$min.rank; max.rank <- par.dim$max.rank
    ncomT <- par.dim$ncomT; ntopT <- par.dim$ntopT; steps <- par.dim$steps
    ncomU <- par.dim$ncomU; ntopU <- par.dim$ntopU; pcs <- par.dim$pcs
    if (is.null(cell)|is.null(sce.var)|is.null(top)) return()
    validate(need(round(min.rank)==min.rank & round(max.rank)==max.rank & max.rank > min.rank, ''))
    validate(need(round(ncomT)==ncomT & round(ntopT)==ntopT, ''))
    validate(need(round(ncomU)==ncomU & round(ntopU)==ntopU & round(pcs)==pcs, ''))
    output$msg.umap <- renderText({
      if(pcs < ncomU) return('In "runUMAP", "Number of PCs" must be >= "Number of dimensions to obtain"!') else return()
    })
    validate(need(pcs >= ncomU, 'In "runUMAP", "Number of PCs" must be >= "Number of dimensions"!'))
    withProgress(message="Reducing dimensions: ", value=0, {
    incProgress(0.3, detail="in progress ...")
    cell <- denoisePCA(cell, technical=sce.var, subset.row=top, min.rank=min.rank, max.rank=max.rank)
    cell <- runTSNE(cell, dimred="PCA", ncomponents=ncomT, ntop=ntopT)
    incProgress(0.3, detail="in progress ...")
    # Avoid warnings due to duplicated column names.
    cna <- colnames(cell)
    colnames(cell) <- seq_len(ncol(cell))
    cell <- runUMAP(cell, dimred="PCA", ncomponents=ncomU, ntop=ntopU, pca=pcs)
    # Row names in colData, reducedDim change accordingly.
    colnames(cell) <- cna
    })
    # covis.man$cell.nor <- NULL
    covis.man$dimred <- cell; cat('Done! \n')
  })
   observeEvent(covis.man$dimred, ignoreInit=FALSE, {
     updateTabsetPanel(session, inputId="tabSetCell", selected='dimred')
  })
  observe({ # Cell data table.
    sce <- covis.man$dimred; if (is.null(sce)) return()
    withProgress(message="Cell data table in co-visualization: ", value=0, {
    incProgress(0.3, detail=" ...")
    output$datCell <- renderDataTable({
      dat_covis_man(sce, nr=1000, nc=100)
    })
   })
  })

  observeEvent(covis.man$dimred, ignoreInit=FALSE, ignoreNULL=FALSE, {
    covis.type <- sce.upl$covis.type
    if (is.null(covis.type)) return()
    if (covis.type %in% c('toBulkAuto', 'toCellAuto')) {
        showTab('tabSetCell', target='autoMatch')
        hideTab('tabSetCell', target='qcTab')
        hideTab('tabSetCell', target='norm')
        hideTab('tabSetCell', target='varMol')
        hideTab('dimredNav', target='dimredPar')
        hideTab('tabSetCell', target='manualMatch')
    } else {
        hideTab('tabSetAuto', target='autoMatch')
        showTab('tabSetCell', target='qcTab')
        showTab('tabSetCell', target='norm')
        showTab('tabSetCell', target='varMol')
        showTab('dimredNav', target='dimredPar')
        showTab('tabSetCell', target='manualMatch')
    }
  })

  sce.dimred <- reactive({ covis.man$dimred })
  # Generates both embedding plots and the colData table in the scell page.
  dim.lis <- reactive({
    dim_server('dim', sce=sce.dimred, sce.upl=sce.upl, section='scell', upl.mod.lis=upl.mod.lis) 
  }) # Avoid endless circle.
  observe({
    covis.man$covisGrp <- dim.lis()$covisGrp$val
  })
   observe({ # The id 'rematchCell' is fixed, since it is recognised internally. 
     if (upl.mod.lis$ipt$fileIn=='customSingleCellData') covis.man$match.mod.lis <- match_server('rematchCell', shm.mod.lis$sam, tab, upl.mod.lis, covis.man=covis.man) else covis.man$match.mod.lis <- NULL
   })
  observe({
   meth <- input$scell.cluster; if (is.null(meth)) return()
   output$clus.par <- renderUI({
     ns <- session$ns
     if (meth=='cluster_walktrap') fluidRow(splitLayout(cellWidths = c('1%', '20%'), '', numericInput(ns('steps'), label='Random walks', value=4, min=1, max=Inf, step=1, width=150))) 
   })
  })
  cat('Module covis_man_server done !\n')
  onBookmark(function(state) { state })
  # return(list(match.mod.res=match.mod.res))
})}
