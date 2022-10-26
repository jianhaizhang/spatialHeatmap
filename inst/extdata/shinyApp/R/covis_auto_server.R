# Module for co-visualization through automatic method.
covis_auto_server <- function(id, sce.upl, upl.mod.lis, shm.mod.lis, tab, covis.auto, session) {
  moduleServer(id, function(input, output, session) {
   ns <- session$ns 
   observe({
     covis.auto$method <- sce.upl$method
     covis.auto$covis.type <- sce.upl$covis.type
   })

  # Normalize bulk and cells.
  norm.coclus <- reactiveValues()
  observeEvent(input$parAutoBut+1, {
    if (is.null(sce.upl$bulk)) return()
    norm.coclus$val <- list(input$normCoclus)
  })
  # eventReactive avoids endless circles.
  dat.nor <- eventReactive(list(norm.coclus$val, sce.upl$bulk, sce.upl$cell), {
    cat('Auto-matching: normalizing ... \n')
    norm.meth <- input$normCoclus
    bulk <- sce.upl$bulk; cell <- sce.upl$cell
    if (is.null(norm.meth)|is.null(bulk)|is.null(cell)) return()
    withProgress(message="Normalizing: ", value=0, {
    incProgress(0.3, detail="in progress ...")
    dat.nor <- norm_cell(sce=cell, bulk=bulk, cpm=ifelse(norm.meth=='cpm', TRUE, FALSE)); cat('Done! \n'); return(dat.nor)
    })
  }) 
 
  # Filtering.
  par.fil.coclus <- reactiveValues()
  observeEvent(input$parAutoBut+1, {
    if (is.null(sce.upl$bulk)) return()
    par.fil.coclus$val <- list(input$filBlkP, input$filBlkA, input$filBlkCV1, input$filBlkCV2, input$filPGen, input$filPCell)
  })
  # eventReactive avoids endless circles.
  dat.fil <- eventReactive(list(par.fil.coclus$val, dat.nor()), {
    cat('Auto-matching: filtering ... \n')
    library(SingleCellExperiment)
    dat.nor <- dat.nor(); if (is.null(dat.nor)) return() 
    blk <- dat.nor$bulk; cell <- dat.nor$cell
    filBlkP <- input$filBlkP; filBlkA <- input$filBlkA
    filBlkCV1 <- input$filBlkCV1
    filBlkCV2 <- input$filBlkCV2
    filPGen <- input$filPGen; filPCell <- input$filPCell
    validate(need((filBlkP >=0 & filBlkP <=1) & (filBlkA >= 0) & (filPGen >=0 & filPGen <=1) & (filPCell >=0 & filPCell <=1) & (filBlkCV1 >= -1000 & filBlkCV1 <= 1000 & filBlkCV2 >= -1000 & filBlkCV2 <= 1000 & filBlkCV1 < filBlkCV2), ''))
    withProgress(message="Filtering: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    blk.aggr <- aggr_rep(data=blk, assay.na='logcounts', sam.factor='sample', aggr='mean')
    blk.fil <- filter_data(data=blk.aggr, pOA=c(filBlkP, filBlkA), CV=c(filBlkCV1, filBlkCV2), verbose=FALSE)
    incProgress(0.3, detail="please wait ...")
    dat.fil <- filter_cell(sce=cell, bulk=blk.fil, gen.rm=NULL, cutoff=1, p.in.cell=filPCell, p.in.gen=filPGen, verbose=FALSE)
    cat('Done! \n'); return(dat.fil)
    })
  })

  observe({
    withProgress(message="Data table in co-visualization: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    dat.fil <- dat.fil()
    bulk <- dat.fil$bulk; cell <- dat.fil$cell
    if (is.null(bulk)|is.null(cell)) return()
    output$datCovisBlk <- renderDataTable({
      dat_covis_man(bulk, nr=1000, nc=100)
    })
    output$datCovisCell <- renderDataTable({
      dat_covis_man(cell, nr=1000, nc=100)
    })
   })
  })

  # Dimension reduction.
  par.dim <- reactiveValues()
  # Button of 0 cannot trigger observeEvent.
  observeEvent(input$parAutoBut+1, {
    if (is.null(dat.fil())) return()
    par.dim$val <- list(input$minRank, input$maxRank)
  })
  dimred <- eventReactive(list(par.dim$val, dat.fil()), {
    cat('Auto-matching: reducing dimensions ... \n')
    dat.fil <- dat.fil()
    bulk <- dat.fil$bulk; cell <- dat.fil$cell
    minRank <- input$minRank; maxRank <- input$maxRank
    if (is.null(cell)|is.null(bulk)|is.null(minRank)|is.null(maxRank)) return()
    validate(need(round(minRank)==minRank & round(maxRank)==maxRank & maxRank > minRank, ''))
    inter <- intersect(rownames(bulk), rownames(cell))
    blk.kp <- bulk[inter, ]; sc.kp <- cell[inter, ]
    com.kp <- cbind(blk.kp, sc.kp)
    com.kp$index <- seq_len(ncol(com.kp))
    com.kp$sample <- colnames(com.kp)
    withProgress(message="Reducing dimensions: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    dimred <- reduce_dim(sce = com.kp, min.dim = minRank, max.dim = maxRank); cat('Done! \n'); return(dimred)
    })
  })
  
  # Bulid graph.
  par.graph <- reactiveValues()
  # Button of 0 cannot trigger observeEvent.
  observeEvent(input$parAutoBut+1, {
    if (is.null(dimred())) return()
    par.graph$val <- list(input$dimSel, input$graphMeth)
  }) 
  # eventReactive avoids endless circles.
  gr <- eventReactive(list(par.graph$val, dimred()), {
    cat('Auto-matching: building graph ... \n')
    dimSel <- input$dimSel; graphMeth <- input$graphMeth
    dimred <- dimred();
    if (is.null(dimred)|is.null(graphMeth)|is.null(dimSel)) return()
    withProgress(message="Buliding graph : ", value=0, {
    incProgress(0.3, detail="please wait ...")
    if (graphMeth == "knn") {
      gr.sc <- do.call(buildKNNGraph, list(x = dimred, use.dimred = dimSel)) 
    } else if (graphMeth == "snn") {
      gr.sc <- do.call(buildSNNGraph, list(x = dimred, use.dimred = dimSel))
    }; cat('Done! \n'); return(gr.sc)
    })
  }) 

  # Detect clusters.
  par.coclus <- reactiveValues()
  observeEvent(input$parAutoBut+1, {
    if (is.null(gr())) return()
    par.coclus$val <- list(input$clusMeth)
  })
  # eventReactive avoids endless circles.
  coclus <- eventReactive(list(par.coclus$val, gr()), {
    cat('Auto-matching: detecting clusters ... \n')
    gr <- gr(); clusMeth <- input$clusMeth; dimred <- dimred()
    if (is.null(gr)|is.null(clusMeth)|is.null(dimred)) return()
    withProgress(message="Detecting clusters: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    clus.all <- detect_cluster(graph = gr, clustering = clusMeth)
    clus <- as.character(clus.all$membership)
    clus <- paste0("clus", clus); cdat.sc <- colData(dimred)
    rna <- rownames(cdat.sc)
    lab.lgc <- "label" %in% make.names(colnames(cdat.sc))
    if (lab.lgc) {
      cdat.sc <- cbind(cluster = clus, cdat.sc)
      idx <- colnames(cdat.sc) %in% c("cluster", "label")
      cdat.sc <- cdat.sc[, c(which(idx), which(!idx))]
    } else cdat.sc <- cbind(cluster = clus, cdat.sc)
    rownames(cdat.sc) <- rna
    colnames(cdat.sc) <- make.names(colnames(cdat.sc))
    colData(dimred) <- cdat.sc
    cat('Done! \n'); return(dimred)
    })
  }) 
  # Assign bulk.
  asg <- eventReactive(list(coclus()), ignoreNULL=FALSE, {
    cat('Auto-matching: assigning bulk tissues ... \n')
    coclus <- coclus(); dimSel <- input$dimSel
    if (is.null(coclus)|is.null(dimSel)) return()
    asg <- com_roc(sce.coclus=coclus, dimred=dimSel, dat.blk = subset(coclus, , bulkCell=='bulk'))
    cat('Done! \n'); return(asg)
  })

  res <- eventReactive(list(asg(), input$asgThr), ignoreNULL=FALSE, {
    cat('Auto-matching: subsetting assignments ... \n')
    # res.coclus() should not be replaced with res.lis$val. The former keeps original coclustering results, and once the desired bulk is reset, the former will be used. If replaced, the reset button is useless, since the original coclustering results already included desired bulk and can not be reverted.
    asg <- asg(); asgThr <- input$asgThr
    if (is.null(asg)|is.null(asgThr)) return()
    validate(need(asgThr >=-Inf & asgThr <=1, ''))
    covis.auto$res <- res <- filter_asg(asg, min.sim=asgThr)
    cat('Done! \n'); return(res)
  })
   observeEvent(res(), ignoreInit=FALSE, {
     updateTabsetPanel(session, inputId="tabSetCellAuto", selected='result')
  })

  tailor.lis <- reactiveValues()
  observe({
    tailor.lis$v <- tailor_match_server('tailor', sce=res, sce.upl=sce.upl, section='scell', upl.mod.lis=upl.mod.lis)
  }) # Avoid endless circle.
  observeEvent(tailor.lis$v$res.tailor$v, {
   # if (!is.null(tailor.lis$v$res.tailor$v)) print(list('tailor.v', colData(tailor.lis$v$res.tailor$v)[c(156, 256), c('assignedBulk', 'similarity')]))
    covis.auto$res <- tailor.lis$v$res.tailor$v
  })
  observe({ covis.auto$tailor.lis <- tailor.lis$v })
  # Generates both embedding plots and the colData table in the scell page.
  dim.lis <- reactive({
    # if (!is.null(covis.auto$res)) print(list('dim', colData(covis.auto$res)[c(156, 256), c('assignedBulk', 'similarity')]))
    dim_server('dim', sce=reactive({covis.auto$res}), sce.upl=sce.upl, section='scell', upl.mod.lis=upl.mod.lis) 
  }) # Avoid endless circle.
  observe({ 
    covis.auto$covisGrp <- dim.lis()$covisGrp$val 
    covis.auto$but.covis <- dim.lis()$but.covis$v 
  })
  onBookmark(function(state) { state })
  # return(list(covis.auto=covis.auto, tailor.lis=tailor.lis))
})}
