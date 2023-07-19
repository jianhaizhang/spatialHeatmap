# Module for co-visualization through automatic method.
covis_auto_server <- function(id, sce.upl, upl.mod.lis, shm.mod.lis, tab, covis.auto, lis.url, parent, session) {
  moduleServer(id, function(input, output, session) {
   ns <- session$ns; cfg <- upl.mod.lis$cfg 
   quick <- reactiveValues(v=0)
   observeEvent(list(parent$input$tabTop, input$tabSetCellAuto), {
     tabTop <- parent$input$tabTop; if (!check_obj(tabTop)) return()
     if (quick$v <= 2 & 'scell' %in% tabTop & 'result' %in% input$tabSetCellAuto & 'auto' %in% covis.auto$method) {
         showModal(modal(title = HTML('<b>Quick start!</b>'), msg = 'Showing 3 times only!', img='coclus_quick.jpg', img.w="100%")) 
         quick$v <- quick$v + 1
     }
   }) 
   cnt.help <- reactiveValues(v=0)
   observeEvent(input$tabSetCellAuto, {
     if (!'tailor' %in% input$tabSetCellAuto | cnt.help$v > 2) return()
      showModal(div(id = 'tailorIns',
      modalDialog(title = HTML('<strong><center>Tailoring Assignment Results (optional)</center></strong>'),
      div(style = 'overflow-y:scroll;overflow-x:scroll',
      HTML('<img src="image/tailoring.jpg">'),
      )))); cnt.help$v <- cnt.help$v+1
    })
   observe({
     covis.auto$method <- sce.upl$method
     covis.auto$covis.type <- sce.upl$covis.type
   })
   observe({
    hideTab(inputId="tabSetCellAuto", target="tailor")
   })
  # Normalize bulk and cells.
  norm.coclus <- reactiveValues()
  observeEvent(input$parAutoBut+1, {
    if (is.null(sce.upl$bulk)) return()
    norm.coclus$val <- list(input$normCoclus)
  })
  # eventReactive avoids endless circles.
  dat.nor <- eventReactive(list(norm.coclus$val, sce.upl$bulk, sce.upl$cell), {
    cat('Co-clustering: normalizing ... \n')
    norm.meth <- input$normCoclus
    bulk <- sce.upl$bulk; cell <- sce.upl$cell
    if (is.null(norm.meth)|is.null(bulk)|is.null(cell)) return()
    withProgress(message="Normalizing: ", value=0, {
    incProgress(0.3, detail="in progress ...")
    dat.nor <- norm_cell_m(sce=cell, bulk=bulk, cpm=ifelse(norm.meth=='cpm', TRUE, FALSE)); cat('Done! \n'); return(dat.nor)
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
    filBlkCV2 <- input$filBlkCV2; cutoff <- input$cutoff
    filPGen <- input$filPGen; filPCell <- input$filPCell
    validate(need((filBlkP >=0 & filBlkP <=1) & (filBlkA >= 0) & (filPGen >=0 & filPGen <=1) & (filPCell >=0 & filPCell <=1) & (filBlkCV1 >= -1000 & filBlkCV1 <= 1000 & filBlkCV2 >= -1000 & filBlkCV2 <= 1000 & filBlkCV1 < filBlkCV2) & is(cutoff,'numeric'), ''))
    withProgress(message="Filtering: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    blk.aggr <- aggr_rep(data=blk, assay.na='logcounts', sam.factor='sample', aggr='mean')
    blk.aggr$spFeature <- NULL
    blk.fil <- filter_data(data=blk.aggr, pOA=c(filBlkP, filBlkA), CV=c(filBlkCV1, filBlkCV2), verbose=FALSE)
    incProgress(0.3, detail="please wait ...")
    dat.fil <- filter_cell(sce=cell, bulk=blk.fil, gen.rm=NULL, cutoff=cutoff, p.in.cell=filPCell, p.in.gen=filPGen, verbose=FALSE)
    cat('Done! \n'); return(dat.fil)
    })
  })
  dat_all_server(id='dat', dat.fil, r2=500, c2=50)

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
    dimred <- reduce_dim_m(sce = com.kp, min.dim = minRank, max.dim = maxRank); cat('Done! \n'); return(dimred)
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
    clus.all <- detect_cluster_m(graph = gr, clustering = clusMeth)
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
     if (!is.null(res()) & 'auto' %in% covis.auto$method) showModal(modal(title = HTML('<b>Quick start!</b>'), msg = NULL, img='coclus_quick.jpg', img.w="100%")) 
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
   output$help <- renderUI({ 
    tags$iframe(seamless="seamless", src= "html/shm_shiny_manual.html#coclus", width='100%', height='100%') 
  })

  observe({
    lis.par <- cfg$lis.par; lis.url
    req(check_obj(list(lis.par, lis.url)))
    updateSelectInput(session, 'normCoclus', selected=url_val('scell-covisAuto-normCoclus', lis.url, def=lis.par$coclus.set['normCoclus', 'default']))  
    updateNumericInput(session, 'filBlkA', value=url_val('scell-covisAuto-filBlkA', lis.url, def=as.numeric(lis.par$coclus.set['A', 'default'])))  
    updateNumericInput(session, 'filBlkP', value=url_val('scell-covisAuto-filBlkP', lis.url, def=as.numeric(lis.par$coclus.set['P', 'default'])))  
    updateNumericInput(session, 'filBlkCV1', value=url_val('scell-covisAuto-filBlkCV1', lis.url, def=as.numeric(lis.par$coclus.set['CV1', 'default'])))  
    updateNumericInput(session, 'filBlkCV2', value=url_val('scell-covisAuto-filBlkCV2', lis.url, def=as.numeric(lis.par$coclus.set['CV2', 'default'])))  
    updateNumericInput(session, 'cutoff', value=url_val('scell-covisAuto-cutoff', lis.url, def=as.numeric(lis.par$coclus.set['cutoff', 'default'])))  
    updateNumericInput(session, 'filPGen', value=url_val('scell-covisAuto-filPGen', lis.url, def=as.numeric(lis.par$coclus.set['filPGen', 'default'])))  
    updateNumericInput(session, 'filPCell', value=url_val('scell-covisAuto-filPCell', lis.url, def=as.numeric(lis.par$coclus.set['filPCell', 'default'])))  
    updateNumericInput(session, 'minRank', value=url_val('scell-covisAuto-minRank', lis.url, def=as.numeric(lis.par$coclus.set['minRank', 'default'])))  
    updateNumericInput(session, 'maxRank', value=url_val('scell-covisAuto-maxRank', lis.url, def=as.numeric(lis.par$coclus.set['maxRank', 'default'])))  
    updateSelectInput(session, 'dimSel', selected=url_val('scell-covisAuto-dimSel', lis.url, def=lis.par$coclus.set['dimSel', 'default']))  
    updateSelectInput(session, 'graphMeth', selected=url_val('scell-covisAuto-graphMeth', lis.url, def=lis.par$coclus.set['graphMeth', 'default']))  
    updateSelectInput(session, 'clusMeth', selected=url_val('scell-covisAuto-clusMeth', lis.url, def=lis.par$coclus.set['clusMeth', 'default']))  
    updateNumericInput(session, 'asgThr', value=url_val('scell-covisAuto-asgThr', lis.url, def=as.numeric(lis.par$coclus.set['asgThr', 'default'])))  
  })
  onBookmark(function(state) { state })
  # return(list(covis.auto=covis.auto, tailor.lis=tailor.lis))
})}
