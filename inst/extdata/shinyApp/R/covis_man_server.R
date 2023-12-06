# Module for co-visualization through annotation/manual methods.
covis_man_server <- function(id, sce.upl, upl.mod.lis, shm.mod.lis, tab, covis.man, lis.url, parent, session) {
  moduleServer(id, function(input, output, session) {
  cat('Module covis_man_server ... \n')
   ns <- session$ns; cfg <- upl.mod.lis$cfg 
   observe({
     # sce.man <- reactiveValues(bulk=sce.upl$bulk): if outside "observe", only runs one time when the app is started.
     # covis.man$bulk <- sce.upl$bulk
     covis.man$method <- sce.upl$method
     covis.man$covis.type <- sce.upl$covis.type
   })
  
  # Normalize bulk and cells.
  norm.par <- reactiveValues()
  observeEvent(input$parManBut+1, {
    if (is.null(sce.upl$bulk)) return()
    norm.par$val <- list(input$norm)
  })
  # eventReactive avoids endless circles.
  dat.nor <- eventReactive(list(norm.par$val, sce.upl$bulk, sce.upl$cell), {
    message('Manual matching: normalizing ...')
    norm.meth <- input$norm
    bulk <- sce.upl$bulk; cell <- sce.upl$cell
    if (is.null(norm.meth)|is.null(bulk)|is.null(cell)) return()
    withProgress(message="Normalizing: ", value=0, {
    incProgress(0.3, detail="in progress ...")
    dat.nor <- norm_cell_m(sce=cell, bulk=bulk, cpm=ifelse(norm.meth=='cpm', TRUE, FALSE))
    incProgress(0.3, detail="in progress ...")
    message('Done!'); return(dat.nor)
    })
  }) 
 
  # Filtering.
  par.fil <- reactiveValues()
  observeEvent(input$parManBut+1, {
    if (is.null(sce.upl$bulk)) return()
    par.fil$val <- list(input$filBlkP, input$filBlkA, input$filBlkCV1, input$filBlkCV2, input$filPGen, input$filPCell)
  })
  # eventReactive avoids endless circles.
  dat.fil <- eventReactive(list(par.fil$val, dat.nor()), {
    cat('Manual matching: filtering ... \n')
    library(SingleCellExperiment)
    dat.nor <- dat.nor(); if (is.null(dat.nor)) return() 
    if (is(dat.nor, 'list')) { 
      blk <- dat.nor$bulk; cell <- dat.nor$cell
    } else { blk <- NULL; cell <- dat.nor }
    
    filBlkP <- input$filBlkP; filBlkA <- input$filBlkA
    filBlkCV1 <- input$filBlkCV1
    filBlkCV2 <- input$filBlkCV2; cutoff <- input$cutoff
    filPGen <- input$filPGen; filPCell <- input$filPCell
    validate(need((filBlkP >=0 & filBlkP <=1) & (filBlkA >= 0) & (filPGen >=0 & filPGen <=1) & (filPCell >=0 & filPCell <=1) & (filBlkCV1 >= -1000 & filBlkCV1 <= 1000 & filBlkCV2 >= -1000 & filBlkCV2 <= 1000 & filBlkCV1 < filBlkCV2) & is(cutoff, 'numeric'), ''))
    withProgress(message="Filtering: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    if (!is.null(blk)) {
      blk.aggr <- aggr_rep(data=blk, assay.na='logcounts', sam.factor='sample', aggr='mean')
      # blk.aggr$spFeature <- NULL
      blk.fil <- filter_data(data=blk.aggr, pOA=c(filBlkP, filBlkA), CV=c(filBlkCV1, filBlkCV2), verbose=FALSE)
    } else blk.fil <- NULL 
    incProgress(0.3, detail="please wait ...")
    dat.fil <- filter_cell(sce=cell, bulk=blk.fil, gen.rm=NULL, cutoff=cutoff, p.in.cell=filPCell, p.in.gen=filPGen, verbose=FALSE)
    incProgress(0.3, detail="please wait ...")
    if (is(dat.nor, 'list')) covis.man$bulk <- dat.fil$bulk
    cat('Done! \n'); return(dat.fil)
    })
  })
  dat_all_server(id='dat', dat.fil, r2=500, c2=50)
  # Dimension reduction.
  par.dim <- reactiveValues()
  # Button of 0 cannot trigger observeEvent.
  observeEvent(input$parManBut+1, {
    if (is.null(dat.fil())) return()
    par.dim$val <- list(input$minRank, input$maxRank)
  })
  dimred <- eventReactive(list(par.dim$val, dat.fil()), {
    cat('Manual matching: reducing dimensions ... \n')
    dat.fil <- dat.fil(); if (is(dat.fil, 'list')) {     
      cell <- dat.fil$cell
    } else cell <- dat.fil 
    minRank <- input$minRank; maxRank <- input$maxRank
    if (is.null(cell)|is.null(minRank)|is.null(maxRank)) return()
    validate(need(round(minRank)==minRank & round(maxRank)==maxRank & maxRank > minRank, ''))
    withProgress(message="Reducing dimensions: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    dimred <- reduce_dim_m(sce = cell, min.dim = minRank, max.dim = maxRank)
    incProgress(0.3, detail="please wait ...")
    cat('Done! \n'); return(dimred)
    })
  })
  observe({ 
    bulk <- covis.man$bulk; cell <- dimred()
    # Keep colData columns consistent between bulk and cell.
    if (!is.null(bulk)) {
      cdat.b <- colData(bulk); cdat.c <- colData(cell)
      int <- intersect(colnames(cdat.b), colnames(cdat.c))
      lgc.na <- length(int) > 0
      if (!lgc.na) showModal(modal(msg = 'No common column names detected between "colData" slots of bulk and single-cell data!')); req(lgc.na)
      colData(bulk) <- cdat.b[, int]; colData(cell) <- cdat.c[, int]
      covis.man$bulk <- bulk; covis.man$dimred <- cell
    } else covis.man$dimred <- cell 
  })
 
   observeEvent(covis.man$dimred, ignoreInit=FALSE, {
     updateTabsetPanel(session, inputId="tabSetCell", selected='dimred')
  })

   cnt.help <- reactiveValues(v=0, notshow=FALSE)
   observeEvent(input$showman, {
     showman <- input$showman; if (!check_obj(showman) | TRUE %in% cnt.help$notshow) return()
     cnt.help$notshow <- showman
   })
   observeEvent(input$tabSetCell, {
     if (!'dimred' %in% input$tabSetCell | cnt.help$v >= 3 | TRUE %in% cnt.help$notshow) return()
      showModal(modal(title='Quick start!', msg = 'Showing 3 times only!', img='ann_quick.jpg', img.w="100%", idshow=ns('showman')))
      cnt.help$v <- cnt.help$v+1
    })
  observeEvent(covis.man$dimred, ignoreInit=FALSE, ignoreNULL=FALSE, {
    covis.type <- sce.upl$covis.type
    if (is.null(covis.type)) return()
    if (covis.type %in% c('toBulkAuto', 'toCellAuto')) {
        showTab('tabSetCell', target='autoMatch')
        hideTab('tabSetCell', target='manualMatch')
    } else {
        showTab('tabSetCell', target='manualMatch')
    }
  })
  output$help <- renderUI({ 
    tags$iframe(seamless="seamless", src= "html/shm_shiny_manual.html#ann", width='100%', height='100%') 
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
     if (grepl(na.sgl, upl.mod.lis$ipt$fileIn)) covis.man$match.mod.lis <- match_server('rematchCell', shm.mod.lis$sam, tab, upl.mod.lis, covis.man=covis.man, col.idp='idp' %in% shm.mod.lis$ipt$profile) else covis.man$match.mod.lis <- NULL
   })
   observe({
     lis.par <- cfg$lis.par; lis.url
     req(check_obj(list(lis.par, lis.url)))
     updateSelectInput(session, 'norm', selected=url_val('scell-covisMan-norm', lis.url, def=lis.par$ann.set['norm', 'default']))
     updateNumericInput(session, 'filBlkA', value=url_val('scell-covisAuto-filBlkA', lis.url, def=as.numeric(lis.par$ann.set['A', 'default']))) 
     updateNumericInput(session, 'filBlkP', value=url_val('scell-covisAuto-filBlkP', lis.url, def=as.numeric(lis.par$ann.set['P', 'default']))) 
     updateNumericInput(session, 'filBlkCV1', value=url_val('scell-covisAuto-filBlkCV1', lis.url, def=as.numeric(lis.par$ann.set['CV1', 'default']))) 
     updateNumericInput(session, 'filBlkCV2', value=url_val('scell-covisAuto-filBlkCV2', lis.url, def=as.numeric(lis.par$ann.set['CV2', 'default']))) 
     updateNumericInput(session, 'cutoff', value=url_val('scell-covisAuto-cutoff', lis.url, def=as.numeric(lis.par$ann.set['cutoff', 'default']))) 
     updateNumericInput(session, 'filPGen', value=url_val('scell-covisAuto-filPGen', lis.url, def=as.numeric(lis.par$ann.set['filPGen', 'default']))) 
     updateNumericInput(session, 'filPCell', value=url_val('scell-covisAuto-filPCell', lis.url, def=as.numeric(lis.par$ann.set['filPCell', 'default']))) 
     updateNumericInput(session, 'minRank', value=url_val('scell-covisAuto-minRank', lis.url, def=as.numeric(lis.par$ann.set['minRank', 'default']))) 
     updateNumericInput(session, 'maxRank', value=url_val('scell-covisAuto-maxRank', lis.url, def=as.numeric(lis.par$ann.set['maxRank', 'default']))) 
  })

  cat('Module covis_man_server done !\n')
  onBookmark(function(state) { state })
  # return(list(match.mod.res=match.mod.res))
})}
