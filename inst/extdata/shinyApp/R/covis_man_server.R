# Module for co-visualization through annotation/manual methods.
covis_man_server <- function(id, sce.upl, upl.mod.lis, shm.mod.lis, tab, covis.man, session) {
  moduleServer(id, function(input, output, session) {
  cat('Module covis_man_server ... \n')
   ns <- session$ns; 
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
      blk.fil <- filter_data(data=blk.aggr, pOA=c(filBlkP, filBlkA), CV=c(filBlkCV1, filBlkCV2), verbose=FALSE)
    } else blk.fil <- NULL 
    incProgress(0.3, detail="please wait ...")
    dat.fil <- filter_cell(sce=cell, bulk=blk.fil, gen.rm=NULL, cutoff=cutoff, p.in.cell=filPCell, p.in.gen=filPGen, verbose=FALSE)
    incProgress(0.3, detail="please wait ...")
    if (is(dat.nor, 'list')) covis.man$bulk <- dat.fil$bulk
    cat('Done! \n'); return(dat.fil)
    })
  })

  observeEvent(input$subdat+1, ignoreInit=FALSE, {
    withProgress(message="Data table in co-visualization: ", value=0, {
    incProgress(0.3, detail="please wait ...")
    dat.fil <- dat.fil() 
    if (is(dat.fil, 'list')) { 
      blk <- dat.fil$bulk; cell <- dat.fil$cell
    } else { blk <- NULL; cell <- dat.fil }
    r1 <- input$r1; r2 <- input$r2
    c1 <- input$c1; c2 <- input$c2
    if (!check_obj(list(r1, r2, c1, c2))) return()
    output$datall <- renderDataTable({
      dat_covis_man(cbind(blk, cell), r1=r1, r2=r2, c1=c1, c2=c2)
    }); incProgress(0.3, detail="please wait ...")
   })
  })

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
    dimred <- reduce_dim(sce = cell, min.dim = minRank, max.dim = maxRank)
    incProgress(0.3, detail="please wait ...")
    cat('Done! \n'); return(dimred)
    })
  }); observe({ covis.man$dimred <- dimred() })
 
   observeEvent(covis.man$dimred, ignoreInit=FALSE, {
     updateTabsetPanel(session, inputId="tabSetCell", selected='dimred')
  })

   cnt.help <- reactiveValues(v=0)
   observeEvent(input$tabSetCell, {
     if (!'dimred' %in% input$tabSetCell | cnt.help$v > 3) return()
      showModal(div(id = 'matchHel',
      modalDialog(title = HTML('<strong><center>Matching Tissues and Cell Groups</center></strong>'),
      div(style = 'overflow-y:scroll;overflow-x:scroll',
      HTML('<img src="image/match.jpg">'),
      )))); cnt.help$v <- cnt.help$v+1
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
     if (grepl(na.sgl, upl.mod.lis$ipt$fileIn)) covis.man$match.mod.lis <- match_server('rematchCell', shm.mod.lis$sam, tab, upl.mod.lis, covis.man=covis.man, col.idp='idp' %in% shm.mod.lis$ipt$profile) else covis.man$match.mod.lis <- NULL
   })
  cat('Module covis_man_server done !\n')
  onBookmark(function(state) { state })
  # return(list(match.mod.res=match.mod.res))
})}
