
# The Shiny modules (e.g. search_server) are temporarily placed in this file only for debugging purpose, and will be moved to independent files in the R folder after the App development is completed.
options(list(stringsAsFactors=FALSE, shiny.fullstacktrace=TRUE))
# options(stringsAsFactors=FALSE) 

# Every variable in every container should be checked at the beginning. E.g. input$fileIn in reactive({}). These checks will avoid amost all errors/warnings.

# enableWGCNAThreads()
# enableBookmarking("url")
server <- function(input, output, session) {
  set.seed(10)
  observe({ hide(id = "absInf") })
  observeEvent(input$tabTop, {
    tabTop <- input$tabTop
    if (!check_obj(tabTop)) { hide(id = "absInf"); return() }
    if (tabTop %in% c('ldg', 'about')) { hide(id = "absInf"); return() }
    showElement(id = "absInf")
  })

  observeEvent(list(input$btnInf), {
    btnInf <- input$btnInf
    if (!check_obj(btnInf)) return()
    if (btnInf >0) {
      updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
      runjs("scrollBott()") 
    }
  })
  observeEvent(input$btnLdg, { 
    updateTabsetPanel(session, "tabTop", selected='landing')
  })

  ldg <- reactiveValues(v=0)
  observeEvent(input$tabTop, { 
    if (ldg$v <= 2 & 'landing' %in% input$tabTop) { showModal(modal(title = HTML('<center><b>Please wait when the App is in progress!</b><center>'), msg = strong('Getting started (showing 3 times only)!'), img='dataset.jpg', img.w="100%"))
    if ('landing' %in% input$tabTop) ldg$v <- ldg$v + 1
    }
  })

  lis.url <- reactiveValues(par=NULL)
  observe({
    # The url on browser is captured only if the url is refreshed or the "Enter" key is pressed, which applies in the cases that shiny app is first launched or users modified parameters in the url. Otherwise the query is null.
    # lis.url$par is an empty list not NULL before refreshing/bookmarking.
    hos.port <- session$clientData
    lis.url$par <- parseQueryString(hos.port$url_search)
    hos.port1 <- paste0(hos.port$url_hostname, ':', hos.port$url_port, hos.port$url_pathname)
    lis.url$hos.port <- hos.port1
    lis.url$url_search <- hos.port$url_search
    # cat('Parameters in URL:', (names(lis.url$par)), '\n')
   })

  observe({
    tabTop <- input$tabTop; if (!check_obj(tabTop)) return()
    if (tabTop %in% c('ldg', 'about')) return()  
    # validate(need(input$start.but>0, ''))
    withProgress(message="Loading dependencies: ", value=0, {
    incProgress(0.2, detail="in progress ...")
    library(spatialHeatmap); library(SummarizedExperiment); library(shiny); library(shinydashboard); library(shinydashboardPlus); library(grImport); library(rsvg); library(ggplot2); library(DT) 
    incProgress(0.2, detail="in progress ...")
    library(gridExtra); library(ggdendro); library(grid); library(xml2); library(plotly); library(data.table); library(genefilter); 
    incProgress(0.2, detail="in progress...")
   library(reshape2); library(igraph); library(animation); library(av); library(shinyWidgets); library(yaml); library(HDF5Array); library(shinyBS); library(shinyjs); library(htmltools)
    # DEG
    library(gplots); library(UpSetR); library(sparkline); library(spsUtil)
    incProgress(0.2, detail="in progress...")
  })
  })

  url.id <- reactiveValues()
  mods <- reactiveValues(upload=NULL, search=NULL, data=NULL, shm=NULL, deg=NULL, scell=NULL)

  # Selected IDs and button in search box.
  ids <- reactiveValues(sel=NULL, but=NULL)
  mods$upload <- upl.mod.lis <- upload_server('upl', lis.url, prt=session)
  ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg 

  sch <- reactiveValues()
  observe({ sch$sch <- input$search; sch$but <- input$search.but })
  #ipt0 <- reactiveValues(url.but=NULL, but=NULL)
  # observe({ ipt0$url.but <- input$url.but })
  mods$data <- dat.mod.lis <- data_server('dat', sch, lis.url, ids, upl.mod.lis, deg.mod.lis=mods$deg, scell.mod.lis=mods$scell, shm.mod=mods$shm, parent=session)

  mods$search <- sch.mod.lis <- search_server('sear', ids, lis.url, url.id, upl.mod.lis, dat.mod.lis)
 
 # Active tab.
 tab <- reactiveValues(); observe({ tab$val <- input$tabTop })

# Variables of "reactiveValues()" are accessible anywhere once created. E.g. created in one module and used instantly in another module without being passed as an argument. So "url.id" can be accessed if not passed as an argument, even in the inner nested modules.

mods$shm <- shm.mod.lis <- shm_server('shmAll', sch, lis.url, url.id, tab, upl.mod.lis, dat.mod.lis, sch.mod.lis, mods$scell, mods$dim, deg.mod=mods$deg, prt=session)

  cnt.ana <- reactiveValues(v=0)
  observeEvent(input$tabTop, {
    tabTop <- input$tabTop
    if (!'ana' %in% tabTop) return()
    if (length(ids$sel) > 0 & cnt.ana$v <=2) showModal(modal(title=HTML(run.msg), msg='Showing 3 times only!', easyClose=TRUE))
    if ('ana' %in% tabTop) cnt.ana$v <- cnt.ana$v + 1
  })

  mods$ana <- ana.mod <- analysis_server('ana', upl.mod.lis, dat.mod.lis, shm.mod.lis, ids)

  svgs.upd <- reactive({ 
    if (is.null(mods$shm$svgs)) return(); mods$shm$svgs() 
  })
  observeEvent(svgs.upd(), ignoreInit=FALSE, ignoreNULL=TRUE, {
    svgs <- svgs.upd()
    if (!is.null(svgs) & (!grepl(na.sgl, mods$upload$ipt$fileIn)) & ldg$v >=2) updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
  })

  cnt.shm <- reactiveValues(v=0)
  observeEvent(input$tabTop, {
    tabTop <- input$tabTop
    if (!'shmPanelAll' %in% tabTop) return()
    lgc <- length(ids$sel)==0 & cnt.shm$v <=2
    if (lgc) { 
      showModal(modal(title = 'Quick start!', msg='Showing 3 times only!', img='select.jpg', img.w="70%")) 
      cnt.shm$v <- cnt.shm$v + 1
    }
  })

  mods$deg <- deg.mod.lis <- deg_server('deg', sch, lis.url, url.id, ids, upl.mod.lis, dat.mod.lis, shm.mod.lis, parent=session)
  output$datInf <- renderUI({
    cfg <- upl.mod.lis$cfg; fileIn <- upl.mod.lis$ipt$fileIn
    se.scl <- mods$data$se.scl(); df.meta <- metadata(se.scl)$df.meta
    if (!check_obj(list(cfg, fileIn, se.scl, df.meta))) return()
    files <- c(cfg$na.cus, cfg$na.def)
    dat <- names(files[files %in% fileIn])
    if (!grepl('^covis_', fileIn)) {
      df.meta <- subset(df.meta, name=='assay')
      sumar <- HTML(paste0('<strong>Dataset: ', dat, '</strong><br/>', '&nbsp;&nbsp;&nbsp;', df.meta$species[1], ', ', df.meta$technology[1], ', ', df.meta$id[1]))
    } else {
      df.m.blk <- subset(df.meta, grepl('^assayBulk$', name))
      df.m.cell <- subset(df.meta, grepl('^assayCell$', name))
      sumar <- HTML(paste0('<strong>Dataset: ', dat, '</strong><br/>', '&nbsp;&nbsp;&nbsp;1. Bulk: ', df.m.blk$species, ', ', df.m.blk$technology, ', ', df.m.blk$id, '<br/>', '&nbsp;&nbsp;&nbsp;2. Cell: ', df.m.cell$species, ', ', df.m.cell$technology, ', ', df.m.cell$id)) 
    }; sumar
  })
  output$covisInf <- renderUI({
    fileIn <- upl.mod.lis$ipt$fileIn
    sce.ipt <- scell.mod.lis$input
    if (!check_obj(list(sce.ipt, fileIn))) return()
    sumar <- NULL; if (grepl('^covis_', fileIn)) {
      meth <- sce.ipt$methCovis; direc <- sce.ipt$direc
      meth.sel <- ifelse('man' %in% meth, 'Annotation (or similar) labels', 'Co-clustering labels') 
      direc.sel <- ifelse('toBulk' %in% direc, 'Cell-to-bulk', 'Bulk-to-cell')
      sumar <- paste0('Cell group labels: ', meth.sel, ', Mapping direction: ', direc.sel)
    }; sumar
  })
  mods$scell <- scell.mod.lis <- scell_server('scell', tab, upl.mod.lis, shm.mod.lis, lis.url=lis.url, parent=session, session)
  observe({
    # observeEvent(mods$scell$sce.upl$cell, ignoreInit=FALSE, ignoreNULL=FALSE, {
    ipt <- mods$upload$ipt; fileIn <- ipt$fileIn; datpa <- ipt$geneInpath
    svgs <- grepl('\\.svg$', c(ipt$svgInpath1, ipt$svgInpath2))
    lgc.upl <- !is.null(mods$scell$sce.upl$cell) & 'customCovisData' %in% fileIn & sum(svgs)>0
    lgc.def <- grepl(na.sgl.def, fileIn)
    cus <- fileIn %in% na.cus
    cus.upl.no <- is.null(datpa)|sum(svgs)==0 
    if (cus) { # Custom data.
      if (cus.upl.no) { message('Select SHM tab ...')
      disable(selector='a[data-value="shmPanelAll"]')
      disable(selector='a[data-value="deg"]') 
      disable(selector='a[data-value="ana"]') 
      disable(selector='a[data-value="scell"]'); message('Done!')
      } else enable(selector='a[data-value="shmPanelAll"]')
    } else enable(selector='a[data-value="shmPanelAll"]')
    if (lgc.upl | lgc.def) { # covis is active.
      message('Select covis tab ...') 
      enable(selector='a[data-value="scell"]')
      disable(selector='a[data-value="deg"]') 
      disable(selector='a[data-value="ana"]') 
      updateTabsetPanel(session, "tabTop", selected='scell')
      ft.match <- mods$scell$covis.man$match.mod.lis$ft.reorder$ft.rematch
      if (length(ft.match)==0 & 'man' %in% scell.mod.lis$input$methCovis) disable(selector='a[data-value="shmPanelAll"]') else enable(selector='a[data-value="shmPanelAll"]'); message('Done!')
    } else if (('customBulkData' %in% fileIn & !cus.upl.no)|(!cus & !grepl(na.sgl.def, fileIn))) {
      message('Select enrichment tab ...')
      enable(selector='a[data-value="deg"]')
      if (check_obj(ids$sel)) enable(selector='a[data-value="ana"]')
      disable(selector='a[data-value="scell"]') 
      message('Done!')
    }
  })
  
  #dimred <- reactive({ 
  #  sce.sub <- mods$scell$sce.rct$sce.sub
  #  sce.clus <- mods$scell$sce.rct$clus 
  #  if (is.null(sce.sub)) sce.clus else NULL
  #})
  # Only generates the colData table in the SHM page.
  #mods$dim <- dim.mod.lis <- dim_server('datDim', sce=dimred, sce.upl=mods$scell$sce.upl, section='data', upl.mod.lis=mods$upload, dat.lis=mods$scell$sce.res)

  # Only after the "Single-cell metadata" is clicked, this "renderUI" is called, so this "ui" is often rendered after "hideElement" in the server (tailor_match_server). As a result, elements are not hidden. Since "hideElement" should be called after the elements are rendered. The solution is to hide the elements through an argument (hide=TRUE) passed to ui. But if "hide=TRUE", "dimCellBut", "coclusPlotBut", "selBlkBut", "selBlkButCan" on ui are NULL, and the server is not executed. As a result, no table is shown on "Single-cell metadata". Thus the ui (tailor_match_ui) should be generated on the general ui part.
  # The downside of "renderUI" is some inputs are NULL before it is called, and the server might not be executed.
 
  # If 'eventExpr' is a list, one slot of NULL will trigger 'observeEvent', though ignoreNULL = TRUE, so these 'observeEvents' are not merged. 
  # Switch to SHM tab.
  #observeEvent(list(mods$deg$but.sgl()), ignoreInit=TRUE, {
  #  but.sgl <- mods$deg$but.sgl()
  #  if (is.null(but.sgl)) return(); if (but.sgl==0) return()
  #  updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
  #})
  #observeEvent(mods$deg$but.mul(), ignoreInit=TRUE, {
  #  but.mul <- mods$deg$but.mul()
  #  if (is.null(but.mul)) return(); if (but.mul==0) return()
  #  updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
  #})

  observeEvent(mods$deg$input$eSHMBut, ignoreInit=TRUE, {
    but <- mods$deg$input$eSHMBut
    if (is.null(but)) return(); if (but==0) return()
    updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
  })

  observeEvent(mods$scell$covis.man$match.mod.lis$but.match$val, ignoreInit=TRUE, {
    ft.match <- mods$scell$covis.man$match.mod.lis$ft.reorder$ft.rematch
    if (length(ft.match)==0) return()
    # Activates diabled tables.
    updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
  })
  observeEvent(mods$scell$covis.auto$but.covis, ignoreInit=TRUE, {
    updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
  })
  tailor.lis.com <- reactive({
    if (is.null(mods$scell)) return()
    if (is.null(mods$scell$covis.auto$res)) return()
    lis <- mods$scell$covis.auto$tailor.lis
    sce.res <- mods$scell$sce.res()
    if (is.null(lis$coclusPlotBut)) return()
    if (is.null(lis)|is.null(sce.res)) return()
    if (lis$selBlkBut()==0 & lis$selBlkCancel()==0 & lis$coclusPlotBut()==0) return()
    list(lis$coclusPlotBut(), lis$selBlkBut(), lis$selBlkCancel())
  })
  observeEvent(tailor.lis.com(), ignoreInit=TRUE, {
    updateTabsetPanel(session, inputId="tabTop", selected='shmPanelAll')
  })
  output$about <- renderUI({
    tags$iframe(seamless="seamless", src= "html/about.html", width='100%', height='100%') 
  })
  output$ldg <- renderUI({
    tags$iframe(seamless="seamless", src= "html/landing.html", width='100%', height='100%') 
  })

  setBookmarkExclude(c("dat-dtSel_rows_all", "dat-dtSel_rows_current", "dat-dtSel_search_columns", "dat-dtAll_rows_all", "dat-dtAll_rows_current", "dat-dtAll_search_columns", "dat-dtAll_state", "cell", "bulk", 'scell-covisMan-rematchCell-matHelp', 'right.bar', 'shmAll-net-dpbMea', 'scell-covisAuto-tailor-selBlkCancel', 'sidebarCollapsed', 'sear-sch.mul.but', 'shmAll-val.lgd', 'scell-covisAuto-tabSetCellAuto', 'dat-dat.all.but', 'scell-covisMan-dimredNav', 'deg-degAll', 'shmAll-shms.in', 'shmAll-col.but', 'sear-sch.mul', 'sear-sch.mode', 'deg-deg-sch.mode', 'dat-dtAll_columns_selected', 'scell-covisAuto-tailor-coclusPlotBut', 'scell-covisAuto-parAutoBut', 'shmAll-rematch-matHelp', 'shmAll-transBut', 'dat-tran.scale.but.sel', 'shmAll-dropdown', 'dat-dtAll_rows_selected', 'upl-tar', 'shmAll-net-gen.sel', 'shmAll-dld.but', 'scell-covisHelp', 'dat-fil.but', 'scell-covisMan-parManBut', 'deg-datDEG-fil.but', 'shmAll-net-col.but.net', 'dat-sig.but', 'dat-tran.scale.but.prof', 'shmAll-net-cpt.nw', 'shmAll-net-mhm.but', 'shmAll-glyBut', 'deg-ssg.update')) 
  observe({
    lis.ipt <- reactiveValuesToList(input); session$doBookmark()
    # lapply(seq_along(lis.ipt), function(i) {if (length(lis.ipt[[i]])<1000) { print(lis.ipt[i]) }})
    # print(length(lis.ipt[['shmAll-dat-dt_rows_all']]))
    # print(getUrlHash()); print(getQueryString())
  })
 onBookmarked(function(url) { updateQueryString(url) })
 # onBookmarked(updateQueryString)
  onStop(function() { ggly_rm(); vdo_rm(); data_mining_rm()
    message("Session stopped! \n") 
  })

}






