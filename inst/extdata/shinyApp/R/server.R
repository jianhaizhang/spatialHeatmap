
# The Shiny modules (e.g. search_server) are temporarily placed in this file only for debugging purpose, and will be moved to independent files in the R folder after the App development is completed.
options(list(stringsAsFactors=FALSE, shiny.fullstacktrace=TRUE))
# options(stringsAsFactors=FALSE) 

# Every variable in every container should be checked at the beginning. E.g. input$fileIn in reactive({}). These checks will avoid amost all errors/warnings.

# enableWGCNAThreads()
# enableBookmarking("url")
server <- function(input, output, session) {
  set.seed(10)
  observeEvent(input$man, {
    man <- input$man; req(man); req(man>=1)
    updateNavbarPage(session, 'navLdg', selected = 'help')
  })
  observe({ hide(id = "absInf") })
  observeEvent(input$tabTop, {
    tabTop <- input$tabTop
    if (!check_obj(tabTop)) { hide(id = "absInf"); return() }
    if (tabTop %in% c('ldg', 'about')) { hide(id = "absInf"); return() }
    showElement(id = "absInf")
  })
  distab <- reactiveValues()
  observeEvent(list(input$btnInf), {
    btnInf <- input$btnInf
    if (!check_obj(btnInf)) return()
    if (btnInf >0) {
      updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
      runjs("scrollBott()") 
    }
  })
  observeEvent(input$btnLdg, { 
    updateTabsetPanel(session, "tabTop", selected='dataset')
  })

  ldg <- reactiveValues(v=0, notshow=FALSE)
  observeEvent(input$showldg, {
    showldg <- input$showldg; if (!check_obj(showldg) | TRUE %in% ldg$notshow) return()
    ldg$notshow <- showldg
  })
  observeEvent(input$tabTop, { 
    if (ldg$v <= 2 & 'dataset' %in% input$tabTop & FALSE %in% ldg$notshow) { showModal(modal(title = HTML('<center><b>Please wait when the App is in progress!</b><center>'), msg = strong('Getting started (showing 3 times only)!'), img='dataset.jpg', img.w="100%", idshow='showldg'))
    if ('dataset' %in% input$tabTop) ldg$v <- ldg$v + 1 
    } 
  }) 

  lis.url <- reactiveValues(par=NULL)
  observeEvent(session$clientData, {
    # The url on browser is captured only if the url is refreshed or the "Enter" key is pressed, which applies in the cases that shiny app is first launched or users modified parameters in the url. Otherwise the query is null.
    # lis.url$par is an empty list not NULL before refreshing/bookmarking.
    hos.port <- session$clientData
    lis.url$par <- parseQueryString(hos.port$url_search)
    # print(list(hos.port$url_search, lis.url$par))
    id.url <- lis.url$par$ids
    if (check_obj(id.url)) {
      id.url <- setdiff(gsub('"', '', id.url), 'null')
      if (length(id.url)>0) lis.url$par$ids <- strsplit(id.url, ',')[[1]]
    }
    hos.port1 <- paste0(hos.port$url_hostname, ':', hos.port$url_port, hos.port$url_pathname)
    lis.url$hos.port <- hos.port1
    lis.url$url_search <- hos.port$url_search
    # cat('Parameters in URL:', (names(lis.url$par)), '\n')
    # setBookmarkExclude(setdiff(names(input), c('upl-fileIn')))
   })
   url.cnt <- reactiveValues(v=0)
   observeEvent(lis.url$par$ids, {
     id.url <- lis.url$par$ids; fileIn <- lis.url$par$'upl-fileIn'
     req(url.cnt$v==0 & check_obj(list(id.url, !dat.no %in% fileIn)))
     Sys.sleep(1); updateTabsetPanel(session, "tabTop", selected='dataset')
     Sys.sleep(0.5); if (!grepl(na.sgl, fileIn)) { 
       updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
       url.cnt$v <- 1
     } else {
       updateTabsetPanel(session, "tabTop", selected='scell'); Sys.sleep(1)
       updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
       url.cnt$v <- 1
     }
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
  # his <- reactiveValues()
  his <- spsUtil::historyStack$new(); his$add('ldg')
  observeEvent(input$tabTop, {
    # save(h, file='h')
    tabTop <- input$tabTop; req(check_obj(list(tabTop, his))) 
    req(!tabTop %in% his$get()$item)
    his$add(tabTop)
  })
  observeEvent(input$backBut, {
    tabTop <- input$tabTop; req(check_obj(list(tabTop, his)))
    first <- his$get()$first
    if (TRUE %in% first) showModal(modal(title='Already at the first tab!', easyClose=TRUE))
    req(!first); his$backward() 
    item <- his$get()$item; lgc.dis <- item %in% distab$v
    if (lgc.dis) showModal(modal(title='The requested tab is disabled at the current status!', easyClose=TRUE))
    req(!lgc.dis)
    updateTabsetPanel(session, "tabTop", selected=item)
  })
  observeEvent(input$forBut, {
    tabTop <- input$tabTop; req(check_obj(list(tabTop, his)))
    last <- his$get()$last
    if (TRUE %in% last) showModal(modal(title='Already at the last tab!', easyClose=TRUE))
    req(!last); his$forward() 
    item <- his$get()$item; lgc.dis <- item %in% distab$v
    if (lgc.dis) showModal(modal(title='The requested tab is disabled at the current status!', easyClose=TRUE))
    req(!lgc.dis)
    updateTabsetPanel(session, "tabTop", selected=item)
  })

  url.id <- reactiveValues()
  mods <- reactiveValues(upload=NULL, search=NULL, data=NULL, shm=NULL, deg=NULL, scell=NULL)

  # Selected IDs and button in search box.
  ids <- reactiveValues(sel=NULL, but=NULL)
  mods$upload <- upl.mod.lis <- upload_server('upl', lis.url, prt=session)
  ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg; covis.pa <- upl.mod.lis$covis.pa 

  sch <- reactiveValues()
  observe({ sch$sch <- input$search; sch$but <- input$search.but })
  #ipt0 <- reactiveValues(url.but=NULL, but=NULL)
  # observe({ ipt0$url.but <- input$url.but })
  mods$data <- dat.mod.lis <- data_server('dat', sch, lis.url, ids, upl.mod.lis, deg.mod.lis=mods$deg, scell.mod.lis=mods$scell, shm.mod=mods$shm, parent=session)

  mods$search <- sch.mod.lis <- search_server('sear', ids, lis.url, upl.mod.lis, dat.mod.lis)
 
 # Active tab.
 tab <- reactiveValues(); observe({ tab$val <- input$tabTop })

# Variables of "reactiveValues()" are accessible anywhere once created. E.g. created in one module and used instantly in another module without being passed as an argument. So "url.id" can be accessed if not passed as an argument, even in the inner nested modules.

mods$shm <- shm.mod.lis <- shm_server('shmAll', sch, lis.url, url.id, tab, upl.mod.lis, dat.mod.lis, sch.mod.lis, mods$scell, mods$dim, deg.mod=mods$deg, prt=session)
  # After plotting SHMs, lis.url$par is set NULL.
  #observeEvent(mods$shm$shmLay$val, {
  #  shm.res <- mods$shm$shmLay$val; req(check_obj(shm.res))
  #  id.url <- lis.url$par$ids
  #  if (check_obj(id.url)) lis.url$par <- NULL
  # })
  cnt.ana <- reactiveValues(v=0, notshow=FALSE)
  observeEvent(input$showana, {
    showana <- input$showana; if (!check_obj(showana) | TRUE %in% cnt.ana$notshow) return()
    cnt.ana$notshow <- showana
  })
  observeEvent(input$tabTop, {
    tabTop <- input$tabTop
    if (!'ana' %in% tabTop) return()
    if (length(ids$sel) > 0 & cnt.ana$v <=2 & FALSE %in% cnt.ana$notshow) showModal(modal(title=HTML(run.msg), msg='Showing 3 times only!', easyClose=TRUE, idshow='showana'))
    if ('ana' %in% tabTop) cnt.ana$v <- cnt.ana$v + 1
  })

  mods$ana <- ana.mod <- analysis_server('ana', upl.mod.lis, dat.mod.lis, shm.mod.lis, ids)

  svgs.upd <- reactive({ 
    if (is.null(mods$shm$svgs)) return(); mods$shm$svgs() 
  })
  observeEvent(svgs.upd(), ignoreInit=FALSE, ignoreNULL=TRUE, {
    svgs <- svgs.upd()
    if (!is.null(svgs) & (!grepl(na.sgl, mods$upload$ipt$fileIn)) & ldg$v >=0) updateTabsetPanel(session, "tabTop", selected='shmPanelAll')
  })
  cnt.shm <- reactiveValues(v=0, notshow=FALSE)
  observeEvent(input$showshm, {
    showshm <- input$showshm; if (!check_obj(showshm) | TRUE %in% cnt.shm$notshow) return()
    cnt.shm$notshow <- showshm
  })
  observeEvent(input$tabTop, {
    tabTop <- input$tabTop
    if (!'shmPanelAll' %in% tabTop) return()
    lgc <- length(ids$sel)==0 & cnt.shm$v <=2 & FALSE %in% cnt.shm$notshow
    if (lgc) { 
      showModal(modal(title = 'Quick start!', msg='Showing 3 times only!', img='select.jpg', img.w="70%", idshow='showshm')) 
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
  # Try to use observeEvent and avoid observe as much as possible, especiallly long code chunks.
  observeEvent(list(mods$upload$ipt$fileIn, mods$upload$ipt$geneInpath, mods$upload$ipt$svgInpath1, mods$upload$ipt$svgInpath2, mods$scell$sce.upl$cell, mods$scell$covis.man$match.mod.lis$ft.reorder$ft.rematch, scell.mod.lis$input$methCovis, ids$sel, covis.pa$dat, covis.pa$svg1, covis.pa$svg2), {
    message('Enable/disable tabs ...')
    # observeEvent(mods$scell$sce.upl$cell, ignoreInit=FALSE, ignoreNULL=FALSE, {
    ipt <- mods$upload$ipt; fileIn <- ipt$fileIn; datpa <- ipt$geneInpath
    svgs <- grepl('\\.svg$', c(ipt$svgInpath1, ipt$svgInpath2))
    if (dat.no %in% fileIn) {
      disable(selector='a[data-value="shmPanelAll"]')
      disable(selector='a[data-value="deg"]') 
      disable(selector='a[data-value="ana"]') 
      disable(selector='a[data-value="scell"]'); return()
    }
    lgc.upl <- !is.null(mods$scell$sce.upl$cell) & 'customCovisData' %in% fileIn & sum(svgs)>0 & check_obj(list(covis.pa$dat, !is.null(covis.pa$svg1)|!is.null(covis.pa$svg2)))
    lgc.def <- grepl(na.sgl.def, fileIn)
    cus <- fileIn %in% na.cus
    cus.upl.no <- is.null(datpa)|sum(svgs)==0 
    if (cus) { # Custom data.
      if (cus.upl.no) { message('Select SHM tab ...')
      disable(selector='a[data-value="shmPanelAll"]')
      disable(selector='a[data-value="deg"]') 
      disable(selector='a[data-value="ana"]') 
      disable(selector='a[data-value="scell"]')
      distab$v <- unique(c(distab$v, 'shmPanelAll', 'deg', 'ana', 'scell'))
      message('Done!')
      } else { 
        enable(selector='a[data-value="shmPanelAll"]')
        distab$v <- setdiff(distab$v, 'shmPanelAll')
      }
    } else { 
      enable(selector='a[data-value="shmPanelAll"]')
      distab$v <- setdiff(distab$v, 'shmPanelAll')
    }
    if (lgc.upl | lgc.def) { # covis is active.
      message('Select covis tab ...') 
      enable(selector='a[data-value="scell"]')
      disable(selector='a[data-value="deg"]') 
      disable(selector='a[data-value="ana"]')
      if (!check_obj(ids$sel)) updateTabsetPanel(session, "tabTop", selected='scell')
      distab$v <- unique(c(distab$v, 'deg', 'ana'))
      distab$v <- setdiff(distab$v, 'scell')
      ft.match <- mods$scell$covis.man$match.mod.lis$ft.reorder$ft.rematch
      if (length(ft.match)==0 & 'man' %in% scell.mod.lis$input$methCovis) { 
        disable(selector='a[data-value="shmPanelAll"]') 
        distab$v <- unique(c(distab$v, 'shmPanelAll'))
        message('Ann labels, not matched, done!')
      } else { 
        enable(selector='a[data-value="shmPanelAll"]')
        distab$v <- setdiff(distab$v, 'shmPanelAll')
        message('Ann labels, matched, done!')
      }
    } else if (('customBulkData' %in% fileIn & !cus.upl.no)|(!cus & !grepl(na.sgl.def, fileIn))) {
      message('Select enrichment tab ...')
      enable(selector='a[data-value="deg"]')
      distab$v <- setdiff(distab$v, 'deg')
      if (check_obj(ids$sel)) { 
        enable(selector='a[data-value="ana"]')
        distab$v <- setdiff(distab$v, 'ana')
      }
      disable(selector='a[data-value="scell"]') 
      distab$v <- unique(c(distab$v, 'scell'))
      message('Done!')
    }; message('Done!')
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
  output$gallery <- renderUI({ 
    tags$iframe(seamless="seamless", src= "html/gallery.html", width='100%', height='100%')
  })
  output$manual <- renderUI({
    tags$iframe(seamless="seamless", src= "html/shm_shiny_manual.html", width='100%', height='100%') 
  })

  showModal(modal(title = HTML("<p>Please note that the forward/backward buttons in the browser are not supported.</p> <p>Instead, kindly use the <b>App's own buttons</b> located at the left of the header.</p>"), show=FALSE))
  observe({
    setBookmarkExclude(setdiff(names(input), c('upl-fileIn', 'dat-normDat', 'dat-log', 'dat-sig.min', 'dat-sig.max', 'dat-scl', 'dat-A', 'deg-A', 'scell-covisAuto-filBlkA', 'scell-covisMan-filBlkA', 'dat-P', 'deg-P', 'scell-covisAuto-filBlkP', 'scell-covisMan-filBlkP', 'dat-CV2', 'deg-CV1', 'scell-covisAuto-filBlkCV1', 'scell-covisMan-filBlkCV1', 'dat-CV2', 'deg-CV2', 'scell-covisAuto-filBlkCV2', 'scell-covisMan-filBlkCV2', 'shmAll-col.n', 'shmAll-genCon', 'shmAll-scale.shm', 'shmAll-scrollH')))
    # nas <- names(input); save(nas, file='nas')
  })
  # Automatically bookmark every time an input changes
  observe({
    lis.ipt <- reactiveValuesToList(input); session$doBookmark()
    # updateQueryString('?ids=null', mode = "push")
    # lapply(seq_along(lis.ipt), function(i) {if (length(lis.ipt[[i]])<1000) { print(lis.ipt[i]) }})
    # print(length(lis.ipt[['shmAll-dat-dt_rows_all']]))
    # print(getUrlHash()); print(getQueryString())
  })
 id.url.cnt <- reactiveValues(v=0)
 onBookmarked(function(url) {
   url1 <- gsub('(http.*\\?_inputs_&)(.*)', '\\1', url)
   url2 <- gsub('(http.*\\?_inputs_&)(.*)', '\\2', url)
   id.url <- lis.url$par$ids; 
   # print(list('sel1', id.url, ids$sel, id.url.cnt$v))
   if (!check_obj(id.url) | id.url.cnt$v!=0) { # If not refreshing the browser or not using URL to start this app in the browser search bar, lis.url$par/id.url will be NULL.
     id.url <- 'ids=null&'; sel <- ids$sel
     if (check_obj(sel)) id.url <- paste0('ids="', paste0(sel, collapse=','), '"&') 
     # print(list('sel2', id.url, ids$sel, id.url.cnt$v))
   } else if (check_obj(id.url) & id.url.cnt$v==0) {
     if (check_obj(id.url)) ids$sel <- id.url
     id.url <- paste0('ids="', paste0(id.url, collapse=','), '"&')
     id.url.cnt$v <- 1
     # print(list('url', id.url, ids$sel, id.url.cnt$v))
   }
   updateQueryString(paste0(url1, id.url, url2))
 })
 # onBookmarked(updateQueryString)
  onStop(function() { ggly_rm(); vdo_rm(); data_mining_rm()
    message("Session stopped! \n") 
  })

}





