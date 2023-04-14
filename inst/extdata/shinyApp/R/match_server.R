# Match spatial features between data and aSVG.
match_server <- function(id, sam, tab, upl.mod.lis, covis.man=NULL, col.idp=FALSE, session) {
  moduleServer(id, function(input, output, session) {
  observeEvent(input$matHelp, {
    showModal(
    div(id = 'matchHel',
    modalDialog(title = HTML('<strong><center>Matching spatial features</center></strong>'),
      div(style = 'overflow-y:scroll;overflow-x:scroll',
      HTML('<img src="image/match.jpg">'),
    ))))
  })
  ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg
  # renderUI: if the tab/page containing uiOutput('svg') is not active/clicked, the input$svg on the server side is NULL. To avoid this, the ui side should have "selectInput".
  output$svgs <- renderUI({
    # When customCovisData is selected, matching is disabled in SHM.
    if(id!='rematchCell' & grepl(na.sgl, ipt$fileIn)) return()
    ns <- session$ns; # nas <- c(names(cfg$pa.svg.reg), names(cfg$svg.def))
    selectInput(ns('svg'), label='Choose an aSVG to match', choices=cfg$na.def, selected=ipt$fileIn)
  })

  output$match.but <- renderUI({
    # When customCovisData is selected, matching is disabled in SHM.
    if(id!='rematchCell' & grepl(na.sgl, ipt$fileIn)) return()
    ns <- session$ns
    actionButton(ns("match"), 'Run', icon=icon("sync"), style=run.col)
  })
  output$match.reset <- renderUI({
    # When customCovisData is selected, matching is disabled in SHM.
    if(id!='rematchCell' & grepl(na.sgl, ipt$fileIn)) return()
    ns <- session$ns
    actionButton(ns("matReset"), 'Reset', icon=icon("sync"))
  })
  ft.reorder <- reactiveValues(ft.dat = NULL, ft.svg = NULL, ft.rematch = NULL)
  # Used to extract coordinates for SHMs.
  svg.na.rematch <- reactiveValues(svg.path=NULL, svg.na=NULL)
  observeEvent(list(ipt$fileIn, ipt$geneInpath, ipt$sglCell), {
    ft.reorder$ft.dat <- ft.reorder$ft.rematch <- NULL
  }) # Setting NULL: should not be merged with ipt$svgInpath1, ipt$svgInpath2.
  observeEvent(list(ipt$fileIn, ipt$svgInpath1, ipt$svgInpath2), {
    ft.reorder$ft.svg <- ft.reorder$ft.rematch <- NULL
    svg.na.rematch$svg.path <- svg.na.rematch$svg.na <- NULL
  })  # Setting NULL: should not be merged with ipt$geneInpath, ipt$sglCell
  # If multiple aSVGs are re-matched to a data, features in all aSVGs are extracted and put together on the user interface. Then each aSVG is rematched to data according to the rematch list sequentially.
  # Extract features in data and aSVG and create user interface to host these features.
  observeEvent(list(input$svg, ipt$svgInpath1, ipt$svgInpath2), {
    cat('Re-matching: features in aSVG ... \n')
    svg.in <- input$svg; svg.def <- cfg$svg.def
    if (is.null(ipt$fileIn)|is.null(svg.def)|is.null(svg.in)) return()
    if (svg.in!='uploaded') {
      svg.path <- svg.def[[svg.in]]
      if ('data_shm.tar' %in% basename(svg.path)) {
        svg.path <- read_hdf5('data/data_shm.tar', svg.in)[[1]]$svg
        validate(need(try(file.exists(svg.path)), svg.path))
      }
      svg.na <- basename(svg.path)
      # Single or multiple svg paths are treated same way.
      # lis <- svg_pa_na(svg.def[[svg.in]], cfg$pa.svg.upl, raster.ext)
      # output$msg.match <- renderText({ validate(need(try(!is.character(lis)), lis)) })
      # validate(need(try(!is.character(lis)), lis))
      # svg.path <- lis$svg.path; svg.na <- lis$svg.na
    } else { # aSVGs uploaded in regular files, not tar.
      svg.path <- cfg$pa.svg.reg[[1]]
      svg.na <- vapply(strsplit(svg.path, '/'), function(x) {x[length(x)]}, character(1))
    }
    cat('Access aSVG path for re-matching ... \n')
    # If multiple svgs, check suffixes.
    lis <- svg_suffix(svg.path, svg.na, raster.ext)
    validate(need(try(!is.character(lis)), lis))
    svg.path <- lis$svg.path; svg.na <- lis$svg.na
    svg.na.rematch$svg.path <- svg.path; svg.na.rematch$svg.na <- svg.na
  withProgress(message="Spatial heatmap re-matching: ", value=0, {  
    incProgress(0.5, detail="parsing aSVGs, please wait ...") 
    sf.all <- NULL
    cat('Extract all spatial features for re-matching ... \n')
    # Whether a single or multiple SVGs, all are returned a coord.
    svg.paths <- grep('\\.svg$', svg.path, value=TRUE)
    svgs <- read_svg_m(svg.path=svg.paths)
    validate(need(!is.character(svgs), svgs))
    sf.all <- unique(unlist(lapply(seq_along(svgs), function(x) { svg_separ(svg.all=svgs[x])$tis.path })))
  })
  # paths and groups are dropped to bottom. 
  # Matching samples are raised to top.
  pas.idx <- grepl('^path\\d+|^g\\d+', sf.all)
  sf.all <- c(sf.all[!pas.idx], sf.all[pas.idx])
  ft.reorder$ft.svg <- sf.all; cat('Done! \n')
  })

  observeEvent(list(ft.reorder$ft.svg, col.idp), {
    ft.svg <- ft.reorder$ft.svg; bulk <- covis.man$bulk
    covisGrp <- covis.man$covisGrp
    if (!check_obj(list(ft.svg, col.idp, bulk))) return()
    if (col.idp==TRUE) {
      # In covis independent coloring, fts abesent in data are excluded, since even if matched with cell groups, they will be transparent.
      ft.dat.blk <- unique(colData(bulk)[, covisGrp][bulk$bulkCell %in% 'bulk'])
      ft.svg <- intersect(ft.svg, ft.dat.blk)
      if (length(ft.svg)==0) ft.svg <- NULL 
      ft.reorder$ft.svg <- ft.svg
    }
  })
  observeEvent(list(sam(), input$svg, ipt$svgInpath1, ipt$svgInpath2), {
    cat('Re-matching: features in data ... \n')
    if (is.null(ipt$fileIn)|is.null(cfg$svg.def)|is.null(input$svg)) return()
    sams <- sam(); if (is.null(sams)) return()
    ft.reorder$ft.dat <- unique(sams); cat('Done! \n')
  })
  clean <- reactiveValues()
  observeEvent(covis.man$covis.type, {  
    covis.type <- covis.man$covis.type
    if (is.null(covis.type)) return()
    if (covis.type %in% c('toBulk','toCell') & 'no' %in% clean$val) {
      # Set NULL to from.ft in last matching.
      for (i in ft.reorder$ft.dat) {
        if (length(input[[i]])>0) {
          runjs(paste0("setipt('", session$ns(i), "', null)"))
          # runjs(paste0("Shiny.onInputChange('", session$ns(i), "', null)"))
        }; clean$val <- 'yes'
      }
    }
  })
  observeEvent(input$matReset, {
    ft.dat <- ft.reorder$ft.dat; if (is.null(ft.dat)) return()
      for (i in ft.reorder$ft.dat) {
        # Set NULL to from.ft in last matching.
        if (length(input[[i]])>0) {
          # Inside this observeEvent, input[[session$ns(i)]] will not change. It is NULL only outside this observeEvent.
          runjs(paste0("setipt('", session$ns(i), "', null)"))
          #runjs(paste0("Shiny.onInputChange('", session$ns(i), "', null)"))  
        }; clean$val <- 'yes'
      }
  })
  # Inside Shiny modules, change input values: always use session$ns; access input values: do not use session$ns, e.g. input$test.
  # observe({ runjs(paste0("setipt('", session$ns('test'), "', 100)")) })

  output$ft.match <- renderUI({
    cat('Re-matching: preparing interface of data/aSVG features ... \n')
    input$matReset
    # When customCovisData is selected, matching is disabled in SHM.
    if(id!='rematchCell' & grepl(na.sgl, ipt$fileIn)) return()
    ns <- session$ns; to.ft <- ft.reorder$ft.svg
    from.ft <- ft.reorder$ft.dat
    if (is.null(to.ft)|is.null(from.ft)) return()
    to.div.id='ftSVG'; to.div.tit='Features in aSVG'
    from.div.tit='Features in data'
    dimred <- covis.man$dimred

    if (!is.null(dimred)) {
      covis.type <- covis.man$covis.type
      if (covis.type %in% 'toBulk') {
        from.div.tit='Cell groups'
      } else if (covis.type %in% 'toCell') {
        covisGrp <- covis.man$covisGrp
        to.ft <- unique(colData(dimred)[, covisGrp])
        to.div.id='cellGroup'; to.div.tit='Cell groups'
        from.div.tit='Bulk tissues'
      }
    }
    frow <- match_interface(to.ft=to.ft, to.div.id=to.div.id, to.div.tit=to.div.tit, from.ft=from.ft, from.div.tit=from.div.tit, ns)
    clean$val <- 'no'; cat('Done! \n'); frow 
  })
  observeEvent(list(input$match), {
    if (is.null(ft.reorder$ft.dat) | is.null(ft.reorder$ft.svg)) return()
    ft.dat <- ft.reorder$ft.dat
    lis0 <- lapply(ft.dat, function(x) input[[x]])
    names(lis0) <- ft.dat 
    lis0 <- lis0[lapply(lis0, function(x) length(x)>0)==TRUE]
    # save(lis0, file='lis0')
    if (length(lis0)==0) showModal(modal(msg ='No spatial features are matched!')); validate(need(try(length(lis0)>0), ''))
    ft.reorder$ft.rematch <- lis0
  })
  observeEvent(list(input$matReset), { ft.reorder$ft.rematch <- list() })

  but.match <- reactiveValues()
  observe({
    match <- input$match; Tab <- tab$val
    if (is.null(match)|is.null(Tab)) return()
    # if (match==0|!Tab %in% c('scell')) return(); 
    but.match$val <- match
  })
  onBookmark(function(state) { state })
  return(list(svg.na.rematch=svg.na.rematch, ft.reorder=ft.reorder, but.match=but.match))

})}
