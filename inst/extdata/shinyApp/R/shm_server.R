
# Module for plotting SHMs.
shm_server <- function(id, sch, lis.url, url.id, tab, upl.mod.lis, dat.mod.lis, sch.mod.lis, scell.mod.lis, dim.mod.lis, deg.mod, prt=NULL) {  
  moduleServer(id, function(input, output, session) {
    message('SHM module starts ... ')
    library(magick); ns <- session$ns 
    ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg
  # The reactive type in and outside module is the same: sear is a reactiveValue in and outside module; geneIn is reactive expression in and outside module. "geneIn()" is accessing the content of a reactive expression, and loses the "reactive" attribute.
  # As long as the content of reactiveValues (col.reorder$col.na.re) is not accessed, the operation does not need to be inside reactive environment (observe).
  se.scl <- dat.mod.lis$se.scl
  se.scl.sel <- dat.mod.lis$se.scl.sel
  con.na <- dat.mod.lis$con.na
  con.na.cell <- dat.mod.lis$con.na.cell
  ipt.dat <- reactiveValues()
  ipt.dat$dat <- dat.mod.lis$ipt.dat; sear <- dat.mod.lis$sear
  geneIn <- dat.mod.lis$geneIn
  scaleDat <- dat.mod.lis$scaleDat
  log <- dat.mod.lis$log; A <- dat.mod.lis$A
  search.but <- dat.mod.lis$search.but
  sig.but <- dat.mod.lis$sig.but
  ids <- sch.mod.lis$ids
  gID <- reactiveValues(geneSel="none", new=NULL, all=NULL)
  observeEvent(scell.mod.lis$sce.upl$covis.type, {
    covis.type <- scell.mod.lis$sce.upl$covis.type
    if (!check_obj(list(covis.type))) return()
    cho <- c('Cell-by-value'='idp', 'Fixed-by-group'='fixed')
    if (covis.type %in% c('toBulk', 'toBulkAuto')) cho <- c(cho, 'Cell-by-group'='cellgrp')
    if (covis.type %in% c('toCell', 'toCellAuto')) cho <- c(cho, 'Feature-by-group'='ftgrp')
    updateSelectizeInput(session, inputId='profile', choices=cho, selected='idp')
  })
  observe({ ipt$geneInpath; ipt$fileIn; gID$geneSel <- "none"
    gID$new <- gID$all <- NULL
  })

  observeEvent(list(dat.mod.lis$sn$input$selRow, deg.mod$input$eSHMBut, scell.mod.lis$covis.man$match.mod.lis$but.match$val, scell.mod.lis$covis.auto$but.covis), {
    updateTabsetPanel(session, inputId="shmMhNet", selected='shm1')
  })

  searbox <- reactiveValues()
  observe({ 
    shinyjs::hide(id = "searshm"); searbox$v <- 'hide'
    hideTab(inputId="shmPar", target="relasize") 
    hideTab(inputId="shmPar", target="rematch") 
  })

  observe({
    if (length(ids$sel)==0) {
      disable(selector='a[data-value="interTab"]')
      disable(selector='a[data-value="vdoTab"]')
    } else { 
      enable(selector='a[data-value="interTab"]')
      enable(selector='a[data-value="vdoTab"]')
    }
  })
  observe({ if (is.null(se.scl())) gID$geneSel <- "none" })
  # To make the "gID$new" and "gID$all" updated with the new "input$fileIn", since the selected row is fixed (3rd row), the "gID$new" is not updated when "input$fileIn" is changed, and the downstream is not updated either. The shoot/root examples use the same data matrix, so the "gID$all" is the same (pre-selected 3rd row) when change from the default "shoot" to others like "organ". As a result, the "gene$new" is null and downstream is not updated. Also the "gene$new" is the same when change from shoot to organ, and downstream is not updated, thus "gene$new" and "gene$all" are both set NULL above upon new "input$fileIn".  

  init <- reactiveValues(but=0, new=0)
  # observeEvent(ids$but, { init$but <- ids$but })
  # observeEvent(session, { init$n <- init$n+1; print(init$n)})
  rna.fil <- reactiveValues(val=NULL)
  observe({ # Filtered data.
    if (length(ids$sel)==0) return(); if (ids$sel[1]=='') return()
    if (!check_obj(list(se.scl()))) return()
    rna.fil$val <- rownames(se.scl())
  })
  # The on-start ID processing is controlled by 0 and 1 states.
  observeEvent(list(session, ipt$fileIn, rna.fil$val), { init$but <- 0 })
  # init$but: triggers id update if the data is same but aSVG is different such as the rice shiny app. 
  observeEvent(list(ids$sel, rna.fil$val, init$but), { # All on-start and on-start similar IDs. Eg. the default first ID after data is filtered.
    cat('New file:', ipt$fileIn, '\n')
    if (length(ids$sel)==0 ) return()
    if ('' %in% ids$sel[1]|init$but>0) return()
    # Avoid selected genes are from last data while new data is used.
    if (!all(ids$sel %in% rna.fil$val)) return()
    # Avoid multiple selected rows from last input$fileIn. Must be behind gID$geneSel. 
    if (length(ids$sel)>1 & is.null(lis.url$par)) return()
    gID$geneSel <- ids$sel; if (!check_obj(gID$geneSel)) gID$geneSel <- 'none'
    gID$all <- gID$new <- NULL
    gID$new <- setdiff(gID$geneSel, gID$all); gID$all <- c(gID$all, gID$new)
    init$new <- 1 # Indicates ids are processed on-start, and no need to re-process in below.
    init$but <- 1
    cat('New ID:', gID$new, 'Selected ID:', gID$geneSel, 'All ID:', gID$all, '\n')
  })

  observeEvent(list(ids$sel), { # Selected IDs after the landing page.
    cat('Confirm selection ... \n')
    # In case init$new is not assigned to 0, since gg_shm processe is not triggered.
    if (!is.null(gID$all) & !is.null(gID$new) & !is.null(ids$sel) & init$new==1) {
      if (!all(ids$sel %in% gID$all)) init$new <- 0
    }
    # Ensure executions only after the landing page.
    if (init$but==0|init$new==1) return()
    if (is.null(ids$but.sgl) & is.null(ids$but.mul) & !'hide' %in% searbox$v) return()
    if (length(ids$sel)==0) return()
    if (ids$sel[1]=='') return()
    # Avoid selected genes are from last data while new data is used.
    if (!all(ids$sel %in% rna.fil$val)) return()
    gID$geneSel <- unique(ids$sel)
    gID$new <- setdiff(gID$geneSel, gID$all); gID$all <- c(gID$all, gID$new) 
    cat('New ID:', gID$new, 'Selected ID:', gID$geneSel, 'All ID:', gID$all, '\n'); cat('Done! \n')
    # if (any(is.na(gID$geneSel))) gID$geneSel <- "none"
  })
  geneV <- reactive({
    cat('All colour key values ... \n')
    se.scl <- se.scl(); se.scl.sel <- se.scl.sel()
    if (!check_obj(list(se.scl, se.scl.sel))) req('')
    assay <- assay(se.scl); assay.sel <- assay(se.scl.sel)
    validate(need(!any(is.na(gID$geneSel)) & gID$geneSel[1]!='', ''))
    # if (any(is.na(gID$geneSel))) return()
    if (sum(gID$geneSel[1]!='none')==0) return(NULL)
    if (input$ckeyV=="Selected rows" & length(ids$sel)==0) return(NULL)
    if (input$ckeyV=="Selected rows") assay <- assay.sel
    if (!all(gID$geneSel %in% rownames(assay))) return()
    bar.v <- seq(min(assay), max(assay), len=1000) # len must be same with that from the function "spatial_hm()". Otherwise the mapping of a gene value to the colour bar is not accurate.
    thr <- c(min(assay), max(assay))
    cat('Done! \n'); return(list(bar.v=bar.v, thr=thr))
  })

  observeEvent(input$colorOp, ignoreNULL=FALSE, {
    colorOp <- input$colorOp
    if ('custom' %in% colorOp) shinyjs::show(id = "colCus") else hide(id = "colCus")
  })

  col.sch <- eventReactive(list(input$colorOp, input$col.but), {
    cat('Color scheme ... \n') 
    colorOp <- input$colorOp; req(check_obj(list(colorOp)))
    but <- input$col.but; color <- input$color
    if ('custom' %in% colorOp) {
      req(check_obj(list(but, color)))
    } else color <- colorOp
    col <- gsub(' |\\.|-|;|,|/', '_', color)
    col <- strsplit(col, '_')[[1]]; col <- col[col!='']
    lgc <- length(col)>=2; if (!lgc) { 
        showModal(modal(msg='At least 2 colors are needed!', easyClose=TRUE))
    }; req(lgc)
    col1 <- col[!col %in% colors()]; lgc <- length(col1)==0
    if (!lgc) { 
        msg <- paste0('Invalide colors: ', paste0(col1, collapse=',')) 
        showModal(modal(msg=msg, easyClose=TRUE))
    }; req(lgc)
    cat('Done! \n'); col
  }); observe({ col.sch() })
  
  color <- reactiveValues(col="none")
  observe({
    cat('Initial color code for color key ... \n')
    session # Avoid color$col is "none", sine new session triggers color <- reactiveValues(col="none")
    lis.par <- cfg$lis.par; req(check_obj(lis.par))
    col0 <- lis.par$shm.img['color', 'default']
    col.but <- input$col.but
    if (is.null(col.but)|is.null(col0)|gID$geneSel[1]=='none') return()
    if(col.but==0) color$col <- colorRampPalette(col_sep(col0))(length(geneV()$bar.v))
    cat('Done! \n')
  })

  # As long as a button is used, observeEvent should be used. All variables inside 'observeEvent' trigger code evaluation, not only 'eventExpr'.  
  observeEvent(list(col.sch(), geneV()), {
    cat('Customized color code for color key ... \n')
    chk <- check_exp(check_obj(list(col.sch(), geneV())))
    req(TRUE %in% chk)
    color$col <- colorRampPalette(col.sch())(length(geneV()$bar.v))
    cat('Done! \n')
  })
  # Should not be the same with profile line graph, since the latter only reflect selected genes not all genes. 
  x.title <- reactiveValues(val='')
  observe({
    thr <- geneV()$thr; if (is.null(thr)) return()
    scale.dat <- scaleDat(); if (is.null(scale.dat)) return()
    if (!is.null(scale.dat)) if (scale.dat=='No') title <- 'No scaling' else if (scale.dat=='Row') title <- 'Scaling by row' else if (scale.dat=='Selected') title <- 'Scaling selected rows together' else if (scale.dat=='All') title <- 'Scaling all rows together' else title <- ''
    x.title$val <- paste0(title, ' (', round(thr[1], 2), '-', round(thr[2], 2), ')')
  })
  shm.bar <- reactive({
    cat('Colour key ... \n')
    if (is.null(gID$all)) return(); se.scl <- se.scl()
    if ((ipt$fileIn %in% cfg$na.def & !is.null(se.scl))|(ipt$fileIn %in% cfg$na.cus & (!is.null(ipt$svgInpath1)|!is.null(ipt$svgInpath2)) & !is.null(se.scl))) {
      bar.v <- geneV()$bar.v
      if (length(color$col=="none")==0|is.null(bar.v)) return(NULL)
      withProgress(message="Color key: ", value = 0, {
        incProgress(0.75, detail="plotting ...")
        cell <- scell.mod.lis$sce.upl$cell
        if (!is.null(cell) & input$profile=='fixed') return(ggplot())
        cs.g <- col_bar(geneV=bar.v, cols=color$col, width=1, x.title=x.title$val, x.title.size=10)
        incProgress(0.1, detail="plotting ...")
        # save(cs.g, file='cs.g')
        cat('Done! \n'); return(cs.g)
      })
    }
  })
  # One output can only be used once in ui.R.
  output$bar1 <- bar2 <- renderPlot({ if (!is.null(shm.bar)) shm.bar() })
  # output$bar2 <- renderPlot({ if (!is.null(shm.bar)) shm.bar() })
  observe({
    glyBut <- input$glyBut
    if (is.null(glyBut)) output$bar2 <- NULL else if (glyBut==0) output$bar2 <- NULL else output$bar2 <- bar2
  }) 

  svg.path <- reactive({ # Organise svg name and path in a nested list.
    message('Access aSVG path ...')
    fileIn <- ipt$fileIn; svg.def <- cfg$svg.def
    req(check_obj(svg.def))
    if (fileIn %in% cfg$na.cus) {
      if (is.null(ipt$svgInpath2)) svgIn.df <- ipt$svgInpath1 else svgIn.df <- ipt$svgInpath2
      svg.path <- svgIn.df$datapath; svg.na <- svgIn.df$name
    } else {
      # Extract svg path and name: single or multiple svg paths are treated same way.
      svg.path <- svg.def[[fileIn]]
      # pa.svg.upl <- cfg$pa.svg.upl
      if ('data_shm.tar' %in% basename(svg.path)) {
        svg.path <- read_hdf5('data/data_shm.tar', fileIn)[[1]]$svg
        validate(need(try(file.exists(svg.path)), svg.path))
      }
      svg.na <- basename(svg.path)
      # Check if SVGs are paired with templates of raster images.
      if (any(!grepl('\\.svg$', svg.na))) svg_raster(svg.na, raster.ext)   
    }
    # If multiple svgs/templates (treated same way), check suffixes.
    lis <- svg_suffix(svg.path, svg.na, raster.ext)
    validate(need(try(!is.character(lis)), lis))
    message('Done!'); return(lis)
  })

  sam <- reactive({
    se.scl <- se.scl(); if (!check_obj(list(se.scl))) req('')
    prof <- input$profile; blk <- scell.mod.lis$sce.upl$bulk
    covis.type <- scell.mod.lis$sce.upl$covis.type
    if ('idp' %in% prof & !is.null(blk) & 'bulkCell' %in% colnames(colData(se.scl))) { 
      if ('toBulk' %in% covis.type) se.scl <- subset(se.scl, , bulkCell=='cell')
      if ('toCell' %in% covis.type) se.scl <- subset(se.scl, , bulkCell=='bulk')
    } 
    cname <- colnames(se.scl); idx <- grep("__", cname)
    c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\1", c.na) else return() 
  })

  svg.na.remat <- reactiveValues(svg.path=NULL, svg.na=NULL)
  ft.reord <- reactiveValues(ft.dat = NULL, ft.svg = NULL, ft.rematch = NULL)

  but.match <- reactiveValues(); match.mod.lis <- NULL
    # Put the code belew in observe below: leads to infinite circles.
    # if (ipt$fileIn!='customCovisData'): if condition cannot supress module execution.
    match.mod.lis <- match_server('rematch', sam, tab, upl.mod.lis)
  observeEvent(list(match.mod.lis$svg.na.rematch$svg.path, match.mod.lis$svg.na.rematch$svg.na, match.mod.lis$ft.reorder$ft.dat, match.mod.lis$ft.reorder$ft.svg, match.mod.lis$ft.reorder$ft.rematch, match.mod.lis$but.match$val), ignoreNULL = FALSE, { # Rematch in bulk data.
    # if (is.null(match.mod.lis$ft.reorder$ft.rematch)) return()
    # svg.na.rematch <- match.mod.lis$svg.na.rematch: does not update svg.na.rematch outside "observe", so svg.path and svg.na are updated separately.
    svg.na.remat$svg.path <- match.mod.lis$svg.na.rematch$svg.path
    svg.na.remat$svg.na <- match.mod.lis$svg.na.rematch$svg.na

    ft.reord$ft.dat <- match.mod.lis$ft.reorder$ft.dat
    ft.reord$ft.svg <- match.mod.lis$ft.reorder$ft.svg
    ft.reord$ft.rematch <- match.mod.lis$ft.reorder$ft.rematch
    but.match$val <- match.mod.lis$but.match$val
  })

  cell.match <- reactiveValues()
  observe({ cell.match$val <- scell.mod.lis$covis.man$match.mod.lis })

  observeEvent(list(cell.match$val$svg.na.rematch$svg.path, cell.match$val$svg.na.rematch$svg.na, cell.match$val$ft.reorder$ft.dat, cell.match$val$ft.reorder$ft.svg, cell.match$val$ft.reorder$ft.rematch, cell.match$val$but.match$val), ignoreNULL = FALSE, {
    # Rematch in single cell data.
    match.lis <- cell.match$val
    # if (is.null(match.lis$ft.reorder$ft.rematch)) return()
    svg.na.remat$svg.path <- match.lis$svg.na.rematch$svg.path
    svg.na.remat$svg.na <- match.lis$svg.na.rematch$svg.na

    ft.reord$ft.dat <- match.lis$ft.reorder$ft.dat
    ft.reord$ft.svg <- match.lis$ft.reorder$ft.svg
    ft.reord$ft.rematch <- match.lis$ft.reorder$ft.rematch
    but.match$val <- match.lis$but.match$val
  })
  observeEvent(list(ipt$fileIn, ipt$geneInpath, ipt$sglCell, ipt$svgInpath1, ipt$svgInpath2), {
    svg.na.remat$svg.path <- svg.na.remat$svg.na <- NULL
  })
  # Reactive object in "observeEvent" is not accessible outside observeEvent. The solution is eventReactive. 
  # svg.path1 stores the final svg path/na after re-matching, and will be used in SHMs.
  svg.path1 <- reactive({
    if (!is.null(svg.na.remat$svg.path) & !is.null(svg.na.remat$svg.na)) { svg.path <- svg.na.remat$svg.path; svg.na <- svg.na.remat$svg.na } else { svg.path <- svg.path()$svg.path; svg.na <- svg.path()$svg.na }
    return(list(svg.path=svg.path, svg.na=svg.na))
  })

  # cna.match <- reactiveValues(cna=NULL)

  svgs <- reactive({
    cat('Reading aSVGs ... \n')
    fileIn <- ipt$fileIn; geneInpath <- ipt$geneInpath; dimName <- ipt$dimName
    svgInpath1 <- ipt$svgInpath1; svgInpath2 <- ipt$svgInpath2
    if (fileIn =='customBulkData') {
      if (is.null(dimName)) return()   
      if (is.null(geneInpath) | dimName == "None") return()
    }
    if ((fileIn %in% cfg$na.cus & 
    (!is.null(svgInpath1)|!is.null(svgInpath2)))|any(fileIn %in% cfg$na.def)) {
      withProgress(message="Spatial heatmap: ", value=0, {
      incProgress(0.5, detail="reading aSVG ...")
      svg.path <- svg.path1()$svg.path; svg.na <- svg.path1()$svg.na
      # Whether a single or multiple SVGs, all are returned a coord.
      svg.paths <- grep('\\.svg$', svg.path, value=TRUE)
      raster.paths <- setdiff(svg.path, svg.paths)
      svgs <- read_svg_m(svg.path=svg.paths, raster.path=raster.paths)
      validate(need(!is.character(svgs), svgs))
      cat('Done! \n'); return(svgs)
      })
    }
  })

  observe({
    if (!is.null(ft.reord$ft.rematch)) return()
    ipt$fileIn; se.scl(); ipt$adj.modInpath; svgs(); input$lgdTog; input$scrollH; svgs <- svgs()
    ft.path.all <- NULL; for (i in seq_along(svgs)) { 
      ft.path.all <- c(ft.path.all, svg_separ(svgs[i])$tis.path)
    }
    # inline=TRUE should not be ignored in update.
    updateSelectizeInput(session, inputId='tis', choices=intersect(unique(sam()), unique(ft.path.all)))
  })
  observe({
    input$svg; svgs <- svgs(); ft.svg.reorder <- ft.reord$ft.svg; but.match$val
    if (is.null(svgs) | is.null(ft.svg.reorder)) return()
    ft.path.all <- NULL; for (i in seq_along(svgs)) { 
      ft.path.all <- c(ft.path.all, svg_separ(svgs[i])$tis.path)
    }
    # inline=TRUE should not be ignored in update.
    updateSelectizeInput(session, inputId='tis', choices=intersect(unique(ft.svg.reorder), unique(ft.path.all)))
  })
  tis.trans <- reactiveValues()
  observeEvent(input$transBut, { tis.trans$v <- input$tis })
  con <- eventReactive(se.scl(), {
    cat('All variables ... \n'); se.scl <- se.scl() 
    if (!check_obj(list(se.scl))) req('')
    cname <- colnames(se.scl); idx <- grep("__", cname)
    c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) cons <- gsub("(.*)(__)(.*$)", "\\3", c.na) else cons <- NULL 
    cat('Done! \n'); cons
  })

  # General selected gene/condition pattern.
  pat.con <- reactive({
   con <- con(); if (!check_obj(list(con))) return()
   con.uni <- unique(con); if (is.null(con.uni)) return()
   paste0(con.uni, collapse='|') 
  })
  pat.gen <- reactive({ if (is.null(gID$geneSel)) return(); if (gID$geneSel[1]=='none') return(NULL);  paste0(gID$geneSel, collapse='|') })

  pat.all <- reactive({ if (is.null(pat.con())|is.null(pat.gen())) return(NULL); paste0('(', pat.gen(), ')_(', pat.con(), ')') })

  # SHM ggplots, grobs legends are stored in gg.all, grob.all, lgd.all respectively for various purposes. grob.gg.all is used in relative scale of multiple SVGs, and the rescaled SVGs are stored in gg.all/grob.all finally. 
  shm <- reactiveValues(grob.all=NULL, grob.all1=NULL, gg.all=NULL, gg.all1=NULL, lgd.all=NULL, grob.gg.all = NULL)
  observeEvent(list(ipt$fileIn, ipt$geneInpath, ipt$sglCell, ipt$svgInpath1, ipt$svgInpath2), { shm$grob.all <- shm$grob.all1 <- shm$gg.all1 <- shm$gg.all <- shm$lgd.all <- shm$lgd.grob.all <- shm$gcol.all <- shm$gcol.lgd.all <- shm$grob.gg.all <- NULL 
  })
  raster.par <- reactiveValues(over='Yes', coal='No', alp=NULL)
  # Use observeEvent: use NULL to replace 'No' to avoid unnecessary trigering of gg_shm, since NULL does not triger observeEvent below.
  observeEvent(list(input$raster, input$coal, input$alpOverBut, svg.path1()$svg.na), {
    svg.na <- svg.path1()$svg.na
    if (is.null(svg.na)|is.null(input$raster)|is.null(input$coal)|is.null(input$alpOver)) return()
    # raster.par$over is NULL or Yes, not No.
    if (input$raster=='Yes' & any(!grepl('\\.svg$', svg.na))) raster.par$over <- 'Yes' else raster.par$over <- NULL
    if (is.null(raster.par$over)) raster.par$coal <- NULL else if (raster.par$over=='Yes') {
      if (input$coal=='Yes') raster.par$coal <- 'Yes' else raster.par$coal <- NULL
      raster.par$alp <- input$alpOver
    }
  })

  # Avoid repetitive computation under input$ckeyV=='All rows'.
  gs.new <- reactive({
    cat('New grob/ggplot: \n ')
    # print(list(is.null(svgs()), is.null(se.scl.sel()), gID$new, gID$all, ids$sel, color$col[1]))
    se.scl.sel <- se.scl.sel(); prof <- input$profile
    lis.par <- cfg$lis.par
    if (!check_obj(list(se.scl.sel, prof, lis.par))|is.null(con.na$v)) return()
    col.idp <- 'idp' %in% prof & grepl(na.sgl, ipt$fileIn)
    validate(
      need(!is.null(svgs()) & !is.null(se.scl.sel) & length(gID$new) > 0 & !is.null(gID$all) & length(ids$sel)>0 & color$col[1]!='none', '')
    )
    scale.shm <- input$scale.shm
    if (!is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    # If color key is build on selected rows, all SHMs should be computed upon selected rows are changed. This action is done through a separate observeEvent triggered by gID$geneSel. So in this "reactive" only one gene is accepted each time.
    # Only works at "Selected rows" and one gene is selected, i.e. when the app is launched.
    # print(list('new', ids$but.sgl, gID$geneSel, gID$new))
    if (length(ids$but.sgl)==0 & length(ids$but.mul)==0 & !'hide' %in% searbox$v) return()
    if (length(url.id$sch.mul)==0|length(url.id$sch.sgl)==0 & !'hide' %in% searbox$v) return()
    urlID <- 'null'
    if (!'hide' %in% searbox$v) {
      if (url.id$sch.sgl[1]!='null') urlID <- url.id$sch.sgl else if (url.id$sch.mul[1]!='null') urlID <- url.id$sch.mul
    }
    # if (length(urlID)==0) return()
    if (input$ckeyV=="Selected rows") ID <- gID$geneSel else if (all(sort(urlID)==sort(gID$geneSel))) ID <- gID$geneSel else if (input$ckeyV=="All rows") ID <- gID$new else return()
    # Works all the time as long as "All rows" selected.
    # if (input$ckeyV!="All rows") return() 
    # ID <- gID$new
    if (is.null(ID)) return()
    # if (length(gID$new)>1|length(ID)>1|ID[1]=='none') return()
    if (ID[1]=='none') return()
    # Avoid repetitive computation.  
    pat.new <- paste0('^(', paste0(ID, collapse='|'), ')_(', pat.con(), ')_\\d+$')
    if (any(grepl(pat.new, names(shm$grob.all)))) return()
    withProgress(message="Spatial heatmap: ", value=0, { 
      incProgress(0.25, detail="preparing data ...")
      gene <- assay(se.scl.sel); 
      # When input$fileIn updates, ID is from last session while gene is from new session.
      if (!all(ID %in% rownames(gene))) return()
      if (is.null(raster.par$coal)) charcoal <- FALSE else if (raster.par$coal=='Yes') charcoal <- TRUE else if (raster.par$coal=='No') charcoal <- FALSE
      alp.over <- 1
      if (!is.null(raster.par$over)) if (raster.par$over=='Yes') alp.over <- raster.par$alp
      svgs <- svgs()
      lis.rematch <- ft.reord$ft.rematch
      ft.trans.shm <- NULL; ft.trans <- tis.trans$v
      covisGrp <- scell.mod.lis$sce.res()$covisGrp
      covis.type <- scell.mod.lis$sce.upl$covis.type
      method <- scell.mod.lis$sce.upl$method
      tar.cell.blk <- input$tarCellBlk
      if (grepl(na.sgl, ipt$fileIn)) { 
        if (is.null(covis.type)|is.null(method)|is.null(tar.cell.blk)) return()
        if ('man' %in% method) {
        dimred <- scell.mod.lis$covis.man$dimred
        bulk <- scell.mod.lis$covis.man$bulk
        if (!check_obj(list(dimred, covisGrp))) return()
        cell.all <- unique(colData(dimred)[, covisGrp]) 
        if (!is.null(bulk)) bulk.all <- unique(colData(bulk)[, covisGrp])
        ft.all <- unique(unlist(lapply(seq_along(svgs), function(i) attribute(svgs[i])[[1]]$feature)))
        covis.trans <- covis_trans(bulk.all=bulk.all, cell.all=cell.all, ft.all=ft.all, tar.bulk=tar.cell.blk, tar.cell=tar.cell.blk, lis.match=lis.rematch, covis.type=covis.type, col.idp=col.idp)
        ft.trans.shm <- covis.trans$ft.trans.shm
        lis.rematch <- covis.trans$lis.match 
        if (col.idp) gene <- gene[, grep('^bulk$', se.scl.sel$bulkCell), drop=FALSE] 
      } else if ('auto' %in% method) {
        res <- scell.mod.lis$covis.auto$res
        cell.all <- setdiff(unique(subset(res, , bulkCell=='cell')$assignedBulk), 'none')
        bulk.all <- unique(subset(res, , bulkCell=='bulk')$sample)
        covis.trans <- covis_trans(bulk.all=bulk.all, cell.all=cell.all, ft.all=NULL, tar.bulk=tar.cell.blk, tar.cell=tar.cell.blk, lis.match=NULL, covis.type=covis.type, col.idp=col.idp)
        ft.trans.shm <- covis.trans$ft.trans.shm
        lis.rematch <- NULL
        if (col.idp) gene <- gene[, grep('^bulk$', res$bulkCell), drop=FALSE] 
      }
      }
      svg.na <- names(svgs[, 'svg']) 
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      grob.all <- gg.all <- lgd.all <- lgd.grob.all <- gcol.all <- gcol.lgd.all <- grob.gg.all <- NULL
      for (i in seq_along(svgs)) {
        cat(ID, ' \n')
        size.key <- as.numeric(lis.par$legend['key.size', 'default']) 
        svg0 <- svgs[i]
        if (is.null(raster.par$over)) raster_pa(svg0)[[1]] <- NULL
        # Cores: the orders in svg.path(), names(svg.df.lis) are same.
        gg.lis <- gg_shm(svg.all=svg0, gene=gene, con.na=con.na$v, geneV=geneV()$bar.v, col.idp=col.idp, charcoal=charcoal, alpha.overlay=alp.over, ID=ID, cols=color$col, covis.type=covis.type, ft.trans=ft.trans, ft.trans.shm=ft.trans.shm, sub.title.size=input$title.size * scale.shm, legend.nrow=as.numeric(lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.width=input$line.size, line.color=input$line.color, lis.rematch = lis.rematch) # Only gID$new is used.
        msg <- paste0(svg.na[i], ': no common spatial features detected between data and aSVG!')
        if (is.null(gg.lis)) {
        showNotification(msg, duration=2, closeButton = TRUE)
        cat(msg, '\n'); 
        }
       validate(need(!is.null(gg.lis), msg)) 
       # Append suffix '_i' for the SHMs of ggplot under SVG[i], and store them in a list.
       ggs <- gg.lis$g.lis.all; names(ggs) <- paste0(names(ggs), '_', i)
       gg.all <- c(gg.all, ggs)
       gcols <- gg.lis$gcol.lis.all; names(gcols) <- paste0(names(gcols), '_', i)
       gcol.all <- c(gcol.all, gcols)
       # Store legend/legend colours of ggplot in a list.
       lgd.all <- c(lgd.all, list(gg.lis$g.lgd))
       gcol.lgd.all <- c(gcol.lgd.all, list(gg.lis$gcol.lgd))
       # Same names with ggs: append suffix '_i' for the SHMs of grob under SVG[i], and store them in a list.
       grob.lis <- grob_shm(ggs, cores=deter_core(1, svg.obj=svgs[, 'svg'][[i]])) 
       grob.all <- c(grob.all, grob.lis)
       lgd.grob.lis <- grob_shm(lgd.all, cores=deter_core(1, svg.obj=svgs[, 'svg'][[i]]), lgd.pos='bottom')
       lgd.grob.all <- c(lgd.grob.all, lgd.grob.lis)
       # All ggplots/grobs are stored in nested lists under each SVG for use in relatice scale. In above, all ggplots/grobs are stored in the same list with suffix '_i' to indicate SVGs.
        lis0 <- list(grob.lis = grob.lis, gg.lis = ggs, lgd.lis = gg.lis$g.lgd, lgd.grob=lgd.grob.lis[[1]], gcol.lis=gcols, gcol.lgd=gg.lis$gcol.lgd)
       grob.gg.all <- c(grob.gg.all, list(lis0)) 
     }
     names(lgd.all) <- names(lgd.grob.all) <- names(grob.gg.all) <- svg.na
     names(gcol.lgd.all) <- paste0('col_', svg.na)
     init$new <- 0 # Terminates gs.new.
     cat('Done! \n'); return(list(gg.all = gg.all, grob.all = grob.all, lgd.all = lgd.all, lgd.grob.all=lgd.grob.all, gcol.all=gcol.all, gcol.lgd.all=gcol.lgd.all, grob.gg.all = grob.gg.all))
    }) # withProgress

  })

  # Extension of 'observeEvent': any of 'input$log; tis.trans$v; input$col.but; input$ckeyV' causes evaluation of all code. 
  # tis.trans$v as an argument in "gg_shm" will not cause evaluation of all code, thus it is listed here.
  # Use "observeEvent" to replace "observe" and list events (input$log, tis.trans$v, ...), since if the events are in "observe", every time a new gene is clicked, "input$dt_rows_selected" causes the evaluation of all code in "observe", and the evaluation is duplicated with "gs.new".
  # Update SHMs, above theme().
  observeEvent(list(log(), tis.trans$v, color$col, sig.but(), input$ckeyV, scaleDat(), but.match$val, ft.reord$ft.rematch, input$line.size, input$line.color, raster.par$over, raster.par$coal, raster.par$alp, input$profile, input$tarCellBlk), {
    shm$grob.all <- shm$gg.all <- shm$lgd.all <- shm$lgd.grob.all <- shm$gcol.all <- shm$gcol.lgd.all <- shm$grob.gg.all <- NULL; gs.all <- reactive({ 
      cat('Updating all SHMs ... \n')
      se.scl.sel <- se.scl.sel(); prof <- input$profile
      lis.par <- cfg$lis.par
      if (!check_obj(list(se.scl.sel, prof, lis.par))|is.null(con.na$v)) req('')
      col.idp <- 'idp' %in% prof & grepl(na.sgl, ipt$fileIn)
      # print(list(is.null(svgs()), is.null(geneIn), ids$sel, color$col[1], gID$geneSel))
      if.con <- is.null(svgs())|is.null(se.scl.sel)| length(ids$sel)==0 |color$col[1]=='none'|gID$geneSel[1]=='none'
      if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
      scale.shm <- input$scale.shm
      if (!is.numeric(scale.shm)) return()
      if (scale.shm <= 0) return()
      withProgress(message="Spatial heatmap: ", value=0, {
      incProgress(0.25, detail="in progress ...")
      #if (input$ckeyV=="Selected rows") gene <- geneIn()[["df.aggr.tran"]][ipt.dat$dat$dt_rows_selected, ]
      #if (input$ckeyV=="All rows") gene <- geneIn()[["df.aggr.tran"]]
      gene <- assay(se.scl.sel)
      alp.over <- 1
      if (!is.null(raster.par$over)) if (raster.par$over=='Yes') alp.over <- raster.par$alp
      svgs <- svgs()
      lis.rematch <- ft.reord$ft.rematch
      ft.trans.shm <- NULL; ft.trans <- tis.trans$v
      covisGrp <- scell.mod.lis$sce.res()$covisGrp
      covis.type <- scell.mod.lis$sce.upl$covis.type
      method <- scell.mod.lis$sce.upl$method
      tar.cell.blk <- input$tarCellBlk
      if (grepl(na.sgl, ipt$fileIn)) { 
        if (is.null(covis.type)|is.null(method)|is.null(tar.cell.blk)) return()
        if ('man' %in% method) {
        dimred <- scell.mod.lis$covis.man$dimred
        bulk <- scell.mod.lis$covis.man$bulk
        if (!check_obj(list(dimred, covisGrp))) return()
        cell.all <- unique(colData(dimred)[, covisGrp]) 
        if (!is.null(bulk)) bulk.all <- unique(colData(bulk)[, covisGrp])
        ft.all <- unique(unlist(lapply(seq_along(svgs), function(i) attribute(svgs[i])[[1]]$feature)))
        covis.trans <- covis_trans(bulk.all=bulk.all, cell.all=cell.all, ft.all=ft.all, tar.bulk=tar.cell.blk, tar.cell=tar.cell.blk, lis.match=lis.rematch, covis.type=covis.type, col.idp=col.idp)
        ft.trans.shm <- covis.trans$ft.trans.shm
        lis.rematch <- covis.trans$lis.match 
        if (col.idp) gene <- gene[, grep('^bulk$', se.scl.sel$bulkCell), drop=FALSE] 
      } else if ('auto' %in% method) {
        res <- scell.mod.lis$covis.auto$res
        cell.all <- setdiff(unique(subset(res, , bulkCell=='cell')$assignedBulk), 'none')
        bulk.all <- unique(subset(res, , bulkCell=='bulk')$sample)
        covis.trans <- covis_trans(bulk.all=bulk.all, cell.all=cell.all, ft.all=NULL, tar.bulk=tar.cell.blk, tar.cell=tar.cell.blk, lis.match=NULL, covis.type=covis.type, col.idp=col.idp)
        ft.trans.shm <- covis.trans$ft.trans.shm
        lis.rematch <- NULL
        if (col.idp) gene <- gene[, grep('^bulk$', res$bulkCell), drop=FALSE] 
      }
      }
      svg.na <- names(svgs[, 'svg'])
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      grob.all <- gg.all <- lgd.all <- lgd.grob.all <- gcol.all <- gcol.lgd.all <- gg.grob.lis <- NULL
      for (i in seq_along(svgs)) { 
        if (is.null(raster.par$coal)) charcoal <- FALSE else if (raster.par$coal=='Yes') charcoal <- TRUE else if (raster.par$coal=='No') charcoal <- FALSE
        cat('All grob/ggplot:', gID$geneSel, ' \n')
        incProgress(0.75, detail=paste0('preparing ', paste0(gID$geneSel, collapse=';')))
        size.key <- as.numeric(lis.par$legend['key.size', 'default'])
        svg0 <- svgs[i]
        if (is.null(raster.par$over)) raster_pa(svg0)[[1]] <- NULL
        gg.lis <- gg_shm(svg.all=svg0, gene=gene, con.na=con.na$v, geneV=geneV()$bar.v, col.idp=col.idp, charcoal=charcoal, alpha.overlay=alp.over, ID=gID$geneSel, cols=color$col, covis.type=covis.type, ft.trans=ft.trans, ft.trans.shm=ft.trans.shm, sub.title.size=input$title.size * scale.shm, legend.nrow=as.numeric(lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.width=input$line.size, line.color=input$line.color, lis.rematch = lis.rematch) # All gene IDs are used.
        msg <- paste0(svg.na[i], ': no common spatial features detected between data and aSVG!')
        if (is.null(gg.lis)) {
        showNotification(msg, duration=2, closeButton = TRUE)
        cat(msg, '\n'); 
        }
        validate(need(!is.null(gg.lis), msg))
       # Append suffix '_i' for the SHMs of ggplot under SVG[i], and store them in a list.
       ggs <- gg.lis$g.lis.all; names(ggs) <- paste0(names(ggs), '_', i)
       gg.all <- c(gg.all, ggs) 
       gcols <- gg.lis$gcol.lis.all; names(gcols) <- paste0(names(gcols), '_', i)
       gcol.all <- c(gcol.all, gcols) 
       # Store legend/colours of ggplot in a list.
       lgd.all <- c(lgd.all, list(gg.lis$g.lgd))
       gcol.lgd.all <- c(gcol.lgd.all, list(gg.lis$gcol.lgd))
       # Same with ggs: append suffix '_i' for the SHMs of grob under SVG[i], and store them in a list.
       grob.lis <- grob_shm(ggs, cores=deter_core(1, svg.obj=svgs[, 'svg'][[i]]))
       grob.all <- c(grob.all, grob.lis)
       lgd.grob.lis <- grob_shm(lgd.all, cores=deter_core(1, svg.obj=svgs[, 'svg'][[i]]), lgd.pos='bottom')
       lgd.grob.all <- c(lgd.grob.all, lgd.grob.lis)
       # All ggplots/grobs are stored in nested lists under each SVG for use in relatice scale.
       lis0 <- list(grob.lis = grob.lis, gg.lis = ggs, lgd.lis = gg.lis$g.lgd, lgd.grob=lgd.grob.lis[[1]], gcol.lis=gcols, gcol.lgd=gg.lis$gcol.lgd)
       gg.grob.lis <- c(gg.grob.lis, list(lis0))
      }
     names(lgd.all) <- names(lgd.grob.all) <- names(gg.grob.lis) <- svg.na
     names(gcol.lgd.all) <- paste0('col_', svg.na)
     init$new <- 0 # Terminates gs.new.
     cat('Done! \n'); return(list(grob.all = grob.all, gg.all = gg.all, lgd.all = lgd.all, lgd.grob.all=lgd.grob.all, gcol.all=gcol.all, gcol.lgd.all=gcol.lgd.all, gg.grob.lis = gg.grob.lis))
     }) # withProgress
    }) # reactive
    shm$grob.all <- gs.all()$grob.all; shm$gg.all <- gs.all()$gg.all
    shm$lgd.all <- gs.all()$lgd.all; shm$lgd.grob.all <- gs.all()$lgd.grob.all;
    shm$gcol.all <- gs.all()$gcol.all
    shm$gcol.lgd.all <- gs.all()$gcol.lgd.all 
    shm$grob.gg.all <- gs.all()$gg.grob.lis
  }) # observeEvent
  # Avoid repetitive computation under input$ckeyV=='All rows'.
  observeEvent(list(gID$geneSel, se.scl.sel()), { 
    cat('Updating all SHMs caused by selected IDs ... \n')
    lis.par <- cfg$lis.par; req(check_obj(lis.par))
    if.con <- is.null(input$ckeyV)|gID$geneSel[1]=='none'|input$ckeyV=='All rows'
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    ID <- gID$geneSel
    shm$grob.all <- shm$gg.all <- shm$lgd.all <- shm$lgd.grob.all <- shm$gcol.all <- shm$gcol.lgd.all <- shm$grob.gg.all <- NULL; gs.all <- reactive({
     # print(list(ID, is.null(svgs()), is.null(geneIn()), ids$sel, color$col[1], class(color$col[1])))
      se.scl.sel <- se.scl.sel(); prof <- input$profile
      if (!check_obj(list(se.scl.sel, prof))|is.null(con.na$v)) req('')
      col.idp <- 'idp' %in% prof & grepl(na.sgl, ipt$fileIn)
      if.con <- is.null(svgs())|is.null(se.scl.sel)| length(ids$sel)==0 |color$col[1]=='none'
      if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
      scale.shm <- input$scale.shm
      if (!is.numeric(scale.shm)) return()
      if (scale.shm <= 0) return()
      withProgress(message="Spatial heatmap: ", value=0, {
      incProgress(0.25, detail="preparing data ...")
      gene <- assay(se.scl.sel); alp.over <- 1
      if (!is.null(raster.par$over)) if (raster.par$over=='Yes') alp.over <- raster.par$alp
      svgs <- svgs()
      lis.rematch <- ft.reord$ft.rematch
      ft.trans.shm <- NULL; ft.trans <- tis.trans$v
      covisGrp <- scell.mod.lis$sce.res()$covisGrp
      covis.type <- scell.mod.lis$sce.upl$covis.type
      method <- scell.mod.lis$sce.upl$method
      tar.cell.blk <- input$tarCellBlk
      if (grepl(na.sgl, ipt$fileIn)) { 
        if (is.null(covis.type)|is.null(method)|is.null(tar.cell.blk)) return()
        if ('man' %in% method) {
        dimred <- scell.mod.lis$covis.man$dimred
        bulk <- scell.mod.lis$covis.man$bulk
        if (!check_obj(list(dimred, covisGrp))) return()
        cell.all <- unique(colData(dimred)[, covisGrp]) 
        if (!is.null(bulk)) bulk.all <- unique(colData(bulk)[, covisGrp])
        ft.all <- unique(unlist(lapply(seq_along(svgs), function(i) attribute(svgs[i])[[1]]$feature)))
        covis.trans <- covis_trans(bulk.all=bulk.all, cell.all=cell.all, ft.all=ft.all, tar.bulk=tar.cell.blk, tar.cell=tar.cell.blk, lis.match=lis.rematch, covis.type=covis.type, col.idp=col.idp)
        ft.trans.shm <- covis.trans$ft.trans.shm
        lis.rematch <- covis.trans$lis.match 
        if (col.idp) gene <- gene[, grep('^bulk$', se.scl.sel$bulkCell), drop=FALSE] 
      } else if ('auto' %in% method) {
        res <- scell.mod.lis$covis.auto$res
        cell.all <- setdiff(unique(subset(res, , bulkCell=='cell')$assignedBulk), 'none')
        bulk.all <- unique(subset(res, , bulkCell=='bulk')$sample)
        covis.trans <- covis_trans(bulk.all=bulk.all, cell.all=cell.all, ft.all=NULL, tar.bulk=tar.cell.blk, tar.cell=tar.cell.blk, lis.match=NULL, covis.type=covis.type, col.idp=col.idp)
        ft.trans.shm <- covis.trans$ft.trans.shm
        lis.rematch <- NULL
        if (col.idp) gene <- gene[, grep('^bulk$', res$bulkCell), drop=FALSE] 
      }
      }
      svg.na <- names(svgs[, 'svg'])
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      grob.all <- gg.all <- gcol.all <- lgd.all <- lgd.grob.all <- gcol.lgd.all <- gg.grob.lis <- NULL
      for (i in seq_along(svgs)) {
        if (is.null(raster.par$coal)) charcoal <- FALSE else if (raster.par$coal=='Yes') charcoal <- TRUE else if (raster.par$coal=='No') charcoal <- FALSE
        cat('All grob/ggplot of row selection:', ID, ' \n')
        incProgress(0.75, detail=paste0('preparing ', paste0(ID, collapse=';')))
        # if (!is.null(cna.match$cna)) { 
		#  if (ncol(gene)!=length(cna.match$cna)) return()
        #  colnames(gene) <- cna.match$cna 
        #}
        size.key <- as.numeric(lis.par$legend['key.size', 'default'])
        svg0 <- svgs[i]
        if (is.null(raster.par$over)) raster_pa(svg0)[[1]] <- NULL
        gg.lis <- gg_shm(svg.all=svg0, gene=gene, con.na=con.na$v, geneV=geneV()$bar.v, col.idp=col.idp, charcoal=charcoal, alpha.overlay=alp.over, ID=ID, cols=color$col, covis.type=covis.type, ft.trans=ft.trans, ft.trans.shm=ft.trans.shm, sub.title.size=input$title.size * scale.shm, legend.nrow=as.numeric(lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.size=input$line.size, line.color=input$line.color, lis.rematch = lis.rematch) # All gene IDs are used.
        msg <- paste0(svg.na[i], ': no common spatial features detected between data and aSVG!')
        if (is.null(gg.lis)) {
        showNotification(msg, duration=2, closeButton = TRUE)
        cat(msg, '\n'); 
        }
       validate(need(!is.null(gg.lis), msg))
       # Append suffix '_i' for the SHMs of ggplot under SVG[i], and store them in a list.
       ggs <- gg.lis$g.lis.all; names(ggs) <- paste0(names(ggs), '_', i)
       gg.all <- c(gg.all, ggs) 
       gcols <- gg.lis$gcol.lis.all; names(gcols) <- paste0(names(gcols), '_', i)
       gcol.all <- c(gcol.all, gcols)
       # Store legend/colours of ggplot in a list.
       lgd.all <- c(lgd.all, list(gg.lis$g.lgd))
       gcol.lgd.all <- c(gcol.lgd.all, list(gg.lis$gcol.lgd))
       # Same with ggs: append suffix '_i' for the SHMs of grob under SVG[i], and store them in a list.
       grob.lis <- grob_shm(ggs, cores=deter_core(1, svg.obj=svgs[, 'svg'][[i]]))
       grob.all <- c(grob.all, grob.lis)
       lgd.grob.lis <- grob_shm(lgd.all, cores=deter_core(1, svg.obj=svgs[, 'svg'][[i]]), lgd.pos='bottom')
       lgd.grob.all <- c(lgd.grob.all, lgd.grob.lis)
       # All ggplots/grobs are stored in nested lists under each SVG for use in relatice scale.
       lis0 <- list(grob.lis = grob.lis, gg.lis = ggs, lgd.lis = gg.lis$g.lgd, lgd.grob=lgd.grob.lis[[1]], gcol.lis=gcols, gcol.lgd=gg.lis$gcol.lgd)
       gg.grob.lis <- c(gg.grob.lis, list(lis0))
      }
     names(lgd.all) <- names(lgd.grob.all) <- names(gg.grob.lis) <- svg.na
     names(gcol.lgd.all) <- paste0('col_', svg.na)
     init$new <- 0 # Terminates gs.new.
     cat('Done! \n'); return(list(grob.all = grob.all, gg.all = gg.all, lgd.all = lgd.all, lgd.grob.all=lgd.grob.all, gcol.all=gcol.all, gcol.lgd.all=gcol.lgd.all, gg.grob.lis = gg.grob.lis))
     }) # withProgress
    }) # reactive
    shm$grob.all <- gs.all()$grob.all; shm$gg.all <- gs.all()$gg.all
    shm$lgd.all <- gs.all()$lgd.all; shm$lgd.grob.all <- gs.all()$lgd.grob.all
    shm$gcol.all <- gs.all()$gcol.all
    shm$gcol.lgd.all <- gs.all()$gcol.lgd.all
    shm$grob.gg.all <- gs.all()$gg.grob.lis
  }) # observeEvent
 
  # when 'color <- reactiveValues(col="none")', upon the app is launched, 'gs.new' is evaluated for 3 time. In the 1st time, 'gID$new'/'gID$all' are NULL, so 'gs.new' is NULL. In the 2nd time, 'color$col[1]=='none'' is TRUE, so NULL is returned to 'gs.new', but 'gID$new'/'gID$all' are 'HRE2'. In the third time, 'color$col[1]=='none'' is FALSE, so 'gs.new' is not NULL, but 'gID$new' is still 'HRE2', so it does not triger evaluation of 'observeEvent' and hence SHMs and legend plot are not returned upon being launched. The solution is to assign colors to 'color$col' in 'observe' upon being launched so that in the 2nd time 'gs.new' is not NULL, and no 3rd time.
  observeEvent(gs.new(), { 
    cat('Updating grobs/ggplots/legends based on new ID ... \n')
    if (is.null(svgs())|is.null(gID$new)|length(gID$new)==0|is.null(gID$all)|is.null(gs.new())) return(NULL)
    grob.gg.lis <- gs.new()
    # Update grobs.
    grobs <- grob.gg.lis[['grob.all']]
    grob.rm <- !names(shm$grob.all) %in% names(grobs)
    shm$grob.all <- c(shm$grob.all[grob.rm], grobs)
    # gs.new() becomes NULL at this step.
    # print(list(0, names(gs.new())))
    # Update ggs.
    ggs <- grob.gg.lis[['gg.all']]
    gg.rm <- !names(shm$gg.all) %in% names(ggs)
    shm$gg.all <- c(shm$gg.all[gg.rm], ggs) 
    # Update colours of ggs.
    gcols <- grob.gg.lis[['gcol.all']]
    gcol.rm <- !names(shm$gcol.all) %in% names(gcols)
    shm$gcol.all <- c(shm$gcol.all[gcol.rm], gcols)
    # Update legend/colours of SVGs.
    lgd0 <- grob.gg.lis[['lgd.all']] 
    shm$lgd.all <- c(shm$lgd.all, lgd0[!names(lgd0) %in% names(shm$lgd.all)])
    lgd.grob0 <- grob.gg.lis[['lgd.grob.all']] 
    shm$lgd.grob.all <- c(shm$lgd.grob.all, lgd.grob0[!names(lgd.grob0) %in% names(shm$lgd.grob.all)])
    gcol.lgd0 <- grob.gg.lis[['gcol.lgd.all']] 
    shm$gcol.lgd.all <- c(shm$gcol.lgd.all, gcol.lgd0[!names(gcol.lgd0) %in% names(shm$gcol.lgd.all)])
    # gs.new() becomes NULL before this step.
    # grob.gg.all <- gs.new()$grob.gg.all
    grob.gg.all <- grob.gg.lis$grob.gg.all
    if (is.null(shm$grob.gg.all)) shm$grob.gg.all <- grob.gg.all else {
      svg.na <- names(grob.gg.all)
      for (i in svg.na) {
        grobs <- grob.gg.all[[i]][['grob.all']]
        grob.rm <- !names(shm$grob.gg.all[[i]]$grob.all) %in% names(grobs)
        shm$grob.gg.all[[i]]$grob.all <- c(shm$grob.gg.all[[i]]$grob.all[grob.rm], grobs) 
        
        ggs <- grob.gg.all[[i]][['gg.all']]
        gg.rm <- !names(shm$grob.gg.all[[i]]$gg.all) %in% names(ggs)
        shm$grob.gg.all[[i]]$gg.all <- c(shm$grob.gg.all[[i]]$gg.all[gg.rm], ggs)
      }
    }; cat('Done! \n')
  })
  
  # Update subtitle size through theme().
  observeEvent(list(input$title.size, input$scale.shm, lis.url), {
    cat('Adjust title size ... \n')
    grob.gg.all <- shm$grob.gg.all; title.size <- input$title.size; scale.shm <- input$scale.shm
    svgs <- svgs()
    if (!is.list(grob.gg.all) | !is.numeric(title.size) | is.null(svgs) | is.null(lay.shm()) | !is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    gg.all <- grob.all <- NULL
    for (i in seq_along(grob.gg.all)) {
      gg.lis <- grob.gg.all[[i]]$gg.lis
      # Also update the central shm$grob.gg.all
      grob.gg.all[[i]]$gg.lis <- gg.lis <- lapply(gg.lis, function(x) { x + theme(plot.title = element_text(hjust = 0.5, size = title.size * scale.shm)) })
    gg.all <- c(gg.all, gg.lis)
    # Also update the central shm$grob.gg.all
    grob.gg.all[[i]]$grob.lis <- grob.lis <- grob_shm(gg.lis, cores = deter_core(2, svg.obj=svgs[, 'svg'][[i]]))
    grob.all <- c(grob.all, grob.lis) 
    }; shm$grob.all <- grob.all; shm$gg.all <- gg.all
    shm$grob.gg.all <- grob.gg.all
    cat('Done!\n')
  })

  observe({
    if (is.null(se.scl.sel())| length(ids$sel)==0 |is.null(svgs())|is.null(shm$grob.all)) return(NULL)
    col.n <- input$col.n; if (!check_obj(list(col.n))) return()
    lgc.nc <- col.n>=1 & as.integer(col.n)==col.n 
    if (!lgc.nc) {
      show_mod(lgc.nc, msg='No. of columns should be a positive integer!')
    }; req(lgc.nc)
  })
  # shm$lgd.all can update itself and lead to endless circles, thus it cannot be used to update the observeEvent below. In addition, when using bookmarked url, shm$lgd.all is first NULL (legend parameters are updating observeEvent below) then real ggplot object (parameters will not update oberverEvent again since they didn't change). Therefore, use lgd.par as an anchor. Only none of shm$lgd.all and legend parameters is NULL, will the observeEvent below be updated. 
  lgd.par <- reactiveValues(par=NULL)
  observe({
    # On the landing page, if the url is taken with all default parameters, after clicking the image/link, the app displays blank page, since input$lgd.key.size, input$lgd.row, input$lgd.label are all NULL, thereby not executing this step. Solution: change at least one of the parameters (e.g. horizontal layout) then take the url. 
    # print(list('adjust', is.null(shm$lgd.all), !is.numeric(input$lgd.key.size), input$lgd.row, input$lgd.label))
    if (is.null(shm$lgd.all)|!is.numeric(input$lgd.key.size)|!is.numeric(input$lgd.row)|is.null(input$lgd.label)) return()
    lgd.par$par <- list(lgd.key.size=input$lgd.key.size, lgd.row=input$lgd.row, lgd.label=input$lgd.label, lgd.lab.size=input$lgd.lab.size, lis.url=lis.url)
  })
  # lis.url is included in lgd.par$par, so it can trigger observeEvent when bookmarked url is used.
  observeEvent(lgd.par$par, {
    cat('Adjust legend size/rows/aspect ratio ... \n')
    lis.par <- lgd.par$par
    lgd.key.size <- lis.par$lgd.key.size; lgd.row <- lis.par$lgd.row
    lgd.label <- lis.par$lgd.label; label.size <- lis.par$lgd.lab.size
    gcol.lgd <- shm$gcol.lgd.all
    if (is.null(shm$lgd.all)|!is.numeric(lgd.key.size)|!is.numeric(lgd.row)|is.null(lgd.label)|is.null(gcol.lgd)) return()
    # Potential endless circles: shm$lgd.all updates itself.
    # gg.all=shm$lgd.all; size.key=lgd.key.size; size.text.key=NULL; row=lgd.row; position.text.key='right'; label=(lgd.label=='Yes'); label.size=label.size
    # save(gg.all, size.key, size.text.key, row, position.text.key, label, label.size, file='gg.lgd.all')
    withProgress(message="Adjusting the legend plot: ", value=0, {
     incProgress(0.25, detail="please wait ...")
     shm$lgd.all <- gg_lgd(gg.all=shm$lgd.all, gcol.lgd=gcol.lgd, size.key=lgd.key.size, size.text.key=NULL, row=lgd.row, position.text.key='right', label=(lgd.label=='Yes'), label.size=label.size); cat('Done! \n')
     incProgress(0.5, detail="please wait ...")
    })
  })

  observeEvent(list(shm$grob.all, input$genCon), {
    cat('Reordering grobs/ggplots ... \n') 
    if (is.null(gID$all)|is.null(shm$grob.all)|is.null(shm$gg.all)) return()
    na.all <- names(shm$grob.all); pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    # gg.all1 <- shm$gg.all1; save(gg.all1, pat.all, na.all, file='gg.all1')
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', na.all)); if (length(con)==0) return()
    na.all <- sort_gen_con(ID.sel=gID$all, na.all=na.all, con.all=con, by=input$genCon)
    # grob1/gg.all1 are used to add/remove 2nd legend.
    shm$grob.all1 <- shm$grob.all[na.all]; shm$gg.all1 <- shm$gg.all[na.all]
    cat('Done! \n')
  })

  # Add value legend to SHMs.
  # 'observeEvent' is able to avoid infinite cycles while 'observe' may cause such cycles. E.g. in the latter, 'is.null(shm$gg.all)' and 'shm$gg.all1 <- gg.all <- gg_2lgd()' would induce each other and form infinit circles.
  observeEvent(list(input$val.lgd, input$val.lgd.row, input$val.lgd.key, input$val.lgd.text, input$val.lgd.feat), {
    
    cat('Adding value legend... \n')
    validate(need(try(as.integer(input$val.lgd.row)==input$val.lgd.row & input$val.lgd.row>0), 'Legend key rows should be a positive integer!'))
    validate(need(try(input$val.lgd.key>0), 'Legend key size should be a positive numeric!'))
    validate(need(try(input$val.lgd.text>0), 'Legend text size should be a positive numeric!'))
    
    if.con <- is.null(shm$gg.all)|is.null(sam())|is.null(input$val.lgd)|is.null(input$val.lgd.feat)|input$val.lgd==0
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    gg.all <- shm$gg.all1
    if ((input$val.lgd %% 2)==1) {

      gg.all <- gg_2lgd(gg.all=gg.all, sam.dat=sam(), ft.trans=tis.trans$v, position.2nd='bottom', legend.nrow.2nd=input$val.lgd.row, legend.key.size.2nd=input$val.lgd.key, legend.text.size.2nd=input$val.lgd.text, add.feature.2nd=(input$val.lgd.feat=='Yes'))
      shm$gg.all1 <- gg.all <- lapply(gg.all, function(x) { x+theme(legend.position="bottom") } )
      png(tmp.file); shm$grob.all1 <- lapply(gg.all, ggplotGrob)
      dev.off(); if (file.exists(tmp.file)) do.call(file.remove, list(tmp.file))
    
    } else if ((input$val.lgd %% 2)==0) { 

      cat('Remove value legend... \n')
      shm$gg.all1 <- gg.all <- lapply(gg.all, function(x) { x+theme(legend.position="none") })
      png(tmp.file); shm$grob.all1 <- lapply(gg.all, ggplotGrob) 
      dev.off(); if (file.exists(tmp.file)) do.call(file.remove, list(tmp.file))

    }; cat('Done! \n')

  })
  # In "observe" and "observeEvent", if one code return (NULL), then all the following code stops. If one code changes, all the code renews.
    lay.shm <- reactive({
    cat('Spatial heatmaps layout ... \n')
    se.scl.sel <- se.scl.sel()
    if (!check_obj(list(se.scl.sel))) req('')
    if.con <- is.null(se.scl.sel)| length(ids$sel)==0 |is.null(svgs())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    if.con <- length(ids$sel)==0 |is.null(svgs())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    if (is.na(color$col[1])|length(color$col=="none")==0) return(NULL)
    grob.na <- names(shm$grob.all1)
    cell <- scell.mod.lis$sce.upl$cell
    profile <- ifelse(input$profile!='fixed', TRUE, FALSE)
    # Select target grobs.
    # Use definite patterns and avoid using '.*' as much as possible. Try to as specific as possible.
    pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    if (profile==TRUE) { 
      grob.lis.p <- shm$grob.all1[grepl(pat.all, grob.na)] 
      con <- unique(gsub(pat.all, '\\2\\3', names(grob.lis.p))); if (length(con)==0) return()
    } else if (profile==FALSE) { 
      grob.lis.p <- dim.shm.grob.all$val; con <- NULL
    }
    # grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    # Indexed cons with '_1', '_2', ... at the end.
    lay <- input$genCon; ID <- gID$geneSel; ncol <- input$col.n
    lay <- lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, lay.mat = TRUE, scell=ifelse(is.null(dim.shm.grob.all$val), FALSE, TRUE), profile=profile)
    # If 'cat' is the last step, NULL is returned.
    cat('Done! \n'); lay
  })

# shm <- shms_rela_size(input, svg.df, shm, lay.shm, svg.path)

  observeEvent(list(input$relaSize), {
    relaSize <- input$relaSize; svgs <- svgs(); grob.gg.all <- shm$grob.gg.all
    cat('Adjust relative plot size in multiple aSVGs ... \n')
    if (!is.numeric(relaSize) | !is(svgs, 'SVG') | length(svgs) == 1 | !is.list(grob.gg.all) | is.null(lay.shm()) | is.null(svg.path())) return()
    if (relaSize < 0) return()
    # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
    w.h.max <- max(unlist(svgs[, 'dimension']))
    gg.all <- grob.all <- NULL
    for (i in seq_along(grob.gg.all)) {
      gg.lis <- grob.gg.all[[i]]$gg.lis 
      # Also update the central shm$grob.gg.all
      grob.gg.all[[i]]$gg.lis <- gg.lis <- rela_size(dimension(svgs[i])[[1]]['height'], w.h.max, relaSize, nrow(lay.shm()), gg.lis)
      gg.all <- c(gg.all, gg.lis)
      # Also update the central shm$grob.gg.all
      grob.gg.all[[i]]$grob.lis <- grob.lis <- grob_shm(gg.lis, cores = deter_core(2, svg.obj=svgs[, 'svg'][[i]]))
      grob.all <- c(grob.all, grob.lis) 
    }; shm$grob.all <- grob.all; shm$gg.all <- gg.all; shm$grob.gg.all <- grob.gg.all
  })
   dim.shm.gg.all <- reactiveValues()
   dim.shm.grob.all <- reactiveValues()
   dim.lgd.lis <- reactiveValues()

   output$tarCellBlk <- renderUI({
     covis.type <- scell.mod.lis$sce.upl$covis.type
     if (is.null(covis.type)) return()
     if ('toCellAuto'%in% covis.type) {
       res <- scell.mod.lis$covis.auto$res
       if (is.null(res)) return()
       cho <- c('all', unique(subset(res, , bulkCell=='bulk')$sample))
       lab <- 'Target bulk'
     } else if ('toBulkAuto'%in% covis.type) {
       res <- scell.mod.lis$covis.auto$res
       if (is.null(res)) return()
       cho <- c('all', unique(subset(res, , bulkCell=='cell')$assignedBulk))
       lab <- 'Target cell'
       cho <- setdiff(cho, 'none')
     } else if ('toCell'%in% covis.type) {
       ft.rematch <- scell.mod.lis$covis.man$match.mod.lis$ft.reorder$ft.rematch
       cho <- c('all', names(ft.rematch)); lab <- 'Target bulk'
     } else if ('toBulk'%in% covis.type) {
       ft.rematch <- scell.mod.lis$covis.man$match.mod.lis$ft.reorder$ft.rematch
       cho <- c('all', names(ft.rematch)); lab <- 'Target cell' 
     }
     ns <- session$ns
     if (!grepl(na.sgl, ipt$fileIn)) return()
     selectInput(ns('tarCellBlk'), label=lab, choices=cho)
   })

   # Ensure the matching button is effective only if shm$grob.all1/shm$gg.all1 is not NULL. Otherwise dim plots are absent.
   #dim.shm.par <- reactiveValues(); observe({
   #  if (is.null(shm$grob.all1)|is.null(shm$gg.all1)) return()
   #  dim.shm.par$val <- scell.mod.lis$, match.lis$val$but.match$val
   # })
   # Single-cell reduced dimensionality for each gene__con. 
   observeEvent(list(shm$grob.all, shm$grob.all1, input$dims, input$dimLgdRows, input$profile), {
     cat('Manual: compiling PCA/TSNE/UMAP and SHMs ... \n')
     if (is.null(ipt$fileIn)) return()
     if (!grepl(na.sgl, ipt$fileIn)) return()
     covis.type <- scell.mod.lis$covis.man$covis.type
     if (all(!c('toCell', 'toBulk') %in% covis.type)) return()
     sce.dimred <- scell.mod.lis$covis.man$dimred
     ft.rematch <- scell.mod.lis$covis.man$match.mod.lis$ft.reorder$ft.rematch
     targ <- tar.bulk <- tar.cell <- input$tarCellBlk; profile <- input$profile
     covisGrp <- scell.mod.lis$sce.res()$covisGrp
     gg.all1 <- shm$gg.all1; grob.all1 <- shm$grob.all1
     lgd.all <- shm$lgd.all; lgd.grob.all <- shm$lgd.grob.all
     gcol.all <- shm$gcol.all; dims <- input$dims
     gcol.lgd.all <- shm$gcol.lgd.all
     con.na <- scell.mod.lis$sce.res()$sce.lis$con.na
     dimLgdRows <- input$dimLgdRows
     if (is.null(sce.dimred)|is.null(ft.rematch)|is.null(grob.all1)|is.null(gg.all1)|is.null(dims)|is.null(covisGrp)|is.null(targ)|is.null(profile)|is.null(dimLgdRows)) return()
     prof <- ifelse(profile!='fixed', TRUE, FALSE)
     if (dims=='TSNE') gg.dim <- plotTSNE(sce.dimred, colour_by=covisGrp)
     if (dims=='PCA') gg.dim <- plotPCA(sce.dimred, colour_by=covisGrp)
     if (dims=='UMAP') gg.dim <- plotUMAP(sce.dimred, colour_by=covisGrp)

     scale.shm <- input$scale.shm
     if (!is.numeric(scale.shm)) return()
     if (scale.shm <= 0) return()
     tit.size <- input$title.size * scale.shm
     if (covis.type %in% 'toBulk') {
       if (targ=='all') { 
         tar.cell <- names(ft.rematch)
         if (length(tar.cell)==0) return()
       }; tar.bulk <- NULL
     } else if (covis.type %in% 'toCell') {
       if (targ=='all') { 
         tar.bulk <- names(ft.rematch)
         if (length(tar.bulk)==0) return()
       }; tar.cell <- NULL
     }
    dim.shm.lis <- NULL
    withProgress(message="Embedding plot: ", value=0, {
     incProgress(0.25, detail="please wait ...")
     if (!'idp' %in% profile) { 
       if (covis.type %in% 'toBulk') {
         # source('~/spatialHeatmap/R/dim_color.R')
         dim.shm.lis <- dim_color(gg.dim=gg.dim, gg.shm.all=gg.all1, grob.shm.all=grob.all1, col.shm.all=gcol.all, cell.group=covisGrp, tar.cell=tar.cell, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=lgd.grob.all, profile=prof, con.na=con.na, lis.match=ft.rematch, sub.title.size=tit.size, dim.lgd.pos='bottom', dim.lgd.nrow=dimLgdRows)
         incProgress(0.25, detail="please wait ...")
        } else if (covis.type %in% 'toCell') {
          # Tar bulk is in SHM already.
          dim.shm.lis <- dim_color2cell(gg.dim=gg.dim, gg.shm.all=gg.all1, grob.shm.all=grob.all1, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=lgd.grob.all, profile=prof, cell.group=covisGrp, con.na=con.na, lis.match=ft.rematch, sub.title.size=tit.size, dim.lgd.pos='bottom', dim.lgd.nrow=2)
          incProgress(0.25, detail="please wait ...")
        } 
    } else {
       ID <- gID$geneSel; se.scl.sel <- se.scl.sel()
       req('bulkCell' %in% colnames(colData(se.scl.sel)))
       geneV <- geneV()$bar.v; cols <- color$col
       if ('none' %in% cols[1]) return()
       if (!check_obj(list(ID, se.scl.sel, geneV, cols))) return()
       vars.cell <- unique(sce.dimred$variable)
       # In data_server, in order to combine bulk and cell for joint scaling, reduced dims in cell are erased.
       gg.dim.lis <- NULL; for (i in vars.cell) {
         sce.dimred0 <- subset(sce.dimred, , variable==i)
         gg.dim <- plot_dim(sce.dimred0, dim=dims, color.by=covisGrp)
         gg.dim.lis <- c(gg.dim.lis, list(gg.dim))
      }; names(gg.dim.lis) <- vars.cell
      blk.n <- sum(se.scl.sel$bulkCell=='bulk')
      cell.n <- sum(se.scl.sel$bulkCell=='cell')
      assay <- assay(se.scl.sel)
      as.cell <- assay[, seq_len(cell.n) + blk.n, drop=FALSE]
      dim.shm.lis <- dim_color_idp(sce=sce.dimred, covis.type=covis.type, ID=ID, gene=as.cell, tar.cell=tar.cell, tar.bulk=tar.bulk, con.na.cell=con.na.cell$v, geneV=geneV, cols=cols, gg.dim=gg.dim.lis, gg.shm.all=gg.all1, grob.shm.all=grob.all1, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=lgd.grob.all, profile=prof, cell.group=covisGrp, lis.match=ft.rematch, sub.title.size=tit.size, dim.lgd.pos='bottom', dim.lgd.nrow=dimLgdRows, lgd.plot.margin=margin(t=0.01, r=0.2, b=0.01, l=0.2, unit="npc"))
      incProgress(0.25, detail="please wait ...")
      dim.lgd.lis$v <- dim.shm.lis$dim.lgd.lis
     }
    })
     dim.shm.grob.all$val <- dim.shm.lis$dim.shm.grob.lis
     dim.shm.gg.all$val <- dim.shm.lis$dim.shm.gg.lis
     cat('Done! \n')
     # save(gg.all1, file='gg.all1'); save(grob.all1, file='grob.all1'); save(gcol.all, file='gcol.all'); save(gg.dim, file='gg.dim'); save(clus, file='clus'); save(ft.rematch, file='ft.rematch')
     # lgd.lis <- shm$lgd.all; save(dim.shm.grob.lis, gg.all1, gg.dim.all, gcol.all, ft.rematch, lgd.lis, file='dgggl')

   })
   # Coclustering: single-cell reduced dimensionality for each gene__con.
   covisAuto <- reactiveValues()
   observe({ covisAuto$v <- scell.mod.lis$covis.auto })
   observeEvent(list(covisAuto$v$res, covisAuto$v$tailor.lis$v$df.sel.cell$val1, input$tarCellBlk, input$profile, shm$grob.all, shm$grob.all1, input$dims, input$dimLgdRows), ignoreNULL=FALSE, {
     cat('Coclustering: compiling PCA/TSNE/UMAP and SHMs ... \n')
     if (is.null(ipt$fileIn)) return()
     if (!grepl(na.sgl, ipt$fileIn)) return()
     covis.type <- scell.mod.lis$covis.auto$covis.type
     if (all(!c('toCellAuto', 'toBulkAuto') %in% covis.type)) return()
     targ <- input$tarCellBlk; profile <- input$profile
     res <- covisAuto$v$res
     # covisGrp <- scell.mod.lis$sce.res()$covisGrp
     # Only cell groups are need here.
     covisGrp <- 'assignedBulk'
     gg.all1 <- shm$gg.all1; grob.all1 <- shm$grob.all1
     gcol.all <- shm$gcol.all; dims <- input$dims
     lgd.all <- shm$lgd.all; lgd.grob.all <- shm$lgd.grob.all
     gcol.lgd.all <- shm$gcol.lgd.all
     dimLgdRows <- input$dimLgdRows
     # tar.cell interace is generated by renderUI, it will not be NULL until the relevant tab is clicked.
     if (is.null(targ)|is.null(grob.all1)|is.null(gg.all1)|is.null(lgd.all)|is.null(dims)|is.null(covisGrp)|is.null(profile)|is.null(res)|is.null(dimLgdRows)) return()
     prof <- !'fixed' %in% profile
     cell <- subset(res, , bulkCell=='cell') 
     if (dims=='TSNE') gg.dim <- plotTSNE(cell, colour_by=covisGrp)
     if (dims=='PCA') gg.dim <- plotPCA(cell, colour_by=covisGrp)
     if (dims=='UMAP') gg.dim <- plotUMAP(cell, colour_by=covisGrp)
     con.na <- scell.mod.lis$sce.res()$sce.lis$con.na
     scale.shm <- input$scale.shm
     if (!is.numeric(scale.shm)) return()
     if (scale.shm <= 0) return()
     tit.size <- input$title.size * scale.shm
     if ('all' %in% targ) {
       if ('toBulkAuto' %in% covis.type) targ <- setdiff(unique(cell$assignedBulk), 'none')
       if ('toCellAuto' %in% covis.type) targ <- unique(subset(res, , bulkCell=='bulk')$sample)
     }
     # source('~/spatialHeatmap/R/dim_color_coclus.R')

    withProgress(message="Embedding plot: ", value=0, {
     incProgress(0.25, detail="please wait ...")
     if (!'idp' %in% profile) {
       message('coloring by group or fixed colors ...')
       dim.shm.lis <- dim_color_coclus(sce=cell, targ=targ, profile=prof, gg.dim = gg.dim, gg.shm.all=gg.all1, grob.shm.all = grob.all1, gg.lgd.all=lgd.all, col.shm.all = gcol.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=lgd.grob.all, con.na=con.na, lis.match=NULL, sub.title.size=tit.size, dim.lgd.pos='bottom', dim.lgd.nrow=dimLgdRows) 
       incProgress(0.25, detail="please wait ...")
     } else {
      message('independent coloring ...')
      ID <- gID$geneSel; se.scl.sel <- se.scl.sel()
      req('bulkCell' %in% colnames(colData(se.scl.sel)))
      geneV <- geneV()$bar.v; cols <- color$col
      if ('none' %in% cols[1]) return()
      if (!check_obj(list(ID, se.scl.sel, geneV, cols))) return()
      gg.dim <- plot_dim(cell, dim=dims, color.by=covisGrp)
      gg.dim.lis <- list(con=gg.dim); assay <- assay(se.scl.sel)
      blk.n <- sum(se.scl.sel$bulkCell=='bulk')
      cell.n <- sum(se.scl.sel$bulkCell=='cell')
      as.cell <- assay[, seq_len(cell.n) + blk.n, drop=FALSE]
      dim.shm.lis <- dim_color_idp(sce=cell, covis.type=covis.type, targ=targ, ID=ID, gene=as.cell, con.na.cell=FALSE, geneV=geneV, cols=cols, gg.dim=gg.dim.lis, gg.shm.all=gg.all1, grob.shm.all=grob.all1, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=lgd.grob.all, profile=prof, cell.group=covisGrp, sub.title.size=tit.size, dim.lgd.pos='bottom', dim.lgd.nrow=dimLgdRows, lgd.plot.margin=margin(t=0.01, r=0.2, b=0.01, l=0.2, unit="npc"))
      incProgress(0.25, detail="please wait ...")
      dim.lgd.lis$v <- dim.shm.lis$dim.lgd.lis
      }
    })
      #grob.all <- dim.shm.lis$dim.shm.grob.lis
      #gg.all <- dim.shm.lis$dim.shm.gg.lis
      #dim.lgd.lis <- dim.shm.lis$dim.lgd.lis
      # save(dim.shm.lis, file='dim.shm.lis')
     dim.shm.grob.all$val <- dim.shm.lis$dim.shm.grob.lis
     dim.shm.gg.all$val <- dim.shm.lis$dim.shm.gg.lis
     #dim.shm.gg.lis <- dim.shm.lis$dim.shm.gg.lis
     #save(dim.shm.gg.lis, file='dim.shm.gg.lis')
     cat('Done! \n')
   }) 
   observeEvent(ipt$fileIn, { 
     dim.shm.gg.all$val <- NULL; dim.shm.grob.all$val <- NULL
     dim.lgd.lis$v <- NULL
   })
   observe({
    raster.ext <- paste0('\\', raster.ext, '$')
    svg.na <- svg.path1()$svg.na; if (!check_obj(svg.na)) return()
    if (!any(unlist(lapply(raster.ext, function(x) grepl(x, svg.na))))) hideTab(inputId="shmPar", target="raster") else showTab(inputId="shmPar", target="raster")
   })
   observe({
   # observeEvent(scell.mod.lis$sce.upl$covis.type, ignoreInit=FALSE, ignoreNULL=FALSE, { 
     covis.type <- scell.mod.lis$sce.upl$covis.type
     if (is.null(covis.type)|!grepl(na.sgl, ipt$fileIn)) { 
       hideTab(inputId="shmPar", target="scellTab") 
       updateTabsetPanel(session, inputId="shmPar", selected='basic')
     } else {
       showTab(inputId="shmPar", target="scellTab")
       updateTabsetPanel(session, inputId="shmPar", selected='scellTab')
     }
  })
   shmLay <- reactiveValues(val=NULL) 
  # Variables in 'observe' are accessible anywhere in the same 'observe'.
  observe({
    lay <- lay.shm(); scale.shm <- input$scale.shm
    dim.shm.grob.lis <- dim.shm.grob.all$val
    if (is.null(lay)|!is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    # subplot: height 300, width 250 
    # Avoid: if one column has all NAs in the layout matrix, the aspect ratio is distroyed. So only take the columns not containing all NAs.
    col.vld <- sum(unlist(lapply(seq_len(ncol(lay)), function(x) !all(is.na(lay[, x])))))
    # width/height relate to scrolling in box.
    shmLay$width <- width <- col.vld * 300 * scale.shm
    shmLay$height <- height <- nrow(lay) * 300 * scale.shm
    output$shm <- renderPlot(width = width, height = height, { 
      cat('Plotting spatial heatmaps ... \n')
      if.con <- length(ids$sel)==0 |is.null(svgs())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)

    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    if (is.na(color$col[1])|length(color$col=="none")==0) return(NULL)
    grob.na <- names(shm$grob.all1)
    # Select target grobs.
    # Use definite patterns and avoid using '.*' as much as possible. Try to as specific as possible.
    pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    grob.lis.p <- shm$grob.all1[grepl(pat.all, grob.na)] # grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', names(grob.lis.p))); if (length(con)==0) return()
    profile <- input$profile; if (is.null(profile)) return()
    profile <- ifelse(profile!='fixed', TRUE, FALSE)
    method <- scell.mod.lis$sce.upl$method
    if (!is.null(dim.shm.grob.lis)) { # Select dimred and SHMs.
      if (is.null(method)) return()
      if ('man' %in% method & profile==TRUE) { 
        pat.all <- paste0('^(dim_|)', pat.all(), '(_\\d+$)')
        grob.lis.p <- dim.shm.grob.lis[grepl(pat.all, names(dim.shm.grob.lis))]
      } else grob.lis.p <- dim.shm.grob.lis # Co-clustering.
    }
    lay <- input$genCon; ID <- gID$geneSel; ncol <- input$col.n
    # This step is plotting.
    shmLay$val <- shm.lay <- lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, scell=ifelse(is.null(dim.shm.grob.lis), FALSE, TRUE), profile=profile, shiny=TRUE)
    cat('Done! \n')
    })
  })
  observeEvent(list(input$ext, input$res, input$lgd.incld, input$lgd.size), {
  output$dldBut <- renderUI({ }) 
  })
  observeEvent(input$dld.but, ignoreInit=TRUE, {
    if (is.null(shmLay$val)) return()
    showNotification(HTML('Please wait till the <strong> "Download" </strong> button shows up!'), closeButton = TRUE)
    shm.arr <- shmLay$val$shm; shm.lay <- shmLay$val$lay
      cat('Downloading SHMs ... \n')
      validate(need(try(input$res>0), 'Resolution should be a positive numeric!'))
      png(paste0(tmp.dir, '/tmp.png'));
      cs.grob <- ggplotGrob(shm.bar()); dev.off()
      cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc"))
      # Legend size in downloaded SHM is reduced.
      lgd.lis <- shm$lgd.all; gcol.lgd <- shm$gcol.lgd.all
      lgd.lis <- gg_lgd(gg.all=lgd.lis, gcol.lgd=gcol.lgd, size.key=input$lgd.key.size*0.5, size.text.key=NULL, label.size=input$lgd.lab.size, row=input$lgd.row, position.text.key='right', label=(input$lgd.label=='Yes'))
      if (input$lgd.incld=='Yes') { 
        png(paste0(tmp.dir, '/tmp.png'));
        grob.lgd.lis <- lapply(lgd.lis, ggplotGrob); dev.off()
        lgd.tr <- lapply(grob.lgd.lis, grobTree)
    # In 'arrangeGrob', if numbers in 'layout_matrix' are more than items in 'grobs', there is no difference. The width/height of each subplot is decided by 'widths' and 'heights'.
      w.lgd <- (1-0.08)/(ncol(shm.lay)+1); shm.w <- 1-0.08-w.lgd
      # If legend.r = 0, legend plot size is a square.
      lgd.size <- input$lgd.size; validate(need(is.numeric(lgd.size), ''))
      lgd.arr <- arrangeGrob(grobs=lgd.tr, layout_matrix=matrix(seq_along(lgd.lis), ncol=1), widths=unit(0.99, "npc"), heights=unit(rep(w.lgd + (0.99 - w.lgd) * lgd.size, length(lgd.lis)), "npc"))
        png(paste0(tmp.dir, '/tmp.png')); shm1 <- grid.arrange(cs.arr, shm.arr, lgd.arr, ncol=3, widths=unit(c(0.08-0.005, shm.w, w.lgd), 'npc')); dev.off() } else { 
        png(paste0(tmp.dir, '/tmp.png')); shm1 <- grid.arrange(cs.arr, shm.arr, ncol=2, widths=unit(c(0.08-0.005, 1-0.08), 'npc')); dev.off() 
      }
      ggsave(paste0(tmp.dir, '/shm.', input$ext), plot=shm1, device=input$ext, width=shmLay$width/72, height=shmLay$height/72, dpi=input$res, unit='in', limitsize = FALSE); cat('Done! \n') 
  output$dldBut <- renderUI({
    ns <- session$ns    
    downloadButton(ns("dld.shm"), "Download", style = "margin-top: 24px;")
  })  
    })

  output$dld.shm <- downloadHandler(
    filename=function() { paste0('shm.', input$ext) },
    content=function(file) { file0 <- paste0(tmp.dir, '/shm.', input$ext); 
    cat("Downloading 'shm' from", tmp.dir, '...\n')
    file.copy(file0, file, overwrite=TRUE) }
  )

  observe({ 
    ipt$fileIn; se.scl.sel(); ipt$adj.modInpath; A(); input$p; input$cv1; input$cv2; ids$sel; tis.trans$v; input$genCon 
    lis.par <- cfg$lis.par; req(check_obj(lis.par))
    url.val <- url_val('shmAll-ext', lis.url)
    updateRadioButtons(session, inputId='ext', selected=ifelse(url.val!='null', url.val, lis.par$shm.img['file.type', 'default']))
  })

  observe({
   input$vdo.key.size; input$vdo.key.row; input$vdo.val.lgd; tis.trans$v; input$vdo.lab.size; input$vdo.res; input$vdo.itvl
   input$vdo.bar.width
  })

  output$lgd1 <- lgd2 <- renderPlot(width='auto', height="auto", { # auto: no need to scroll. 
    cat('Plotting legend plot ... \n')
    lgd.row <- input$lgd.row; lgd.key.size <- input$lgd.key.size
    validate(need(try(as.integer(lgd.row)==lgd.row & lgd.row>0), ''))
    validate(need(try(lgd.key.size>0 & lgd.key.size<1), 'Legend key size should be between 0 and 1!'))
    svg.path <- svg.path1()
    if (is.null(svg.path1())|is.null(shm$lgd.all)|(length(svg.path$svg.na)>1 & is.null(input$shms.in))) return(ggplot())
    if (grepl(na.sgl, ipt$fileIn) & 'fixed' %in% input$profile) return(ggplot())
      # Width and height in original SVG.
    if (length(svg.path$svg.na)>1) svg.na <- input$shms.in else svg.na <- 1
    g.lgd <- shm$lgd.all[[svg.na]]
    if (!grepl(na.sgl, ipt$fileIn) | !'idp' %in% input$profile) {
      message('Done! \n'); return(g.lgd)
    } else {
     dim.lgd <- dim.lgd.lis$v; se.scl <- se.scl()
     if (!check_obj(list(se.scl))|is.null(dim.lgd)) return()
     # In covis, "variable" is added by users or upstream in the app.
     vars.cell <- unique(se.scl$variable); lgd.lis <- shm$lgd.all
     for (i in vars.cell) {
       dim.lgd0 <- dim.lgd[[i]]$dim.lgd
       lgd.lis <- c(lgd.lis, setNames(list(dim.lgd0), paste0(i, '_dim.lgd')))
     }; na.lgd <- names(lgd.lis); len <- length(na.lgd)
     grob.lgd.lis <- lapply(lgd.lis, ggplotGrob)
     res <- grid.arrange(grobs=grob.lgd.lis, layout_matrix=matrix(seq_along(na.lgd), ncol=1), newpage=TRUE, widths=unit(c(0.99), 'npc'), heights=unit(rep(c(0.99/len), len), 'npc')); message('Done!'); res
    }
  })
  observe({
    glyBut <- input$glyBut
    if (is.null(glyBut)) output$lgd2 <- NULL else if (glyBut==0) output$lgd2 <- NULL else output$lgd2 <- lgd2
  }) 

  output$lgd.ui <- renderUI({ 
    ns <- session$ns
    lis.par <- cfg$lis.par; req(check_obj(lis.par))
    if (is.null(input$lgdTog)) return(NULL) 
    if (input$lgdTog %% 2 == 1) return(NULL)
    url.lgd.row <- url_val('shmAll-lgd.row', lis.url)
    url.lgd.key.size <- url_val('shmAll-lgd.key.size', lis.url)
    url.lgd.label <- url_val('shmAll-lgd.label', lis.url)
    url.lgd.lab.size <- url_val('shmAll-lgd.lab.size', lis.url)
    box(title="Legend Plot", status="primary", solidHeader=TRUE, collapsible=TRUE, width = 3, 
    navbarPage('Settings:',
    tabPanel("Basic",
    splitLayout(cellWidths=c("32%", "1%", '32%', '1%', '35%'),
    numericInput(inputId=ns('lgd.row'), label='Key rows', value=ifelse(url.lgd.row!='null', url.lgd.row, as.numeric(lis.par$legend['key.row', 'default'])), min=1, max=Inf, step=1, width=150), '',
    numericInput(inputId=ns('lgd.key.size'), label='Key size', value=ifelse(url.lgd.key.size!='null', url.lgd.key.size, as.numeric(lis.par$legend['key.size', 'default'])), min=0, max=1, step=0.02, width=150), ''
    # numericInput(inputId=ns('lgd.ratio1'), label='Aspect ratio', value=as.numeric(lis.par$legend['aspect.ratio', 'default']), min=0.01, max=Inf, step=0.01, width=150)
    )), # tabPanel

    tabPanel("Feature labels",
    splitLayout(cellWidths=c("30%", "1%", '30%'),
    radioButtons(inputId=ns("lgd.label"), label="Feature labels", choices=c("Yes", "No"), selected=ifelse(url.lgd.label!='null', url.lgd.label, lis.par$legend['label', 'default']), inline=TRUE), '',
    numericInput(inputId=ns('lgd.lab.size'), label='Label size', value=ifelse(url.lgd.lab.size!='null', url.lgd.lab.size, as.numeric(lis.par$legend['label.size', 'default'])), min=0, max=Inf, step=0.5, width=150)
    )) # tabPanel
    ), # navbarPage
    uiOutput(ns('lgds.sel')), splitLayout(cellWidths=c("99%", "1%"), plotOutput(ns("lgd")), "")) # box

  })

  observeEvent(list(ipt$fileIn, log(), tis.trans$v, input$col.but, input$sig.but, input$ckeyV, input$preScale), { ggly_rm(); vdo_rm() })

  gly.par <- reactiveValues()
  observeEvent(input$glyBut, {
    if (is.null(input$glyBut)) return()
    gly.par$val <- list(input$aspr, input$scale.ly, shm$gg.all1, shm.bar(), dim.lgd.lis$v, dim.shm.gg.all$val)
  })
  # eventReactive avoids endless circles.
  res.gly <- eventReactive(list(gly.par$val), {
    aspr <- input$aspr; scl <- input$scale.ly
    pat.all <- pat.all(); gg.all <- shm$gg.all1
    cs.g <- shm.bar(); na <- names(gg.all); fileIn <- ipt$fileIn
    if (!check_obj(list(aspr, scl, pat.all, gg.all, cs.g, fileIn))) return()
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all, '_\\d+$'), na)]
    if (length(na) == 0) return(); gg.all <- gg.all[na]

    if (grepl(na.sgl, fileIn)) {
      prof <- input$profile; dim.shm <- dim.shm.gg.all$val
      if (!check_obj(list(prof, dim.shm))) return()
      dim.gg <- dim.shm[paste0('dim_', na)]
    } else dim.gg <- NULL

    out.dir <- 'www'; html.dir <- file.path(out.dir, 'html_shm')
    if (dir.exists(html.dir)) unlink(html.dir, recursive=TRUE)
    withProgress(message="Animation: ", value=0, {
      incProgress(0.25, detail="preparing frames ...") 
      lis <- html_ly(gg.all=c(dim.gg, gg.all), cs.g=cs.g, aspr=aspr, anm.scale=scl, selfcontained=FALSE, out.dir=out.dir)
      if (!is(lis, 'list')) return()
      incProgress(0.25, detail="preparing frames ...") 
    }); return(c(lis, list(na=na)))
  }); observe({ res.gly() })

  output$ggly <- renderUI({
    message('Animation: accessing HTML files ...')
    lis <- res.gly(); if (!is(lis, 'list')) return()
    fm <- input$fm; if (is.null(fm)) return()
    na <- lis$na; na.sel <- na[as.integer(fm)]
    na.sel <- list.files('www/html_shm', pattern=na.sel)
    if (!check_obj(list(na.sel))) return(); 
    message(na.sel, ' Done!')
    tags$iframe(src=file.path('html_shm', na.sel), height = lis$height, width=lis$width, scrolling='yes') 
  })
  
  output$sld.fm <- renderUI({
    if (input$glyBut==0) return(); ns <- NS(id) 
    if (is.null(shm$gg.all)|is.null(pat.all())|is.null(gID$geneSel)) return(NULL) 
    gen.con.pat <- paste0('^', pat.all(), '_\\d+$') 
    sliderInput(inputId=ns('fm'), 'Frames', min=1, max=sum(grepl(gen.con.pat, names(shm$gg.all1))), step=1, value=1, animate=animationOptions(interval=input$t*10^3, loop=FALSE, playButton=icon('play'), pauseButton=icon('pause'))) 
  })

  anm.dld <- reactive({
    scale.ly <- input$scale.ly; gly.url <- gly.url()
    if (input$glyBut==0|is.null(gly.url)) return()
    if (is.null(svgs())|is.null(se.scl.sel())| length(ids$sel)==0 |color$col[1]=='none') return(NULL) 
    withProgress(message="Downloading animation: ", value=0, {
    incProgress(0.1, detail="in progress ...")
    gg.all <- shm$gg.all1; na <- names(gg.all)
    gg.na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]
    gg <- gg.all[gg.na]
    html_ly(gg.all=gg, cs.g=shm.bar(), anm.scale=1, out.dir=out.dir)
    incProgress(0.5, detail=paste0('preparing HTML files ...'))
   })
  })

  # This step leaves 'fil.na' in 'output$dld.anm' being a global variable.
  output$dld.anm <- downloadHandler( 
    filename=function(){ "html_shm.zip" },
    fil.na <- file.path(tmp.dir, 'html_shm.zip'),
    content=function(fil.na){ cat('Downloading animation... \n'); zip(fil.na, 'www/html_shm/'); message('Done!') }
  )

  observe({
    lis.par <- cfg$lis.par; req(check_obj(lis.par))
    updateSelectInput(session, "vdo.dim", selected=lis.par$shm.video['dimension', 'default'])
  })
  observe({
   fileIn  <- ipt$fileIn; if(!check_obj(list(fileIn))) return()
   if (!grepl(na.sgl, fileIn)) { 
     shinyjs::hide(id = "lgdVdo"); shinyjs::show(id = "shmVdo")
   } else { 
     shinyjs::show(id = "lgdVdo"); shinyjs::hide(id = "shmVdo")
   }
  })

  vdo.par <- reactiveValues(); observeEvent(input$vdo.but, {
    if (is.null(input$vdo.but)) return()
    ffm <- check_exp(test_ffm())
    idx <- 'w' %in% ffm | 'e' %in% ffm
    if (idx) {
      msg <- "'ffmpeg' is required to make videos!"
      show_mod(!idx, msg=msg); return()
    }
    vdo.par$val <- list(input$vdo.key.row, input$vdo.key.size, input$vdo.val.lgd, input$vdo.label, input$vdo.lab.size, input$vdoLgdDimRow, input$vdoLgdDimText, input$vdoLgdDimkey, input$vdoLgdKeyRow, input$vdoLgdText, input$vdoLgdkey, input$vdoH, input$vdo.bar.width, input$lgdR, input$vdo.dim, input$vdo.itvl, input$vdo.res, shm$gg.all1, dim.lgd.lis$v, dim.shm.gg.all$val)
  })
  observeEvent(list(vdo.par$val), {
    message('Making video ... \n')
    vdo.key.row <- input$vdo.key.row
    vdo.key.size <- input$vdo.key.size
    vdo.val.lgd <- input$vdo.val.lgd
    vdoText2 <- input$vdoText2; prof <- input$profile
    vdo.label <- input$vdo.label
    vdo.lab.size <- input$vdo.lab.size
    vdoLgdDimRow <- input$vdoLgdDimRow
    vdoLgdDimText <- input$vdoLgdDimText;
    vdoLgdDimkey <- input$vdoLgdDimkey
    vdoLgdKeyRow <- input$vdoLgdKeyRow
    vdoLgdText <- input$vdoLgdText
    vdoLgdkey <- input$vdoLgdkey; vdo.res <- input$vdo.res
    vdoH <- input$vdoH; lgdR <- input$lgdR
    vdo.bar.width <- input$vdo.bar.width
    vdo.dim <- input$vdo.dim; vdo.itvl <- input$vdo.itvl
    pat.all <- pat.all(); svgs <- svgs()
    se.scl.sel <- se.scl.sel(); gg.all <- shm$gg.all1
    shm.bar <- shm.bar(); fileIn <- ipt$fileIn
 
    if (!check_obj(list(vdo.key.row, vdo.key.size, vdo.val.lgd, vdoText2, vdo.label, vdo.lab.size, vdo.res, vdoH, vdo.bar.width, vdo.dim, vdo.itvl, pat.all, svgs, se.scl.sel, gg.all, shm.bar, fileIn))) return()
    if (length(ids$sel)==0 |color$col[1]=='none') return()
    idx <- vdo.res>=1 & vdo.res<=700 
    if (!idx) {
      msg <- "Resolution should be between 1 and 700!"
      show_mod(idx, msg=msg); return()
    }
    withProgress(message="Video: ", value=0, {
    incProgress(0.75, detail="in progress ...")
    na <- names(gg.all); pat <- paste0('^', pat.all, '_\\d+$')
    na <- na[grepl(pat, na)]; gg.sel <- gg.all[na]
    res <- vdo.res; dim <- input$vdo.dim
    if (vdo.dim %in% c('1280x800', '1280x1024', '1280x720') & vdo.res > 450) vdo.res <- 450
    if (vdo.dim=='1920x1080' & vdo.res > 300) vdo.res <- 300
  
    dim.lgd <- dim.gg <- NULL
    if (!grepl(na.sgl, fileIn)) type <- 'shm' else {
      dim.lgd <- dim.lgd.lis$v; se.scl <- se.scl()
      dim.shm <- dim.shm.gg.all$val; lgd.lis <- shm$lgd.all
      gcol.lgd <- shm$gcol.lgd.all
      if (!check_obj(list(prof, dim.lgd, se.scl, dim.shm, lgd.lis, gcol.lgd, lgdR, vdoLgdDimRow, vdoLgdDimText, vdoLgdDimkey, vdoLgdKeyRow, vdoLgdText, vdoLgdkey))) return()
      if ('idp' %in% prof) type <- 'col.idp' else type <- 'col.grp'
      dim.gg <- dim.shm[paste0('dim_', na)]
      vars.cell <- unique(se.scl$variable) 
      for (i in vars.cell) { # Combine SHM and dim legends.
        dim.lgd0 <- dim.lgd[[i]]$dim.lgd
        lgd.lis <- c(lgd.lis, setNames(list(dim.lgd0), paste0(i, '_dim.lgd')))
      }
    }
    vdo <- video(gg=c(dim.gg, gg.sel), cs.g=shm.bar, lgd=lgd.lis, gcol.lgd=gcol.lgd, lgd.r=lgdR, lgd.title='Legend', h=vdoH, type=type, sub.title.size=7, bar.width=vdo.bar.width, bar.value.size=4, lgd.key.size=vdoLgdkey, lgd.text.size=vdoLgdText, lgd.key.size.2nd=vdo.key.size, lgd.text.size.2nd=vdoText2, lgd.row=vdoLgdKeyRow, lgd.row.2nd=vdo.key.row, legend.value.vdo=('Yes' %in% vdo.val.lgd), label=('Yes' %in% vdo.label), label.size=vdo.lab.size, dim.lgd.text.size=vdoLgdDimText, dim.lgd.key.size=vdoLgdDimkey, dim.lgd.nrow=vdoLgdDimRow, video.dim=vdo.dim, res=vdo.res, interval=vdo.itvl, out.dir='www/video')
    if (is.null(vdo)) return()
    cat('Presenting video ... \n')
    incProgress(0.95, detail="Presenting video ...")
    w.h <- as.numeric(strsplit(vdo.dim, 'x')[[1]])
    output$video <-renderUI({ tags$video(id="video", type="video/mp4", src="video/shm.mp4", width=w.h[1], height=w.h[2], controls="controls") }); cat('Done! \n')
    })
  })
    scroll.h <- reactiveValues()
    observe({ h <- input$scrollH; scroll.h$h <- ifelse(is.null(h), 450, h) })
 output$shm.ui <- renderUI({
    ns <- session$ns; if (is.null(input$togSld)) return()
    lis.par <- cfg$lis.par; req(check_obj(lis.par))
    url.lgd.row <- url_val('shmAll-lgd.row', lis.url)
    url.lgd.key.size <- url_val('shmAll-lgd.key.size', lis.url)
    url.lgd.label <- url_val('shmAll-lgd.label', lis.url)
    url.lgd.lab.size <- url_val('shmAll-lgd.lab.size', lis.url)
    column(12, 
    fluidRow(splitLayout(id='barSHM', cellWidths=c("10px", "70px", paste0(input$togSld*92, '%'), paste0((1-input$togSld)*92, '%')), "",  
    plotOutput(ns("bar1")),
    if (input$togSld!=0) div(id='divSHM', style=paste0('overflow-y:scroll;height:', scroll.h$h, 'px;overflow-x:scroll'), plotOutput(ns("shm"), height='100%', width='100%')),

    if (input$togSld!=1) navbarPage('',
    tabPanel('Legend', list(uiOutput(ns('lgds.sel')), plotOutput(ns("lgd1")))),
    tabPanel("Settings",
    div(id=ns('setLgd'), 
    splitLayout(cellWidths=c("32%", "1%", '32%', '1%', '35%'),
    numericInput(inputId=ns('lgd.row'), label='Key rows', value=ifelse(url.lgd.row!='null', url.lgd.row, as.numeric(lis.par$legend['key.row', 'default'])), min=1, max=Inf, step=1, width=150), '',
    numericInput(inputId=ns('lgd.key.size'), label='Key size', value=ifelse(url.lgd.key.size!='null', url.lgd.key.size, as.numeric(lis.par$legend['key.size', 'default'])), min=0, max=1, step=0.02, width=150), ''
    ),
    splitLayout(cellWidths=c("30%", "1%", '30%'),
    radioButtons(inputId=ns("lgd.label"), label="Feature label", choices=c("Yes", "No"), selected=ifelse(url.lgd.label!='null', url.lgd.label, lis.par$legend['label', 'default']), inline=TRUE), '',
    numericInput(inputId=ns('lgd.lab.size'), label='Label size', value=ifelse(url.lgd.lab.size!='null', url.lgd.lab.size, as.numeric(lis.par$legend['label.size', 'default'])), min=0, max=Inf, step=0.5, width=150)
    )),
    bsTooltip(ns('setLgd'), title="Adjust legend keys and text in the spatial heatmap legend plot. <br/> Feature label: label spatial features with text or not.", placement = "top", trigger = "hover")
    ) # tabPanel
    ) # navbarPage
  )) # splitLayout(cellWidths
  ) # column
  })

 # addTooltip(session=session, id=ns('setLgd'), title="test", placement = "top", trigger = "hover", options = NULL)  
  
 output$lgds.sel <- renderUI({
    ns <- session$ns 
    if (is.null(svg.path1())) return(NULL)
    if (length(svg.path1()$svg.na)==1) return(NULL)
    svg.na <- svg.path1()[['svg.na']]
    svg.na <- svg.na[grepl('\\.svg$', svg.na)]
    selectInput(ns('shms.in'), label='Select plots', choices=svg.na, selected=svg.na[1])
  })
  observe({
    shmMhNet <- input$shmMhNet; interNav <- input$interNav
    if (is.null(shmMhNet)|is.null(interNav)) return()
    tab.inter <- ifelse(shmMhNet=='interTab' & interNav=='interPlot', 'yes', 'no')
    if (input$glyBut==0 & tab.inter=='yes') showModal(modal(msg=HTML(run.msg), easyClose=TRUE))
  })
  observe({
    shmMhNet <- input$shmMhNet; vdoNav <- input$vdoNav
    if (is.null(shmMhNet)|is.null(vdoNav)) return()
    tab.vdo <- ifelse(shmMhNet=='vdoTab' & vdoNav=='video', 'yes', 'no')
    if (input$vdo.but==0 & tab.vdo=='yes') showModal(modal(msg=HTML(run.msg), easyClose=TRUE))
  })
  # analysis_server('net', upl.mod.lis, dat.mod.lis, shm.mod.lis=list(gID=gID, tab.act.lis=tab.act.lis), sch.mod.lis)
  output$helpStatic <- renderUI({ 
    tags$iframe(seamless="seamless", src= "html/shm_shiny_manual.html#21_Static_image", width='100%', height='100%')
  }) 
  output$helpInter <- renderUI({ 
    tags$iframe(seamless="seamless", src= "html/shm_shiny_manual.html#22_Interactive_image", width='100%', height='100%')
  }) 
  output$helpVdo <- renderUI({ 
    tags$iframe(seamless="seamless", src= "html/shm_shiny_manual.html#23_Video", width='100%', height='100%')
  })
  observeEvent(list(prt$input$btnInf), {
    btnInf <- prt$input$btnInf
    if (!check_obj(btnInf)) return()
    if (btnInf > 0) updateTabsetPanel(session, inputId="shmMhNet", selected='shm1')
  }) 

  observe({
    ipt$fileIn; ipt$geneInpath; lis.par <- cfg$lis.par
    lis.par <- cfg$lis.par; req(check_obj(lis.par))
    url.val <- url_val('shmAll-ckeyV', lis.url)
    updateRadioButtons(session, inputId='ckeyV', selected=ifelse(url.val!='null', url.val, lis.par$shm.img['color.scale', 'default']))
    updateSliderInput(session, inputId='col.n', value=url_val('shmAll-col.n', lis.url, def=as.numeric(lis.par$shm.img['columns','default'])))
    updateRadioButtons(session, inputId="genCon", selected = url_val('shmAll-genCon', lis.url, def=lis.par$shm.img['layout.by', 'default']))
  # addPopover(session, "genCon", title="Data column: by the column order in data matrix.", placement="bottom", trigger='hover')
    url.val <- url_val('shmAll-scale.shm', lis.url)
  updateSliderInput(session, inputId='scale.shm', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['scale.plots', 'default'])))
    url.val <- url_val('shmAll-title.size', lis.url)
  updateSliderInput(session, inputId='title.size', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['title.size', 'default'])))
    url.val <- url_val('shmAll-color', lis.url)
    col.def <- ifelse(url.val!='null', url.val, lis.par$shm.img['color', 'default']) 
    updateSelectInput(session, 'colorOp', choices=unique(c(col.def, 'yellow,orange,red', 'green,yellow,orange', 'custom')), selected=col.def)
  updateTextInput(session, "color", value=col.def, placeholder=paste0('Eg: ', lis.par$shm.img['color', 'default']))
  url.val <- url_val('shmAll-val.lgd.row', lis.url)
  updateNumericInput(session, inputId='val.lgd.row', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['value.legend.rows', 'default'])))
  url.val <- url_val('shmAll-val.lgd.key', lis.url)
  updateNumericInput(session, inputId='val.lgd.key', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['value.legend.key', 'default'])))
  url.val <- url_val('shmAll-val.lgd.text', lis.url)
  updateNumericInput(session, inputId='val.lgd.text', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['value.legend.text', 'default'])))
  url.val <- url_val('shmAll-val.lgd.feat', lis.url)
  updateRadioButtons(session, inputId='val.lgd.feat', selected=ifelse(url.val!='null', url.val, lis.par$shm.img['include.feature', 'default']), inline=TRUE)
  url.val <- url_val('shmAll-line.color', lis.url)
  updateSelectInput(session, 'line.color', selected=ifelse(url.val!='null', url.val, lis.par$shm.img['line.color', 'default']))
  url.val <- url_val('shmAll-line.size', lis.url)
  updateNumericInput(session, inputId='line.size', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['line.size', 'default']))) 
  url.val <- url_val('shmAll-ext', lis.url)
  updateRadioButtons(session, inputId='ext', selected=ifelse(url.val!='null', url.val, lis.par$shm.img['file.type', 'default']))
  url.val <- url_val('shmAll-res', lis.url)
  updateNumericInput(session, inputId='res', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['dpi', 'default'])))
  url.val <- url_val('shmAll-lgd.incld', lis.url)
  updateRadioButtons(session, inputId='lgd.incld', choices=c('Yes', 'No'), selected=ifelse(url.val!='null', url.val, lis.par$shm.img['include.legend.plot', 'default']))
  url.val <- url_val('shmAll-lgd.size', lis.url) 
  updateNumericInput(session, inputId='lgd.size', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['legend.plot.size', 'default'])))
  url.val <- url_val('shmAll-relaSize', lis.url)
  updateNumericInput(session, inputId='relaSize', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.img['relative.size', 'default'])))
  url.val <- url_val('shmAll-vdo.key.row', lis.url)
  updateNumericInput(session, inputId='vdo.key.row', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.video['key.rows', 'default'])))
  url.val <- url_val('shmAll-vdo.key.size', lis.url)
  updateNumericInput(session, inputId='vdo.key.size', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.video['key.size', 'default'])))
  url.val <- url_val('shmAll-vdo.val.lgd', lis.url)
  updateSelectInput(session, inputId="vdo.val.lgd", selected=ifelse(url.val!='null', url.val, lis.par$shm.video['value.legend', 'default']))
  url.val <- url_val('shmAll-vdo.label', lis.url)
  updateSelectInput(session, inputId="vdo.label", selected=ifelse(url.val!='null', url.val, lis.par$shm.video['feature.label', 'default']))
  url.val <- url_val('shmAll-vdo.lab.size', lis.url)
  updateNumericInput(session, inputId='vdo.lab.size', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.video['label.size', 'default'])))
  url.val <- url_val('shmAll-vdo.bar.width', lis.url)
  updateNumericInput(session, inputId='vdo.bar.width', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.video['bar.width.video', 'default'])))
  url.val <- url_val('shmAll-vdo.itvl', lis.url)
  updateNumericInput(session, inputId='vdo.itvl', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.video['transition', 'default'])))
  url.val <- url_val('shmAll-vdo.res', lis.url)
  updateNumericInput(session, inputId='vdo.res', value=ifelse(url.val!='null', url.val, as.numeric(lis.par$shm.video['dpi', 'default'])))

  })
  onBookmark(function(state) { state })
  return(list(gID=gID, sam=sam, svgs=svgs, shmLay=shmLay, ipt=input))
})} # shm_server
