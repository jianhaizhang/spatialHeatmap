
options(stringsAsFactors=FALSE) 

# Every variable in every container should be checked at the beginning. E.g. input$fileIn in reactive({}). These checks will avoid amost all errors/warnings.

# enableWGCNAThreads()
# enableBookmarking("url")
server <- function(input, output, session) {
  tmp.dir <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
  lis.url <- reactiveValues(par=NULL)
  observe({
    # The url on browser is captured only if the url is refreshed or the "Enter" key is pressed, which applies in the cases shiny app is first launched or users modified parameters in the url. Otherwise the query is null.
    # lis.url$par is an empty list not NULL before refreshing/bookmarking.
    hos.port <- session$clientData
    lis.url$par <- parseQueryString(hos.port$url_search)
    hos.port1 <- paste0(hos.port$url_hostname, ':', hos.port$url_port, hos.port$url_pathname)
    lis.url$hos.port <- hos.port1
    lis.url$url_search <- hos.port$url_search
    # cat('Parameters in URL:', (names(lis.url$par)), '\n')
   })

    #observeEvent(input$show, {
    observe({
    lis.cfg <- yaml.load_file('config/config.yaml')
    lis.cfg <- lis.cfg[!vapply(lis.cfg, is.null, logical(1))]
    # Separate data sets, download files, and parameters.
    lis.dat <- lis.cfg[grepl('^dataset\\d+', names(lis.cfg))]
    lis.dld <- lis.cfg[grepl('download_single|download_multiple|download_spatial_temporal|download_batched_data_aSVGs', names(lis.cfg))]
    lis.par <- lis.cfg[!grepl('^dataset\\d+|download_single|download_multiple|download_spatial_temporal|download_batched_data_aSVGs', names(lis.cfg))]

    na.ipt <- NULL; na.ipt <- unlist(lapply(lis.dat, function(x) na.ipt <- c(na.ipt, x$name)))
    names(na.ipt) <- NULL
    if (FALSE) showModal(modalDialog(
        selectInput("fileIn", 'Select a start data set', na.ipt, lis.par$default.dataset),
        footer=tagList(actionButton("start.but", label='Confirm'))
    ))
    
    })
  observe({
    # validate(need(input$start.but>0, ''))
    withProgress(message="Loading dependencies: ", value=0, {
    incProgress(0.3, detail="in progress ...")
    library(spatialHeatmap); library(SummarizedExperiment); library(shiny); library(shinydashboard); library(shinydashboardPlus); library(grImport); library(rsvg); library(ggplot2); library(DT) 
    incProgress(0.6, detail="in progress ...")
    library(gridExtra); library(ggdendro); library(WGCNA); library(grid); library(xml2); library(plotly); library(data.table); library(genefilter); library(flashClust); library(visNetwork); 
    showModal(modal(title = HTML('<center><b>Welcome to spatialHeatmap!</b><center>'), msg = strong('Please wait for the landing page to finish!')))
    incProgress(0.9, detail="in progress...")
   library(reshape2); library(igraph); library(animation); library(av); library(shinyWidgets); library(yaml); library(HDF5Array); library(sortable); library(shinyBS); library(shinyjs); library(htmltools)
    # DEG
    library(gplots); library(UpSetR)
  })
  })

  url.id <- reactiveValues()
  mods <- reactiveValues(upload=NULL, search=NULL, data=NULL, shm=NULL)

search_server <- function(id, ids, cfg, lis.url, url.id, dat.mod.lis, session) {
  moduleServer(id, function(input, output, session) {
    cat('ID search box ... \n')
    observeEvent(session$clientData$url_search, {
      id.sch.sgl <- url_val('sear-sch.sgl', lis.url)
      # lis <- reactiveValuesToList(lis.url)
      # In bookmarked url of single search box, if there is a comma in gene ID and annotation combination, the combined ID and annotation will be broken at the comma, resulting two strings.
      if (id.sch.sgl[1]!='null') {
        # Multiple ID-discs are concatenated and separated by comma. E.g [ID1 disc1,ID2 disc2].
        id.sch.sgl <- unlist(strsplit(gsub('^\\[|\\]$', '', id.sch.sgl), ','))
      }
      url.id$sch.sgl <- lis.url$par[['sear-sch.sgl']] <- id.sch.sgl
      url.id$sch.mul <- lis.url$par[['sear-sch.mul']] <- unlist(strsplit(gsub(' |,', '__', url_val('sear-sch.mul', lis.url)), '__'))
      cat('Parsing IDs from URL ... done! \n')
    })
    geneIn <- reactiveValues(lis=NULL)
    observe({
      if (is.null(dat.mod.lis)) return()
      gene.all <- dat.mod.lis$geneIn; geneIn$lis <- gene.all()
    })
    pre.id <- reactiveValues(id=NULL)
    observe({ # Pre-selected ids in config file.
      if (is.null(geneIn$lis)) return()
      # geneIn <- dat.mod.lis$geneIn
      df.aggr.tran <- geneIn$lis$df.aggr.tran 
      id <- cfg$lis.par$data.matrix['selected.id', 'default']
      id <- make.names(strsplit(id, ',')[[1]])
      rna <- rownames(df.aggr.tran)
      id <- id[id %in% rownames(df.aggr.tran)]
      if (length(id)==0) id <- rna[1]; pre.id$id <- id
    })
    sch.sgl <- reactiveValues(cho=NULL, sel=NULL)
    observe({ # Accounts for IDs in URL for single search box.
      if (is.null(geneIn$lis)) return()
      # geneIn <- dat.mod.lis$geneIn; gen.lis <- geneIn()
      df.aggr.tran <- geneIn$lis$df.aggr.tran
      rna <- rownames(df.aggr.tran); df.met <- geneIn$lis$df.met
      cho <- paste0(rna, ' ', ifelse(rep('metadata' %in% colnames(df.met), length(rna)), df.met[, 'metadata'], ''))
      id <- pre.id$id
      if (!is.null(url.id$sch.sgl)) if (url.id$sch.sgl[1]!='null') {
        # Only obtain the ID, no annotation.
        sel <- sub(' .*', '', url.id$sch.sgl)
        id <- sel[sel %in% rna]
      }
      sel <- paste0(id, ' ', df.met[id, 'metadata'])
      sch.sgl$cho <- cho; sch.sgl$sel <- sel
      # updateSelectizeInput(session, 'sch.sgl', choices=cho, selected=sel, server=TRUE)
    })
  output$sch.box <- renderUI({
    ns <- session$ns; sch.mod <- input$sch.mode
    mul.val <- pre.id$id
    if (!is.null(url.id$sch.mul)) if (url.id$sch.mul[1]!='null') mul.val <- url.id$sch.mul
    mul.val <- paste0(mul.val, collapse=',')
    # At the very beginning, two search buttons and two search boxes are all NULL, since they are not rendered. When toggling between search modes, one of two buttons is reverted to 0 not NULL, while the search content is reverted to NULL then strings.
    if (sch.mod=='Single') selectizeInput(ns('sch.sgl'), p(lab.sgl, actionButton(ns("sch.sgl.but"), "Confirm selection")), choices='none', selected='none', multiple = TRUE, options=list(placeholder = 'Supports partial matching.')) else if (sch.mod=='Multiple') textInput(inputId=ns('sch.mul'), p(lab.mul, actionButton(ns("sch.mul.but"), "Confirm selection")), value=mul.val, placeholder='Muliple IDs must ONLY be separated by space or comma.', width='100%')
  })
  rna <- reactiveValues(val=NULL)
  observe({ # rna$val: accounts for filtered data.
    if (is.null(geneIn$lis)) return()
    rna$val <- rownames(geneIn$lis$df.aggr.tran) 
  })
  upd.sel <- reactiveValues()
  observe({
    if (is.null(input$sch.sgl)) return()
    pars <- list(sch.sgl$cho, sch.sgl$sel, input$sch.mode, rna$val)
    upd.sel$pars <- pars
  })
  observeEvent(upd.sel$pars, {
    updateSelectizeInput(session, 'sch.sgl', choices=sch.sgl$cho, selected=sch.sgl$sel, server=TRUE)
  }) 
  # The button controls all IDs in the downstream SHMs.
  observeEvent(list(input$sch.sgl.but, rna$val), { 
    cat('Single search IDs ... \n')
    if (is.null(input$sch.sgl)) return()
    if (is.null(geneIn$lis)|'none' %in% input$sch.sgl) return()
    sel <- sub(' .*', '', input$sch.sgl)
    # validate: holds on til condition is met, can address issues in execution orders.
    validate(need(all(sel %in% rna$val), ''))
    ids$sel <- sel; ids$but.sgl <- input$sch.sgl.but
    cat('Done! \n')
  })

  observe({ # On-start IDs in single search mode.
    if (length(ids$sel %in% rna$val)>0|is.null(input$sch.sgl)) return()
    if (input$sch.sgl=='none') return()
    ids$sel <- sub(' .*', '', input$sch.sgl); ids$but.sgl <- input$sch.sgl.but
  })
  # The button controls all IDs in the downstream SHMs.
  observeEvent(list(input$sch.mul.but, rna$val), {
    cat('Multiple search IDs ... \n')
    if (is.null(rna$val)|!is.character(input$sch.mul)) return()
    if (input$sch.mul=='') sel <- pre.id$id else {
      gens <- strsplit(gsub(' |,', '__', input$sch.mul), '__')[[1]]
      sel <- gens[tolower(gens) %in% tolower(rna$val)]
      if (length(sel)==0) sel <- pre.id$id
    }; ids$sel <- sel; ids$but.mul <- input$sch.mul.but
    cat('Done! \n') 
  })
 
  observe({ # On-start IDs in multiple search mode.
    sch.mul <- input$sch.mul
    if (length(ids$sel %in% rna$val)>0|is.null(sch.mul)) return()
    if (sch.mul=='') return()
    gens <- strsplit(gsub(' |,', '__', sch.mul), '__')[[1]]
    ids$sel <- gens[tolower(gens) %in% tolower(rna$val)]
    ids$but.mul <- input$sch.mul.but
  })

  return(list(ids=ids))
  onBookmark(function(state) { state })
  })
}

upload_server <- function(id, lis.url=NULL, session) {
  moduleServer(id, function(input, output, session) {

  observeEvent(input$cusHelp, {
    showModal(
    div(id = 'datFormat', modalDialog(title = HTML('<strong><center>Data Formats</center></strong>'),
      div(style = 'overflow-y:scroll;overflow-x:scroll',
      HTML('<img src="image/data_format_shiny.jpg">'),
      tags$a(href="https://bioconductor.org/packages/devel/bioc/vignettes/spatialHeatmap/inst/doc/spatialHeatmap.html", target="_blank", "Package vignette")
      ))))
    })

  cfg <- reactiveValues(lis.dat=NULL, lis.dld=NULL, lis.par=NULL, na.def=NULL, dat.def=NULL, svg.def=NULL, pa.upl=NULL, pa.dat.upl=NULL, pa.svg.upl=NULL, na.cus=NULL)
  observe({
    lis.cfg <- yaml.load_file('config/config.yaml')
    lis.cfg <- lis.cfg[!vapply(lis.cfg, is.null, logical(1))]
    # Separate data sets, download files, and parameters.
    lis.dat <- lis.cfg[grepl('^dataset\\d+', names(lis.cfg))]
    lis.dld <- lis.cfg[grepl('download_single|download_multiple|download_spatial_temporal|download_batched_data_aSVGs', names(lis.cfg))]
    if (is.null(input$config)) lis.par <- lis.cfg[!grepl('^dataset\\d+|download_single|download_multiple|download_spatial_temporal|download_batched_data_aSVGs', names(lis.cfg))] else lis.par <- yaml.load_file(input$config$datapath[1])
    upl.size <- toupper(lis.par$max.upload.size)
    num <- as.numeric(gsub('(\\d+)(G|M)', '\\1', upl.size))
    if (grepl('\\d+G$', upl.size)) max.size <- num*1024^3
    if (grepl('\\d+M$', upl.size)) max.size <- num*1024^2 
    options(shiny.maxRequestSize=max.size) 
    # if (!any(lis.par$hide.legend %in% c('Yes', 'No'))) lis.par$hide.legend <- ifelse(lis.par$hide.legend==TRUE, 'Yes', 'No')
    # Organise configuration parameters in a data frame.
    for (i in seq_along(lis.par)) {
      lis0 <- lis.par[[i]]; if (length(lis0)>1) {
        name <- default <- NULL; for (j in seq_along(lis0)) {
          pair <- strsplit(lis0[j], ':')[[1]]
          name <- c(name, pair[1]); default <- c(default, pair[2])
        }; df0 <- data.frame(name=name, default=default)
        rownames(df0) <- df0$name; lis.par[[i]] <- df0
      }
    }

    # Separate data, svg.
    na.ipt <- dat.ipt <- svg.ipt <- NULL; for (i in lis.dat) {  
      na.ipt <- c(na.ipt, i$name); dat.ipt <- c(dat.ipt, i$data)
      svg.ipt <- c(svg.ipt, list(i$svg))
    }; names(dat.ipt) <- names(svg.ipt) <- na.ipt

    df.tar <- input$tar; dat.upl <- svg.upl <- NULL
    tar.num <- grepl('\\.tar$', df.tar$datapath)
    if (!is.null(df.tar)) validate(need(try(sum(tar.num)==2), 'Two separate tar files of data and aSVGs respectively are expected!'))
    if (sum(tar.num)==2) {

      cat('Processing uploaded tar files ... \n')
      p <- df.tar$datapath[1]; strs <- strsplit(p, '/')[[1]]
      cfg$pa.upl <- pa.svg <- paste0(strs[grep('\\.tar$', strs, invert=TRUE)], collapse='/')
      dat.idx <- grepl('data_shm.tar$', df.tar$name) 
      cfg$pa.svg.upl <- df.tar$datapath[!dat.idx]
      # system(paste0('tar -xf', ' ', svg.tar, ' -C ', pa.svg))
      cfg$pa.dat.upl <- dat.pa <- df.tar$datapath[dat.idx]
      df.pair.upl <- read_hdf5(dat.pa, 'df_pair')[[1]]
      pair.na <- df.pair.upl$name; dat.upl <- df.pair.upl$data
      svg.upl <- as.list(df.pair.upl$aSVG); names(dat.upl) <- names(svg.upl) <- pair.na
      # Process multiple aSVGs under the same data.
      for (i in seq_along(svg.upl)) {
        svg0 <- svg.upl[[i]]; if (grepl(';| |,', svg0)) {
          strs <- strsplit(svg0, ';| |,')[[1]]; svg.upl[[i]] <- strs[strs!='']
        }
      }; cat('Done! \n')
    }
    # Separate data, svg of default and customization. 
    na.def <- na.ipt[!grepl('^none$|^customData$|^customComputedData$', na.ipt)]
    #if (na.cus) 
    na.cus <- c('customData', 'customComputedData')

    dat.def <- c(dat.upl, dat.ipt[na.def]); svg.def <- c(svg.upl, svg.ipt[na.def])
    # If data/svg are duplicated between the server and upload, the data/svg on server side is removed.
    dat.def <- dat.def[unique(names(dat.def))]; svg.def <- svg.def[unique(names(svg.def))]
    cfg$lis.dat <- lis.dat; cfg$lis.dld <- lis.dld; cfg$lis.par <- lis.par; cfg$na.def <- names(dat.def); cfg$svg.def <- svg.def; cfg$dat.def <- dat.def; cfg$na.cus <- na.cus

    dat.nas <- c('none', 'customData', names(dat.def))
    url.val <- url_val('upl-fileIn', lis.url)
    updateSelectInput(session, 'fileIn', NULL, dat.nas, ifelse(url.val=='null', lis.par$default.dataset, url.val))
    updateRadioButtons(session, inputId='dimName', label='2B: is column or row gene?', choices=c("None", "Row", "Column"), selected=lis.par$col.row.gene, inline=TRUE)
    updateRadioButtons(session, inputId='measure', choices=c('Correlation', 'Distance'), selected=lis.par$mhm['measure', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId="cor.abs", choices=c('No', 'Yes'), selected=lis.par$mhm['cor.absolute', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId="thr", selected=lis.par$mhm['select.by', 'default'], inline=TRUE)
    updateNumericInput(session, inputId='mhm.v', value=as.numeric(lis.par$mhm['cutoff', 'default']))
    updateRadioButtons(session, inputId="mat.scale", selected=lis.par$mhm['scale', 'default'], inline=TRUE)
    updateSelectInput(session, inputId="net.type", label="Network type:", choices=c('signed', 'unsigned', 'signed hybrid', 'distance'), selected=lis.par$network['net.type', 'default'])
    updateNumericInput(session, "min.size", "Minmum module size:", value=as.numeric(lis.par$network['min.size', 'default']), min=15, max=5000)
    updateSelectInput(session, "ds","Module splitting sensitivity level", 3:2, selected=lis.par$network['ds', 'default'])
    updateNumericInput(session, "max.edg", "Maximun edges (too many edges may crash the app):", value=cfg$lis.par$network['max.edges', 'default'], min=1, max=500)

  })
  
  observe({
    input$fileIn; input$geneInpath
    updateRadioButtons(session, inputId="dimName", label="Step 2B: is column or row gene?", 
    inline=TRUE, choices=c("None", "Row", "Column"), selected="None")
  })
 observe({ 
    dld.exp <- reactiveValues(sgl=NULL, mul=NULL, st=NULL, bat = NULL)
    dld.exp$sgl <- cfg$lis.dld$download_single
    dld.exp$mul <- cfg$lis.dld$download_multiple
    dld.exp$st <- cfg$lis.dld$download_spatial_temporal
    dld.exp$bat <- cfg$lis.dld$download_batched_data_aSVGs

    output$dld.cfg <- downloadHandler(
      filename=function(){ "config_par.yaml" },
 content=function(file=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/config_par.yaml')){  
        lis.cfg <- yaml.load_file('config/config.yaml')
        par.na <- c("max.upload.size", "default.dataset", "col.row.gene", "separator", "data.matrix", "shm.img", "shm.anm", "shm.video", "legend", "mhm", "network")
        par.na <- par.na[par.na %in% names(lis.cfg)]
        lis.par <- lis.cfg[par.na]; write_yaml(lis.par, file)
      }
    )
    output$dld.sgl <- downloadHandler(
      filename=function(){ "single_aSVG_data.zip" },  content=function(file=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/single_aSVG_data.zip')){ zip(file, c(dld.exp$sgl$data, dld.exp$sgl$svg)) }
    )
    output$dld.mul <- downloadHandler(   
      filename=function(){ "multiple_aSVG_data.zip" }, content=function(file=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/multiple_aSVG_data.zip')){ zip(file, c(dld.exp$mul$data, dld.exp$mul$svg)) }
  )
    output$dld.st <- downloadHandler(   
      filename=function(){ "spatiotemporal_aSVG_data.zip" },
content=function(file=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/spatiotemporal_aSVG_data.zip')){ zip(file, c(dld.exp$st$data, dld.exp$st$svg)) }
  )
    output$dld.bat <- downloadHandler(   
      filename=function(){ "batched_data_aSVGs.zip" },
content=function(file=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/batched_data_aSVGs.zip')){ zip(file, c(dld.exp$bat$data, dld.exp$bat$svg)) }
  )
  })

  # URLs on the landing page.
  output$brain.hum <-renderUI({
  tagList(
    p('Human brain', style='font-size:18px'),
  a(img(width='97%', src="image/brain_hum.png"), href=paste0('http://', lis.url$hos.port, brain.hum.url), target='_blank')
    )
  })
  output$mouse <-renderUI({
  tagList(
    p('Mouse organ', style='font-size:18px'),
    a(img(width='97%', src="image/mouse.png"), href=paste0('http://', lis.url$hos.port, mouse.url), target="_blank")
  )
  })
  output$chicken <-renderUI({
  tagList(
    p('Chicken organ', style='font-size:18px'),
    a(img(width='97%', src="image/chicken.png"), href=paste0('http://', lis.url$hos.port, chicken.url), target="_blank")
    )
  })
  output$organ.arab <-renderUI({
  tagList(
    p('Organ', style='font-size:18px'),
    a(img(width='97%', src="image/organ_arab.png"), href=paste0('http://', lis.url$hos.port, organ.arab.url), target="_blank")
    )
  })
  output$shoot.arab <-renderUI({
  tagList(
    p('Shoot tissue', style='font-size:18px'),
    a(img(width='97%', src="image/shoot_arab.png"), href=paste0('http://', lis.url$hos.port, shoot.arab.url), target="_blank")
  )
  })
  output$root.arab <-renderUI({
  tagList(
    p('Root tissue', style='font-size:18px'),
    a(img(width='97%', src="image/root_arab.png"), href=paste0('http://', lis.url$hos.port, root.arab.url), target="_blank")
    )
  })
  output$stage.arab <-renderUI({
  tagList(
    p('Developmental stage', style='font-size:18px'),
    a(img(width='97%', src="image/stage_arab.png"), href=paste0('http://', lis.url$hos.port, stage.arab.url), target="_blank")
    )
  })
  output$clp.rice <-renderUI({
  tagList(
    p('Rice spatiotemporal coleoptile', style='font-size:18px'),
    a(img(width='97%', src="image/clp_rice.png"), href=paste0('http://', lis.url$hos.port, clp.rice.url), target="_blank")
    )
  })

  # Instruction.
  output$dld <-renderUI({ includeHTML("instruction/download.html") })
  output$sum <-renderUI({ includeHTML("instruction/summary.html") })
  output$input <-renderUI({ includeHTML("instruction/input.html") })
  output$matrix <-renderUI({ includeHTML("instruction/data_matrix.html") })
  output$shm.ins <-renderUI({ includeHTML("instruction/spatial_heatmap.html") })
  output$mhm.ins <-renderUI({ includeHTML("instruction/matrix_heatmap.html") })
  output$net.ins <-renderUI({ includeHTML("instruction/network.html") })
  # Acknowledgement.
  output$ack <-renderUI({ includeHTML("instruction/acknowledgement.html") })

  observe({
    toggleState(id = "geneInpath", condition = input$fileIn == 'customData')
    toggleState(id = "dimName", condition = input$fileIn == 'customData')
    toggleState(id = "target", condition = input$fileIn == 'customData')
    toggleState(id = "met", condition = input$fileIn == 'customData')
    toggleState(id = "svgInpath1", condition = input$fileIn == 'customData')
    toggleState(id = "svgInpath2", condition = input$fileIn == 'customData')
  })
  onBookmark(function(state) { state })

  return(list(ipt = input, cfg = cfg))
})}
  # Selected IDs and button in search box.
  ids <- reactiveValues(sel=NULL, but=NULL)
  mods$upload <- ipt.cfg <- upload_server('upl', lis.url)
  ipt <- ipt.cfg$ipt; cfg <- ipt.cfg$cfg 

data_server <- function(id, ipt, cfg, sch, lis.url, ids, deg = FALSE, session) {
  moduleServer(id, function(input, output, session) {

  # Filter parameters.
  fil <- reactiveValues(P=0, A=0, CV1=-Inf, CV2=Inf)
  observe({
    if (!is.null(lis.url)) return()
    ipt$fileIn; ipt$geneInpath; input$log 
    fil$P <- 0; fil$A <- 0; fil$CV1 <- -Inf; fil$CV2 <- Inf
  })

  # By default, observeEvent is trigerred on start.
  observeEvent(list(input$fil.but, lis.url), {
    # if (ipt$fileIn=="none") return(NULL)
    fil$P <- input$P; fil$A <- input$A; fil$CV1 <- input$CV1; fil$CV2 <- input$CV2
  })

  output$fil.par <- renderText({    
    if (ipt$fileIn=="none") return(NULL)  
    P <- input$P
    validate(need(try(P<=1 & P>=0), 'P should be between 0 to 1 !'))
  })

  # Import data, row metadata, targets file, aggregate replicates.
  geneIn0 <- reactive({
    if (ipt$fileIn=="none") return(NULL)
    withProgress(message="Loading data: ", value = 0, {
    if (any(ipt$fileIn %in% cfg$na.def)) {

      incProgress(0.5, detail="loading matrix, please wait ...")
      dat.na <- cfg$dat.def[ipt$fileIn]
      # All default example data have genes in rows.
      if ('example' %in% strsplit(dat.na, '/')[[1]]) df.te <- fread_df(input=dat.na, isRowGene=TRUE) else { 
        # Exrtact data from uploaded tar.
        dat <- NULL; if (!is.null(cfg$pa.dat.upl)) if (file.exists(cfg$pa.dat.upl)) {
          cat('Extracting uploaded data... \n')
          # The prefix is not input$fileIn. The returned value of read_hdf5 is either data.frame or SE.
          dat <- read_hdf5(cfg$pa.dat.upl, prefix=dat.na)[[1]]
        }
        if (is.null(dat)|is.character(dat)|is.null(ipt$tar)) {
          cat('Extracting internal data... \n')
          dat <- read_hdf5('example/data_shm.tar', dat.na)[[1]]
        }
        validate(need(try(is(dat, 'data.frame')|is(dat, 'SummarizedExperiment')), 'The selected data is empty! Solution: 1) select another data, then rematch its samples to the selected aSVG; 2) include a data for the selected aSVG file in the backend database or the uploaded database.'))
        if (is(dat, 'SummarizedExperiment')) {
          dat <- se_from_db(dat)
          validate(need(try(!is.character(dat)), dat))
        }; df.te <- fread_df(input=dat, isRowGene=TRUE)
      }; return(df.te)
    }; dimNa <- ipt$dimName
    if (ipt$fileIn %in% cfg$na.cus & 
    !is.null(ipt$geneInpath) & dimNa!="None") {
      incProgress(0.25, detail="importing matrix, please wait ...")
      geneInpath <- ipt$geneInpath$datapath; targetInpath <- ipt$target$datapath; metInpath <- ipt$met$datapath
      # Keep replicates unchaged, and compared with targets/metadata files.
      df.upl <- fread_df(read_fr(geneInpath), isRowGene=(dimNa=='Row'), rep.aggr=NULL)
      df.rep <- df.upl$df.rep; df.met <- df.upl$df.met
      if (!is.null(targetInpath)) {
        df.tar <- read_fr(targetInpath)
        # If errors detected, modalDialog terminates the app, while validate only stops the current step. 
        if (nrow(df.tar) != ncol(df.rep)) showModal(modal(msg = 'Ensure "columns" in the data matrix corresponds with "rows" in the targets file respectively!'))
        validate(need(try(nrow(df.tar) == ncol(df.rep)), 'Ensure "columns" in the data matrix corresponds with "rows" in the targets file respectively!'))
        # Check feature/factor columns in targets file.
        cna.tar <- colnames(df.tar) <- tolower(colnames(df.tar))
        if (all(c('feature', 'factor') %in% cna.tar)) cna <- paste0(df.tar$feature, '__', df.tar$factor) else if ('feature' %in% cna.tar) cna <- paste0(df.tar$feature, '__', 'con')
        colnames(df.rep) <- cna
      }     

      if (!is.null(metInpath)) df.met <- read_fr(metInpath)
      if (nrow(df.met) != nrow(df.rep)) showModal(modal(msg = 'Ensure "rows" in the data matrix corresponds with "rows" in the row metadata file respectively!'))
      validate(need(try(nrow(df.met) == nrow(df.rep)), 'Ensure "rows" in the data matrix corresponds with "rows" in the row metadata file respectively!'))
      rownames(df.met) <- rownames(df.rep)
      # Aggregate replicates after targets file is processed.
      df.upl1 <- fread_df(df.rep, isRowGene=(dimNa=='Row'), rep.aggr = 'mean')
      df.upl$df.aggr <- df.upl1$df.aggr
      df.upl$df.rep <- df.rep; df.upl$df.met <- df.met 
      df.upl$con.na <- df.upl1$con.na; cat('Done! \n')
      if ((!is.null(ipt$svgInpath1) | !is.null(ipt$svgInpath2))) return(df.upl) 
    }

    })

  })
  # Transform data.
  geneIn1 <- reactive({
    cat('Scale/tranform data ... \n')
    norm.meth <- input$normDat; scale.dat <- input$scaleDat
    if (deg == FALSE) {
      if (is.null(geneIn0())|is.null(norm.meth)|is.null(scale.dat)|is.null(input$log)) return()
    df.aggr <- df.aggr.tran <- geneIn0()[['df.aggr']]
    # Normalization
    if (grepl('^CNF-', norm.meth)) {
      norm.fun <- 'CNF'
      par.lis <- list(method=sub('.*-', '', norm.meth))
    } else { norm.fun <- norm.meth; par.lis <- NULL }
    if (all(df.aggr.tran==round(df.aggr.tran))) df.aggr.tran <- norm_data(data=df.aggr.tran, norm.fun=norm.fun, parameter.list=par.lis, log2.trans=FALSE)

    if (input$log=='log2') { 
      g.min <- min(df.aggr.tran)
      if (g.min<0) df.aggr.tran <- df.aggr.tran-g.min+1; if (g.min==0) df.aggr.tran <- df.aggr.tran+1; df.aggr.tran <- log2(df.aggr.tran)
    } else if (input$log=='exp2') df.aggr.tran <- 2^df.aggr.tran
    # Scale by row/column
    if (scale.dat=='Row') { df.aggr.tran <- t(scale(t(df.aggr.tran))) } else if (scale.dat=='Selected') {
      # validate(need(!is.null(ids$sel), ''))
      idx.sel <- rownames(df.aggr.tran) %in% ids$sel
      df.sel <- df.aggr.tran[idx.sel, ]
      df.aggr.tran[idx.sel, ] <- scale_all(df.sel) 
    } else if (scale.dat=='All') {
      df.aggr.tran <- scale_all(df.aggr.tran)
    }
    } else if (deg == TRUE) {   
      if (is.null(geneIn0())) return()
      df.aggr <- df.aggr.tran <- geneIn0()[['df.aggr']] 
    }
    df.met <- geneIn0()[['df.met']]; df.rep <- geneIn0()[['df.rep']]; cat('Done! \n')
    return(list(df.aggr=df.aggr, df.aggr.tran=df.aggr.tran, df.met=df.met, df.rep=df.rep))

  })

  output$col.order <- renderUI({
    if (is.null(geneIn1())) return()
    col.nas <- colnames(geneIn1()[['df.aggr.tran']])
    column(12, actionButton("col.cfm", "Confirm", icon=icon("refresh")), 
    selectizeInput(inputId="col.na", label='', choices=col.nas, selected=col.nas, multiple=TRUE, options= list(plugins=list('drag_drop')))
    )
  })
  sear <- reactiveValues(id=NULL)
  observeEvent(ipt$fileIn, { sear$id <- NULL })
  observeEvent(sch$but, {
    if (is.null(geneIn1())) return()
    if (sch$sch=='') sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) else {
      gens <- strsplit(gsub(' |,', '_', sch$sch), '_')[[1]]
      pat <- paste0('^', gens, '$', collapse='|')
      sel <- which(grepl(pat, x=rownames(geneIn1()[['df.aggr.tran']]), ignore.case=TRUE, perl=TRUE))
      if (length(sel)==0) sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default'])
     }; sear$id <- sel
  })

  col.reorder <- reactiveValues(col.re='Y', col.na.re = NULL)
  observe({
    if (is.null(geneIn1())) return()
    # input$col.na is NULL on loading.
    if (is.null(input$col.na)) col.reorder$col.na.re <- colnames(geneIn1()$df.aggr.tran) else col.reorder$col.na.re <- input$col.na
  })

  geneIn <- reactive({
    cat('Filtering data ... \n')
    gene.lis <- geneIn1(); if (is.null(gene.lis)) return()
    df.aggr <- gene.lis[['df.aggr']]; df.aggr.tran <- gene.lis[['df.aggr.tran']]
    df.met <- gene.lis[['df.met']]; df.rep <- gene.lis[['df.rep']]; input$fil.but; lis.url
    if (deg == FALSE) df2fil <- df.aggr.tran else df2fil <- df.rep
    # When ipt$fileIn changes, col.reorder$col.na.re of last session still exists and causes errors.
    if (deg == FALSE) if (!identical(sort(col.reorder$col.na.re), sort(colnames(df2fil)))) return()
    # Input variables in "isolate" will not triger re-excution, but if the whole reactive object is trigered by "input$fil.but" then code inside "isolate" will re-excute.
    scale.dat <- input$scaleDat
    isolate({ 
      se <- SummarizedExperiment(assays=list(expr=as.matrix(df2fil)), rowData=df.met)
      if (ncol(df.met)>0 & 'metadata' %in% colnames(df.met)) ann.col <- 'metadata' else ann.col <- NULL
      # If scaled by row, sd is 1, mean is 0, cv is Inf.
      if (deg == FALSE) CVs <- c(ifelse(scale.dat %in% c('Row', 'Selected'), -Inf, fil$CV1), ifelse(scale.dat %in% c('Row', 'Selected'), Inf, fil$CV2)) else CVs <- c(fil$CV1, fil$CV2)
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV = CVs, dir=NULL, verbose=FALSE)
      if (nrow(se)==0) { validate(need(try(nrow(se)>0), 'All rows are filtered out!')); return() }
      # In case all rows are filtered, the app continues to work without refreshing after the filter parameters are reduced.
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=CVs, dir=NULL, verbose=FALSE)
      df2fil <- as.data.frame(assay(se), stringsAsfactors=FALSE)
      colnames(df2fil) <- make.names(colnames(df2fil))
      df.aggr <- df.aggr[rownames(df2fil), , drop=FALSE]
      df.met <- as.data.frame(rowData(se))[, , drop=FALSE]
    })
    cat('Preparing data matrix ... \n')
    if (deg == FALSE) {
    if (is.null(sear$id)) rows <- seq_len(nrow(df2fil)) else rows <- sear$id
    if (length(rows)==1 & rows[1]==as.numeric(cfg$lis.par$data.matrix['row.selected', 'default'])) rows <- seq_len(nrow(df2fil))
    df.aggr.tran.order <- df2fil[rows, col.reorder$col.na.re, drop=FALSE]
    df.aggr.tran <- df2fil
    df.aggr.tran <- df.aggr.tran[rows, , drop=FALSE]
    df.aggr <- df.aggr[rows, , drop=FALSE]
    df.met <- df.met[rows, , drop = FALSE]
    } else {
      # Only "df.rep" is used.
      df.aggr.tran.order <- df.aggr.tran
      df.rep <- df2fil 
    }
    # df.aggr.tran is used for SHMs.
    cat('Done! \n'); return(list(df.aggr = as.data.frame(df.aggr), df.aggr.tran = as.data.frame(df.aggr.tran), df.aggr.tran.order = as.data.frame(df.aggr.tran.order), df.met=df.met, df.rep = as.data.frame(df.rep)))

  })
  dt.shm <- reactive({
    cat('Preparing data matrix ... \n')
    if (is.null(geneIn())|length(ids$sel)==0) return()
    if ((ipt$fileIn %in% cfg$na.cus & is.null(geneIn()))|ipt$fileIn=="none") return()
    withProgress(message="Data table: ", value = 0, {
      incProgress(0.5, detail="Preparing data matrix, please wait ...")
      if (ipt$fileIn!="none") {
      df.all <- geneIn()
      if (deg == FALSE) gene.dt <- cbind.data.frame(df.all[["df.met"]][, , drop=FALSE], df.all[["df.aggr.tran"]][, , drop=FALSE], stringsAsFactors=FALSE) else gene.dt <- cbind.data.frame(df.all[["df.met"]][, , drop=FALSE], df.all[["df.rep"]][, , drop=FALSE], stringsAsFactors=FALSE)
   } 
   # if (is.null(sear$id)) sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) else sel <- sear$id
   # if (length(sel)==1 & sel[1]==as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) & nrow(gene.dt)>1) sel <- sel else if (nrow(gene.dt)==1)  sel <- 1 else if (length(sel)>1) sel <- seq_along(sel)
   # if (deg == FALSE) selection <- list(mode="multiple", target="row", selected=sel) else selection <- 'none'
   cat('Done!\n'); gene.dt
    })
  })

  if (deg==FALSE) output$selProf <- renderPlot(height=800, {
    cat('Profile of selected genes ... \n')
    gene.dt <- dt.shm(); df.all <- geneIn()
    if (is.null(gene.dt)|length(ids$sel)==0) return()
    validate(need(length(ids$sel)<=100, 'Due to space limitation, profiles of 100+ genes are not plotted!'))
    dt.sel <- gene.dt[ids$sel, , drop=FALSE]
    dt.sel <- dt.sel[, !colnames(dt.sel) %in% c('metadata', 'link'), drop=FALSE]
    scale.dat <- input$scaleDat; if (is.null(scale.dat)) return()
    if (!is.null(scale.dat)) if (scale.dat=='No') title <- 'No scaling' else if (scale.dat=='Row') title <- 'Scaled by row' else if (scale.dat=='Selected') title <- 'Scaled across selected genes' else if (scale.dat=='All') title <- 'Scaled across all genes' else title <- ''
    g1 <- profile_gene(dt.sel, y.title=title); 
    g2 <- profile_gene(df.all$df.aggr[ids$sel, , drop=FALSE], y.title='No scaling')
    cat('Done! \n'); grid.arrange(g1, g2, nrow=2)
  })
  if (deg==FALSE) output$dtSel <- renderDataTable({
    cat('Preparing selected data matrix ... \n')
    gene.dt <- dt.shm()
    if (is.null(gene.dt)|length(ids$sel)==0) return()
    # Decimals.
    df.num <- gene.dt[, !colnames(gene.dt) %in% c('metadata', 'link'), drop=FALSE]
    if (all(df.num==round(df.num))) deci <- 0 else deci <- 2
    dt.sel <- gene.dt[ids$sel, , drop=FALSE]
    # Tooltip on metadata.
    col1 <- list(list(targets = c(1), render = JS("$.fn.dataTable.render.ellipsis(5, false)")))
    # In case no metadata column.
    if (colnames(dt.sel)[1]!='metadata') col1 <- NULL
    dtab <- datatable(dt.sel, selection='none', escape=FALSE, filter="top", extensions=c('Scroller'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE, searchHighlight=FALSE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs=col1), 
   class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
   formatRound(colnames(geneIn()[["df.aggr.tran"]]), deci)
    cat('Done! \n'); dtab
  })
  dt.sel <- reactiveValues(val='none')
  if (deg==FALSE) output$dtAll <- renderDataTable({
    cat('Preparing complete data matrix ... \n')
    gene.dt <- dt.shm(); if (is.null(gene.dt)) return()
    # Decimals.
    df.num <- gene.dt[, !colnames(gene.dt) %in% c('metadata', 'link'), drop=FALSE]
    if (all(df.num==round(df.num))) deci <- 0 else deci <- 2
    # Tooltip on metadata.
    col1 <- list(list(targets = c(1), render = JS("$.fn.dataTable.render.ellipsis(5, false)")))
    if (colnames(gene.dt)[1]!='metadata') col1 <- NULL
    dtab <- datatable(gene.dt, selection=list(mode="multiple", target="row", selected=dt.sel$val), escape=FALSE, filter="top", extensions=c('Scroller'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE, searchHighlight=FALSE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs=col1), 
   class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
   formatRound(colnames(geneIn()[["df.aggr.tran"]]), deci)
    cat('Done! \n'); dtab
  })

  if (deg==TRUE) output$dtRep <- renderDataTable({
    cat('Preparing complete data matrix ... \n')
    gene.dt <- dt.shm(); if (is.null(gene.dt)) return()
    col1 <- list(list(targets = c(1), render = JS("$.fn.dataTable.render.ellipsis(5, false)")))
    # In case no metadata column.
    if (colnames(gene.dt)[1]!='metadata') col1 <- NULL
    dtab <- datatable(gene.dt, selection='none', escape=FALSE, filter="top", extensions=c('Scroller'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE, searchHighlight=FALSE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs=col1), 
   class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
   formatRound(colnames(geneIn()[["df.aggr.tran"]]), 0)
    cat('Done! \n'); dtab
  })
  observeEvent(input$dat.all.but, { # Select genes in table.
    gene.dt <- dt.shm(); if (is.null(gene.dt)) return()
    ids$sel <- rownames(gene.dt)[input$dtAll_rows_selected]
    dt.sel$val <- input$dtAll_rows_selected
  })
  observe({
    ipt$fileIn; ipt$geneInpath; lis.par <- cfg$lis.par; lis.url
    url.val <- url_val('dat-log', lis.url)
    updateSelectInput(session, inputId='log', label='Log/exp transform', selected=ifelse(url.val!='null', url.val, cfg$lis.par$data.matrix['log.exp', 'default']), )
    url.val <- url_val('dat-scale', lis.url)
    updateSelectInput(session, 'scaleDat', label='Scale by', selected=ifelse(url.val!='null', url.val, cfg$lis.par$data.matrix['scale', 'default']))
  })

  observe({
    ipt$fileIn; ipt$geneInpath; input$log; lis.url
    url.val <- url_val('dat-A', lis.url)
    updateNumericInput(session, inputId="A", label="Threshold (A) to exceed", value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$data.matrix['A', 'default']))) 
    url.val <- url_val('dat-P', lis.url)
    updateNumericInput(session, inputId="P", label="Proportion (P) of samples with values >= A", value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$data.matrix['P', 'default'])), min=0, max=1)
    url.val <- url_val('dat-CV1', lis.url)
    updateNumericInput(session, inputId="CV1", label="Min coefficient of variation (CV1)", value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$data.matrix['CV1', 'default'])))
    url.val <- url_val('dat-CV2', lis.url)
    updateNumericInput(session, inputId="CV2", label="Max coefficient of variation (CV2)", value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$data.matrix['CV2', 'default']))) 
  })
  
  ipt.dat <- reactiveValues()
  observe({ ipt.dat$dt_rows_selected <- input$dt_rows_selected })
  col.cfm <- reactive({ input$col.cfm })
  col.na <- reactive({ input$na })
  scaleDat <- reactive({ input$scaleDat })
  log <- reactive({ input$log }); A <- reactive({ input$A })
  CV1 <- reactive({ input$CV1 }); CV2 <- reactive({ input$CV2 })
  P <- reactive({ input$P })
  search.but <- reactive({ sch$but })
  onBookmark(function(state) { state })
  return(list(geneIn0 = geneIn0, geneIn1 = geneIn1, geneIn = geneIn, sear = sear, ipt.dat = ipt.dat, col.reorder = col.reorder, col.cfm = col.cfm, col.na = col.na, scaleDat=scaleDat, log = log, A = A, P = P, CV1 = CV1, CV2 = CV2, search.but = search.but))
  })

}

  sch <- reactiveValues()
  observe({ sch$sch <- input$search; sch$but <- input$search.but })
  #ipt0 <- reactiveValues(url.but=NULL, but=NULL)
  # observe({ ipt0$url.but <- input$url.but })
  mods$data <- dat.mod.lis <- data_server('dat', ipt, cfg, sch, lis.url, ids, deg=FALSE)

mods$search <- sch.mod.lis <- search_server('sear', ids, cfg, lis.url, url.id, dat.mod.lis)
 
network_server <- function(id, ipt, cfg, dat.mod.lis, shm.mod.lis, sch.mod.lis, session) {
  
  moduleServer(id, function(input, output, session) {

  ipt.dat <- reactiveValues()
  ipt.dat$dat <- dat.mod.lis$ipt.dat; sear <- dat.mod.lis$sear
  col.reorder <- reactiveValues(); col.reorder <- dat.mod.lis$col.reorder
  geneIn0 <- dat.mod.lis$geneIn0; geneIn1 <- dat.mod.lis$geneIn1
  geneIn <- dat.mod.lis$geneIn
  ipt.dat$dat <- dat.mod.lis$ipt.dat; A <- dat.mod.lis$A
  P <- dat.mod.lis$P; CV1 <- dat.mod.lis$CV1
  CV2 <- dat.mod.lis$CV2
  gID <- shm.mod.lis$gID; geneIn <- dat.mod.lis$geneIn
  ids <- sch.mod.lis$ids
  observe({
    geneIn(); ipt$adj.modInpath; input$A; input$P; input$CV1
    input$CV2; input$min.size; input$net.type
    input$measure; input$cor.abs; input$thr; input$mhm.v
    updateRadioButtons(session, "mat.scale", choices=c("No", "Column", "Row"), selected="Row", inline=TRUE)
  })

  observe({  
    ipt$fileIn; geneIn(); input$adj.modInpath; input$A; input$P; input$CV1; input$CV2; ids$sel
    updateActionButton(session, inputId='mhm.but', icon=icon("refresh"))
    #updateRadioButtons(session, inputId="mhm.but", label="Show plot:", choices=c("Yes", "No"), selected=cfg$lis.par$mhm['show', 'default'], inline=TRUE)
  })
# Avoid unnecessay correlation/distance computation if geneIn() updates due to column reordering. For example, correlation/distance, extracting coordinates, which only depend on the df.aggr.
  df.net <- reactive({
    if (is.null(geneIn())) return()
    return(list(df.aggr=geneIn()[['df.aggr']], df.met=geneIn()[['df.met']]))
  })
  tab.act.lis <- shm.mod.lis$tab.act.lis
  tab.mhm <- reactiveValues(val='no')
  observe({
    shmMhNet <- tab.act.lis$shmMhNet; clusNav <- input$clusNav
    if (is.null(shmMhNet)|is.null(clusNav)) return()
    tab.mhm$val <- ifelse(shmMhNet=='clus' & clusNav=='mhmPlot', 'yes', 'no')
    if (tab.mhm$val=='yes' & input$mhm.but==0) showModal(modal(msg=HTML('To see the latest matrix heatmap, always click the button <strong>"Click to show/update"</strong>!'), easyClose=TRUE))
  })
  # Combine all relevant parameters to pars. After all parameters are adjusted, they only are controlled by buttons, which avoids execution after each parameter is adjusted.
  cor.dis.par <- reactiveValues()
  observeEvent(list(input$mhm.but, input$cpt.nw), ignoreInit=TRUE, {
    pars <- list(cor.abs=input$cor.abs, measure=input$measure)
    cor.dis.par$pars <- pars
  })
  # cor.dis <- reactiveValues(val=NULL)
  # Calculate whole correlation or distance matrix.
  cor.dis <- eventReactive(list(cor.dis.par$pars), ignoreInit=TRUE, valueExpr={
    cat('Correlation/distance matrix ... \n')
    # if (input$mhm.but==0 & input$cpt.nw==0) return()
    pars <- cor.dis.par$pars 
    if (is.null(ipt$fileIn)) return()
    if (ipt$fileIn %in% cfg$na.cus & is.null(ipt$svgInpath1) & is.null(ipt$svgInpath2)) return()
    # if (is.null(input$mhm.but)) return() 
    gene <- df.net()$df.aggr; if (is.null(gene)) return()
    if ((ipt$fileIn %in% cfg$na.cus & is.null(gene))|ipt$fileIn=="none") return(NULL)
    withProgress(message="Compute similarity/distance matrix: ", value = 0, {
      incProgress(0.5, detail="please wait ...")
      # Too many genes may crash the app.
      if (nrow(gene)>15000) showModal(modal(msg=strong('Too many genes (e.g. 15,000+) may crash the app!'), easyClose=TRUE))
      # if (nrow(gene)>15000 & input$mhm.but==0) return()
      if (input$measure=='Correlation') {
        m <- cor(x=t(gene)); cat('Done! \n')
        if (input$cor.abs==TRUE) { m <- abs(m) }; return(m)
      } else if (input$measure=='Distance') { cat('Done! \n'); return(-as.matrix(dist(x=gene))) }
    })
  })
  # Subset nearest neighbours for target genes based on correlation or distance matrix.
  # observe({
    # The submat reactive expression is only accessible inside the same obeserve environment.
    # submat <- reactive({})
   #  if (tab.mhm$val!='yes') { return() }
  submat.par <- reactiveValues()
  observeEvent(list(input$mhm.but, input$cpt.nw), ignoreInit=TRUE, {
    pars <- list(input$thr, input$mhm.v, cor.dis.par$pars)
    submat.par$pars <- pars
  })
  # Avoid calling eventReactive with ignoreInit=TRUE in another eventReactive also with ignoreInit=TRUE, since the latter will be slient even though the former is triggered.
  submat <- eventReactive(list(submat.par$pars), { 
    cat('Subsetting nearest neighbors ... \n')
    if (ipt$fileIn=="none") return()
    # if (is.null(input$mhm.but)) return()
    corDis <- cor.dis()
    if (is.null(corDis)) return()
    gene <- df.net()$df.aggr; rna <- rownames(gene)
    gen.tar<- gID$geneSel; 
    # mat <- cor.dis$val()
    # Validate filtering parameters in matrix heatmap. 
    measure <- input$measure; cor.abs <- input$cor.abs; mhm.v <- input$mhm.v; thr <- input$thr
    if (input$thr=='p') {
      validate(need(try(mhm.v>0 & mhm.v<=1), 'Proportion should be between 0 to 1 !'))
    } else if (input$thr=='n') {
      validate(need(try(mhm.v>=1 & as.integer(mhm.v)==mhm.v & !is.na(mhm.v)), 'Number should be a positive integer !'))
    } else if (input$thr=='v' & measure=='correlation') {
      validate(need(try(mhm.v>-1 & mhm.v <1), 'Correlation value should be between -1 to 1 !'))
    } else if (input$thr=='v' & measure=='distance') {
      validate(need(try(mhm.v>=0), 'Distance value should be non-negative !'))
    }
    withProgress(message="Selecting nearest neighbours: ", value = 0, {
      incProgress(0.5, detail="please wait ...")
      arg <- list(p=NULL, n=NULL, v=NULL)
      arg[names(arg) %in% input$thr] <- input$mhm.v
      if (input$measure=='Distance' & input$thr=='v') arg['v'] <- -arg[['v']]
      if (!all(gen.tar %in% rownames(corDis))) return()    
      validate(need(try(ncol(gene)>4), 'The "sample__condition" variables in the Data Matrix are less than 5, so no coexpression analysis is applied!'))
      gen.na <- do.call(sub_na, c(mat=list(corDis), ID=list(gen.tar), arg))
      if (any(is.na(gen.na))) return() 
      validate(need(try(length(gen.na)>=2), paste0('Only ', gen.na, ' is remaining!'))); cat('Done! \n'); return(gene[gen.na, ])
    })
  })
 # })

  mhm.par <- reactiveValues()
  observeEvent(input$mhm.but, ignoreInit=TRUE, {
    # Include all upstream parameter changes. Since correlation/distance paremeter changes are included in submat.par$pars, only submat.par$pars represents upstream parameter change.
    pars <- list(mat.scale=input$mat.scale, submat.par$pars)
    mhm.par$pars <- pars
  })
  # mhm <- reactiveValues(hm=NULL)
  # Plot matrix heatmap.
  mhm <- eventReactive(list(mhm.par$pars), {
    cat('Initial matrix heatmap ... \n')
    if (is.null(mhm.par$pars)) return()
    #if (is.null(input$mhm.but)) return() # Matrix heatmap sections is removed.
    # if (input$mhm.but!=0) return();
    # mhm$hm <- NULL
    # In submat, there is a validate expression. Only if submat is included in reactive({}) not observe({}), can the message in validate be propagated to the user interface.
    # sub.mat <- submat$val; if (is.null(sub.mat)) return()
    # Reactive expression of sub.mat() is NULL, but sub.mat is not NULL.
    # if (is.null(sub.mat())) return()
    sub.mat <- submat(); if (is.null(sub.mat)) return()
    gene <- df.net()$df.aggr; rna <- rownames(gene)
    gen.tar <- gID$geneSel; # if (length(gen.tar)>1) return()
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="Plotting ...")
      if (input$mat.scale=="Column") scale.hm <- 'column' else if (input$mat.scale=="Row") scale.hm <- 'row' else scale.hm <- 'no'

      hm <- matrix_hm(ID=gen.tar, data=sub.mat, col=c('yellow', 'red'), scale=scale.hm, main='Target Genes and Their Nearest Neighbours', title.size=10, static=FALSE)
      cat('Done! \n'); hm
    })
  })
  if (0) hmly <- eventReactive(input$mhm.but, {
    cat('Matrix heatmap ... \n')
    #if (is.null(submat())|input$mhm.but=='No') return()
    if (is.null(submat())) return()
    gene <- df.net()$df.aggr; rna <- rownames(gene)
    gen.tar<- gID$geneSel
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="preparing ...")
      if (input$mat.scale=="Column") scale.hm <- 'column' else if (input$mat.scale=="Row") scale.hm <- 'row' else scale.hm <- 'no'  
      cat('Done!')
      matrix_hm(ID=gen.tar, data=submat(), scale=scale.hm, main='Target Genes and Their Nearest Neighbours', title.size=10, static=FALSE)
    })
  })

  output$HMly <- renderPlotly({
    # if (tab.mhm$val!='yes') return()
    if (length(ids$sel)==0) return()
    if (is.null(gID$geneSel)|is.null(submat())) return()
    if (gID$geneSel[1]=='none'|is.na(gID$geneSel[1])) return()
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="plotting ...")
    # if (input$mhm.but!=0) hmly() else if (input$mhm.but==0) mhm$hm else return()
      mhm()
    })
  })

  tab.net <- reactiveValues(val='no')
  observe({
    shmMhNet <- tab.act.lis$shmMhNet; clusNav <- input$clusNav
    if (is.null(shmMhNet)|is.null(clusNav)) return()
    tab.net$val <- ifelse(shmMhNet=='clus' & clusNav=='netPlot', 'yes', 'no')
    if (input$cpt.nw==0 & tab.net$val=='yes') showModal(modal(msg=HTML('To see the latest network, always click the button <strong>"Click to show/update"</strong>!'), easyClose=TRUE))
  })
    #gene <- geneIn()[["df.aggr.tran"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    if (0) observe({
      cat('Initial adjacency matrix and modules ...\n')
      if (ipt$fileIn=="none"|tab.net$val!='yes') return()
      if (is.null(input$cpt.nw)) return() # Network section is removed.
      if (is.null(submat())|input$cpt.nw!=0|length(gID$geneSel)>1) return()
    if (ipt$fileIn=="customData"|any(ipt$fileIn %in% cfg$na.def)) {
      gene <- df.net()$df.aggr; if (is.null(gene)) return()
      type <- input$net.type; sft <- if (type=='distance') 1 else 6
      withProgress(message="Computing: ", value = 0, {
        incProgress(0.3, detail="adjacency matrix ...")
        incProgress(0.5, detail="topological overlap matrix ...")
        incProgress(0.1, detail="dynamic tree cutting ...")    
        adj.mods$lis <- reactive({ 
          lis <- adj_mod(data=submat(), type=type, minSize=input$min.size, dir=NULL); cat('Done! \n'); lis
        })
      })
    }
    })

  # er <- eventReactive(exp, {}). If its reactive value "er()" is called before eventReactive is triggered, the code execution stops where "er()" is called. 
  adj.mod.par <- reactiveValues()
  observeEvent(input$cpt.nw, ignoreInit=TRUE, {
    pars <- list(gen.sel=input$gen.sel, ds=input$ds, type=input$net.type, min.size=input$min.size, submat.par$pars)
    adj.mod.par$pars <- pars
  })
  adj.mods <- eventReactive(adj.mod.par$pars, {
    cat('Adjacency and modules ... \n')
    if (ipt$fileIn=="none") return()
    if (ipt$fileIn %in% cfg$na.cus & is.null(ipt$svgInpath1) & is.null(ipt$svgInpath2)) return()
    #gene <- geneIn()[["df.aggr.tran"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    # if (is.null(submat$val())|input$cpt.nw==0) return()
    # sub.mat <- submat$val; if (is.null(sub.mat)) return()
    # if (is.null(sub.mat())) return()
    sub.mat <- submat(); if (is.null(sub.mat)) return()
    # gene <- pars$df.aggr; if (is.null(gene)) return()
    withProgress(message="Computing: ", value = 0, {
      incProgress(0.3, detail="adjacency matrix ...")
      incProgress(0.5, detail="topological overlap matrix ...")
      incProgress(0.1, detail="dynamic tree cutting ...")
        lis <- adj_mod(data=sub.mat, type=input$net.type, minSize=input$min.size); cat('Done! \n'); lis
      })
    })

  observe({
    if (is.null(geneIn())) return(NULL)
    if (length(ids$sel)==0) return()
    gens.sel <- ids$sel
    updateSelectInput(session, inputId="gen.sel", label="", choices=c("None", gens.sel), selected=gens.sel[1])
  })
  observe({ 
    input$gen.sel; input$measure; input$cor.abs; input$thr; input$mhm.v
    updateSelectInput(session, 'ds', choices=3:2, selected=cfg$lis.par$network['ds', 'default'])
  })

  #mcol <- reactive({

   # if ((input$cpt.nw!=0|!is.null(adj.mods$lis)) & (input$fileIn=="customData"|any(input$fileIn %in% cfg$na.def))) { 

    #withProgress(message="Computing dendrogram:", value=0, {
     # incProgress(0.7, detail="hierarchical clustering.")
      #if (input$cpt.nw!=0) mod4 <- adj.mods$lis[['mod']] else mod4 <- adj.mods$lis[['mod']]
 
    #}); return(mod4) 
    
 # }

 #})


  #observe({
  
   # geneIn(); gID$geneSel; input$adj.in; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    #updateRadioButtons(session, inputId="cpt.nw", label="Show plot:", choices=c("Yes", "No"), selected=ifelse(nrow(visNet()[["link"]])<120, "Yes", "No"), inline=TRUE)
  
  #})

  col.sch.net <- reactive({ 

    if(input$color.net=="") return(NULL) 
    col <- gsub(' |\\.|-|;|,|/', '_', input$color.net)
    col <- strsplit(col, '_')[[1]]
    col <- col[col!='']; col1 <- col[!col %in% colors()]
    if (length(col1>0)) validate(need(try(col1 %in% colors()), paste0('Colors not valid: ', paste0(col1, collapse=', '), '!'))); col

  }); color.net <- reactiveValues(col.net="none")

  len.cs.net <- 350
  observeEvent(input$col.but.net, {
    if (is.null(col.sch.net())) return (NULL)
    color.net$col.net <- colorRampPalette(col.sch.net())(len.cs.net)
  })
  nod.edg.par <- reactiveValues()
  observeEvent(input$cpt.nw, ignoreInit=TRUE, {
    pars <- list(input$gen.sel, input$ds, input$net.type, input$min.size, input$max.edg, input$adj.in, adj.mod.par$pars)
    nod.edg.par$pars <- pars
  })
  # visNet <- reactiveValues(val=NULL)
  nod.edg <- eventReactive(nod.edg.par$pars, {
    cat('Network nodes and edges ... \n')
    # input$cpt.nw; 
    adj.mods.lis <- adj.mods()
    if (ipt$fileIn=="none"|is.null(adj.mods.lis)) return()
    if (is.null(input$gen.sel)) return() # Matrix heatmap section is removed.
    # if (ipt$fileIn=='customComputedData' & is.null(df.net())) return()
    # if (input$adj.in=="None") return(NULL)
    # if (ipt$fileIn=="customComputedData") { adj <- adj.mod()[['adj']]; mods <- adj.mod()[['mcol']] } else 
    if (ipt$fileIn=="customData"|ipt$fileIn %in% cfg$na.def) { 
      adj <- adj.mods.lis[['adj']]; mods <- adj.mods.lis[['mod']]
    }
    # if (ipt$fileIn=='customComputedData') gene <- df.net()$df.aggr else {
    gene <- submat(); if (is.null(gene)) return()
    if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    lab <- mods[, input$ds][rownames(gene)==input$gen.sel]
    # validate(need(try(length(lab)==1 & !is.na(lab) & nrow(mods)==nrow(gene)), 'Click "Update" to display new network!'))
    if (length(lab)==0) return()
    if (length(lab)>1|is.na(lab)) return() # When input$fileIn is changed, gene is changed also, but mods is not since it is controled by observeEvent.
    validate(need(try(lab!='0'), 'Warning: the selected gene is not assigned to any module. Please select a different one or adjust the "Minmum module size"!'))
    idx.m <- mods[, input$ds]==lab; adj.m <- adj[idx.m, idx.m]; gen.na <- colnames(adj.m) 
    idx.sel <- grep(paste0("^", input$gen.sel, "$"), gen.na); gen.na[idx.sel] <- paste0(input$gen.sel, "_target")
    colnames(adj.m) <- rownames(adj.m) <- gen.na
    withProgress(message="Computing network:", value=0, { 
      incProgress(0.8, detail="making network data frame ...")
      cat('Extracting nodes and edges... \n')
      # Identify adjcency threshold with edges < max number (e.g. 300) 
      ID <- input$gen.sel; adjs <- 1; lin <- 0; adj.lin.vec <- NULL
      validate(need(try(as.integer(input$max.edg)==input$max.edg), 'The number of edges should be an integer!'))
      # Compute the min adj.
      while (lin<input$max.edg) {
          
          adjs <- adjs-0.002; if (adjs<=10^-15) adjs <- 0
          nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adjs)
          lin <- nrow(nod.lin[['link']])
          vec0 <- adjs; names(vec0) <- lin
          adj.lin.vec <- c(adj.lin.vec, vec0)
          if (adjs==0) break

      }; cat('Adjacency-edge pairs done! \n')
      # The first version of links computed from the min adj or the input adj, which is subjected to the following check.
      nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=ifelse(input$adj.in==1, adjs, input$adj.in))
      link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'
      # If the links are 0 due to the input adj, change the "adjs" to the value bringing 1 or 2 links.
      lins <- NULL; if (nrow(link1)==0) {

        adjs <- adj.lin.vec[names(adj.lin.vec)>=1][1]
        nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adjs)
        link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)

      } else if (nrow(link1)>input$max.edg) {
       
        # If the links are larger than max links due to the input adj, change the "adjs" to the value producing max links.
        adjs <- adjs+0.002
        nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adjs)
        link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)

      } else if (nrow(link1)<=input$max.edg & input$adj.in!=1) {
        
        # If 0<link total<max links, use the input adj.
        adjs <- input$adj.in
        nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adjs)
        link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)

      }
      node <- nod.lin[['node']]; colnames(node) <- c('id', 'value') 
      if (nrow(link1)!=0) { 
        
        link1$title <- link1$value # 'length' is not well indicative of adjacency value, so replaced by 'value'.
        link1$color <- 'lightblue'
        
      }; df.met <- df.net()$df.met
      if ('metadata' %in% colnames(df.met)) meta <- df.met[, 'metadata', drop=FALSE] else meta <- data.frame()
      if (ncol(meta)>0) node <- cbind(node, title=meta[node$id, ], borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      if (ncol(meta)==0) node <- cbind(node, borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      net.lis <- list(node=node, link=link1, adjs=adjs, lins=lins)

    }); cat('Done! \n'); net.lis

  })
  # The order of reactive expression matters so "updateSelectInput" of "adj.in" should be after visNet().
  observe({
 
    if (ipt$fileIn=="none") return()
    geneIn(); gID$geneSel; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$P; input$CV1; input$CV2; input$min.size; input$net.type
    input$gen.sel; input$measure; input$cor.abs; input$thr; input$mhm.v; input$cpt.nw
     #if ((input$adj.in==1 & is.null(visNet()[["adjs1"]]))|(input$cpt.nw!=cfg$lis.par$network['max.edges', 'default'] & is.null(visNet()[["adjs1"]]))) { updateSelectInput(session, "adj.in", "Adjacency threshold:", sort(seq(0, 1, 0.002), decreasing=TRUE), visNet()[["adjs"]]) } else if (!is.null(visNet()[["adjs1"]])) updateSelectInput(session, "adj.in", "Adjacency threshold:", sort(seq(0, 1, 0.002), decreasing=TRUE), visNet()[["adjs1"]])
    # if (is.null(visNet$val)) return(); if (is.null(visNet$val())) return()
    nd.eg <- nod.edg(); if (is.null(nd.eg)) return()
    lins <- nd.eg[["lins"]]
    if (is.null(input$adj.in)) return() # Network section is removed. 
    if (input$adj.in==1|is.null(lins)|is.numeric(lins)) updateSelectInput(session, "adj.in", choices=sort(seq(0, 1, 0.002), decreasing=TRUE), selected=as.numeric(nd.eg[["adjs"]])) 
  
  })
  output$bar.net <- renderPlot({ 
    cat('Network bar ... \n')
    nd.eg <- nod.edg(); if (is.null(nd.eg)) return()
    #if (input$adj.in=="None"|input$cpt.nw=="No") return(NULL)
    if (input$adj.in=="None") return(NULL)
    if (length(color.net$col.net=="none")==0) return(NULL)
    gene <- df.net()$df.aggr; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    if(input$col.but.net==0) color.net$col.net <- colorRampPalette(col_sep(cfg$lis.par$network['color', 'default']))(len.cs.net) # color.net$col.net is changed alse outside renderPlot, since it is a reactive value.
      withProgress(message="Color scale: ", value = 0, {
      incProgress(0.25, detail="preparing data, please wait ...")
      incProgress(0.75, detail="plotting, please wait ...")
      node <- nd.eg[["node"]]; if (is.null(node)) return()
      node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
      cs.net <- col_bar(geneV=v.net, cols=color.net$col.net, width=1)
      cat('Done! \n'); return(cs.net) # '((max(v.net)-min(v.net))/len.cs.net)*0.7' avoids bar overlap.

      })

  })

  observeEvent(nod.edg(), ignoreInit=TRUE, {
    nd.eg <- nod.edg(); if (is.null(nd.eg)) return()
    output$edge <- renderUI({ 
      cat('Remaining edges ... \n')
      if (input$adj.in=="None"|is.null(nd.eg)) return(NULL)
      if (ipt$fileIn=="none"|(ipt$fileIn=="customData" & is.null(geneIn()))|
      input$gen.sel=="None") return(NULL)
      span(style = "color:black;font-weight:NULL;", HTML(paste0("Remaining edges: ", dim((nd.eg[["link"]]))[1])))
      cat('Done! \n')
    })

  })
  vis.net <- reactive({
    cat('Network ... \n')
    #if (input$adj.in=="None"|input$cpt.nw=="No") return(NULL)
    if (input$adj.in=="None") return(NULL)
    nd.eg <- nod.edg(); if (is.null(nd.eg)) return()
    gene <- df.net()$df.aggr; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    withProgress(message="Network:", value=0.5, {
    incProgress(0.3, detail="prepare for plotting ...")
    # Match colours with gene connectivity by approximation.
    node <- nd.eg[["node"]]; if (is.null(node)) return() 
    node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
    col.nod <- NULL; for (i in node$value) {
      ab <- abs(i-v.net); col.nod <- c(col.nod, color.net$col.net[which(ab==min(ab))[1]])
    }; node$color <- col.nod
    cat('Done! \n')
    visNetwork(node, nd.eg[["link"]], height="300px", width="100%", background="", main=paste0("Network Module Containing ", input$gen.sel), submain="", footer= "") %>% visIgraphLayout(physics=FALSE, smooth=TRUE) %>% visOptions(highlightNearest=list(enabled=TRUE, hover=TRUE), nodesIdSelection=TRUE)
    })
    
  })

  output$vis <- renderVisNetwork({

    if (length(ids$sel)==0) return()
    if (ipt$fileIn=="none"|is.null(vis.net())) return(NULL)
    # if (input$cpt.nw=="No") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="plotting ...")
      cat('Rendering network...\n'); vis.net()

    })

  })
  onBookmark(function(state) { state })

})

}
# network_server('net', ipt, cfg, dat.mod.lis, shm.mod.lis)

shm_server <- function(id, ipt, cfg, sch, lis.url, url.id, dat.mod.lis, sch.mod.lis, session) {  
  moduleServer(id, function(input, output, session) {
  # The reactive type in and outside module is the same: sear is a reactiveValue in and outside module; geneIn is reactive expression in and outside module. "geneIn()" is accessing the content of a reactive expression, and loses the "reactive" attribute.
  # As long as the content of reactiveValues (col.reorder$col.na.re) is not accessed, the operation does not need to be inside reactive environment (observe).
  ipt.dat <- reactiveValues()
  ipt.dat$dat <- dat.mod.lis$ipt.dat; sear <- dat.mod.lis$sear
  col.reorder <- reactiveValues(); col.reorder <- dat.mod.lis$col.reorder
  geneIn0 <- dat.mod.lis$geneIn0; geneIn1 <- dat.mod.lis$geneIn1
  geneIn <- dat.mod.lis$geneIn
  col.na <- dat.mod.lis$col.na; col.cfm <- dat.mod.lis$col.cfm 
  scaleDat <- dat.mod.lis$scaleDat
  log <- dat.mod.lis$log; A <- dat.mod.lis$A
  search.but <- dat.mod.lis$search.but
  ids <- sch.mod.lis$ids
  gID <- reactiveValues(geneSel="none", new=NULL, all=NULL)
  observe({ ipt$geneInpath; ipt$fileIn; gID$geneSel <- "none" })
  observe({ if (is.null(geneIn())) gID$geneSel <- "none" })
  # To make the "gID$new" and "gID$all" updated with the new "input$fileIn", since the selected row is fixed (3rd row), the "gID$new" is not updated when "input$fileIn" is changed, and the downstream is not updated either. The shoot/root examples use the same data matrix, so the "gID$all" is the same (pre-selected 3rd row) when change from the default "shoot" to others like "organ". As a result, the "gene$new" is null and downstream is not updated. Also the "gene$new" is the same when change from shoot to organ, and downstream is not updated, thus "gene$new" and "gene$all" are both set NULL above upon new "input$fileIn".  

  init <- reactiveValues(but=0, new=0)
  # observeEvent(ids$but, { init$but <- ids$but })
  # observeEvent(session, { init$n <- init$n+1; print(init$n)})
  rna.fil <- reactiveValues(val=NULL)
  observe({ # Filtered data.
    if (length(ids$sel)==0) return(); if (ids$sel[1]=='') return()
    rna.fil$val <- rownames(geneIn()$df.aggr.tran)
  })
  # The on-start ID processing is controlled by 0 and 1 states.
  observeEvent(list(session, ipt$fileIn, rna.fil$val), { init$but <- 0 })
  observeEvent(list(ids$sel, rna.fil$val), { # All on-start and on-start similar IDs. Eg. the default first ID after data is filtered.
    cat('New file:', ipt$fileIn, '\n')
    if (length(ids$sel)==0) return()
    if (ids$sel[1]==''|init$but>0) return()
    if (!all(ids$sel %in% rna.fil$val)) return()
    # Avoid multiple selected rows from last input$fileIn. Must be behind gID$geneSel. 
    if (length(ids$sel)>1 & is.null(lis.url$par)) return()
    gID$geneSel <- ids$sel; gID$all <- gID$new <- NULL
    gID$new <- setdiff(gID$geneSel, gID$all); gID$all <- c(gID$all, gID$new)
    init$new <- 1 # Indicates ids are processed on-start, and no need to re-process in below.
    init$but <- 1
    cat('New ID:', gID$new, 'Selected ID:', gID$geneSel, 'All ID:', gID$all, '\n')
  })
  observeEvent(list(ids$sel), { # Selected IDs after the landing page.
    cat('Confirm selection ... \n')
    # Ensure executions only after the landing page.
    if (init$but==0|init$new==1) return()
    if (is.null(ids$but.sgl) & is.null(ids$but.mul)) return()
    if (length(ids$sel)==0|ids$sel[1]=='') return()
    gID$geneSel <- unique(ids$sel)
    gID$new <- setdiff(gID$geneSel, gID$all); gID$all <- c(gID$all, gID$new) 
    cat('New ID:', gID$new, 'Selected ID:', gID$geneSel, 'All ID:', gID$all, '\n'); cat('Done! \n')
    # if (any(is.na(gID$geneSel))) gID$geneSel <- "none"
  })

  geneV <- reactive({
    cat('All colour key values ... \n')
    validate(need(!any(is.na(gID$geneSel)) & gID$geneSel[1]!='', ''))
    # if (any(is.na(gID$geneSel))) return()
    if (is.null(geneIn())|sum(gID$geneSel[1]!='none')==0) return(NULL)
    if (input$cs.v=="Selected rows" & length(ids$sel)==0) return(NULL)
    if (ipt$fileIn!="none") { if (input$cs.v=="Selected rows") gene <- geneIn()[["df.aggr.tran"]][gID$geneSel, ]
    if (input$cs.v=="All rows") gene <- geneIn()[["df.aggr.tran"]] }
    if (!all(gID$geneSel %in% rownames(gene))) return()
    cat('Done!\n'); seq(min(gene), max(gene), len=1000) # len must be same with that from the function "spatial_hm()". Otherwise the mapping of a gene value to the colour bar is not accurate. 

  })

  col.sch <- reactive({
    cat('Color scheme ... \n') 
    if(input$color=="") return(NULL)
    col <- gsub(' |\\.|-|;|,|/', '_', input$color)
    col <- strsplit(col, '_')[[1]]
    col <- col[col!='']; col1 <- col[!col %in% colors()]
    if (length(col1>0)) validate(need(try(col1 %in% colors()), paste0('Colors not valid: ', col1, ' !')))
    cat('Done! \n'); col
  })
  
  color <- reactiveValues(col="none")
  observe({
    cat('Initial color code for color key ... \n')
    session # Avoid color$col is "none", sine new session triggers color <- reactiveValues(col="none")
    col0 <- cfg$lis.par$shm.img['color', 'default']
    col.but <- input$col.but 
    if (is.null(col.but)|is.null(col0)|gID$geneSel[1]=='none') return()
    if(col.but==0) color$col <- colorRampPalette(col_sep(col0))(length(geneV()))
    cat('Done! \n')
  })
  # As long as a button is used, observeEvent should be used. All variables inside 'observeEvent' trigger code evaluation, not only 'eventExpr'.  
  observeEvent(input$col.but, {
    cat('Customized color code for color key ... \n') 
    validate(need(col.sch(), ''))
    if (is.null(col.sch())) return (NULL)
    if (ipt$fileIn!="none") { color$col <- colorRampPalette(col.sch())(length(geneV())) }
    cat('Done! \n')
  })
  x.title <- reactiveValues(val='')
  observe({
    scale.dat <- scaleDat(); if (is.null(scale.dat)) return()
    if (!is.null(scale.dat)) if (scale.dat=='No') title <- 'No scaling' else if (scale.dat=='Row') title <- 'Scaled by row' else if (scale.dat=='Selected') title <- 'Scaled across selected genes' else if (scale.dat=='All') title <- 'Scaled across all genes' else title <- ''
    x.title$val <- title
  })
  shm.bar <- reactive({
    cat('Colour key ... \n')
    if (is.null(gID$all)) return(NULL)
    if ((any(ipt$fileIn %in% cfg$na.def) & !is.null(geneIn()))|(ipt$fileIn %in% cfg$na.cus & (!is.null(ipt$svgInpath1)|!is.null(ipt$svgInpath2)) & !is.null(geneIn()))) {
      if (length(color$col=="none")==0|input$color==""|is.null(geneV())) return(NULL)
      withProgress(message="Color scale: ", value = 0, {
        incProgress(0.75, detail="plotting, please wait ...")
        cs.g <- col_bar(geneV=geneV(), cols=color$col, width=1, x.title=x.title$val, x.title.size=10)
        cat('Done! \n'); return(cs.g)
      })
    }
  })
  # One output can only be used once in ui.R.
  output$bar1 <- bar2 <- renderPlot({ if (!is.null(shm.bar)) shm.bar() })
  # output$bar2 <- renderPlot({ if (!is.null(shm.bar)) shm.bar() })
  observe({
    ggly.but <- input$ggly.but
    if (is.null(ggly.but)) output$bar2 <- NULL else if (ggly.but==0) output$bar2 <- NULL else output$bar2 <- bar2
  }) 

  svg.na.rematch <- reactiveValues(svg.path=NULL, svg.na=NULL)
  svg.path <- reactive({
    if (ipt$fileIn=='none') return()
    if (ipt$fileIn %in% cfg$na.cus) {
      if (is.null(ipt$svgInpath2)) svgIn.df <- ipt$svgInpath1 else svgIn.df <- ipt$svgInpath2
      svg.path <- svgIn.df$datapath; svg.na <- svgIn.df$name
    } else {
      # Extract svg path and name: single or multiple svg paths are treated same way.
      lis <- svg_pa_na(cfg$svg.def[[ipt$fileIn]], cfg$pa.svg.upl)
      validate(need(try(!is.character(lis)), lis))
      svg.path <- lis$svg.path; svg.na <- lis$svg.na
    }; cat('Access aSVG path... \n')
    # If multiple svgs, check suffixes.
    lis <- svg_suffix(svg.path, svg.na)
    validate(need(try(!is.character(lis)), lis)); return(lis)
  })

  sam <- reactive({ 
    if (is.null(geneIn())) return() 
    cname <- colnames(geneIn()$df.aggr); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\1", c.na) else return(NULL) 
  })
  # renderUI: if the tab/page containing uiOutput('svg') is active clicked, the input$svg on the server side is NULL. To avoid this, the ui side should have "selectInput".
  output$svg <- renderUI({
    ns <- session$ns
    nas <- names(cfg$svg.def)
    selectInput(ns('svg'), label='Choose an aSVG to re-match', choices=nas, selected=ipt$fileIn)
  })

  # The object scope in reactive expression is restricted to the reactive evironment, thus sf.sep and sam.sep are moved outside 'observe'.
  sf.sep <- 'only_above_features_are_re-matched'
  sam.sep <- 'only_above_samples_are_re-matched'
  ft.reorder <- reactiveValues(ft.dat = NULL, ft.svg = NULL, ft.rematch = NULL)
  observeEvent(list(input$svg), {
    if (ipt$fileIn=='none'|is.null(cfg$svg.def)|is.null(input$svg)) return()
    if (!(any(ipt$fileIn %in% cfg$na.def) & is.null(ipt$svgInpath1))) return()
    # Single or multiple svg paths are treated same way.
    lis <- svg_pa_na(cfg$svg.def[[input$svg]], cfg$pa.svg.upl)
    output$msg.match <- renderText({ validate(need(try(!is.character(lis)), lis)) })
    validate(need(try(!is.character(lis)), lis))
    svg.path <- lis$svg.path; svg.na <- lis$svg.na
    cat('Access aSVG path for re-matching ... \n')
    # If multiple svgs, check suffixes.
    lis <- svg_suffix(svg.path, svg.na)
    validate(need(try(!is.character(lis)), lis))
    svg.path <- lis$svg.path; svg.na <- lis$svg.na
 
  withProgress(message="Spatial heatmap re-matching: ", value=0, {  
    incProgress(0.5, detail="extracting coordinates, please wait ...") 
    # Whether a single or multiple SVGs, all are returned in a list.
    sf.all <- NULL; for (i in seq_along(svg.na)) { 
      cat('Extract all spatial features for re-matching:', svg.na[i], '\n')
      df_tis <- svg_df(svg.path=svg.path[i], feature=sam(), cores=deter_core(2, svg.path[i]))
      validate(need(!is.character(df_tis), paste0(svg.na[i], ': ', df_tis)))
      sf.all <- c(sf.all, df_tis$tis.path)
    }
  })
  # paths and gs are dropped to bottom. 
  sf.all <- unique(sf.all); pas.idx <- grepl('^path\\d+|^g\\d+', sf.all)
  sf.all <- c(sf.all[!pas.idx], sf.all[pas.idx])
  # Matching samples are raised to top.
  sam.all <- unique(sam()); inter <- intersect(sam.all, sf.all)
  cat('Adding separators to samples and features for re-matching ... \n')
  if (length(inter)!=0) { 
    int.idx1 <- sam.all %in% inter; inter1 <- sam.all[int.idx1]
    sam.all <- c(inter1, sam.sep, sam.all[!int.idx1])
    int.idx2 <- sf.all %in% inter; sf.all <- c(inter1, sf.sep, sf.all[!int.idx2])
  } else { sf.all <- c(sf.sep, sf.all); sam.all <- c(sam.sep, sam.all) }
  ft.reorder$ft.dat <- sam.all1 <- sam.all[sam.all != sam.sep]
  ft.reorder$ft.svg <- sf.all1 <- sf.all[sf.all != sf.sep]

  cat('Preparing interface of data/aSVG features ... \n')
  output$ft.match <- renderUI({
    ns <- session$ns
      fluidRow(
      span(class = "panel panel-default", style = 'margin-left:0px',
        div(class = "panel-heading", strong("Features in aSVG")),
        div(class = "panel-body", id = ns("ftSVG"), ft2tag(sf.all1))
      ),
      div(class = "panel panel-default",
        div(class = "panel-heading", strong("Features in data")),
        lapply(sam.all1, ft_dat, ns = ns)
      ), lapply(c('ftSVG', sam.all1), ft_js, ns = ns) # Items are interchangeable across ftSVG and sam.all1.
      )
  }) 


  output$sam <- renderUI({
    rank_list(text="Drag features in data in any desired order to match right", labels=sam.all, input_id="sam")
  })

  output$sf <- renderUI({
    rank_list(text="Drag features in aSVG in any desired order to match left", labels=sf.all, input_id="sf"
    )
  }); svg.na.rematch$svg.path <- svg.path; svg.na.rematch$svg.na <- svg.na

  })

#  ft.rematch <- reactive({
  observeEvent(input$match, {
    if (is.null(ft.reorder$ft.dat) | is.null(ft.reorder$ft.svg)) return()
    ft.dat <- ft.reorder$ft.dat
    lis0 <- lapply(ft.dat, function(x) input[[x]])
    names(lis0) <- ft.dat; ft.reorder$ft.rematch <- lis0
  })  
  observeEvent(ipt$fileIn, {
    ft.reorder$ft.dat <- ft.reorder$ft.svg <- ft.reorder$ft.rematch <- NULL
    cna.match$cna <- NULL
    svg.na.rematch$svg.path <- svg.na.rematch$svg.na <- NULL
  })
  # Reactive object in "observeEvent" is not accessible outside observeEvent. The solution is eventReactive. 
  # svg.path1 stores the final svg path/na after re-matching, and will be used in SHMs.
  svg.path1 <- reactive({
    if (!is.null(svg.na.rematch$svg.path) & !is.null(svg.na.rematch$svg.na)) { svg.path <- svg.na.rematch$svg.path; svg.na <- svg.na.rematch$svg.na } else { svg.path <- svg.path()$svg.path; svg.na <- svg.path()$svg.na }
    return(list(svg.path=svg.path, svg.na=svg.na))
  })

  cna.match <- reactiveValues(cna=NULL)
  observeEvent(input$match1, {
  
    sam <- input$sam; sf <- input$sf
    cna <- colnames(geneIn()[["df.aggr.tran"]])
    w.sf <- which(sf %in% sf.sep); w.sam <- which(sam %in% sam.sep)
    # If no "renderText", the validate message is not assigned to any reactive values and thus not seen on the user interface.
    output$msg.match <- renderText({
    validate(need(try(w.sf==w.sam), 'Samples and features should be one-to-one re-matched!')); NULL })
    sf <- sf[-w.sf]; w.sf <- w.sf-1; sam <- sam[-w.sam]; w.sam <- w.sam-1
    output$msg.match <- renderText({
    validate(need(try(w.sf!=0|w.sam!=0), 'No samples or features are re-matched!')); NULL })
    if (w.sf>0) {
      cat('Re-match samples and features ... \n')
      sf.mat <- sf[seq_len(w.sf)]
      cna <- colnames(geneIn()[["df.aggr.tran"]]);
      # If sf1 is re-matched to sam1, and sf1 is included in all sams, sf1 in all sams is replaced by sam1. 
      sam.idx1 <- sam.idx2 <- list(); for (i in seq_len(w.sf)){
        if (sam[i]==sf[i]) next
        idx1 <- list(grep(paste0('^', sam[i], '__'), cna))
        names(idx1) <- sf.mat[i]; sam.idx1 <- c(sam.idx1, idx1)
        if (sf[i] %in% sam) {
          idx2 <- list(grep(paste0('^', sf[i], '__'), cna))
          names(idx2) <- sam[i]; sam.idx2 <- c(sam.idx2, idx2)
        }
      }
      # sf1 is re-matched to sam1.
      if (length(sam.idx1)>0) for (i in seq_along(sam.idx1)) {
        cna[sam.idx1[[i]]] <- sub('.*__', paste0(names(sam.idx1)[i], '__'), cna[sam.idx1[[i]]])
      }
      # sf1 in sams is replaced by sam1.
      if (length(sam.idx2)>0) for (i in seq_along(sam.idx2)) {
        cna[sam.idx2[[i]]] <- sub('.*__', paste0(names(sam.idx2)[i], '__'), cna[sam.idx2[[i]]])
      }; cna.match$cna <- cna 
    }
  })

  svg.df <- reactive({
    if (ipt$fileIn %in% cfg$na.cus & is.null(ipt$geneInpath) & ipt$dimName == "None") return()
    if ((ipt$fileIn %in% cfg$na.cus & 
    (!is.null(ipt$svgInpath1)|!is.null(ipt$svgInpath2)))|(any(ipt$fileIn %in% cfg$na.def) & is.null(ipt$svgInpath1))) {
      withProgress(message="Spatial heatmap: ", value=0, {
        incProgress(0.5, detail="extracting coordinates, please wait ...")
          svg.path <- svg.path1()$svg.path; svg.na <- svg.path1()$svg.na
          # Whether a single or multiple SVGs, all are returned in a list.
         svg.df.lis <- NULL; for (i in seq_along(svg.na)) {
            cat('Coordinate:', svg.na[i], '\n')
            df_tis <- svg_df(svg.path=svg.path[i], feature=sam(), cores=deter_core(2, svg.path[i]))
            validate(need(!is.character(df_tis), paste0(svg.na[i], ': ', df_tis)))
            svg.df.lis <- c(svg.df.lis, list(df_tis))
          }; names(svg.df.lis) <- svg.na; return(svg.df.lis)
      })
    }
  })


  observe({
    if (!is.null(ft.reorder$ft.rematch)) return()
    ipt$fileIn; geneIn(); ipt$adj.modInpath; svg.df(); input$lgdTog; input$scrollH 
    ft.path.all <- NULL; for (i in seq_along(svg.df())) { ft.path.all <- c(ft.path.all, svg.df()[[i]][['tis.path']]) }
    updateCheckboxGroupInput(session, inputId="tis", label='Select features to be transparent', choices=intersect(unique(sam()), unique(ft.path.all)), selected='', inline=TRUE)
  })
  observe({
    input$svg; svg.df <- svg.df(); ft.svg.reorder <- ft.reorder$ft.svg; input$match
    if (is.null(svg.df) | is.null(ft.svg.reorder)) return()
    ft.path.all <- NULL; for (i in seq_along(svg.df())) { ft.path.all <- c(ft.path.all, svg.df()[[i]][['tis.path']]) }
    updateCheckboxGroupInput(session, inputId="tis", label='Select features to be transparent', choices=intersect(unique(ft.svg.reorder), unique(ft.path.all)), selected='', inline=TRUE)
  })

  con <- reactive({
    if (is.null(geneIn())) return() 
    cname <- colnames(geneIn()$df.aggr); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\3", c.na) else return(NULL) 
  })

  # General selected gene/condition pattern.
  pat.con <- reactive({ con.uni <- unique(con()); if (is.null(con.uni)) return(NULL); paste0(con.uni, collapse='|') })
  pat.gen <- reactive({ if (is.null(gID$geneSel)) return(); if (gID$geneSel[1]=='none') return(NULL);  paste0(gID$geneSel, collapse='|') })
  pat.all <- reactive({ if (is.null(pat.con())|is.null(pat.gen())) return(NULL); paste0('(', pat.gen(), ')_(', pat.con(), ')') })

  # SHM ggplots, grobs legends are stored in gg.all, grob.all, lgd.all respectively for various purposes. grob.gg.all is used in relative scale of multiple SVGs, and the rescaled SVGs are stored in gg.all/grob.all finally. 
  shm <- reactiveValues(grob.all=NULL, grob.all1=NULL, gg.all=NULL, gg.all1=NULL, lgd.all=NULL, grob.gg.all = NULL)
  observeEvent(ipt$fileIn, { shm$grob.all <- shm$grob.all1 <- shm$gg.all1 <- shm$gg.all <- shm$lgd.all <- shm$grob.gg.all <- NULL })
  
  # url.id <- reactiveValues(id=NULL)
  # observe({ url.id$id <- sub(' .*', '', lis.url$par[['sear-sch.sgl']]) })
  msg.shm <- reactiveValues(msg=NULL)
  # Avoid repetitive computation under input$cs.v=='All rows'.
  gs.new <- reactive({
     cat('New grob/ggplot: ')
     # print(list(is.null(svg.df()), is.null(geneIn()), gID$new, gID$all, ids$sel, color$col[1]))
    validate(
      need(!is.null(svg.df()) & !is.null(geneIn()) & length(gID$new) > 0 & !is.null(gID$all) & length(ids$sel)>0 & color$col[1]!='none', '')
    )
    scale.shm <- input$scale.shm
    if (!is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    # If color key is build on selected rows, all SHMs should be computed upon selected rows are changed. This action is done through a separate observeEvent triggered by gID$geneSel. So in this "reactive" only one gene is accepted each time.
    # Only works at "Selected rows" and one gene is selected, i.e. when the app is launched.
    # print(list('new', ids$but.sgl, gID$geneSel, url.id$id, gID$new))
    if (length(ids$but.sgl)==0 & length(ids$but.mul)==0) return()
    if (length(url.id$sch.mul)==0|length(url.id$sch.sgl)==0) return()
    if (url.id$sch.sgl[1]!='null') urlID <- url.id$sch.sgl else if (url.id$sch.mul[1]!='null') urlID <- url.id$sch.mul else urlID <- 'null'
    if (length(url.id$id)==0) return()
    if (input$cs.v=="Selected rows" & ids$but.sgl==0) ID <- gID$geneSel else if (all(sort(url.id$id)==sort(gID$geneSel))) ID <- gID$geneSel else if (input$cs.v=="All rows") ID <- gID$new else return()
    # Works all the time as long as "All rows" selected.
    # if (input$cs.v!="All rows") return() 
    # ID <- gID$new
    if (is.null(ID)) return()
    # if (length(gID$new)>1|length(ID)>1|ID[1]=='none') return()
    if (ID[1]=='none') return()
    # Avoid repetitive computation.  
    pat.new <- paste0('^(', paste0(ID, collapse='|'), ')_(', pat.con(), ')_\\d+$')
    # print(grepl(pat.new, names(shm$grob.all)))
    if (any(grepl(pat.new, names(shm$grob.all)))) return()
    withProgress(message="Spatial heatmap: ", value=0, { 
      incProgress(0.25, detail="preparing data ...")
      gene <- geneIn()[["df.aggr.tran"]]
      # When input$fileIn updates, ID is from last session while gene is from new session.
      if (!all(ID %in% rownames(gene))) return()
      svg.df.lis <- svg.df(); 
      #w.h.all <- NULL
      # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
      # for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      svg.na <- names(svg.df.lis); grob.all <- gg.all <- lgd.all <- grob.gg.all <- NULL
      for (i in seq_along(svg.df.lis)) {

        svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]; w.h <- svg.df[['w.h']]
        tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']]
        # if (input$preScale=='Yes') mar <- (1-w.h/w.h.max*0.99)/2 else mar <- NULL
        cat(ID, ' \n')
        if (!is.null(cna.match$cna)) { 
		  if (ncol(gene)==length(cna.match$cna)) colnames(gene) <- cna.match$cna 
        }
        size.key <- as.numeric(cfg$lis.par$legend['key.size', 'default'])
        # Cores: the orders in svg.path(), names(svg.df.lis) are same.
        gg.lis <- gg_shm(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=ID, legend.col=fil.cols, cols=color$col, tis.path=tis.path, ft.trans=input$tis, sub.title.size=input$title.size * scale.shm, aspect.ratio = svg.df$aspect.r, legend.nrow=as.numeric(cfg$lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.size=input$line.size, line.color=input$line.color, lis.rematch = ft.reorder$ft.rematch) # Only gID$new is used.
        msg <- paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')
        if (is.null(gg.lis)) {
          cat(msg, '\n'); msg.shm$msg <- msg
        }
       validate(need(!is.null(gg.lis), msg)) 
       # Append suffix '_i' for the SHMs of ggplot under SVG[i], and store them in a list.
       ggs <- gg.lis$g.lis.all; names(ggs) <- paste0(names(ggs), '_', i)
       gg.all <- c(gg.all, ggs)
       # Store legend of ggplot in a list.
       lgd.all <- c(lgd.all, list(gg.lis$g.lgd))
       # Same names with ggs: append suffix '_i' for the SHMs of grob under SVG[i], and store them in a list.
       grob.lis <- grob_shm(ggs, cores=deter_core(2, svg.path()$svg.path[i])) 
       grob.all <- c(grob.all, grob.lis)
       # All ggplots/grobs are stored in nested lists under each SVG for use in relatice scale.
        lis0 <- list(grob.lis = grob.lis, gg.lis = ggs, lgd.lis = gg.lis$g.lgd)
       grob.gg.all <- c(grob.gg.all, list(lis0)) 
     }; names(lgd.all) <- names(grob.gg.all) <- svg.na
     init$new <- 0 # Terminates gs.new.
     cat('Done! \n'); return(list(gg.all = gg.all, grob.all = grob.all, lgd.all = lgd.all, grob.gg.all = grob.gg.all))
    }) # withProgress

  })
  output$msgSHM <- renderText({ 
    msg <- msg.shm$msg; validate(need(is.null(msg), msg)) 
  })
  # Extension of 'observeEvent': any of 'input$log; input$tis; input$col.but; input$cs.v' causes evaluation of all code. 
  # input$tis as an argument in "gg_shm" will not cause evaluation of all code, thus it is listed here.
  # Use "observeEvent" to replace "observe" and list events (input$log, input$tis, ...), since if the events are in "observe", every time a new gene is clicked, "input$dt_rows_selected" causes the evaluation of all code in "observe", and the evaluation is duplicated with "gs.new".
  observeEvent(col.na(), { if (col.cfm()>0) col.reorder$col.re <- 'N' })
  # Update SHMs, above theme().
  observeEvent(list(log(), input$tis, color$col, input$cs.v, col.cfm(), scaleDat(), input$match, input$line.size, input$line.color), {
    shm$grob.all <- shm$gg.all <- shm$lgd.all <- shm$grob.gg.all <- NULL; gs.all <- reactive({ 
      cat('Updating all SHMs ... \n')
      if.con <- is.null(svg.df())|is.null(geneIn())|length(ids$sel)==0|color$col[1]=='none'|gID$geneSel[1]=='none'
      if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
      scale.shm <- input$scale.shm
      if (!is.numeric(scale.shm)) return()
      if (scale.shm <= 0) return()
      withProgress(message="Spatial heatmap: ", value=0, {
      incProgress(0.25, detail="preparing data ...")
      #if (input$cs.v=="Selected rows") gene <- geneIn()[["df.aggr.tran"]][ipt.dat$dat$dt_rows_selected, ]
      #if (input$cs.v=="All rows") gene <- geneIn()[["df.aggr.tran"]]
      gene <- geneIn()[["df.aggr.tran"]][gID$geneSel, ]
      svg.df.lis <- svg.df()
      # w.h.all <- NULL
      # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
      # for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      svg.na <- names(svg.df.lis); grob.all <- gg.all <- lgd.all <- gg.grob.lis <- NULL
      for (i in seq_along(svg.df.lis)) {
        
        svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]
        tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']]
        # w.h <- svg.df[['w.h']]
        # if (input$preScale=='Yes') mar <- (1-w.h/w.h.max*0.99)/2 else mar <- NULL
        cat('All grob/ggplot:', gID$geneSel, ' \n')
        incProgress(0.75, detail=paste0('preparing ', paste0(gID$geneSel, collapse=';')))
        if (!is.null(cna.match$cna)) { 
		  if (ncol(gene)!=length(cna.match$cna)) return()
          colnames(gene) <- cna.match$cna 
        }
        size.key <- as.numeric(cfg$lis.par$legend['key.size', 'default'])
        gg.lis <- gg_shm(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=gID$geneSel, legend.col=fil.cols, cols=color$col, tis.path=tis.path, ft.trans=input$tis, sub.title.size=input$title.size * scale.shm, aspect.ratio = svg.df$aspect.r, legend.nrow=as.numeric(cfg$lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.size=input$line.size, line.color=input$line.color, lis.rematch = ft.reorder$ft.rematch) # All gene IDs are used.
        msg <- paste0(svg.na[i], ': no matching features are detected between data and aSVG!')
        if (is.null(gg.lis)) { cat(msg, '\n') }
        validate(need(!is.null(gg.lis), msg)) 
       # Append suffix '_i' for the SHMs of ggplot under SVG[i], and store them in a list.
       ggs <- gg.lis$g.lis.all; names(ggs) <- paste0(names(ggs), '_', i)
       gg.all <- c(gg.all, ggs)
       # Store legend of ggplot in a list.
       lgd.all <- c(lgd.all, list(gg.lis$g.lgd))
       # Same with ggs: append suffix '_i' for the SHMs of grob under SVG[i], and store them in a list.
       grob.lis <- grob_shm(ggs, cores=deter_core(2, svg.path()$svg.path[i]))
       grob.all <- c(grob.all, grob.lis)
       # All ggplots/grobs are stored in nested lists under each SVG for use in relatice scale.
       lis0 <- list(grob.lis = grob.lis, gg.lis = ggs, lgd.lis = gg.lis$g.lgd)
       gg.grob.lis <- c(gg.grob.lis, list(lis0))
      }; names(lgd.all) <- names(gg.grob.lis) <- svg.na
     init$new <- 0 # Terminates gs.new.
     cat('Done! \n'); return(list(grob.all = grob.all, gg.all = gg.all, lgd.all = lgd.all, gg.grob.lis = gg.grob.lis))
     }) # withProgress
    }) # reactive
    shm$grob.all <- gs.all()$grob.all; shm$gg.all <- gs.all()$gg.all
    shm$lgd.all <- gs.all()$lgd.all; shm$grob.gg.all <- gs.all()$gg.grob.lis
  }) # observeEvent
  # Avoid repetitive computation under input$cs.v=='All rows'.
  observeEvent(list(gID$geneSel), { 
    cat('Updating all SHMs caused by selected IDs ... \n')
    if.con <- is.null(input$cs.v)|gID$geneSel[1]=='none'|input$cs.v=='All rows'
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    ID <- gID$geneSel
    shm$grob.all <- shm$gg.all <- shm$lgd.all <- shm$grob.gg.all <- NULL; gs.all <- reactive({
     # print(list(ID, is.null(svg.df()), is.null(geneIn()), ids$sel, color$col[1], class(color$col[1])))
      if.con <- is.null(svg.df())|is.null(geneIn())|length(ids$sel)==0|color$col[1]=='none'
      if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
      scale.shm <- input$scale.shm
      if (!is.numeric(scale.shm)) return()
      if (scale.shm <= 0) return()
      withProgress(message="Spatial heatmap: ", value=0, {
      incProgress(0.25, detail="preparing data ...")
      gene <- geneIn()[["df.aggr.tran"]][gID$geneSel, ]
      svg.df.lis <- svg.df()
      # w.h.all <- NULL
      # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
      # for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      svg.na <- names(svg.df.lis); grob.all <- gg.all <- lgd.all <- gg.grob.lis <- NULL
      for (i in seq_along(svg.df.lis)) {

        svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]
        tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']];
        # w.h <- svg.df[['w.h']]
        # if (input$preScale=='Yes') mar <- (1-w.h/w.h.max*0.99)/2 else mar <- NULL
        cat('All grob/ggplot of row selection:', ID, ' \n')
        incProgress(0.75, detail=paste0('preparing ', paste0(ID, collapse=';')))
        if (!is.null(cna.match$cna)) { 
		  if (ncol(gene)!=length(cna.match$cna)) return()
          colnames(gene) <- cna.match$cna 
        }
        size.key <- as.numeric(cfg$lis.par$legend['key.size', 'default'])
        gg.lis <- gg_shm(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=ID, legend.col=fil.cols, cols=color$col, tis.path=tis.path, ft.trans=input$tis, sub.title.size=input$title.size * scale.shm, aspect.ratio = svg.df$aspect.r, legend.nrow=as.numeric(cfg$lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.size=input$line.size, line.color=input$line.color, lis.rematch = ft.reorder$ft.rematch) # All gene IDs are used.
        msg <- paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')
        if (is.null(gg.lis)) cat(msg, '\n')
        validate(need(!is.null(gg.lis), msg)) 
       # Append suffix '_i' for the SHMs of ggplot under SVG[i], and store them in a list.
       ggs <- gg.lis$g.lis.all; names(ggs) <- paste0(names(ggs), '_', i)
       gg.all <- c(gg.all, ggs)
       # Store legend of ggplot in a list.
       lgd.all <- c(lgd.all, list(gg.lis$g.lgd))
       # Same with ggs: append suffix '_i' for the SHMs of grob under SVG[i], and store them in a list.
       grob.lis <- grob_shm(ggs, cores=deter_core(2, svg.path()$svg.path[i]))
       grob.all <- c(grob.all, grob.lis)
       # All ggplots/grobs are stored in nested lists under each SVG for use in relatice scale.
       lis0 <- list(grob.lis = grob.lis, gg.lis = ggs, lgd.lis = gg.lis$g.lgd)
       gg.grob.lis <- c(gg.grob.lis, list(lis0))
      }; names(lgd.all) <- names(gg.grob.lis) <- svg.na
     init$new <- 0 # Terminates gs.new.
     cat('Done! \n'); return(list(grob.all = grob.all, gg.all = gg.all, lgd.all = lgd.all, gg.grob.lis = gg.grob.lis))
     }) # withProgress
    }) # reactive
    shm$grob.all <- gs.all()$grob.all; shm$gg.all <- gs.all()$gg.all
    shm$lgd.all <- gs.all()$lgd.all; shm$grob.gg.all <- gs.all()$gg.grob.lis
  }) # observeEvent
  
  # when 'color <- reactiveValues(col="none")', upon the app is launched, 'gs.new' is evaluated for 3 time. In the 1st time, 'gID$new'/'gID$all' are NULL, so 'gs.new' is NULL. In the 2nd time, 'color$col[1]=='none'' is TRUE, so NULL is returned to 'gs.new', but 'gID$new'/'gID$all' are 'HRE2'. In the third time, 'color$col[1]=='none'' is FALSE, so 'gs.new' is not NULL, but 'gID$new' is still 'HRE2', so it does not triger evaluation of 'observeEvent' and hence SHMs and legend plot are not returned upon being launched. The solution is to assign colors to 'color$col' in 'observe' upon being launched so that in the 2nd time 'gs.new' is not NULL, and no 3rd time.
  observeEvent(gs.new(), { 
    cat('Updating grobs/ggplots/legends based on new ID ... \n')
    if (is.null(svg.df())|is.null(gID$new)|length(gID$new)==0|is.null(gID$all)|is.null(gs.new())) return(NULL)
    grob.gg.lis <- gs.new()
    grobs <- grob.gg.lis[['grob.all']]
    grob.rm <- !names(shm$grob.all) %in% names(grobs)
    shm$grob.all <- c(shm$grob.all[grob.rm], grobs)
    # gs.new() becomes NULL at this step.
    # print(list(0, names(gs.new())))
    
    ggs <- grob.gg.lis[['gg.all']]
    gg.rm <- !names(shm$gg.all) %in% names(ggs)
    shm$gg.all <- c(shm$gg.all[gg.rm], ggs) 
    lgd0 <- grob.gg.lis[['lgd.all']]
    
    shm$lgd.all <- c(shm$lgd.all, lgd0[!names(lgd0) %in% names(shm$lgd.all)])
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
    if (!is.list(grob.gg.all) | !is.numeric(title.size) | is.null(svg.path()) | is.null(lay.shm()) | !is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    gg.all <- grob.all <- NULL
    for (i in seq_along(grob.gg.all)) {
      gg.lis <- grob.gg.all[[i]]$gg.lis
      # Also update the central shm$grob.gg.all
      grob.gg.all[[i]]$gg.lis <- gg.lis <- lapply(gg.lis, function(x) { x + theme(plot.title = element_text(hjust = 0.5, size = title.size * scale.shm)) })
    gg.all <- c(gg.all, gg.lis)
    # Also update the central shm$grob.gg.all
    grob.gg.all[[i]]$grob.lis <- grob.lis <- grob_shm(gg.lis, cores = deter_core(2, svg.path()$svg.path[i]))
    grob.all <- c(grob.all, grob.lis) 
    }; shm$grob.all <- grob.all; shm$gg.all <- gg.all
    shm$grob.gg.all <- grob.gg.all
    cat('Adjust title size done ... \n')
  })
 
  output$h.w.c <- renderText({
    
    if (is.null(geneIn())|length(ids$sel)==0|is.null(svg.df())|is.null(shm$grob.all)) return(NULL)

    height <- input$height; width <- input$width
    col.n <- input$col.n;
    # validate(need(height>=0.1 & !is.na(height), 'Height should be a positive numeric !'))
    # validate(need(width>=0.1 & !is.na(width), 'Width should be a positive numeric !'))
    validate(need(col.n>=1 & as.integer(col.n)==col.n & !is.na(col.n), 'No. of columns should be a positive integer !'))

  })
  # shm$lgd.all can update itself and lead to endless circles, thus it cannot be used to update the observeEvent below. In addition, when using bookmarked url, shm$lgd.all is first NULL (legend parameters are updating observeEvent below) then real ggplot object (parameters will not update oberverEvent again since they didn't change). Therefore, use lgd.par as an anchor. Only none of shm$lgd.all and legend parameters is NULL, will the observeEvent below be updated. 
  lgd.par <- reactiveValues(par=NULL)
  observe({
    # On the landing page, if the url is taken with all default parameters, after clicking the image/link, the app displays blank page, since input$lgd.key.size, input$lgd.row, input$lgd.label are all NULL, thereby not executing this step. Solution: change at least one of the parameters (e.g. horizontal layout) then take the url. 
    # print(list('adjust', is.null(shm$lgd.all), !is.numeric(input$lgd.key.size), input$lgd.row, input$lgd.label))
    if (is.null(shm$lgd.all)|!is.numeric(input$lgd.key.size)|!is.numeric(input$lgd.row)|is.null(input$lgd.label)) return()
    lgd.par$par <- list(lgd.key.size=input$lgd.key.size, lgd.row=input$lgd.row, tis=input$tis, lgd.label=input$lgd.label, lgd.lab.size=input$lgd.lab.size, lis.url=lis.url)
  })
  # lis.url is included in lgd.par$par, so it can trigger observeEvent when bookmarked url is used.
  observeEvent(lgd.par$par, {
    cat('Adjust legend size/rows/aspect ratio ... \n')
    lis.par <- lgd.par$par
    lgd.key.size <- lis.par$lgd.key.size; lgd.row <- lis.par$lgd.row
    lgd.label <- lis.par$lgd.label; label.size <- lis.par$lgd.lab.size
    if (is.null(shm$lgd.all)|!is.numeric(lgd.key.size)|!is.numeric(lgd.row)|is.null(lgd.label)) return()
    # Potential endless circles: shm$lgd.all updates itself.
    shm$lgd.all <- gg_lgd(gg.all=shm$lgd.all, size.key=lgd.key.size, size.text.key=NULL, row=lgd.row, sam.dat=sam(), ft.trans=input$tis, position.text.key='right', label=(lgd.label=='Yes'), label.size=label.size); cat('Done! \n') 
  })
  observeEvent(list(grob.all=shm$grob.all, gen.con=input$genCon), {
    cat('Reordering grobs/ggplots ... \n') 
    if (is.null(gID$all)|is.null(shm$grob.all)|is.null(shm$gg.all)) return()
    na.all <- names(shm$grob.all); pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', na.all)); if (length(con)==0) return()
    na.all <- sort_gen_con(ID.sel=gID$all, na.all=na.all, con.all=con, by=input$genCon)
    # grob1/gg.all1 are used to add/remove 2nd legend.
    shm$grob.all1 <- shm$grob.all[na.all]; shm$gg.all1 <- shm$gg.all[na.all]
    cat('Done! \n')
  })
  # Add value legend to SHMs.
  # 'observeEvent' is able to avoid infinite cycles while 'observe' may cause such cycles. E.g. in the latter, 'is.null(shm$gg.all)' and 'shm$gg.all1 <- gg.all <- gg_2lgd()' would induce each other and form infinit circles.
  observeEvent(list(val.lgd=input$val.lgd, row=input$val.lgd.row, key=input$val.lgd.key, text=input$val.lgd.text, feat=input$val.lgd.feat), {
    
    validate(need(try(as.integer(input$val.lgd.row)==input$val.lgd.row & input$val.lgd.row>0), 'Legend key rows should be a positive integer!'))
    validate(need(try(input$val.lgd.key>0), 'Legend key size should be a positive numeric!'))
    validate(need(try(input$val.lgd.text>0), 'Legend text size should be a positive numeric!'))
    
    cat('Adding value legend... \n')
    if.con <- is.null(shm$gg.all)|is.null(sam())|is.null(input$val.lgd)|is.null(input$val.lgd.feat)|input$val.lgd==0
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    gg.all <- shm$gg.all1
    if ((input$val.lgd %% 2)==1) {

      gg.all <- gg_2lgd(gg.all=gg.all, sam.dat=sam(), ft.trans=input$tis, position.2nd='bottom', legend.nrow.2nd=input$val.lgd.row, legend.key.size.2nd=input$val.lgd.key, legend.text.size.2nd=input$val.lgd.text, add.feature.2nd=(input$val.lgd.feat=='Yes'))
      shm$gg.all1 <- gg.all <- lapply(gg.all, function(x) { x+theme(legend.position="bottom") } )
      tmp <- normalizePath(tempfile(), winslash='/', mustWork=FALSE)
      png(tmp); shm$grob.all1 <- lapply(gg.all, ggplotGrob)
      dev.off(); if (file.exists(tmp)) do.call(file.remove, list(tmp))
    
    } else if ((input$val.lgd %% 2)==0) { 

      cat('Remove value legend... \n')
      shm$gg.all1 <- gg.all <- lapply(gg.all, function(x) { x+theme(legend.position="none") })
      tmp <- normalizePath(tempfile(), winslash='/', mustWork=FALSE); png(tmp); shm$grob.all1 <- lapply(gg.all, ggplotGrob) 
      dev.off(); if (file.exists(tmp)) do.call(file.remove, list(tmp))

    }

  })
  observeEvent(col.cfm(), { col.reorder$col.re <- 'Y' })
  # In "observe" and "observeEvent", if one code return (NULL), then all the following code stops. If one code changes, all the code renews.
    lay.shm <- reactive({
    cat('Spatial heatmaps layout ... \n')
    if.con <- is.null(geneIn())|length(ids$sel)==0|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
  if (col.reorder$col.re=='N') return()
    if.con <-  length(ids$sel)==0|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    if (is.na(color$col[1])|length(color$col=="none")==0|input$color=="") return(NULL)
    r.na <- rownames(geneIn()[["df.aggr.tran"]])
    grob.na <- names(shm$grob.all1)
    # Select target grobs.
    # Use definite patterns and avoid using '.*' as much as possible. Try to as specific as possible.
    pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    grob.lis.p <- shm$grob.all1[grepl(pat.all, grob.na)] # grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', names(grob.lis.p))); if (length(con)==0) return()
    lay <- input$genCon; ID <- gID$geneSel; ncol <- input$col.n
    lay <- lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, lay.mat = TRUE)
    # If 'cat' is the last step, NULL is returned.
    cat('Done! \n'); lay
  })

# shm <- shms_rela_size(input, svg.df, shm, lay.shm, svg.path)

  observeEvent(list(input$relaSize), {
    relaSize <- input$relaSize; svg.df <- svg.df(); grob.gg.all <- shm$grob.gg.all
    cat('Adjust relative plot size in multiple aSVGs ... \n')
    if (!is.numeric(relaSize) | !is.list(svg.df) | length(svg.df) == 1 | !is.list(grob.gg.all) | is.null(lay.shm()) | is.null(svg.path())) return()
    if (relaSize < 0) return()
    # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
    w.h.all <- NULL; for (i in seq_along(svg.df)) { w.h.all <- c(w.h.all, svg.df[[i]][['w.h']]); w.h.max <- max(w.h.all) }
    gg.all <- grob.all <- NULL
    for (i in seq_along(grob.gg.all)) {
      gg.lis <- grob.gg.all[[i]]$gg.lis 
      # Also update the central shm$grob.gg.all
      grob.gg.all[[i]]$gg.lis <- gg.lis <- rela_size(svg.df[[i]]$w.h['height'], w.h.max, relaSize, nrow(lay.shm()), gg.lis)
      gg.all <- c(gg.all, gg.lis)
      # Also update the central shm$grob.gg.all
      grob.gg.all[[i]]$grob.lis <- grob.lis <- grob_shm(gg.lis, cores = deter_core(2, svg.path()$svg.path[i]))
      grob.all <- c(grob.all, grob.lis) 
    }; shm$grob.all <- grob.all; shm$gg.all <- gg.all; shm$grob.gg.all <- grob.gg.all
  })
  shmLay <- reactiveValues(val=NULL) 
  # Variables in 'observe' are accessible anywhere in the same 'observe'.
  observe({
    lay <- lay.shm(); scale.shm <- input$scale.shm
    if (is.null(lay)|!is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    # subplot: height 300, width 250 
    # Avoid: if one column has all NAs in the layout matrix, the aspect ratio is distroyed. So only take the columns not containing all NAs.
    col.vld <- sum(unlist(lapply(seq_len(ncol(lay)), function(x) !all(is.na(lay[, x])))))
    # width/height relate to scrolling in box.
    shmLay$width <- width <- col.vld * 250 * scale.shm
    shmLay$height <- height <- nrow(lay) * 300 * scale.shm
    output$shm <- renderPlot(width = width, height = height, { 
    cat('Plotting spatial heatmaps ... \n')
    if (col.reorder$col.re=='N') return()
    if.con <- length(ids$sel)==0|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)

    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    if (is.na(color$col[1])|length(color$col=="none")==0|input$color=="") return(NULL)
    grob.na <- names(shm$grob.all1)
    # Select target grobs.
    # Use definite patterns and avoid using '.*' as much as possible. Try to as specific as possible.
    pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    grob.lis.p <- shm$grob.all1[grepl(pat.all, grob.na)] # grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', names(grob.lis.p))); if (length(con)==0) return()
    lay <- input$genCon; ID <- gID$geneSel; ncol <- input$col.n
    # This step is plotting.
    shmLay$val <- shm.lay <- lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, shiny=TRUE)
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
      lgd.lis <- shm$lgd.all
      lgd.lis <- gg_lgd(gg.all=lgd.lis, size.key=input$lgd.key.size*0.5, size.text.key=NULL, label.size=input$lgd.lab.size, row=input$lgd.row, sam.dat=sam(), ft.trans=input$tis, position.text.key='right', label=(input$lgd.label=='Yes'))
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
      ggsave(paste0(tmp.dir, '/shm.', input$ext), plot=shm1, device=input$ext, width=shmLay$width/72, height=shmLay$height/72, dpi=input$res, unit='in'); cat('Done! \n') 
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
  
    ipt$fileIn; geneIn(); ipt$adj.modInpath; A(); input$p; input$cv1; input$cv2; ids$sel; input$tis; input$genCon  
    url.val <- url_val('shmAll-ext', lis.url)
    updateRadioButtons(session, inputId='ext', selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['file.type', 'default']))
    url.val <- url_val('shmAll-ggly.but', lis.url)
    # updateRadioButtons(session, inputId="ggly.but", label="Show animation", choices=c("Yes", "No"), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.anm['show', 'default']), inline=TRUE)
    url.val <- url_val('shmAll-vdo.but', lis.url)
    # updateRadioButtons(session, inputId="vdo.but", label="Show/update video", choices=c("Yes", "No"), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.video['show', 'default']), inline=TRUE)

  })

  observe({
   input$vdo.key.size; input$vdo.key.row; input$vdo.val.lgd; input$tis; input$vdo.lab.size; input$vdo.res; input$vdo.itvl
   url.val <- url_val('shmAll-vdo.but', lis.url)
   # updateRadioButtons(session, inputId="vdo.but", label="Show/update video", choices=c("Yes", "No"), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.video['show', 'default']), inline=TRUE)
  })

  output$lgd1 <- lgd2 <- renderPlot(width='auto', height="auto", { # auto: no need to scroll. 
    cat('Plotting legend plot ... \n')
    lgd.row <- input$lgd.row; lgd.key.size <- input$lgd.key.size
    validate(need(try(as.integer(lgd.row)==lgd.row & lgd.row>0), ''))
    validate(need(try(lgd.key.size>0 & lgd.key.size<1), 'Legend key size should be between 0 and 1!'))
    svg.path <- svg.path1()
    if (is.null(svg.path1())|is.null(shm$lgd.all)|(length(svg.path$svg.na)>1 & is.null(input$shms.in))) return(ggplot())
      # Width and height in original SVG.
      if (length(svg.path$svg.na)>1) svg.na <- input$shms.in else svg.na <- 1
      # w.h <- svg.df()[[svg.na]][['w.h']]
      # w.h <- as.numeric(gsub("^(\\d+\\.\\d+|\\d+).*", "\\1", w.h)); r <- w.h[1]/w.h[2]
      # if (is.na(r)) return(); 
      g.lgd <- shm$lgd.all[[svg.na]]
      # g.lgd <- g.lgd+coord_fixed(ratio=r); # Aspect.ratio is fixed allready through theme(aspect.ratio). 
      cat('Done! \n'); return(g.lgd)
  })
  observe({
    ggly.but <- input$ggly.but
    if (is.null(ggly.but)) output$lgd2 <- NULL else if (ggly.but==0) output$lgd2 <- NULL else output$lgd2 <- lgd2
  }) 

  output$lgd.ui <- renderUI({ 
    ns <- session$ns    
    if (is.null(input$lgdTog)) return(NULL) 
    if (input$lgdTog %% 2 == 1) return(NULL)
    url.lgd.row <- url_val('shmAll-lgd.row', lis.url)
    url.lgd.key.size <- url_val('shmAll-lgd.key.size', lis.url)
    url.lgd.label <- url_val('shmAll-lgd.label', lis.url)
    url.lgd.lab.size <- url_val('shmAll-lgd.lab.size', lis.url)
    box(title="Legend Plot", status="primary", solidHeader=TRUE, collapsible=TRUE, width = 3, 
    navbarPage('Parameters:',
    tabPanel("Basic",
    splitLayout(cellWidths=c("32%", "1%", '32%', '1%', '35%'),
    numericInput(inputId=ns('lgd.row'), label='Key rows', value=ifelse(url.lgd.row!='null', url.lgd.row, as.numeric(cfg$lis.par$legend['key.row', 'default'])), min=1, max=Inf, step=1, width=150), '',
    numericInput(inputId=ns('lgd.key.size'), label='Key size', value=ifelse(url.lgd.key.size!='null', url.lgd.key.size, as.numeric(cfg$lis.par$legend['key.size', 'default'])), min=0, max=1, step=0.02, width=150), ''
    # numericInput(inputId=ns('lgd.ratio1'), label='Aspect ratio', value=as.numeric(cfg$lis.par$legend['aspect.ratio', 'default']), min=0.01, max=Inf, step=0.01, width=150)
    )), # tabPanel

    tabPanel("Feature labels",
    splitLayout(cellWidths=c("30%", "1%", '30%'),
    radioButtons(inputId=ns("lgd.label"), label="Feature labels", choices=c("Yes", "No"), selected=ifelse(url.lgd.label!='null', url.lgd.label, cfg$lis.par$legend['label', 'default']), inline=TRUE), '',
    numericInput(inputId=ns('lgd.lab.size'), label='Label size', value=ifelse(url.lgd.lab.size!='null', url.lgd.lab.size, as.numeric(cfg$lis.par$legend['label.size', 'default'])), min=0, max=Inf, step=0.5, width=150)
    )) # tabPanel
    ), # navbarPage
    uiOutput(ns('lgds.sel')), splitLayout(cellWidths=c("99%", "1%"), plotOutput(ns("lgd")), "")) # box

  })


  output$tran <- renderText({
    if (is.null(geneIn())|length(ids$sel)==0|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$grob.all)) return(NULL)
    if (!is.null(input$t)) validate(need(try(input$t>=0.1), 'Transition time should be at least 0.1 second!'))
  })

  observeEvent(list(fineIn=ipt$fileIn, log=log(), tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, preScale=input$preScale), { ggly_rm(); vdo_rm() })

  # Once dimension of each frame is changed, delete previous frames.
  observeEvent(list(input$scale.ly), {
    if (dir.exists('html_shm/')) { unlink('html_shm/lib', recursive=TRUE)
      file.remove(list.files('html_shm/', '*.html$', full.names=TRUE))
    } else dir.create('html_shm/')
  })

  observeEvent(list(log=log(), tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, preScale=input$preScale, ggly.but=input$ggly.but, gID.new=gID$new), {
    cat('Preparing animation frames ... \n')
    scale.ly <- input$scale.ly; ggly.but <- input$ggly.but
    if (is.null(scale.ly)|is.null(ggly.but)) return()
    if (ggly.but==0) return()
    if (is.null(geneIn())|is.null(gID$new)|length(ids$sel)==0|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$gg.all1)) return(NULL)
    if (length(color$col=="none")==0|input$color=="") return(NULL)

    withProgress(message="Animation: ", value=0, {
    incProgress(0.25, detail="preparing frames ...") 
    gg.all <- shm$gg.all1; na <- names(gg.all)
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]; gg.all <- gg.all[na]
    for (i in seq_along(gg.all)) {
      na0 <- paste0(na[i], ".html")
      if (length(list.files('www/ggly/', na0))>0) next
      # Aspect ratio is not accepted in 'ggplotly'.
      gg.all[[i]]$theme$aspect.ratio <- NULL
      gg0 <- gg.all[[i]]
      # tit.size <- gg0$theme$plot.title$size
      # This step is invalid due to the "next" above.
      # gg0$theme$plot.title$size <- tit.size*scale.ly
      gly <- ggplotly(gg0, tooltip='text') %>% layout(showlegend=FALSE)
      gly$sizingPolicy$padding <- 0
      cat('Animation: saving', na0, '\n')
      saveWidget(gly, na0, selfcontained=FALSE, libdir="lib")
      file.rename(na0, paste0('www/ggly/', na0))

    }
    if (!dir.exists('www/ggly/lib')) file.rename('lib', 'www/ggly/lib') else if (dir.exists('lib/')) unlink('lib', recursive=TRUE); cat('Done! \n')
    })

  })

  output$sld.fm <- renderUI({
    if (input$ggly.but==0) return()
    ns <- NS(id) 
    if (is.null(shm$gg.all)|is.null(pat.all())|is.null(gID$geneSel)) return(NULL) 
    gen.con.pat <- paste0('^', pat.all(), '_\\d+$') 
    sliderInput(inputId=ns('fm'), 'Frames', min=1, max=sum(grepl(gen.con.pat, names(shm$gg.all1))), step=1, value=1, animate=animationOptions(interval=input$t*10^3, loop=FALSE, playButton=icon('play'), pauseButton=icon('pause')))
  
  })

  # As long as the variable of 'reactive' is used in the 'ui.R', changes of elements in 'reactive' would cause chain change all the way to 'ui.R'. E.g. the change in "input$ggly.but=='No'" leads to changes in 'output$ggly' and 'ui.R', not necessarily changes in 'output$ggly' call changes in 'gly.url'.
  gly.url <- reactive({
    ggly.but <- input$ggly.but; fm <- input$fm 
    if (is.null(ggly.but)|is.null(fm)) return() 
    if (is.null(shm$gg.all1)|ggly.but==0|gID$geneSel[1]=='none'|is.null(pat.all())) return()
    gg.all <- shm$gg.all1; na <- names(gg.all)
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]
    if (length(na) == 0) return()
    na1 <- na[as.integer(fm)]
    asp.r <- gg.all[[na1]]$theme$aspect.ratio
    na2 <- list.files('www/ggly', pattern=na1)
    if (length(na2) == 0) return(); if (is.na(na2)) return()
    cat('Animation: access', na2, 'path \n')
    return(list(url = paste0('ggly/', na2), asp.r = asp.r))
  })
  observe({
    cfg
    updateNumericInput(session, 't', label='Transition time (s)', value=as.numeric(cfg$lis.par$shm.anm['transition', 'default']), min=0.1, max=Inf, step=0.5)
    updateNumericInput(session, 'scale.ly', label='Scale plot', value = as.numeric(cfg$lis.par$shm.anm['scale.plot', 'default']), min=0.1, max=Inf, step=0.1)
  })
  observeEvent(list(log=log(), tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, preScale=input$preScale, ggly.but=input$ggly.but, fm=input$fm), {
  
  output$ggly <- renderUI({
    scale.ly <- input$scale.ly; gly.url <- gly.url()
    ggly.but <- input$ggly.but
    cat('Animation: plotting', gly.url$url, '\n')
    if (is.null(ggly.but)) return()
    if (ggly.but==0|is.null(gly.url) | is.null(scale.ly)) return()
    if (is.null(svg.df())|is.null(geneIn())|length(ids$sel)==0|color$col[1]=='none') return(NULL)
    withProgress(message="Animation: ", value=0, {
    incProgress(0.75, detail="plotting ...")
    width <- 700*scale.ly; cat('Done! \n')
    tags$iframe(src=gly.url$url, height = width*gly.url$asp.r, width=width, scrolling='yes')  
    })
  })
  })

  anm.dld <- reactive({
    scale.ly <- input$scale.ly; gly.url <- gly.url()
    if (input$ggly.but==0|is.null(gly.url)) return()
    if (is.null(svg.df())|is.null(geneIn())|length(ids$sel)==0|color$col[1]=='none') return(NULL) 
    withProgress(message="Downloading animation: ", value=0, {
    incProgress(0.1, detail="in progress ...")
    gg.all <- shm$gg.all1; na <- names(gg.all)
    gg.na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]
    gg <- gg.all[gg.na]
    pro <- 0.1; for (i in seq_along(gg.na)) {
    incProgress(pro+0.2, detail=paste0('preparing ', gg.na[i], '.html...'))
    anm.width <- 550 * scale.ly / gly.url$asp.r
    html_ly(gg=gg[i], cs.g=shm.bar(), ft.trans=input$tis, sam.uni=sam(), anm.width = anm.width, anm.height = 550 * scale.ly, out.dir='.') }
   })
  })

  # This step leaves 'fil.na' in 'output$dld.anm' being a global variable.
  output$dld.anm <- downloadHandler( 
    # The rest code will run only after 'anm.dld()' is done.
    filename=function(){ anm.dld(); "html_shm.zip" },
    fil.na <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/html_shm.zip'),
    content=function(fil.na){ cat('Downloading animation... \n'); zip(fil.na, 'html_shm/') }
  )

  observe({
    cfg; updateSelectInput(session, "vdo.dim", label="Fixed dimension", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected=cfg$lis.par$shm.video['dimension', 'default'])
  })

  output$ffm <- renderText({
    ffm <- tryCatch({ test_ffm() }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
    if (grepl('error|warning', ffm)) paste("<span style=\"color:red\">Error: \"ffmpeg\" is required to make videos!\"</span>")
  })

  observeEvent(list(input$vdo.but), {
    cat('Making video ... \n')
    vdo.itvl <- input$vdo.itvl; vdo.res <- input$vdo.res
    vdo.but <- input$vdo.but
    if (is.null(vdo.but)|!is.numeric(vdo.itvl)|!is.numeric(vdo.res)) return(NULL)
    if (vdo.but==0|is.null(pat.all())) return(NULL)
    if (is.null(svg.df())|is.null(geneIn())|length(ids$sel)==0|color$col[1]=='none') return(NULL)
    validate(need(try(!is.na(vdo.itvl) & vdo.itvl>0), 'Transition time should be a positive numeric!'))
    validate(need(try(!is.na(vdo.res) & vdo.res>=1 & vdo.res<=700), 'Resolution should be between 1 and 700!'))
 
    withProgress(message="Video: ", value=0, {
    incProgress(0.75, detail="in progress ...")
    gg.all <- shm$gg.all1; na <- names(gg.all)
    pat <- paste0('^', pat.all(), '_\\d+$'); na <- na[grepl(pat, na)]
    gg.all1 <- gg.all[na]
    res <- vdo.res; dim <- input$vdo.dim
    if (dim %in% c('1280x800', '1280x1024', '1280x720')&res>450) res <- 450
    if (dim=='1920x1080'&res>300) res <- 300
    # selectInput("vdo.dim", label="Fixed dimension:", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected='640x480', width=110)
    vdo <- video(gg=gg.all1, cs.g=shm.bar(), sam.uni=sam(), ft.trans=input$tis, lgd.key.size=input$vdo.key.size, lgd.text.size=NULL, position.text.key='right', legend.value.vdo=(input$'vdo.val.lgd'=='Yes'), label=(input$vdo.label=='Yes'), label.size=input$vdo.lab.size, sub.title.size=8, bar.value.size=6, lgd.row=input$vdo.key.row, video.dim=dim, interval=vdo.itvl, res=res, out.dir='./www/video')
    if (is.null(vdo)) return()
    cat('Presenting video ... \n')
    incProgress(0.95, detail="Presenting video...")
    w.h <- as.numeric(strsplit(input$vdo.dim, 'x')[[1]])
    output$video <-renderUI({ tags$video(id="video", type="video/mp4", src="video/shm.mp4", width=w.h[1], height=w.h[2], controls="controls") }); cat('Done! \n')
    })

  })
    scroll.h <- reactiveValues()
    observe({ h <- input$scrollH; scroll.h$h <- ifelse(is.null(h), 450, h) })
 output$shm.ui <- renderUI({
    ns <- session$ns; if (is.null(input$togSld)) return()
    url.lgd.row <- url_val('shmAll-lgd.row', lis.url)
    url.lgd.key.size <- url_val('shmAll-lgd.key.size', lis.url)
    url.lgd.label <- url_val('shmAll-lgd.label', lis.url)
    url.lgd.lab.size <- url_val('shmAll-lgd.lab.size', lis.url)
    column(12,
 
    fluidRow(splitLayout(id='barSHM', cellWidths=c("0.5%", "6%", paste0(input$togSld*92, '%'), paste0((1-input$togSld)*92, '%')), "",  
    plotOutput(ns("bar1")),
    if (input$togSld!=0) div(id='divSHM', style=paste0('overflow-y:scroll;height:', scroll.h$h, 'px;overflow-x:scroll'), plotOutput(ns("shm"), height='100%', width='100%')),

    if (input$togSld!=1) navbarPage('',
    tabPanel('Legend', list(uiOutput(ns('lgds.sel')), plotOutput(ns("lgd1")))),
    tabPanel("Parameters",
    splitLayout(cellWidths=c("32%", "1%", '32%', '1%', '35%'),
    numericInput(inputId=ns('lgd.row'), label='Key rows', value=ifelse(url.lgd.row!='null', url.lgd.row, as.numeric(cfg$lis.par$legend['key.row', 'default'])), min=1, max=Inf, step=1, width=150), '',
    numericInput(inputId=ns('lgd.key.size'), label='Key size', value=ifelse(url.lgd.key.size!='null', url.lgd.key.size, as.numeric(cfg$lis.par$legend['key.size', 'default'])), min=0, max=1, step=0.02, width=150), ''
    ),
    splitLayout(cellWidths=c("30%", "1%", '30%'),
    radioButtons(inputId=ns("lgd.label"), label="Feature labels", choices=c("Yes", "No"), selected=ifelse(url.lgd.label!='null', url.lgd.label, cfg$lis.par$legend['label', 'default']), inline=TRUE), '',
    numericInput(inputId=ns('lgd.lab.size'), label='Label size', value=ifelse(url.lgd.lab.size!='null', url.lgd.lab.size, as.numeric(cfg$lis.par$legend['label.size', 'default'])), min=0, max=Inf, step=0.5, width=150)
    )) # tabPanel
    ) # navbarPage
  )) # splitLayout(cellWidths
  ) # column
  })

# addPopover(session=session, id="height", title="", content="Check 'Yes' to preserve the aspect ratio defined in the aSVG file.", placement = "bottom", trigger = "hover", options = NULL)  
  
 output$lgds.sel <- renderUI({
    ns <- session$ns 
    if (is.null(svg.path1())) return(NULL)
    if (length(svg.path1()$svg.na)==1) return(NULL)
    svg.pa <- svg.path1()[['svg.na']]
    selectInput(ns('shms.in'), label='Select plots', choices=svg.pa, selected=svg.pa[1])
  })
  # If tab.act.lis is defined inside observe, then it is not accessible outside observe.
  tab.act.lis <- reactiveValues()
  observe({
    tab.act.lis$shmMhNet <- input$shmMhNet
  })
  observe({
    shmMhNet <- input$shmMhNet; interNav <- input$interNav
    if (is.null(shmMhNet)|is.null(interNav)) return()
    tab.inter <- ifelse(shmMhNet=='interTab' & interNav=='interPlot', 'yes', 'no')
    if (input$ggly.but==0 & tab.inter=='yes') showModal(modal(msg=HTML('To see the latest interactive image, always click the button <strong>"Click to show/update"</strong>!'), easyClose=TRUE))
  })
  observe({
    shmMhNet <- input$shmMhNet; vdoNav <- input$vdoNav
    if (is.null(shmMhNet)|is.null(vdoNav)) return()
    tab.vdo <- ifelse(shmMhNet=='vdoTab' & vdoNav=='video', 'yes', 'no')
    if (input$vdo.but==0 & tab.vdo=='yes') showModal(modal(msg=HTML('To see the latest video, always click the button <strong>"Click to show/update"</strong>!'), easyClose=TRUE))
  })
  network_server('net', ipt, cfg, dat.mod.lis, shm.mod.lis=list(gID=gID, tab.act.lis=tab.act.lis), sch.mod.lis)

  observe({
    ipt$fileIn; ipt$geneInpath; lis.par <- cfg$lis.par
    url.val <- url_val('shmAll-cs.v', lis.url)
    updateRadioButtons(session, inputId='cs.v', label='Color key based on', choices=c("Selected rows", "All rows"), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['color.scale', 'default']), inline=TRUE)
    url.val <- url_val('shmAll-col.n', lis.url)
    updateSliderInput(session, inputId='col.n', label='', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['columns', 'default'])), min=1, max=50, step=1)
    url.val <- url_val('shmAll-genCon', lis.url)
    updateRadioButtons(session, inputId="genCon", label="", choices = c("Gene"="gene", "Condition"="con"), selected = ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['display.by', 'default']), inline=FALSE)
  # addPopover(session, "genCon", title="Data column: by the column order in data matrix.", placement="bottom", trigger='hover')
    url.val <- url_val('shmAll-scale.shm', lis.url)
  updateSliderInput(session, inputId='scale.shm', label='', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['scale.plots', 'default'])), min=0.1, max=10, step=0.1)
    url.val <- url_val('shmAll-title.size', lis.url)
  updateSliderInput(session, inputId='title.size', label='', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['title.size', 'default'])), min=0, max=100, step=0.5)
    url.val <- url_val('shmAll-color', lis.url)
  updateTextInput(session, "color", "Color scheme", ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['color', 'default']), placeholder=paste0('Eg: ', cfg$lis.par$shm.img['color', 'default']))
  url.val <- url_val('shmAll-cs.v', lis.url) 
  updateRadioButtons(session, inputId='cs.v', label='Color key based on', choices=c("Selected rows", "All rows"), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['color.scale', 'default']), inline=TRUE)
  url.val <- url_val('shmAll-val.lgd.row', lis.url)
  updateNumericInput(session, inputId='val.lgd.row', label='Rows', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['value.legend.rows', 'default'])), min=1, max=Inf, step=1)
  url.val <- url_val('shmAll-val.lgd.key', lis.url)
  updateNumericInput(session, inputId='val.lgd.key', label='Key size', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['value.legend.key', 'default'])), min=0.0001, max=1, step=0.01)
  url.val <- url_val('shmAll-val.lgd.text', lis.url)
  updateNumericInput(session, inputId='val.lgd.text', label='Text size', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['value.legend.text', 'default'])), min=0.0001, max=Inf, step=1)
  url.val <- url_val('shmAll-val.lgd.feat', lis.url)
  updateRadioButtons(session, inputId='val.lgd.feat', label='Include features', choices=c('No', 'Yes'), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['include.feature', 'default']), inline=TRUE)
  url.val <- url_val('shmAll-line.color', lis.url)
  updateSelectInput(session, 'line.color', label='Line color', choices=c('grey70', 'black', 'red', 'green', 'blue'), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['line.color', 'default']))
  url.val <- url_val('shmAll-line.size', lis.url)
  updateNumericInput(session, inputId='line.size', label='Line size', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['line.size', 'default'])), min=0.05, max=Inf, step=0.05) 
  url.val <- url_val('shmAll-ext', lis.url)
  updateRadioButtons(session, inputId='ext', selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['file.type', 'default']))
  url.val <- url_val('shmAll-res', lis.url)
  updateNumericInput(session, inputId='res', label='Resolution (dpi)', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['dpi', 'default'])), min=10, max=Inf, step=10)
  url.val <- url_val('shmAll-lgd.incld', lis.url)
  updateRadioButtons(session, inputId='lgd.incld', label='Include legend plot', choices=c('Yes', 'No'), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.img['include.legend.plot', 'default']), inline=TRUE)
  url.val <- url_val('shmAll-lgd.size', lis.url) 
  updateNumericInput(session, inputId='lgd.size', label='Legend plot size', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['legend.plot.size', 'default'])), min=-1, max=Inf, step=0.1)
  url.val <- url_val('shmAll-relaSize', lis.url)
  updateNumericInput(session, inputId='relaSize', label='Relative sizes', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.img['relative.size', 'default'])))
  url.val <- url_val('shmAll-vdo.key.row', lis.url)
  updateNumericInput(session, inputId='vdo.key.row', label='Key rows', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.video['key.rows', 'default'])), min=1, max=Inf, step=1)
  url.val <- url_val('shmAll-vdo.key.size', lis.url)
  updateNumericInput(session, inputId='vdo.key.size', label='Key size', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.video['key.size', 'default'])), min=0.01, max=Inf, step=0.1)
  url.val <- url_val('shmAll-vdo.val.lgd', lis.url)
  updateRadioButtons(session, inputId="vdo.val.lgd", label="Key value", choices=c("Yes", "No"), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.video['value.legend', 'default']), inline=TRUE)
  url.val <- url_val('shmAll-vdo.label', lis.url)
  updateRadioButtons(session, inputId="vdo.label", label="Feature label", choices=c("Yes", "No"), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.video['feature.label', 'default']), inline=TRUE)
  url.val <- url_val('shmAll-vdo.lab.size', lis.url)
  updateNumericInput(session, inputId='vdo.lab.size', label='Label size', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.video['label.size', 'default'])), min=0, max=Inf, step=0.5)
  url.val <- url_val('shmAll-vdo.itvl', lis.url)
  updateNumericInput(session, inputId='vdo.itvl', label='Transition time (s)', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.video['transition', 'default'])), min=0.1, max=Inf, step=1)
  url.val <- url_val('shmAll-vdo.res', lis.url)
  updateNumericInput(session, inputId='vdo.res', label='Resolution (dpi)', value=ifelse(url.val!='null', url.val, as.numeric(cfg$lis.par$shm.video['dpi', 'default'])), min=1, max=1000, step=5)
  url.val <- url_val('shmAll-vdo.but', lis.url)
  # updateRadioButtons(session, inputId="vdo.but", label="Show/update video", choices=c("Yes", "No"), selected=ifelse(url.val!='null', url.val, cfg$lis.par$shm.video['show', 'default']), inline=TRUE)

  })
  onBookmark(function(state) { state })
  return(list(gID=gID))
})} # shm_server


mods$shm <- shm.mod.lis <- shm_server('shmAll', ipt, cfg, sch, lis.url, url.id, dat.mod.lis, sch.mod.lis)


deg_server <- function(id, ipt, cfg, sch, lis.url, ids, dat.mod.lis, shm.mod.lis, session) {

  moduleServer(id, function(input, output, session) {

  gID <- shm.mod.lis$gID
  cat('Presenting data matrix (DEG) ... \n')
  dat.deg.mod.lis <- data_server('datDEG', ipt, cfg, sch, lis.url, ids, deg=TRUE)
  cat('Done! \n')
  geneIn1 <- dat.mod.lis$geneIn1 # Take transformed values in DEG table.
  geneIn <- dat.deg.mod.lis$geneIn # Take filted matrix with replicates. 
  ssg.dat <- reactive({
    # if (input$fileIn!="customData") return(NULL)
    cat('DEG: SummarizedExperiment, features, factors ... \n')
    if (is.null(geneIn())) return(NULL)
    gen.lis <- geneIn()
    df.rep <- as.matrix(gen.lis[['df.rep']])
    int <- all(round(df.rep) == df.rep)
    if (!int) showModal(modal(msg = strong('Only count matrix is accepted!'))); validate(need(int, ''))
    rows <- nrow(df.rep) >= 50
    if (!rows) showModal(modal(msg = strong('Make sure count matrix includes at least 50 genes!'))); validate(need(rows, ''))
    cna <- colnames(df.rep)
    sams <- gsub('^(.*)__(.*)$', '\\1', cna)
    cons <- gsub('^(.*)__(.*)$', '\\2', cna)
    validate(need(try(length(unique(sams))>1 & length(unique(cons))>1), 'At least 2 samples and 2 conditions are required!')) 
    se <- SummarizedExperiment(assays=list(expr=df.rep), colData=data.frame(sample=sams, condition=cons, samCon=cna, rep=cna)); cat('Done! \n')
    return(list(se=se, sams=sams, cons=cons))
  })

  output$ssg.sam <- renderUI({
    ns <- session$ns
    ssg.dat <- ssg.dat(); if (is.null(ssg.dat)) return()
    cho <- c('all', unique(ssg.dat$sams))
    selectInput(ns('ssg.sam'), label='Select features', choices=cho, selected=cho[2:3], multiple=TRUE)
  })
  output$ssg.con <- renderUI({
    ns <- session$ns
    ssg.dat <- ssg.dat(); if (is.null(ssg.dat)) return()
    cho <- c('all', unique(ssg.dat$cons))
    selectInput(ns('ssg.con'), label='Select factors', choices=cho, selected=cho[2:3], multiple=TRUE)
  })

  # Subset SE.
  sub_se <- function(se, sams, cons) {
    cna.sel <- NULL; for (i in cons) { 
      cna.sel <- c(cna.sel, paste0(sams, '__', i)) }
    se <- se[, colnames(se) %in% cna.sel]; return(se)
  }

  se <- reactive({
    cat('Subsetting SE with input features/factors ... \n')
    # if (input$fileIn!="customData") return(NULL)
    dat <- ssg.dat(); sam.con <- input$sam.con 
    if (is.null(geneIn())|is.null(dat)|is.null(sam.con)) return(NULL)
    sam <- input$ssg.sam; con <- input$ssg.con; se <- dat$se
    if (is.null(sam)|is.null(con)) return() 
    if ('all' %in% sam) sam <- unique(dat$sams)
    if ('all' %in% con) con <- unique(dat$cons)
    se <- sub_se(se, sam, con)
    if (sam.con=='feature') fct <- 'sample' else if (sam.con=='factor') fct <- 'condition' else if (sam.con=='feature__factor') fct <- 'samCon'
    fct.tab <- table(colData(se)[, fct])
    fct.na <- names(fct.tab)[fct.tab==1] 
    validate(need(length(fct.na)==0, paste0('At least 2 replicates are required: ', paste0(fct.na, collapse=', '))))
    cat('Done! \n'); return(list(se=se, fct=fct))
  })

  se.nor.rok.dis <- reactive({
    cat('Normalizing data for ROKU/distinct ... \n')
    # if (input$fileIn!="customData") return(NULL)
    dat <- ssg.dat() 
    if (is.null(geneIn())|is.null(dat)) return(NULL)
    nor.meth <- input$rok.dis.nor
    if (grepl('^CNF-', nor.meth)) {
      norm.fun <- 'CNF'
      par.lis <- list(method=sub('.*-', '', nor.meth))
    } else if (nor.meth=='none') return(list(se.nor=dat$se, fct=se()$fct)) else { norm.fun <- nor.meth; par.lis <- NULL }
    se.nor <- norm_data(data=dat$se, norm.fun=norm.fun, parameter.list=par.lis, log2.trans=TRUE); cat('Done! \n')
    return(list(se.nor=se.nor))
  })

  output$ssg.tar <- renderUI({
    ns <- session$ns
    # if (input$fileIn!="customData") return(NULL)
    if (is.null(se())) return()
    lis <- se(); se <- lis$se; fct <- lis$fct
    if ('ROKU' %in% input$ssg.meth) {
      validate(need(try(length(unique(colData(se)[, fct]))>2), 'The number of targets (samples/conditions) to compare should be at least 3 for ROKU!')) 
    }
    selectInput(ns('ssg.tar'), label='Target', choices=unique(colData(se)[, fct]), selected=NULL, multiple=FALSE)
  }) 
  # Pairwise comparison coefficients.
  output$dt.vs1 <- output$dt.vs2 <- output$dt.vs3 <- renderDataTable({
    cat('Pairwise comparison coefficients ... \n')
    dat <- ssg.dat(); sam <- input$ssg.sam; con <- input$ssg.con
    sam.con <- input$sam.con; tar <- input$ssg.tar
    if (is.null(dat)|!is.character(sam)|!is.character(con)|!is.character(sam.con)|!is.character(tar)) return()
    if ('all' %in% sam) sam <- unique(dat$sams)
    if ('all' %in% con) con <- unique(dat$cons)
    vs <- data.frame()
    # save(dat, sam, con, sam.con, tar, file='dscstf')
    # Reference, query.
    if (sam.con %in% c('feature', 'factor')) {
      if (sam.con == 'feature') {  
        under <- con; ref <- setdiff(sam, tar)
        validate(need(length(ref)>0, 'If compare by "feature", at least 2 features should be selected!'))
      } else if (sam.con == 'factor') {
        under <- sam; ref <- setdiff(con, tar)
        validate(need(length(ref)>0, 'If compare by "factor", at least 2 factors should be selected!'))
      }
      tar.all <- paste0(tar, '__', under)
      ref.all <- paste0(unlist(lapply(ref, function(i) paste0(i, '__', under))))
      vs <- rbind(vs, c(tar.all, 'VS', ref.all))
      colnames(vs) <- c(paste0('target', seq_along(tar.all)), 'VS', paste0('reference', seq_along(ref.all))) 
    } else if (sam.con == 'feature__factor') {
      coms <- unlist(lapply(sam, function(i) {paste0(i, '__', con)} ))
      ref <- setdiff(coms, tar)
      vs <- rbind(vs, c(tar, 'VS', ref))
      colnames(vs) <- c('target', 'VS', paste0('reference', seq_along(ref))) 
    }
    d.tab <- datatable(vs, selection='none', extensions='Scroller', plugins = "ellipsis", class='cell-border strip hover', options = list(dom = 't', scrollX = TRUE)); cat('Done! \n')
    d.tab
  })
  edg0 <- reactive({
    cat('edgeR all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (is.null(geneIn())|is.null(se())) return(NULL)
    withProgress(message="edgeR: ", value=0, {
    incProgress(0.5, detail="in progress ...")
    edg <- edgeR(se=se()$se, method.norm=sub('.*-', '', input$edg.lim.nor), com.factor=se()$fct, method.adjust='BH', return.all=TRUE); cat('Done! \n'); edg
    })
  })
  edg <- eventReactive(input$ssg.update, {
    cat('edgeR log2/fc ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (!'edgeR' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(edg0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct]) 
    up.dn <- up_dn(sam.all=lvl, df.all=edg0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='logFC', fdr.na='FDR'); cat('Done! \n')
    up.dn
  })

  dsq0 <- reactive({
    cat('DESeq2 all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (is.null(geneIn())|is.null(se())) return(NULL)
    withProgress(message="DESeq2: ", value=0, {
    incProgress(0.5, detail="in progress ...")
    dsq <- deseq2(se=se()$se, com.factor=se()$fct, method.adjust='BH', return.all=TRUE); cat('Done! \n')
    dsq
    })
  })
  dsq <- eventReactive(input$ssg.update, {
    cat('DESeq2 log2/fc ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (!'DESeq2' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(dsq0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct])
    up.dn <- up_dn(sam.all=lvl, df.all=dsq0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='log2FoldChange', fdr.na='padj'); cat('Done! \n'); up.dn
  })

  lim0 <- reactive({
    cat('limma all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (is.null(geneIn())|is.null(se())) return(NULL)
    withProgress(message="limma: ", value=0, {
    incProgress(0.5, detail="in progress ...") 
    lim <- limma(se=se()$se, method.norm=sub('.*-', '', input$edg.lim.nor), m.array=FALSE, com.factor=se()$fct, method.adjust='BH', return.all=TRUE); cat('Done! \n'); lim
    })
  })
  lim <- eventReactive(input$ssg.update, {
    cat('limma log2/fc ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (!'limma' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(lim0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct])
    up.dn <- up_dn(sam.all=lvl, df.all=lim0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='logFC', fdr.na='adj.P.Val'); cat('Done! \n'); up.dn
  })

  rok0 <- reactive({
    cat('ROKU all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    lis <- se.nor.rok.dis()
    if (is.null(geneIn())|is.null(se())|is.null(lis)) return(NULL)
    withProgress(message="ROKU: ", value=0, {
      incProgress(0.5, detail="in progress ...")
      # Column names and factors are same between se and se.nor.
      se.nor <- lis$se.nor[, colnames(se()$se)] 
      fct <- se()$fct
      validate(need(try(length(unique(colData(se.nor)[, fct]))>2), 'The number of targets (samples/conditions) to compare should be at least 3 for ROKU!')) 
      rok <- roku(data=se.nor, com.factor=fct, aggr='mean', n=1, log2.trans=FALSE, return.all=TRUE); cat('Done! \n'); rok
    })
  })
  rok <- eventReactive(input$ssg.update, {
    # if (input$fileIn!="customData") return(NULL)
    if (!'ROKU' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(rok0())) return(NULL)
    cat('ROKU up/down ... \n')
    up_dn_roku(rok0())
  })

  dis0 <- reactive({
    cat('distinct all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    lis <- se.nor.rok.dis()
    if (is.null(geneIn())|is.null(se())|is.null(lis)) return(NULL)
    withProgress(message="distinct: ", value=0, {
      incProgress(0.5, detail="in progress, this method takes longer time ...")
      se.nor <- lis$se.nor[, colnames(se()$se)] 
      fct <- se()$fct
      dis <- distt(se.nor, fct, return.all=TRUE); cat('Done! \n')
      dis
    })
  })
  dis <- eventReactive(input$ssg.update, {
    cat('distinct log2/fc ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (!'distinct' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(dis0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct])
    up.dn <- up_dn(sam.all=lvl, df.all=dis0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='log2FC', fdr.na='FDR'); cat('Done! \n'); up.dn
  })
  
  deg.lis <- reactive({
    cat('DEG list across methods ... \n')
    if (is.null(geneIn())|is.null(input$ssg.tar)) return(NULL)
    meth.sel <- c(!is.null(edg()), !is.null(lim()), !is.null(dsq()), !is.null(rok()), !is.null(dis()))
    validate(need(try(sum(meth.sel) >= 1), 'At least 1 method should be selected!'))
    lis <- list(edgeR=edg(), limma=lim(), DESeq2=dsq(), ROKU=rok(), distinct=dis())[meth.sel]
    deg.lis1 <- deg_lis(lis, sam=input$ssg.tar, 'up')
    deg.lis2 <- deg_lis(lis, sam=input$ssg.tar, 'down')
    if (length(deg.lis1) == 1) { # Only one method is selected.
      deg.lis1 <- c(deg.lis1[1], deg.lis1[1])
      names(deg.lis1) <- paste0(names(deg.lis1), c('', '.1'))
    }
    if (length(deg.lis2) == 1) { # Only one method is selected.
      deg.lis2 <- c(deg.lis2[1], deg.lis2[1])
      names(deg.lis2) <- paste0(names(deg.lis2), c('', '.1'))
    }
    cat('Done! \n'); return(list(up.lis = deg.lis1, down.lis = deg.lis2))
  })
  output$upset1 <- renderPlot({
    cat('Upset up ... \n'); deg.lis <- deg.lis()
    if (is.null(geneIn())|is.null(input$ssg.tar) | is.null(deg.lis)) return(NULL)
    up.lis <- deg.lis$up.lis
    validate(need(!is.null(up.lis), 'No over-expressed genes are detected!'))
    # up <- upset(fromList(up.lis), order.by="degree", nintersects=40, point.size=3, line.size=1, mb.ratio=c(0.6, 0.4), text.scale=1.5); cat('Done! \n'); up
    up <- deg_ovl(deg.lis, type='up', plot='upset'); cat('Done! \n'); up
  })
  output$upset2 <- renderPlot({
    cat('Upset down ... \n'); deg.lis <- deg.lis()
    if (is.null(geneIn())|is.null(input$ssg.tar) | is.null(deg.lis)) return(NULL)
    down.lis <- deg.lis$down.lis
    validate(need(!is.null(down.lis), 'No under-expressed genes are detected!'))
    dn <- deg_ovl(deg.lis, type='down', plot='upset'); cat('Done! \n'); dn
    # dn <- upset(fromList(down.lis), order.by="degree", nintersects=40, point.size=3, line.size=1, mb.ratio=c(0.6, 0.4), text.scale=1.5); cat('Done! \n'); up
  })


  output$ovl1 <- renderPlot({
    cat('Overlap up DEG ... \n'); degLis <- deg.lis()
    up.lis <- degLis$up.lis
    if (is.null(geneIn())|is.null(input$ssg.tar) | is.null(up.lis)) return(NULL)
    g <- deg_ovl(degLis, type='up', plot='matrix'); cat('Done! \n'); g
  })

  output$ovl2 <- renderPlot({
    cat('Overlap down DEG ... \n'); degLis <- deg.lis()
    down.lis <- degLis$down.lis
    if (is.null(geneIn())|is.null(input$ssg.tar) | is.null(down.lis)) return(NULL)
    g <- deg_ovl(degLis, type='down', plot='matrix'); cat('Done! \n'); g
  })

  dt.deg <- reactive({
    cat('DEG data frame ... \n')
    geneIn1 <- geneIn1()
    if (is.null(geneIn())|is.null(geneIn1)|is.null(input$ssg.tar)) return(NULL)
    meth.sel <- c(!is.null(edg()), !is.null(lim()), !is.null(dsq()), !is.null(rok()), !is.null(dis()))
    validate(need(try(sum(meth.sel) >= 1), 'At least 1 method should be selected!'))
    lis <- list(edgeR=edg(), limma=lim(), DESeq2=dsq(), ROKU=rok(), distinct=dis())[meth.sel]
    deg.lis1 <- deg_lis(lis, sam=input$ssg.tar, 'up')
    deg.lis2 <- deg_lis(lis, sam=input$ssg.tar, 'down')
    # venn_inter returns matrix, since matrix accepts duplicate row names while data frame not. Some genes might be up in one method while down in other methods, so there could be duplicated rows in df.deg.
    df.deg <- rbind(venn_inter(deg.lis1), venn_inter(deg.lis2))
    validate(need(nrow(df.deg)>0, 'No up or down regulated genes are detected!'))
    df.met <- geneIn1[["df.met"]][, , drop=FALSE]
    # The data (geneIn1) before filtering is used, so even though all genes are filtered out, the DEG SHMs can still be displayed.
    na.deg <- rownames(df.deg); df.aggr.tran <- round(geneIn1$df.aggr.tran, 2)
    # Switch data sets.
    if (!all(na.deg %in% rownames(df.aggr.tran))) return()
    if (ncol(df.met) > 0) df.deg <- cbind.data.frame(df.met[na.deg, , drop = FALSE], df.deg, stringsAsFactors=FALSE)
    df.deg <- cbind.data.frame(df.deg, df.aggr.tran[na.deg, , drop = FALSE], stringsAsFactors=FALSE)
    cat('Done! \n'); df.deg 
  })
    output$dld.ssg.tab <- downloadHandler(
      filename=function(){ "tissue_specific_genes.txt" },  content=function(file=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/tissue_specific_genes.txt')){
      write.table(dt.deg(), file, sep='\t', col.names=TRUE, row.names=TRUE) }
    )
    
  output$dt.deg <- renderDataTable({
    cat('DEG summary table ... \n'); gene.dt <- dt.deg()
    if (is.null(gene.dt)) return()
    col1 <- list(list(targets=c(1), render=JS("$.fn.dataTable.render.ellipsis(5, false)")))
    # In case no metadata column.
    if (colnames(gene.dt)[1]!='metadata') col1 <- NULL
    d.tab <- datatable(gene.dt, selection='none', escape=FALSE, filter="top", extensions='Scroller', plugins = "ellipsis",
      options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE, searchHighlight = TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs=col1), 
      class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer'); cat('Done! \n')
      d.tab
  }) 
  dat.mod.lis.deg <- reactiveValues()
  observe({
    cat('Preparing search box in DEG section ... \n')
    dat.deg <- dt.deg(); if (is.null(dat.deg)) return()
    dat.lis <- dat.mod.lis; if (is.null(dat.lis)) return()
    deg.rna <- rownames(dat.deg); gen.lis <- dat.lis$geneIn()
    dat.mod.lis.deg$geneIn <- reactive({  
    return(list(df.aggr=gen.lis$df.aggr[deg.rna, , drop=FALSE], df.aggr.tran=gen.lis$df.aggr.tran[deg.rna, , drop=FALSE], df.aggr.tran.order=gen.lis$df.aggr.tran.order[deg.rna, , drop=FALSE], df.met=gen.lis$df.met[deg.rna, , drop=FALSE], df.rep=gen.lis$df.rep[deg.rna, , drop=FALSE]))
    })
   sch.mod.lis <- search_server('deg', ids, cfg, lis.url, url.id, dat.mod.lis.deg)
   cat('Done! \n')
  })

  # gID is a reactiveValue. It is shared between modules, so no need to pass the updated gID to shm_server.
  observeEvent(list(input$dt.deg_rows_selected), {
    cat('Selected genes in DEG summary table ... \n')
    if (is.null(input$dt.deg_rows_selected)) return()
    r.na <- rownames(dt.deg()); gID$geneSel <- unique(r.na[input$dt.deg_rows_selected])
    if (any(is.na(gID$geneSel))) gID$geneSel <- "none"
    gID$new <- setdiff(gID$geneSel, gID$all); gID$all <- c(gID$all, gID$new); cat('Done! \n')
  })

  output$w.table1 <- renderDataTable({
    
    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
    if (is.null(geneIn0())) return(NULL)
    sam.con <- unique(colnames(geneIn()))[1:3]
    lis.all <- list(edg=edg(), dsq=dsq()); sam.all <- names(lis.all[[1]])
    sam <- sam.all[1:3]; sam.tar <- sam.all[2]; sam.vs <- sam.all[c(1, 3)]
    meth <- 'edg'
    df.up <- lis.all[[meth]][[sam[2]]][[1]]
    df.dn <- lis.all[[meth]][[sam[2]]][[2]]
    df.up.dn <- rbind(df.up, df.dn)
    if (nrow(df.up.dn)==0) return("No SSGs detected.")
    # print(df.up.dn)
    pat <- c(paste0('^', sam.tar, '_VS_', sam.vs, '_'), paste0('^', sam.vs, '_VS_', sam.tar, '_'))
    df.up.dn <- df.up.dn[, c(1, grep(paste0(pat, collapse='|'), colnames(df.up.dn)))]
    datatable(df.up.dn, selection=list(mode="multiple", target="row", selected=NULL), filter="top", extensions='Scroller', options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer')
 
  }) 

  output$ssg.sum <- renderPlot({ 

    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
    # if (input$fileIn!="customData") return(NULL)
    if (is.null(geneIn0())) return(NULL)
    per=0.5; width=0.85
    lis.all <- list(edg(), dsq()); sam.all <- names(lis.all[[1]])
    sam <- sam.all[1:3]
    lis.aggr <- sig_frq(lis.all=lis.all, per=per, sam=sam)
    ssg_sum(lis.aggr=lis.aggr, width=width)

  })


  output$a.table <- renderDataTable({
    
    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
    if (is.null(geneIn0())) return(NULL)
    per=0.5
    sam.con <- unique(colnames(geneIn()))[1:3]
    lis.all <- list(edg=edg(), dsq=dsq()); sam.all <- names(lis.all[[1]])
    lis.aggr.r <- sig_frq(lis.all=lis.all, per=per)
    df.sum <- as.data.frame(lis.aggr.r[[1]][[1]][[2]])
    datatable(df.sum, selection=list(mode="multiple", target="row", selected=NULL), filter="top", extensions='Scroller', options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer')

  })

    expr.nor <- reactive({

    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
      if (is.null(geneIn0())) return(NULL)
      df.rep <- geneIn0()[['df.rep']]
      se <- SummarizedExperiment(assays=list(expr=as.matrix(df.rep)), colData=data.frame(fct=colnames(df.rep)))
      norm_aggr(se=ssg.dat()$se, method.norm='TMM', log2.trans=TRUE, sample.factor='fct', rep.aggr='mean')

    })
  onBookmark(function(state) { state })

})}

  deg.mod.lis <- deg_server('deg', ipt, cfg, sch, lis.url, ids, dat.mod.lis, shm.mod.lis)
  setBookmarkExclude(c("dat-dtSel_rows_all", "dat-dtSel_rows_current", "dat-dtSel_search_columns", "dat-dtAll_rows_all", "dat-dtAll_rows_current", "dat-dtAll_search_columns")) 
  observe({
    lis.ipt <- reactiveValuesToList(input); session$doBookmark()
    # lapply(seq_along(lis.ipt), function(i) {if (length(lis.ipt[[i]])<1000) { print(lis.ipt[i]) }})
    # print(length(lis.ipt[['shmAll-dat-dt_rows_all']]))
    # print(getUrlHash()); print(getQueryString())
  })
 onBookmarked(function(url) { updateQueryString(url) })
# onBookmarked(updateQueryString)

  onStop(function() { ggly_rm(); vdo_rm(); cat("Session stopped! \n") })

}



