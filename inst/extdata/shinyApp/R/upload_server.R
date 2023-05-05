# Module for uploading data.

upload_server <- function(id, lis.url=NULL, session) {
  moduleServer(id, function(input, output, session) {
  message('Upload module starts ... ')
  output$bulk.sce <- renderUI({
    ns <- session$ns; fileIn <- input$fileIn
    if (fileIn=='customBulkData') {
    list(
    fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '10%'), '', h4(strong("Step 2: upload custom data")), '',
      actionButton(ns("cusHelp"), "Help", icon = icon('question-circle')))),
      fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '18%', '1%', '25%', '1%', '25%'), '',
      fileInput(ns("geneInpath"), "2A: upload formatted data matrix", accept=c(".txt", ".csv"), multiple=FALSE), '',
      radioButtons(inputId=ns('dimName'), label='2B: is column or row gene?', choices=c("None", "Row", "Column"), selected='None', inline=TRUE), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "columns in the data matrix corresponds with "rows" in the targets file respectively.'),
      fileInput(ns("target"), "2C (optional): upload targets file for columns", accept=c(".txt", ".csv"), multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "rows" in the data matrix corresponds with "rows" in the row metadata file respectively.'),
      fileInput(ns("met"), "2D (optional): upload metadata file for rows", accept=c(".txt", ".csv"), multiple=FALSE))
      ))
    )
    } else if (fileIn=='customCovisData') {
      list(
      h4(strong('Step 2: single-cell/bulk data')),
      fileInput(ns("sglCell"), "", accept=c(".rds"), multiple=FALSE)
      )
   }
  })
  output$svg.upl <- renderUI({
    ns <- session$ns; fileIn <- input$fileIn
    if (fileIn %in% na.cus) {
    list(
      h5(strong("Step 3: upload custom aSVG(s)")),
      fluidRow(splitLayout(cellWidths=c('1%', '27%', '1%', '28%'), '',
      tags$div(class='tp', span(class='tpt', 'The data is matched with a single aSVG file.'),
      fileInput(ns("svgInpath1"), "3A: upload one aSVG file", accept=c('.svg', raster.ext), multiple=TRUE)), '',
      tags$div(class='tp', span(class='tpt', 'The data is matched with multiple aSVG files (e.g. developmental stages).'),
      fileInput(ns("svgInpath2"), "3B (optional): upload multiple aSVG files", accept=c('.svg', raster.ext), multiple=TRUE))
      ))
    )}
  })
  observeEvent(input$cusHelp, {
    showModal(
    div(id = 'datFormat', modalDialog(title = HTML('<strong><center>Data Formats</center></strong>'),
      div(style = 'overflow-y:scroll;overflow-x:scroll',
      HTML('<img src="image/data_format_shiny.jpg">'),
      tags$a(href="https://bioconductor.org/packages/devel/bioc/vignettes/spatialHeatmap/inst/doc/spatialHeatmap.html", target="_blank", "Package vignette")
      ))))
    })

  cfg <- reactiveValues(lis.dat=NULL, lis.dld=NULL, lis.par=NULL, na.def=NULL, dat.def=NULL, svg.def=NULL, pa.upl=NULL, pa.dat.upl=NULL, pa.svg.upl=NULL, na.cus=NULL, pa.svg.reg=NULL)
  lis.cfg <- yaml.load_file('config/config.yaml')
  lis.cfg <- lis.cfg[!vapply(lis.cfg, is.null, logical(1))]
  # Separate default datasets, downloadable datasets, and parameters.
  lis.dat <- lis.cfg[grepl('^dataset\\d+', names(lis.cfg))]
  db.pa <- 'data/data_shm.tar'
  # Merge separate data sets and data base.
  if (file.exists(db.pa)) { 
    lis.dat.db <- ovl_dat_db(data=lis.dat, db=db.pa)
    lis.dat <- lis.dat.db$data; db.dat <- lis.dat.db$dat.db
    lis.dat <- c(lis.dat, db.dat)
    names(lis.dat) <- paste0('dataset', seq_along(lis.dat))
  }
  dld.na <- c('download_single', 'download_multiple', 'download_multiple_variables', 'download_batched_data_aSVGs', 'download_covisualization')
  lis.dld <- lis.cfg[grepl(paste0(dld.na, collapse='|'), names(lis.cfg))]
  observe({
    if (is.null(input$config)) lis.par <- lis.cfg[!grepl(paste0(c('^dataset\\d+', dld.na), collapse='|'), names(lis.cfg))] else lis.par <- yaml.load_file(input$config$datapath[1])
    upl.size <- toupper(lis.par$max.upload.size)
    num <- as.numeric(gsub('(\\d+)(G|M)', '\\1', upl.size))
    if (grepl('\\d+G$', upl.size)) max.size <- num*1024^3
    if (grepl('\\d+M$', upl.size)) max.size <- num*1024^2 
    options(shiny.maxRequestSize=max.size) 
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
    na.ipt <- dis.ipt <- dat.ipt <- svg.ipt <- NULL; for (i in lis.dat) {  
      na.ipt <- c(na.ipt, i$name); dat.ipt <- c(dat.ipt, i$data)
      svg.ipt <- c(svg.ipt, list(i$svg)); dis.ipt <- c(dis.ipt, i$display);
    }; names(dat.ipt) <- names(svg.ipt) <- na.ipt
    # Uploaded tar files.
    df.tar <- input$tar; dat.upl <- svg.upl <- NULL
    tar.num <- grepl('\\.tar$', df.tar$datapath)
    if (!is.null(df.tar)) validate(need(try(sum(tar.num)==2), 'Two separate tar files of respective data and aSVGs are expected!'))
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
    idx.def <- !na.ipt %in% na.cus; na.def <- na.ipt[idx.def]
    dis.def <- dis.ipt[idx.def] 
    # Data in uploaded tar files are also included in default.
    dat.def <- c(dat.upl, dat.ipt[na.def]); svg.def <- c(svg.upl, svg.ipt[na.def])
    dis.def <- unique(c(names(dat.upl), dis.def))
    # If data/svg are duplicated between the server and upload, the data/svg on server side is removed.
    dat.def <- dat.def[unique(names(dat.def))]; svg.def <- svg.def[unique(names(svg.def))]
    cfg$lis.dat <- lis.dat; cfg$lis.dld <- lis.dld; cfg$lis.par <- lis.par; cfg$na.def <- setNames(names(dat.def), dis.def)
    cfg$svg.def <- svg.def; cfg$dat.def <- dat.def; cfg$na.cus <- setNames(na.cus, na.cus.dis)
    dat.nas <- c(na.cus, names(dat.def)); names(dat.nas) <- c(na.cus.dis, dis.def) 
    url.val <- url_val('upl-fileIn', lis.url)
    updateSelectInput(session, 'fileIn', choices=dat.nas, selected=ifelse(url.val=='null', lis.par$default.dataset, url.val))
    updateRadioButtons(session, inputId='dimName', selected=lis.par$col.row.gene, inline=TRUE)

  })
  observe({ # aSVG uploaded in regular files (not tar), used in re-matching.
    if (is.null(input$svgInpath2)) svgIn.df <- input$svgInpath1 else svgIn.df <- input$svgInpath2
    if (is.null(svgIn.df)) return()
    if (input$fileIn %in% na.cus) {
      svg.path <- svgIn.df$datapath; svg.na <- svgIn.df$name
      # Raster images uploaded.
      if (any(!grepl('\\.svg$', svg.na))) svg_raster(svg.na, raster.ext)
      # SVG and template paths are processed in the same way, and both are placed in the same list.
      pa.svg.reg <- list(unlist(lapply(seq_along(svg.na), function(x) {
          strs <- strsplit(svg.path[x], '/')[[1]]
          strs <- strs[-length(strs)]
          paste0(c(strs, svg.na[x]), collapse='/')
        }
      ))); names(pa.svg.reg) <- 'uploaded'
    # Original copy used for regular SHMs, the 2nd copy used in rematching.
    file.copy(svg.path, pa.svg.reg[[1]])
    cfg$pa.svg.reg <- pa.svg.reg
    }
  })
  observe({
    input$fileIn; input$geneInpath
    updateRadioButtons(session, inputId="dimName", selected="None")
  })
 observe({ 
    dld.exp <- reactiveValues(sgl=NULL, mul=NULL, st=NULL, bat = NULL)
    dld.exp$sgl <- cfg$lis.dld$download_single
    dld.exp$mul <- cfg$lis.dld$download_multiple
    dld.exp$st <- cfg$lis.dld$download_multiple_variables
    dld.exp$bat <- cfg$lis.dld$download_batched_data_aSVGs
    dld.exp$covis <- cfg$lis.dld$download_covisualization

    output$dld.cfg <- downloadHandler(
      filename=function(){ "config_par.yaml" },
 content=function(file=paste0(tmp.dir, '/config_par.yaml')){  
        lis.cfg <- yaml.load_file('config/config.yaml')
        par.na <- c("max.upload.size", "default.dataset", "col.row.gene", "separator", "data.matrix", "shm.img", "shm.anm", "shm.video", "legend", "mhm", "network")
        par.na <- par.na[par.na %in% names(lis.cfg)]
        lis.par <- lis.cfg[par.na]; write_yaml(lis.par, file)
      }
    )
    output$dld.sgl <- downloadHandler(
      filename=function(){ "single_aSVG_data.zip" },  content=function(file=paste0(tmp.dir, '/single_aSVG_data.zip')){ zip(file, c(dld.exp$sgl$data, dld.exp$sgl$svg)) }
    )
    output$dld.mul <- downloadHandler(   
      filename=function(){ "multiple_aSVG_data.zip" }, content=function(file=paste0(tmp.dir, '/multiple_aSVG_data.zip')){ zip(file, c(dld.exp$mul$data, dld.exp$mul$svg)) }
  ) 
    output$dld.st <- downloadHandler(   
      filename=function(){ "multiVariables_aSVG_data.zip" },
content=function(file=paste0(tmp.dir, '/multiVariables_aSVG_data.zip')){ zip(file, c(dld.exp$st$data, dld.exp$st$svg)) }
  )
    output$dld.covis <- downloadHandler(   
      filename=function(){ "covisualization_aSVG_data.zip" },
content=function(file=paste0(tmp.dir, '/covisualization_aSVG_data.zip')){ zip(file, c(dld.exp$covis$data, dld.exp$covis$svg)) }
  )
    output$dld.bat <- downloadHandler(   
      filename=function(){ "batched_data_aSVGs.zip" },
content=function(file=paste0(tmp.dir, '/batched_data_aSVGs.zip')){ zip(file, c(dld.exp$bat$data, dld.exp$bat$svg)) }
  )
  })

  # URLs on the landing page.
  output$brain.hum <-renderUI({
  tagList(
    p('Human brain', style='font-size:18px'),
  a(img(width='97%', src="image/brain_hum.png"), href='')
    )
  })
  output$mouse <-renderUI({
  tagList(
    p('Mouse organ', style='font-size:18px'),
    a(img(width='97%', src="image/mouse.png"), href='')
  )
  })
  output$chicken <-renderUI({
  tagList(
    p('Chicken organ', style='font-size:18px'),
    a(img(width='97%', src="image/chicken.png"), href='')
    )
  })
  output$organ.arab <-renderUI({
  tagList(
    p('Organ', style='font-size:18px'),
    a(img(width='97%', src="image/organ_arab.png"), href='')
    )
  })
  output$shoot.arab <-renderUI({
  tagList(
    p('Shoot tissue', style='font-size:18px'),
    a(img(width='97%', src="image/shoot_arab.png"), href='')
  )
  })
  output$root.arab <-renderUI({
  tagList(
    p('Root tissue', style='font-size:18px'),
    a(img(width='97%', src="image/root_arab.png"), href='')
    )
  })
  output$stage.arab <-renderUI({
  tagList(
    p('Developmental stage', style='font-size:18px'),
    a(img(width='97%', src="image/stage_arab.png"), href='')
    )
  })
  output$clp.rice <-renderUI({
  tagList(
    p('Mouse brain multi-variable data', style='font-size:18px'),
    a(img(width='97%', src="image/mus_multi_dim.png"), href='')
    )
  })

  observe({
    toggleState(id = "geneInpath", condition = input$fileIn %in% cfg$na.cus)
    toggleState(id = "dimName", condition = input$fileIn %in% cfg$na.cus)
    toggleState(id = "target", condition = input$fileIn %in% cfg$na.cus)
    toggleState(id = "met", condition = input$fileIn %in% cfg$na.cus)
    toggleState(id = "svgInpath1", condition = input$fileIn %in% cfg$na.cus)
    toggleState(id = "svgInpath2", condition = input$fileIn %in% cfg$na.cus)
  })
  # Switch to avoid files uploaded previously. E.g. 1. upload 'sce.rds' under 'customSingleCell'. 2. select 'brain_Prudencio'. 3. re-select 'customSingleCell', and 'sce.rds' in step 1 is avoided.
  sce.pa <- reactiveValues(val=TRUE)
  observeEvent(list(input$sglCell), { sce.pa$val <- TRUE })
  observeEvent(list(input$fileIn), { sce.pa$val <- FALSE })
  # observeEvent(input$fileIn, { })
  sce <- reactiveValues(); observe({
  library(SingleCellExperiment)
  library(scater); library(scran); library(BiocSingular)
  sgl.cell.ipt <- input$sglCell; pa <- NULL
  # save(sgl.cell.ipt, file='sgl.cell.ipt')
  if (is.null(sgl.cell.ipt) | sce.pa$val==FALSE) { sce$val <- NULL } else pa <- sgl.cell.ipt$datapath
  if (!is.null(input$fileIn)) if (grepl(na.sgl.def, input$fileIn)) pa <- cfg$dat.def[input$fileIn] 
  if (is.null(pa)) return()
  if (grepl('\\.rds$', pa)) sce$val <- readRDS(pa)
  })
  onBookmark(function(state) { state })
  return(list(ipt = input, cfg = cfg, sce=sce))
})}
