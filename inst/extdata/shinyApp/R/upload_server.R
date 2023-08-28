# Module for uploading data.

upload_server <- function(id, lis.url=NULL, prt=NULL, session) {
  moduleServer(id, function(input, output, session) {
  message('Upload module starts ... '); ns <- session$ns
  observeEvent(input$dathelp, {  
    showModal(modal(title='Quick start!', msg = NULL, img='dataset.jpg', img.w="100%"))
  })
  output$bulk.sce <- renderUI({
    ns <- session$ns; fileIn <- input$fileIn
    if (fileIn=='customBulkData') {
    list(
    fluidRow(splitLayout(cellWidths=c('10px', '280px', '1px'), '', h4(strong("Step2: upload custom assay data")), '')),
      fluidRow(splitLayout(cellWidths=c('10px', '180px', '1px', '200px', '1px', '230px', '1px', '190px', '1px', '210px'), '',
      tags$div(class='tp', span(class='tpt', 'A tabular file or ".rds" file ("SummarizedExperiment" saved with "saveRDS")'), fileInput(ns("geneInpath"), "2A: formatted assay data", accept=c(".txt", ".csv", '.rds'), multiple=FALSE)), '',
      selectInput(ns('dimName'), label='2B: genes in column or row?', choices=c("Row", "Column"), selected='Row'), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "columns in the data matrix corresponds with "rows" in the targets file respectively.'),
      fileInput(ns("target"), "2C (optional): sample targets file", accept=c(".txt", ".csv"), multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "rows" in the data matrix corresponds with "rows" in the row metadata file respectively.'),
      fileInput(ns("met"), "2D (optional): row metadata", accept=c(".txt", ".csv"), multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'Assay metadata in a tabular file.'),
      fileInput(ns("asymet"), "2E (optional): assay metadata", accept=c(".txt", ".csv"), multiple=FALSE)))
      )
    )
    } else if (fileIn=='customCovisData') {
      list(
      fluidRow(splitLayout(cellWidths=c('10px', '280px', '1px'), '', h4(strong("Step2: upload custom assay data")), '')),
      fluidRow(splitLayout(cellWidths=c('10px', '300px', '1px', '200px', '1px', '230px', '1px', '190px', '1px', '210px'), '',
      div(class='tp', span(class='tpt', 'A tabular file or ".rds" file ("SingleCellExperiment" that combines bulk & single-cell data, saved with "saveRDS")'),
      fileInput(ns("sglCell"), "2A: formatted bulk & single-cell assay data", accept=c(".txt", ".csv", ".rds"), multiple=FALSE)
      ), '', 
      selectInput(ns('dimNaCovis'), label='2B: genes in column or row?', choices=c("Row", "Column"), selected='Row'), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "columns in the data matrix corresponds with "rows" in the targets file, respectively.'),
      fileInput(ns("tarCovis"), "2C: sample targets file", accept=c(".txt", ".csv"), multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "rows" in the data matrix corresponds with "rows" in the row metadata file, respectively.'),
      fileInput(ns("rmetCovis"), "2D (optional): row metadata", accept=c(".txt", ".csv"), multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'Assay metadata in a tabular file.'),
      fileInput(ns("metCovis"), "2E (optional): assay metadata", accept=c(".txt", ".csv"), multiple=FALSE)))
      )
   )
   }
  })
  output$svg.upl <- renderUI({
    ns <- session$ns; fileIn <- input$fileIn
    if (fileIn %in% na.cus) {
    list(
      h4(strong("Step3: upload custom aSVG(s)")),
      fluidRow(splitLayout(cellWidths=c('10px', '500px', '1px', '500px'), '',
      tags$div(class='tp', span(class='tpt', 'The assay data is matched with a single aSVG file.'),
      fileInput(ns("svgInpath1"), "3A: one aSVG file", accept=c('.svg', raster.ext), multiple=TRUE)), '',
      tags$div(class='tp', span(class='tpt', 'The assay data is matched with multiple aSVG files (e.g. developmental stages).'),
      fileInput(ns("svgInpath2"), "3B (optional): multiple aSVG files", accept=c('.svg', raster.ext), multiple=TRUE))
      ))
    )}
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
  cnt.ldg <- reactiveValues(v=0)
  observe({
    message('Config file ...')
    tabTop <- prt$input$tabTop; if (!check_obj(tabTop)) return()
    if (tabTop %in% c('ldg', 'about') & cnt.ldg$v==0) return()
    cnt.ldg$v <- 1
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
    dat.nas <- c(dat.no, na.cus, names(dat.def)); names(dat.nas) <- c(dat.no.dis, na.cus.dis, dis.def)
    cfg$dat.nas <- dat.nas
  })
  observeEvent(cfg$dat.nas, {
    tabTop <- prt$input$tabTop; lis.par <- cfg$lis.par; dat.nas <- cfg$dat.nas
    req(check_obj(list(tabTop, lis.par, dat.nas)))
    if (tabTop %in% c('ldg', 'about') & cnt.ldg$v==0) return()
    url.val <- url_val('upl-fileIn', lis.url)
    updateSelectInput(session, 'fileIn', choices=dat.nas, selected=ifelse(url.val=='null', lis.par$default.dataset, url.val))
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
      ))); names(pa.svg.reg) <- 'uploaded'; #save(pa.svg.reg, file='pa.svg.reg')
    # Original copy used for regular SHMs, the 2nd copy used in rematching.
    file.copy(svg.path, pa.svg.reg[[1]])
    cfg$pa.svg.reg <- pa.svg.reg
    }
  })
  observe({
    input$fileIn; input$geneInpath
    #updateRadioButtons(session, inputId="dimName", selected="None")
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

  )
    output$dld.bat <- downloadHandler(   
      filename=function(){ "batched_data_aSVGs.zip" },
content=function(file=paste0(tmp.dir, '/batched_data_aSVGs.zip')){ zip(file, c(dld.exp$bat$data, dld.exp$bat$svg)) }
  )
  })

  observeEvent(list(input$geneInpath), {
    pa <- input$geneInpath$datapath; 
    if (sum(grepl('\\.rds$', pa))==1) { 
      hide(id='dimName'); hide(id='target')
      hide(id='met'); hide(id='asymet')
    } else {
      shinyjs::show(id='dimName')
      shinyjs::show(id='target')
      shinyjs::show(id='met')
      shinyjs::show(id='asymet')
    }
  })
  observeEvent(list(input$sglCell), {
    pa <- input$sglCell$datapath; 
    if (sum(grepl('\\.rds$', pa))==1) { 
      hide(id='dimNaCovis'); hide(id='tarCovis')
      hide(id='rmetCovis'); hide(id='metCovis')
    } else {
      shinyjs::show(id='dimNaCovis')
      shinyjs::show(id='tarCovis')
      shinyjs::show(id='rmetCovis')
      shinyjs::show(id='metCovis')
    }
  })
  output$help <- renderUI({
    tags$iframe(seamless="seamless", src= "html/shm_shiny_manual.html#1_Datasets", width='100%', height='100%')
  })
  # Switch to avoid files uploaded previously. E.g. 1. upload 'sce.rds' under 'customSingleCell'. 2. select 'human_brain'. 3. re-select 'customSingleCell', and 'sce.rds' in step 1 is avoided.
  # fileInput cannot be set NULL with Shiny.setInputValue, since the uploaded file is cached.
  covis.pa <- reactiveValues(val=TRUE)
  observeEvent(list(input$fileIn), { 
    sce$val <- NULL
    covis.pa$dat <- covis.pa$svg1 <- covis.pa$svg2 <- NULL
  })
  observeEvent(list(input$sglCell), { covis.pa$dat <- input$sglCell$datapath })
  observeEvent(list(input$svgInpath1), { covis.pa$svg1 <- input$svgInpath1$datapath })
  observeEvent(list(input$svgInpath2), { covis.pa$svg2 <- input$svgInpath2$datapath })
  sce <- reactiveValues()
  observeEvent(list(covis.pa$dat, input$fileIn, cfg$dat.def, cfg$na.cus, input$sglCell, input$dimNaCovis, input$rmetCovis, input$metCovis, input$svgInpath1, input$svgInpath2), {
  library(SingleCellExperiment)
  library(scater); library(scran); library(BiocSingular)
  fileIn <- input$fileIn; pa <- NULL
  req(check_obj(list(fileIn, cfg$na.cus, cfg$dat.def))) 
  # Uploaded data path.
  if (is.null(covis.pa$dat)) { sce$val <- NULL } else {
    pa <- covis.pa$dat
    svgInpa1 <- covis.pa$svg1; svgInpa2 <- covis.pa$svg2
    req(check_obj(list(pa, !is.null(svgInpa1)|!is.null(svgInpa2))))
    if (grepl('\\.rds$', pa)) sce$val <- readRDS(pa) else { # Tabular files uploaded.
      message('Importing covis data from tabular files ... ')
      if (fileIn %in% na.cus.covis) {
        withProgress(message="Loading covis data from tabular files: ", value = 0, {
        dimNa <- input$dimNaCovis; req(check_obj(list(dimNa)))
        incProgress(0.25, detail="importing data matrix ...")
        message('Importing covis assay data ... ')
        dat.cus <- fread_df(read_fr(pa), isRowGene=(dimNa=='Row'), rep.aggr=NULL)
        tarCovis <- input$tarCovis$datapath; lgc.tar <- is.null(tarCovis)      
        if (lgc.tar) showModal(modal(msg = 'When uploading tabular files for covisualization, a targets file is required!')); req(!lgc.tar)
        incProgress(0.3, detail="importing targets file ...")
        sce.rep <- as(dat.cus$se.rep, "SingleCellExperiment")
        message('Importing covis targets file ... ')
        df.tar <- read_fr(tarCovis)
        lgc.tar <- nrow(df.tar) == ncol(sce.rep)
        if (!lgc.tar) showModal(modal(msg = 'Ensure "columns" in the assay matrix corresponds with "rows" in the targets file, respectively!')); req(lgc.tar)
        colData(sce.rep) <- DataFrame(df.tar) 
        rmetCovis <- input$rmetCovis$datapath
        if (!is.null(rmetCovis)) { 
          incProgress(0.3, detail="importing row metadata ...")
          message('Importing covis row metadata ... ')
          df.rmet <- read_fr(rmetCovis); lgc.met <- nrow(df.rmet) == nrow(sce.rep)
          if (!lgc.met) showModal(modal(msg = 'Ensure "rows" in the assay matrix corresponds with "rows" in the row metadata file, respectively!'))
          req(lgc.met); rowData(sce.rep) <- DataFrame(df.rmet)
        }
        metCovis <- input$metCovis$datapath
        if (!is.null(metCovis)) { 
          incProgress(0.3, detail="importing assay metadata ...")
          message('Importing covis metadata ... ')
          df.meta <- read_fr(metCovis)
          if (!is.null(df.meta)) {
            lgc.met <- ncol(data.frame(df.meta))<2
            if (lgc.met) { 
              msg <- 'The assay metadata should be a "data.frame" with at least two columns!';  
              if (lgc.met) showModal(modal(msg = msg)); req(!lgc.met)
            }; metadata(sce.rep)$df.meta <- df.meta
          }
        }; incProgress(0.3, detail="done!")
          sce$val <- sce.rep; message('Done! \n')
        })
      }
    }
   }
  # Default data path of covis data.
  if (grepl(na.sgl.def, fileIn)) { 
    pa <- cfg$dat.def[fileIn]; req(check_obj(list(pa)))
    if (grepl('\\.rds$', pa)) sce$val <- readRDS(pa)
  }
  if (!is.null(sce$val)) {
    lgc.na <- length(assayNames(sce$val))>1
    if (lgc.na) showModal(modal(msg = 'Only one count matrix is expected in "assay(<SingleCellExperiment>)"!')); req(!lgc.na)
    assayNames(sce$val) <- 'counts' 
  }
  })
  onBookmark(function(state) { state })
  return(list(ipt = input, cfg = cfg, sce=sce, covis.pa=covis.pa))
})}


