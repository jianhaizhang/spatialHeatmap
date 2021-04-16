
options(stringsAsFactors=FALSE) 

# source('~/tissue_specific_gene/function/fun.R')
# Every variable in every container should be checked at the beginning. E.g. input$fileIn in reactive({}). These checks will avoid amost all errors/warnings.
# Right before submit the package the following functions will be deleted, and they will be imported as above. They are listed here now for the convenience of functionality development.


# enableWGCNAThreads()
shinyServer(function(input, output, session) {
  observe({
    withProgress(message="Loading dependencies: ", value=0, {
    incProgress(0.3, detail="in progress ...")
    library(SummarizedExperiment); library(shiny); library(shinydashboard); library(shinydashboardPlus); library(grImport); library(rsvg); library(ggplot2); library(DT) 
    incProgress(0.6, detail="in progress ...")
    library(gridExtra); library(ggdendro); library(WGCNA); library(grid); library(xml2); library(plotly); library(data.table); library(genefilter); library(flashClust); library(visNetwork); 
    showModal(modal(title = HTML('<center><b>Welcome to spatialHeatmap!</b><center>'), msg = strong('Please wait for the landing page to finish!')))
    incProgress(0.9, detail="in progress...")
   library(reshape2); library(igraph); library(animation); library(av); library(shinyWidgets); library(yaml); library(HDF5Array); library(sortable); library(shinyBS); library(shinyjs); library(htmltools)
    # DEG
    library(UpSetR)
  })
  })

upload_server <- function(id) {
  moduleServer(id, function(input, output, session) {
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

      cat('Processing uploaded tar... \n')
      p <- df.tar$datapath[1]; strs <- strsplit(p, '/')[[1]]
      cfg$pa.upl <- pa.svg <- paste0(strs[grep('\\.tar$', strs, invert=TRUE)], collapse='/')
      dat.idx <- grepl('data_shm.tar$', df.tar$name) 
      cfg$pa.svg.upl <- df.tar$datapath[!dat.idx]
      # system(paste0('tar -xf', ' ', svg.tar, ' -C ', pa.svg))
      cfg$pa.dat.upl <- dat.pa <- df.tar$datapath[dat.idx]
      df.pair.upl <- read_hdf5(dat.pa, 'df_pair')[[1]]
      pair.na <- df.pair.upl$name; dat.upl <- df.pair.upl$data
      svg.upl <- as.list(df.pair.upl$aSVG); names(dat.upl) <- names(svg.upl) <- pair.na

      for (i in seq_along(svg.upl)) {
        svg0 <- svg.upl[[i]]; if (grepl(';| |,', svg0)) {
          strs <- strsplit(svg0, ';| |,')[[1]]; svg.upl[[i]] <- strs[strs!='']
        }
      }
    }
    # Separate data, svg of default and customization. 
    na.def <- na.ipt[!grepl('^none$|^customData$|^customComputedData$', na.ipt)]
    #if (na.cus) 
    na.cus <- c('customData', 'customComputedData')

    dat.def <- c(dat.upl, dat.ipt[na.def]); svg.def <- c(svg.upl, svg.ipt[na.def])
    # If data/svg are duplicated between the server and upload, the data/svg on server side is removed.
    dat.def <- dat.def[unique(names(dat.def))]; svg.def <- svg.def[unique(names(svg.def))]
    cfg$lis.dat <- lis.dat; cfg$lis.dld <- lis.dld; cfg$lis.par <- lis.par; cfg$na.def <- names(dat.def); cfg$svg.def <- svg.def; cfg$dat.def <- dat.def; cfg$na.cus <- na.cus
    # output$spatialHeatmap <- renderText({ lis.par$title['title', 'default'] })
    # output$title.w <- renderText({ lis.par$title['width', 'default'] })

    dat.nas <- c('none', 'customData', names(dat.def))
    updateSelectInput(session, 'fileIn', NULL, dat.nas, lis.par$default.dataset)
    updateRadioButtons(session, inputId='dimName', label='Step 2B: is column or row gene?', choices=c("None", "Row", "Column"), selected=lis.par$col.row.gene, inline=TRUE)
    updateNumericInput(session, inputId="A", label="Value (A) to exceed", value=as.numeric(lis.par$data.matrix['A', 'default']))
    updateNumericInput(session, inputId="P", label="Proportion (P) of samples with values >= A", value=as.numeric(lis.par$data.matrix['P', 'default']))
    updateNumericInput(session, inputId="CV1", label="Min coefficient of variation (CV1)", value=as.numeric(lis.par$data.matrix['CV1', 'default']))
    updateNumericInput(session, inputId="CV2", label="Max coefficient of variation (CV2)", value=as.numeric(lis.par$data.matrix['CV2', 'default']))
    updateRadioButtons(session, inputId='log', label='Log/exp transform', choices=c("No", "log2", "exp2"), selected=lis.par$log.exp, inline=TRUE)
    updateRadioButtons(session, inputId='scale', label='Scale by', choices=c('No', 'Row', 'Column'), selected=lis.par$data.matrix.scale, inline=TRUE)
    updateRadioButtons(session, inputId='measure', label="Measure:", choices=c('correlation', 'distance'), selected=lis.par$mhm['measure', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId="cor.abs", label="Cor.absolute:", choices=c('No', 'Yes'), selected=lis.par$mhm['cor.absolute', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId="thr", label="Select by:", choices=c('proportion'='p', 'number'='n', 'value'='v'), selected=lis.par$mhm['select.by', 'default'], inline=TRUE)
    updateNumericInput(session, inputId='mhm.v', label='Cutoff: ', value=as.numeric(lis.par$mhm['cutoff', 'default']), min=-Inf, max=Inf, step=NA)
    updateRadioButtons(session, inputId="mat.scale", label="Scale by:", choices=c("No", "Column", "Row"), selected=lis.par$mhm['scale', 'default'], inline=TRUE)
    updateSelectInput(session, inputId="net.type", label="Network type:", choices=c('signed', 'unsigned', 'signed hybrid', 'distance'), selected=lis.par$network['net.type', 'default'])
    updateNumericInput(session, "min.size", "Minmum module size:", value=as.numeric(lis.par$network['min.size', 'default']), min=15, max=5000)
    updateSelectInput(session, "ds","Module splitting sensitivity level:", 3:2, selected=lis.par$network['ds', 'default'])
    updateTextInput(session, "color.net", "Color scheme:", lis.par$network['color', 'default'], placeholder=paste0('Eg: ', lis.par$network['color', 'default']))
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
        lis.par <- lis.cfg[c("max.upload.size", "default.dataset", "col.row.gene", "separator", "hide.legend", "data.matrix", "shm.img", "shm.anm", "shm.video", "legend", "mhm", "network")]
        write_yaml(lis.par, file)
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
  return(list(ipt = input, cfg = cfg))
})}
  ipt.cfg <- upload_server('upl')
  ipt <- ipt.cfg$ipt; cfg <- ipt.cfg$cfg 

  ipt.fileIn <- ipt.geneInpath <- reactiveValues()
  observe({ ipt.fileIn <- ipt$fileIn; ipt.geneInpath <- ipt$ipt.geneInpath })


data_server <- function(id, ipt, cfg) {
  moduleServer(id, function(input, output, session) {
  # Filter parameters.
  fil <- reactiveValues(P=0, A=0, CV1=-Inf, CV2=Inf)

  observe({
    ipt$fileIn; ipt$geneInpath; input$log 
    fil$P <- 0; fil$A <- 0; fil$CV1 <- -Inf; fil$CV2 <- Inf
  })

  observeEvent(input$fil.but, {
    if (ipt$fileIn=="none") return(NULL)  
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
        }; df.te <- fread_df(input=dat)
      }; return(df.te)
    }
    if (ipt$fileIn %in% cfg$na.cus & 
    !is.null(ipt$geneInpath) & ipt$dimName!="None") {
      incProgress(0.25, detail="importing matrix, please wait ...")
      geneInpath <- ipt$geneInpath$datapath; targetInpath <- ipt$target$datapath; metInpath <- ipt$met$datapath
      # Keep replicates unchaged, and compared with targets/metadata files.
      df.upl <- fread_df(read_fr(geneInpath), rep.aggr = NULL)
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
      df.upl1 <- fread_df(df.rep, rep.aggr = 'mean')
      df.upl$df.aggr <- df.upl1$df.aggr
      df.upl$df.rep <- df.rep; df.upl$df.met <- df.met 
      df.upl$con.na <- df.upl1$con.na
      if ((!is.null(ipt$svgInpath1) | !is.null(ipt$svgInpath2))) return(df.upl) 
    }

    })

  })

  # Transform data.
  geneIn1 <- reactive({
    if (is.null(geneIn0())|is.null(input$scale) | is.null(input$log)) return(NULL)
    df.aggr <- df.aggr.tran <- geneIn0()[['df.aggr']] 
    if (input$log=='log2') { 
      g.min <- min(df.aggr.tran)
      if (g.min<0) df.aggr.tran <- df.aggr.tran-g.min+1; if (g.min==0) df.aggr.tran <- df.aggr.tran+1; df.aggr.tran <- log2(df.aggr.tran)
    } else if (input$log=='exp2') df.aggr.tran <- 2^df.aggr.tran
    # Scale by row/column
    if (input$scale=='Row') { df.aggr.tran <- t(scale(t(df.aggr.tran))) } else if (input$scale=='Column') { df.aggr.tran <- scale(df.aggr.tran) }
    df.met <- geneIn0()[['df.met']]; df.rep <- geneIn0()[['df.rep']]
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
  observeEvent(input$search.but, {
    if (is.null(geneIn1())) return()
    if (input$search=='') sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) else {
      gens <- strsplit(gsub(' |,', '_', input$search), '_')[[1]]
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
    if (is.null(geneIn1())) return(NULL)    
    df.aggr <- geneIn1()[['df.aggr']]; df.aggr.tran <- geneIn1()[['df.aggr.tran']]; df.met <- geneIn1()[['df.met']]; input$fil.but
    # When ipt$fileIn changes, col.reorder$col.na.re of last session still exists and causes errors.
    if (!identical(sort(col.reorder$col.na.re), sort(colnames(df.aggr.tran)))) return()
    # Input variables in "isolate" will not triger re-excution, but if the whole reactive object is trigered by "input$fil.but" then code inside "isolate" will re-excute.
    isolate({ 
      se <- SummarizedExperiment(assays=list(expr=as.matrix(df.aggr.tran)), rowData=df.met)
      if (ncol(df.met)>0) ann.col <- colnames(df.met)[1] else ann.col <- NULL
      # If scaled by row, sd is 1, mean is 0, cv is Inf.
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=c(ifelse(input$scale=='Row', -Inf, fil$CV1), ifelse(input$scale=='Row', Inf, fil$CV2)), dir=NULL, verbose=FALSE)
      if (nrow(se)==0) { validate(need(try(nrow(se)>0), 'All rows are filtered out!')); return() }
      # In case all rows are filtered, the app continues to work without refreshing after the filter parameters are reduced.
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=c(ifelse(input$scale=='Row', -Inf, fil$CV1), ifelse(input$scale=='Row', Inf, fil$CV2)), dir=NULL, verbose=FALSE)
      df.aggr.tran <- as.data.frame(assay(se), stringsAsfactors=FALSE)
      colnames(df.aggr.tran) <- make.names(colnames(df.aggr.tran))
       df.aggr <- df.aggr[rownames(df.aggr.tran), ]
       df.met <- as.data.frame(rowData(se))[, , drop=FALSE]
    })
    cat('Preparing data matrix... \n')
    if (is.null(sear$id)) rows <- seq_len(nrow(df.aggr.tran)) else rows <- sear$id
    if (length(rows)==1 & rows[1]==as.numeric(cfg$lis.par$data.matrix['row.selected', 'default'])) rows <- seq_len(nrow(df.aggr.tran))
    df.aggr.tran.order <- df.aggr.tran[rows, col.reorder$col.na.re, drop=FALSE]
    return(list(df.aggr=df.aggr, df.aggr.tran=df.aggr.tran, df.aggr.tran.order = df.aggr.tran.order, df.met=df.met[rows, , drop=FALSE]))

  })
  
  output$dt <- renderDataTable({

    if (is.null(geneIn())) return()
    if ((ipt$fileIn %in% cfg$na.cus & is.null(geneIn()))|ipt$fileIn=="none") return(NULL)

    withProgress(message="Data table: ", value = 0, {
      incProgress(0.5, detail="displaying, please wait ...")
      if (ipt$fileIn!="none") {
      df.all <- geneIn(); gene.dt <- cbind.data.frame(df.all[["df.met"]][, , drop=FALSE], df.all[["df.aggr.tran"]][, , drop=FALSE], stringsAsFactors=FALSE)
   }; cat('Presenting data matrix... \n')
   if (is.null(sear$id)) sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) else sel <- sear$id
   if (length(sel)==1 & sel[1]==as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) & nrow(gene.dt)>1) sel <- sel else if (nrow(gene.dt)==1)  sel <- 1 else if (length(sel)>1) sel <- seq_along(sel)
   datatable(gene.dt, selection=list(mode="multiple", target="row", selected=sel),
   filter="top", extensions=c('Scroller'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE, searchHighlight=FALSE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs = list(list(targets = c(1), render = JS("$.fn.dataTable.render.ellipsis(5, false)")))), 
   class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
   formatRound(colnames(geneIn()[["df.aggr.tran"]]), 2)
    })
  })

  observe({
    ipt$fileIn; ipt$geneInpath; lis.par <- cfg$lis.par
    updateRadioButtons(session, inputId='log', label='Log/exp transform', choices=c("No", "log2", "exp2"), selected=cfg$lis.par$data.matrix['log.exp', 'default'], inline=TRUE)
    updateRadioButtons(session, 'scale', label='Scale by', choices=c('No', 'Row', 'Column'), selected=cfg$lis.par$data.matrix['scale', 'default'], inline=TRUE)
  })

  observe({
    ipt$fileIn; ipt$geneInpath; input$log
    updateNumericInput(session, inputId="A", label="Value (A) to exceed:", value=as.numeric(cfg$lis.par$data.matrix['A', 'default'])) 
    updateNumericInput(session, inputId="P", label="Proportion (P) of samples with values >= A:", value=as.numeric(cfg$lis.par$data.matrix['P', 'default']), min=0, max=1)
    updateNumericInput(session, inputId="CV1", label="Min coefficient of variation (CV1):", value=as.numeric(cfg$lis.par$data.matrix['CV1', 'default']))
    updateNumericInput(session, inputId="CV2", label="Max coefficient of variation (CV2):", value=as.numeric(cfg$lis.par$data.matrix['CV2', 'default'])) 
  })
  
  ipt.dat <- reactiveValues()
  observe({ ipt.dat$dt_rows_selected <- input$dt_rows_selected })
  col.cfm <- reactive({ input$col.cfm })
  col.na <- reactive({ input$na })
  log <- reactive({ input$log }); A <- reactive({ input$A })
  CV1 <- reactive({ input$CV1 }); CV2 <- reactive({ input$CV2 })
  P <- reactive({ input$P })
  search.but <- reactive({ input$search.but })
  return(list(geneIn0 = geneIn0, geneIn1 = geneIn1, geneIn = geneIn, sear = sear, ipt.dat = ipt.dat, col.reorder = col.reorder, col.cfm = col.cfm, col.na = col.na, log = log, A = A, P = P, CV1 = CV1, CV2 = CV2, search.but = search.but))
  })

}

  dat.mod.lis <- data_server('dat', ipt, cfg)
  # The reactive type in and outside module is the same: sear is a reactiveValue in and outside module; geneIn is reactive expression in and outside module. "geneIn()" is accessing the content of a reactive expression, and loses the "reactive" attribute.
  # As long as the content of reactiveValues (col.reorder$col.na.re) is not accessed, the operation does not need to be inside reactive environment (observe).
  ipt.dat <- reactiveValues()
  ipt.dat$dat <- dat.mod.lis$ipt.dat; sear <- dat.mod.lis$sear
  col.reorder <- reactiveValues(); col.reorder <- dat.mod.lis$col.reorder
  geneIn0 <- dat.mod.lis$geneIn0; geneIn1 <- dat.mod.lis$geneIn1
  geneIn <- dat.mod.lis$geneIn

shm_server <- function(id, ipt, cfg, dat.mod.lis) {
  
  ipt.dat$dat <- dat.mod.lis$ipt.dat; sear <- dat.mod.lis$sear
  col.reorder <- reactiveValues(); col.reorder <- dat.mod.lis$col.reorder
  geneIn0 <- dat.mod.lis$geneIn0; geneIn1 <- dat.mod.lis$geneIn1
  geneIn <- dat.mod.lis$geneIn
  col.na <- dat.mod.lis$col.na; col.cfm <- dat.mod.lis$col.cfm 
  log <- dat.mod.lis$log; A <- dat.mod.lis$A
  search.but <- dat.mod.lis$search.but
  moduleServer(id, function(input, output, session) {

  gID <- reactiveValues(geneSel="none", new=NULL, all=NULL)
  observe({ ipt$geneInpath; ipt$fileIn; gID$geneSel <- "none" })
  observe({ if (is.null(geneIn())) gID$geneSel <- "none" })
  # To make the "gID$new" and "gID$all" updated with the new "input$fileIn", since the selected row is fixed (3rd row), the "gID$new" is not updated when "input$fileIn" is changed, and the downstream is not updated either. The shoot/root examples use the same data matrix, so the "gID$all" is the same (pre-selected 3rd row) when change from the default "shoot" to others like "organ". As a result, the "gene$new" is null and downstream is not updated. Also the "gene$new" is the same when change from shoot to organ, and downstream is not updated, thus "gene$new" and "gene$all" are both set NULL above upon new "input$fileIn".  
  observeEvent(ipt$fileIn, {

    if (is.null(ipt.dat$dat$dt_rows_selected)) return()
    gID$all <- gID$new <- NULL
    r.na <- rownames(geneIn()[["df.aggr.tran"]]); gID$geneSel <- r.na[ipt.dat$dat$dt_rows_selected]
    # Avoid multiple selected rows from last input$fileIn. Must be behind gID$geneSel. 
    if (length(ipt.dat$dat$dt_rows_selected)>1) return()
    gID$new <- setdiff(gID$geneSel, gID$all); gID$all <- c(gID$all, gID$new)
    if (is.null(r.na)) gID$geneSel <- "none"
    cat('New file:', ipt$fileIn, '\n')
    cat('New ID:', gID$new, 'Selected ID:', gID$geneSel, 'All ID:', gID$all, '\n')
    })
  observeEvent(list(ipt.dat$dat$dt_rows_selected, geneIn()), {
    if (is.null(ipt.dat$dat$dt_rows_selected)) return()
    r.na <- rownames(geneIn()[["df.aggr.tran"]]); gID$geneSel <- r.na[ipt.dat$dat$dt_rows_selected]
    if (any(is.na(gID$geneSel))) gID$geneSel <- "none"
    gID$new <- setdiff(gID$geneSel, gID$all); gID$all <- c(gID$all, gID$new) 
    cat('Change in data matrix/selected rows:', '\n')
    cat('New ID:', gID$new, 'Selected ID:', gID$geneSel, 'All ID:', gID$all, '\n')
  })
  observeEvent(list(search.but()), {
  
    if (is.null(search.but())|is.null(sear$id)|is.null(geneIn())) return()
    gID$geneSel <- rownames(geneIn()[["df.aggr.tran"]])[sear$id]
    if (any(is.na(gID$geneSel))) gID$geneSel <- "none"
 
  })

  geneV <- reactive({

    if (any(is.na(gID$geneSel))) return()
    if (is.null(geneIn())|sum(gID$geneSel[1]!='none')==0) return(NULL)
    if (input$cs.v=="Selected rows" & is.null(ipt.dat$dat$dt_rows_selected)) return(NULL)
    if (ipt$fileIn!="none") { if (input$cs.v=="Selected rows") gene <- geneIn()[["df.aggr.tran"]][gID$geneSel, ]
    if (input$cs.v=="All rows") gene <- geneIn()[["df.aggr.tran"]] }
    seq(min(gene), max(gene), len=1000) # len must be same with that from the function "spatial_hm()". Otherwise the mapping of a gene value to the colour bar is not accurate. 

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
    col0 <- cfg$lis.par$shm.img['color', 'default']
    if (is.null(input$col.but) | is.null(col0)) return()
    if(input$col.but==0) color$col <- colorRampPalette(col_sep(col0))(length(geneV()))
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

  shm.bar <- reactive({
    cat('Colour key ... \n')
    if (is.null(gID$all)) return(NULL)
    if ((any(ipt$fileIn %in% cfg$na.def) & !is.null(geneIn()))|(ipt$fileIn %in% cfg$na.cus & (!is.null(ipt$svgInpath1)|!is.null(ipt$svgInpath2)) & !is.null(geneIn()))) {
      if (length(color$col=="none")==0|input$color==""|is.null(geneV())) return(NULL)
      withProgress(message="Color scale: ", value = 0, {
        incProgress(0.75, detail="plotting, please wait ...")
        cs.g <- col_bar(geneV=geneV(), cols=color$col, width=1)
        cat('Done! \n'); return(cs.g)
      })
    }
  })
  # One output can only be used once in ui.R.
  output$bar1 <- output$bar2 <- renderPlot({ if (!is.null(shm.bar)) shm.bar() })
  # output$bar2 <- renderPlot({ if (!is.null(shm.bar)) shm.bar() })

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

  output$svg <- renderUI({
    nas <- names(cfg$svg.def)
    selectInput('svg', label='Choose an aSVG to re-match', choices=nas, selected=ipt$fileIn)
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
      fluidRow(
      span(class = "panel panel-default", style = 'margin-left:0px',
        div(class = "panel-heading", strong("Features in aSVG")),
        div(class = "panel-body", id = "ftSVG", ft2tag(sf.all1))
      ),
      div(class = "panel panel-default",
        div(class = "panel-heading", strong("Features in data")),
        lapply(sam.all1, ft_dat)
      ), lapply(c('ftSVG', sam.all1), ft_js) # Items are interchangeable across ftSVG and sam.all1.
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
      withProgress(message="Tissue heatmap: ", value=0, {
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
    ipt$fileIn; geneIn(); ipt$adj.modInpath; svg.df(); input$lgdB; 
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
  
  # Avoid repetitive computation under input$cs.v=='All rows'.
  gs.new <- reactive({
     # if.con <- is.null(svg.df())|is.null(geneIn())|is.null(gID$new)|length(gID$new)==0|is.null(gID$all)|is.null(ipt.dat$dat$dt_rows_selected)|color$col[1]=='none'
    # if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    validate(
      need(!is.null(svg.df()) & !is.null(geneIn()) & !is.null(geneIn()) & length(gID$new) > 0 & !is.null(gID$all) & !is.null(ipt.dat$dat$dt_rows_selected) & color$col[1]!='none', '')
    )
    scale.shm <- input$scale.shm 
    if (!is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    # If color key is build on selected rows, all SHMs should be computed upon selected rows are changed. This action is done through a separate observeEvent triggered by gID$geneSel. So in this "reactive" only one gene is accepted each time.
    if (input$cs.v=="Selected rows") ID <- gID$geneSel
    if (input$cs.v=="All rows") ID <- gID$new
    if (is.null(ID)) return()
    if (length(gID$new)>1|length(ID)>1|ID[1]=='none') return()
    # Avoid repetitive computation.  
    pat.new <- paste0('^', gID$new, '_(', pat.con(), ')_\\d+$')
    if (any(grepl(pat.new, names(shm$grob.all)))) return()
    withProgress(message="Tissue heatmap: ", value=0, { 
      incProgress(0.25, detail="preparing data ...")
      gene <- geneIn()[["df.aggr.tran"]]
      # When input$fileIn updates, ID is from last session while gene is from new session.
      if (!ID %in% rownames(gene)) return()
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
        cat('New grob/ggplot:', ID, ' \n')
        if (!is.null(cna.match$cna)) { 
		  if (ncol(gene)==length(cna.match$cna)) colnames(gene) <- cna.match$cna 
        }
        size.key <- as.numeric(cfg$lis.par$legend['key.size', 'default'])
        # Cores: the orders in svg.path(), names(svg.df.lis) are same.
        gg.lis <- gg_shm(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=ID, legend.col=fil.cols, cols=color$col, tis.path=tis.path, ft.trans=input$tis, sub.title.size=15 * scale.shm, aspect.ratio = svg.df$aspect.r, legend.nrow=as.numeric(cfg$lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.size=input$line.size, line.color=input$line.color, lis.rematch = ft.reorder$ft.rematch) # Only gID$new is used.
        msg <- paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')
        if (is.null(gg.lis)) cat(msg, '\n')
        output$msg.shm <- ({ validate(need(!is.null(gg.lis), msg)) })
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
     return(list(gg.all = gg.all, grob.all = grob.all, lgd.all = lgd.all, grob.gg.all = grob.gg.all))
    }) # withProgress

  })

  # Extension of 'observeEvent': any of 'input$log; input$tis; input$col.but; input$cs.v' causes evaluation of all code. 
  # input$tis as an argument in "gg_shm" will not cause evaluation of all code, thus it is listed here.
  # Use "observeEvent" to replace "observe" and list events (input$log, input$tis, ...), since if the events are in "observe", every time a new gene is clicked, "input$dt_rows_selected" causes the evaluation of all code in "observe", and the evaluation is duplicated with "gs.new".
  observeEvent(col.na(), { if (col.cfm()>0) col.reorder$col.re <- 'N' })
  # Update SHMs, above theme().
  observeEvent(list(log(), input$tis, input$col.but, input$cs.v, col.cfm(), input$scale, input$match, input$line.size, input$line.color), {
    shm$grob.all <- shm$gg.all <- shm$lgd.all <- shm$grob.gg.all <- NULL; gs.all <- reactive({ 
      cat('Updating all SHMs ... \n')
      if.con <- is.null(svg.df())|is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|color$col[1]=='none'|gID$geneSel[1]=='none'
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
        gg.lis <- gg_shm(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=gID$geneSel, legend.col=fil.cols, cols=color$col, tis.path=tis.path, ft.trans=input$tis, sub.title.size=15 * scale.shm, aspect.ratio = svg.df$aspect.r, legend.nrow=as.numeric(cfg$lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.size=input$line.size, line.color=input$line.color, lis.rematch = ft.reorder$ft.rematch) # All gene IDs are used.
        msg <- paste0(svg.na[i], ': no matching features are detected between data and aSVG!')
        if (is.null(gg.lis)) cat(msg, '\n')
        output$msg.shm <- ({ validate(need(!is.null(gg.lis), msg)) }) 
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
      return(list(grob.all = grob.all, gg.all = gg.all, lgd.all = lgd.all, gg.grob.lis = gg.grob.lis))
     }) # withProgress
    }) # reactive
    shm$grob.all <- gs.all()$grob.all; shm$gg.all <- gs.all()$gg.all
    shm$lgd.all <- gs.all()$lgd.all; shm$grob.gg.all <- gs.all()$gg.grob.lis
  }) # observeEvent
  # Avoid repetitive computation under input$cs.v=='gen.sel'.
  observeEvent(list(gID$geneSel), {
    
    cat('Updating all SHMs caused by selected rows ... \n')
    if.con <-  is.null(input$cs.v)|gID$geneSel[1]=='none'|input$cs.v=='All rows'
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    ID <- gID$geneSel
    shm$grob.all <- shm$gg.all <- shm$lgd.all <- shm$grob.gg.all <- NULL; gs.all <- reactive({
      if.con <- is.null(svg.df())|is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|color$col[1]=='none'
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
        gg.lis <- gg_shm(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=ID, legend.col=fil.cols, cols=color$col, tis.path=tis.path, ft.trans=input$tis, sub.title.size=15 * scale.shm, aspect.ratio = svg.df$aspect.r, legend.nrow=as.numeric(cfg$lis.par$legend['key.row', 'default']), legend.key.size=size.key, legend.text.size=8*size.key*33, line.size=input$line.size, line.color=input$line.color, lis.rematch = ft.reorder$ft.rematch) # All gene IDs are used.
        msg <- paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')
        if (is.null(gg.lis)) cat(msg, '\n')
        output$msg.shm <- ({ validate(need(!is.null(gg.lis), msg)) })
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
      return(list(grob.all = grob.all, gg.all = gg.all, lgd.all = lgd.all, gg.grob.lis = gg.grob.lis))
     }) # withProgress
    }) # reactive
    shm$grob.all <- gs.all()$grob.all; shm$gg.all <- gs.all()$gg.all
    shm$lgd.all <- gs.all()$lgd.all; shm$grob.gg.all <- gs.all()$gg.grob.lis
  }) # observeEvent
  
  # when 'color <- reactiveValues(col="none")', upon the app is launched, 'gs.new' is evaluated for 3 time. In the 1st time, 'gID$new'/'gID$all' are NULL, so 'gs.new' is NULL. In the 2nd time, 'color$col[1]=='none'' is TRUE, so NULL is returned to 'gs.new', but 'gID$new'/'gID$all' are 'HRE2'. In the third time, 'color$col[1]=='none'' is FALSE, so 'gs.new' is not NULL, but 'gID$new' is still 'HRE2', so it does not triger evaluation of 'observeEvent' and hence SHMs and legend plot are not returned upon being launched. The solution is to assign colors to 'color$col' in 'observe' upon being launched so that in the 2nd time 'gs.new' is not NULL, and no 3rd time.
  observeEvent(gs.new(), { 
    if (is.null(svg.df())|is.null(gID$new)|length(gID$new)==0|is.null(gID$all)|is.null(gs.new())) return(NULL)
    cat('Updating grobs/ggplots/legends based on new ID ... \n')
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

    }

  })
  
  # Update subtitle size through theme().
  observeEvent(list(input$title.size, input$scale.shm), {
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
    
    if (is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|is.null(svg.df())|is.null(shm$grob.all)) return(NULL)

    height <- input$height; width <- input$width
    col.n <- input$col.n;
    # validate(need(height>=0.1 & !is.na(height), 'Height should be a positive numeric !'))
    # validate(need(width>=0.1 & !is.na(width), 'Width should be a positive numeric !'))
    validate(need(col.n>=1 & as.integer(col.n)==col.n & !is.na(col.n), 'No. of columns should be a positive integer !'))

  })
  observeEvent(list(input$lgd.key.size, input$lgd.row, input$tis, input$lgd.label, input$lgd.lab.size, input$lgd.ratio1), {
    lgd.key.size <- input$lgd.key.size; lgd.row <- input$lgd.row
    lgd.label <- input$lgd.label; label.size <- input$lgd.lab.size; lgd.ratio1 <- input$lgd.ratio1
    if (is.null(shm$lgd.all)|!is.numeric(lgd.key.size)|!is.numeric(lgd.row)|is.null(lgd.label)|!is.numeric(lgd.ratio1)) return()
    cat('Adjust legend size/rows/aspect ratio ... \n')
    shm$lgd.all <- gg_lgd(gg.all=shm$lgd.all, size.key=lgd.key.size, size.text.key=NULL, row=lgd.row, sam.dat=sam(), ft.trans=input$tis, position.text.key='right', label=(lgd.label=='Yes'), label.size=label.size, aspect.ratio = lgd.ratio1) 
  })
  observeEvent(list(grob.all=shm$grob.all, gen.con=input$gen.con), {
  
    if (is.null(gID$all)|is.null(shm$grob.all)|is.null(shm$gg.all)) return()
    cat('Reordering grobs/ggplots... \n')
    na.all <- names(shm$grob.all); pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', na.all)); if (length(con)==0) return()
    na.all <- sort_gen_con(ID.sel=gID$all, na.all=na.all, con.all=con, by=input$gen.con)
    # grob1/gg.all1 are used to add/remove 2nd legend.
    shm$grob.all1 <- shm$grob.all[na.all]; shm$gg.all1 <- shm$gg.all[na.all]

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
    if.con <- is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)

  if (col.reorder$col.re=='N') return()
    if.con <-  is.null(ipt.dat$dat$dt_rows_selected)|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)
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
    lay <- input$gen.con; ID <- gID$geneSel; ncol <- input$col.n
    cat('Done! \n') # If 'cat' is the last step, NULL is returned.
    lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, lay.mat = TRUE)
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

  observe({
    lay <- lay.shm(); scale.shm <- input$scale.shm
    if (is.null(lay)|!is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    # subplot: height 300, width 250 
    # Avoid: if one column has all NAs in the layout matrix, the aspect ratio is distroyed. So only take the columns not containing all NAs.
    col.vld <- sum(unlist(lapply(seq_len(ncol(lay)), function(x) !all(is.na(lay[, x])))))
    # width/height relate to scrolling in box. 
    output$shm <- renderPlot(width = col.vld * 250 * scale.shm, height = nrow(lay) * 300 * scale.shm, { 
    if (col.reorder$col.re=='N') return()
    if.con <-  is.null(ipt.dat$dat$dt_rows_selected)|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$grob.all1)
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
    cat('Plotting spatial heatmaps... \n')
    lay <- input$gen.con; ID <- gID$geneSel; ncol <- input$col.n
    # This step is plotting.
    shm.lay <- lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, shiny=TRUE); shm <- shm.lay$shm
    # Adjust the dimension in chicken example. 
    #svg.df <- svg.df(); w.all <- h.all <- 0
    # for (i in svg.df) { w.all <- w.all+i$w.h['width']; h.all <- h.all+i$w.h['height'] }
    #shm.row <- nrow(shm.lay$lay)
    #if (length(svg.df)==1) updateNumericInput(session, inputId="height", label="Overall height:", value=as.numeric(h.all/w.all*shm.row*input$width), min=0.1, max=Inf, step=NA)
    if (input$ext!='NA') {
      validate(need(try(input$res>0), 'Resolution should be a positive numeric!'))
      validate(need(try(input$lgd.w>=0 & input$lgd.w <1), 'Legend width should be between 0 to 1!'))
      validate(need(try(input$lgd.ratio>0), 'Legend aspect ratio should be a positive numeric!'))
      cs.grob <- ggplotGrob(shm.bar())
      cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc"))
      # Legend size in downloaded SHM is reduced.
      lgd.lis <- shm$lgd.all; lgd.lis <- gg_lgd(gg.all=lgd.lis, sam.dat=sam(), ft.trans=input$tis, label=FALSE)
      lgd.lis <- gg_lgd(gg.all=lgd.lis, size.key=input$lgd.key.size*0.5, size.text.key=NULL, label.size=input$lgd.lab.size, row=input$lgd.row, sam.dat=sam(), ft.trans=input$tis, position.text.key='right', label=(input$lgd.label=='Yes'))
      if (input$lgd.w>0) {
  
        grob.lgd.lis <- lapply(lgd.lis, ggplotGrob)
        lgd.tr <- lapply(grob.lgd.lis, grobTree)
    # In 'arrangeGrob', if numbers in 'layout_matrix' are more than items in 'grobs', there is no difference. The width/height of each subplot is decided by 'widths' and 'heights'.
    lgd.arr <- arrangeGrob(grobs=lgd.tr, layout_matrix=matrix(seq_along(lgd.lis), ncol=1), widths=unit(1, "npc"), heights=unit(rep(1/length(lgd.lis)/input$lgd.ratio, length(lgd.lis)), "npc"))
    w.lgd <- (1-0.08)/(ncol+1)*input$lgd.w # Legend is reduced.
    png(paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/tmp.png')); shm1 <- grid.arrange(cs.arr, shm, lgd.arr, ncol=3, widths=unit(c(0.08-0.005, 1-0.08-w.lgd, w.lgd), 'npc')); dev.off() } else { png(paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/tmp.png')); shm1 <- grid.arrange(cs.arr, shm, ncol=2, widths=unit(c(0.08-0.005, 1-0.08), 'npc')); dev.off() }
    ggsave(paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/shm.', input$ext), plot=shm1, device=input$ext, width=input$width/72, height=input$height/72, dpi=input$res, unit='in') 
    }
    
    })

  })

#  shm_server('shm.img', ipt$all, lay.shm(), svg.df(), geneIn(), pat.all(), shm.bar(), sam(), col.reorder, gID, color)

  output$dld.shm <- downloadHandler(
    filename=function() { paste0('shm.', input$ext) },
    content=function(file) { file0 <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/shm.', input$ext); 
    cat("Downloading 'shm' from", normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '...\n')
    file.copy(file0, file, overwrite=TRUE) }
  )

  observe({
  
    ipt$fileIn; geneIn(); ipt$adj.modInpath; A(); input$p; input$cv1; input$cv2; ipt.dat$dat$dt_rows_selected; input$tis; input$gen.con  
    updateRadioButtons(session, inputId='ext', label='File type', choices=c('NA', "png", "jpg", "pdf"), selected=cfg$lis.par$shm.img['file.type', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId="ggly.but", label="Show animation", choices=c("Yes", "No"), selected=cfg$lis.par$shm.anm['show', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId="vdo.but", label="Show/update video", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['show', 'default'], inline=TRUE)

  })

  observe({
   input$vdo.key.size; input$vdo.key.row; input$vdo.val.lgd; input$tis; input$vdo.lab.size; input$vdo.res; input$vdo.itvl
   updateRadioButtons(session, inputId="vdo.but", label="Show/update video", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['show', 'default'], inline=TRUE)
  })


  output$lgd <- renderPlot(width='auto', height="auto", { # auto: no need to scroll. 
    cat('Plotting legend plot... \n')
    validate(need(try(as.integer(input$lgd.row)==input$lgd.row & input$lgd.row>0), 'Legend key rows should be a positive integer!'))
    validate(need(try(input$lgd.key.size>0&input$lgd.key.size<1), 'Legend key size should be between 0 and 1!'))
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

  output$lgd.ui <- renderUI({ 
    ns <- session$ns    
    if (is.null(input$lgdB)) return(NULL) 
    if (input$lgdB %% 2 == 1) return(NULL)
    box(title="Legend Plot", status="primary", solidHeader=TRUE, collapsible=TRUE, width = 3, 
    navbarPage('Parameters:',
    tabPanel("Basic",
    splitLayout(cellWidths=c("32%", "1%", '32%', '1%', '35%'),
    numericInput(inputId=ns('lgd.row'), label='Key rows', value=as.numeric(cfg$lis.par$legend['key.row', 'default']), min=1, max=Inf, step=1, width=150), '',
    numericInput(inputId=ns('lgd.key.size'), label='Key size', value=as.numeric(cfg$lis.par$legend['key.size', 'default']), min=0, max=1, step=0.02, width=150), ''
    # numericInput(inputId=ns('lgd.ratio1'), label='Aspect ratio', value=as.numeric(cfg$lis.par$legend['aspect.ratio', 'default']), min=0.01, max=Inf, step=0.01, width=150)
    )), # tabPanel

    tabPanel("Feature labels",
    splitLayout(cellWidths=c("43%", "8%", '43%'),
    radioButtons(inputId=ns("lgd.label"), label="Feature labels", choices=c("Yes", "No"), selected=cfg$lis.par$legend['label', 'default'], inline=TRUE), '',
    numericInput(inputId=ns('lgd.lab.size'), label='Label size', value=as.numeric(cfg$lis.par$legend['label.size', 'default']), min=0, max=Inf, step=0.5, width=150)
    )) # tabPanel
    ), # navbarPage
    uiOutput(ns('shms.o')), splitLayout(cellWidths=c("99%", "1%"), plotOutput(ns("lgd")), "")) # box

  })


  output$tran <- renderText({
    if (is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$grob.all)) return(NULL)
    if (!is.null(input$t)) validate(need(try(input$t>=0.1), 'Transition time should be at least 0.1 second!'))
  })


  observeEvent(list(fineIn=ipt$fileIn, log=log(), tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, preScale=input$preScale), { ggly_rm(); vdo_rm() })

  observeEvent(list(width.ly=input$width.ly, height.ly=input$height.ly), {
    if (dir.exists('html_shm/')) { unlink('html_shm/lib', recursive=TRUE)
      file.remove(list.files('html_shm/', '*.html$', full.names=TRUE))
    } else dir.create('html_shm/')
  })
  observeEvent(list(log=log(), tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, preScale=input$preScale, ggly.but=input$ggly.but, gID.new=gID$new), {

    if (is.null(input$ggly.but)) return() 
    if (input$ggly.but=='No') return()
    if (is.null(geneIn())|is.null(gID$new)|is.null(ipt.dat$dat$dt_rows_selected)|is.null(svg.df())|gID$geneSel[1]=="none"|is.null(shm$gg.all1)|input$ggly.but=='No') return(NULL)
    if (length(color$col=="none")==0|input$color=="") return(NULL)

    withProgress(message="Animation: ", value=0, {
    incProgress(0.25, detail="preparing frames ...") 
    gg.all <- shm$gg.all1; na <- names(gg.all)
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]; gg.all <- gg.all[na]
    for (i in seq_along(gg.all)) {

      na0 <- paste0(na[i], ".html")
      if (length(list.files('www/ggly/', na0))>0) next
      gly <- ggplotly(gg.all[[i]], tooltip='text') %>% layout(showlegend=FALSE)
      gly$sizingPolicy$padding <- 0
      cat('Animation: saving', na0, '\n')
      saveWidget(gly, na0, selfcontained=FALSE, libdir="lib")
      file.rename(na0, paste0('www/ggly/', na0))

    }
    if (!dir.exists('www/ggly/lib')) file.rename('lib', 'www/ggly/lib') else if (dir.exists('lib/')) unlink('lib', recursive=TRUE)
    })

  })

  output$sld.fm <- renderUI({
    ns <- NS(id) 
    if (is.null(shm$gg.all)|is.null(pat.all())|is.null(gID$geneSel)) return(NULL) 
    gen.con.pat <- paste0('^', pat.all(), '_\\d+$') 
    sliderInput(inputId=ns('fm'), 'Frames', min=1, max=sum(grepl(gen.con.pat, names(shm$gg.all1))), step=1, value=1, animate=animationOptions(interval=input$t*10^3, loop=FALSE, playButton=icon('play'), pauseButton=icon('pause')))
  
  })

  # As long as the variable of 'reactive' is used in the 'ui.R', changes of elements in 'reactive' would cause chain change all the way to 'ui.R'. E.g. the change in "input$ggly.but=='No'" leads to changes in 'output$ggly' and 'ui.R', not necessarily changes in 'output$ggly' call changes in 'gly.url'.
  gly.url <- reactive({ 
    
    if (is.null(input$ggly.but)) return() 
    if (is.null(shm$gg.all1)|input$ggly.but=='No'|gID$geneSel[1]=='none'|is.null(pat.all())) return(NULL)
    gg.all <- shm$gg.all1; na <- names(gg.all)
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]; na1 <- na[as.integer(input$fm)]
    na2 <- list.files('www/ggly', pattern=na1)
    if (length(na2)==0|is.na(na2)) return(NULL)
    cat('Animation: access', na2, 'path \n')
    paste0('ggly/', na2)
  
  })

  # Variables in 'observe' are accessible anywhere in the same 'observe'.
  observe({

    if (is.null(input$ggly.but)|is.null(input$fm)) return() 
    if (input$ggly.but=='No'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|color$col[1]=='none') return(NULL)
    
    gg.all <- shm$gg.all1; na <- names(gg.all)
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]; na1 <- na[as.integer(input$fm)]
    na2 <- list.files('www/ggly', pattern=paste0(na1, '\\.html$')); if (length(na2)==0) return(NULL)
    gg <- gg.all[[na1]]
    dat <- layer_data(gg); x.max <- max(dat$x); y.max <- max(dat$y)
    w <- cfg$lis.par$shm.anm['width', 'default']
    if (w!='NA') w <- as.numeric(w) else w <- NA
    h <- cfg$lis.par$shm.anm['height', 'default']
    if (h!='NA') h <- as.numeric(h) else h <- NA
    if (!is.na(w)) {
      h <- y.max/x.max*w; if (h>550) { h <- 550; w <- x.max/y.max*h }
    } else if (is.na(w) & !is.na(h)) { w <- x.max/y.max*h }
    output$tran.t <- renderUI({
      ns <- NS(id)
      numericInput(inputId=ns('t'), label='Transition time (s)', value=as.numeric(cfg$lis.par$shm.anm['transition', 'default']), min=0.1, max=Inf, step=NA, width=270)
    }) 
    output$anm.w <- renderUI({
      ns <- session$ns
      numericInput(inputId=ns('width.ly'), label='Width', value=w, min=1, max=Inf, step=NA, width=170)
    })
    output$anm.h <- renderUI({
      ns <- session$ns
      numericInput(inputId=ns('height.ly'), label='Height', value=h, min=1, max=Inf, step=NA, width=170)
    })
  output$dld.anm.but <- renderUI({ 
    ns <- session$ns
    downloadButton(ns("dld.anm"), "Download") 
  })

  })

  
  observeEvent(list(log=log(), tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, preScale=input$preScale, ggly.but=input$ggly.but, fm=input$fm), {
  
  output$ggly <- renderUI({
    if (input$ggly.but=='No'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|color$col[1]=='none') return(NULL)
    withProgress(message="Animation: ", value=0, {
    incProgress(0.75, detail="plotting ...")
    gly.url <- gly.url(); cat('Animation: plotting', gly.url, '\n')
    tags$iframe(src=gly.url, height=input$height.ly, width=input$width.ly, scrolling='auto')
  
    })

  })

  })

  anm.dld <- reactive({
    
    if (input$ggly.but=='No'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|color$col[1]=='none') return(NULL) 
    withProgress(message="Downloading animation: ", value=0, {
    incProgress(0.1, detail="in progress ...")
    gg.all <- shm$gg.all1; na <- names(gg.all)
    gg.na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]
    gg <- gg.all[gg.na]
    pro <- 0.1; for (i in seq_along(gg.na)) {
    incProgress(pro+0.2, detail=paste0('preparing ', gg.na[i], '.html...'))
    html_ly(gg=gg[i], cs.g=shm.bar(), ft.trans=input$tis, sam.uni=sam(), anm.width=input$width.ly, anm.height=input$height.ly, out.dir='.') }

   })

  })
  
  # This step leaves 'fil.na' in 'output$dld.anm' being a global variable.
  output$dld.anm <- downloadHandler( 
    # The rest code will run only after 'anm.dld()' is done.
    filename=function(){ anm.dld(); "html_shm.zip" },
    fil.na <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/html_shm.zip'),
    content=function(fil.na){ cat('Downloading animation... \n'); zip(fil.na, 'html_shm/') }
  )

  output$video.dim <- renderUI({
    ns <- session$ns
    selectInput(ns("vdo.dim"), label="Fixed dimension", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected=cfg$lis.par$shm.video['dimension', 'default'], width=110)
  })

  output$ffm <- renderText({
    ffm <- tryCatch({ test_ffm() }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
    if (grepl('error|warning', ffm)) paste("<span style=\"color:red\">Error: \"ffmpeg\" is not detected!\"</span>")
  })

  observeEvent(list(log=log(), tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, preScale=input$preScale, vdo.but=input$vdo.but, vdo.dim=input$vdo.dim, vdo.itvl=input$vdo.itvl, vdo.height=input$vdo.height, vdo.width=input$vdo.width, vdo.res=input$vdo.res, vdo.val.lgd=input$'vdo.val.lgd', input$vdo.label), {

    cat('Making video ... \n')
    if (is.null(input$vdo.but)) return(NULL) 
    if (input$vdo.but=='No'|is.null(pat.all())) return(NULL)
    if (is.null(svg.df())|is.null(geneIn())|is.null(ipt.dat$dat$dt_rows_selected)|color$col[1]=='none') return(NULL)
    validate(need(try(!is.na(input$vdo.itvl)&input$vdo.itvl>0), 'Transition time should be a positive numeric!'))
    validate(need(try(!is.na(input$vdo.height)&input$vdo.height>0&input$vdo.height<=0.99), 'Height should be between 0.1 and 0.99!'))
    validate(need(try(!is.na(input$vdo.width)&input$vdo.width>0&input$vdo.width<=0.92), 'Width should be between 0.1 and 0.92!'))
    validate(need(try(!is.na(input$vdo.res)&input$vdo.res>=1&input$vdo.res<=700), 'Resolution should be between 1 and 700!'))
    
    withProgress(message="Video: ", value=0, {
    incProgress(0.75, detail="in progress ...")
    gg.all <- shm$gg.all1; na <- names(gg.all)
    pat <- paste0('^', pat.all(), '_\\d+$'); na <- na[grepl(pat, na)]
    gg.all1 <- gg.all[na]
    res <- input$vdo.res; dim <- input$vdo.dim
    if (dim %in% c('1280x800', '1280x1024', '1280x720')&res>450) res <- 450
    if (dim=='1920x1080'&res>300) res <- 300
    # selectInput("vdo.dim", label="Fixed dimension:", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected='640x480', width=110)
    vdo <- video(gg=gg.all1, cs.g=shm.bar(), sam.uni=sam(), ft.trans=input$tis, lgd.key.size=input$vdo.key.size, lgd.text.size=NULL, position.text.key='right', legend.value.vdo=(input$'vdo.val.lgd'=='Yes'), label=(input$vdo.label=='Yes'), label.size=input$vdo.lab.size, sub.title.size=8, bar.value.size=6, lgd.row=input$vdo.key.row, width=input$vdo.width, height=input$vdo.height, video.dim=dim, interval=input$vdo.itvl, res=res, out.dir='./www/video')
    cat('Done! \n'); if (is.null(vdo)) return()
    cat('Presenting video ... \n')
    incProgress(0.95, detail="Presenting video...")
    w.h <- as.numeric(strsplit(input$vdo.dim, 'x')[[1]])
    output$video <-renderUI({ tags$video(id="video", type="video/mp4", src="video/shm.mp4", width=w.h[1], height=w.h[2], controls="controls") })
    cat('Done! \n')
    })

  })

 output$shm.ui <- renderUI({
    ns <- session$ns
    column(width = ifelse(input$lgdB %% 2 == 0, 9, 12), 
      # boxPad(color = NULL, title = NULL, solidHeader = FALSE, 
    
    tabsetPanel(type = "pills", id=NULL, selected="shm1",
  
      tabPanel(title="Image", value='shm1',  
      navbarPage('Parameters:',
      tabPanel("Basic",
      fluidRow(splitLayout(cellWidths=c('11%', '1%', '8%', '1%', '10%', '1%', '10%', '1%', '10%', '1%', '13%'), 
      div(style = "margin-top: 23px;",
      actionButton(ns("fs"), "Full screen", onclick = "openFullscreen(document.getElementById('shmAll-shm'))")), '',
      div(title = 'Number of columns for the subplots.', 
      numericInput(inputId=ns('col.n'), label='Columns', value=as.numeric(cfg$lis.par$shm.img['columns', 'default']), min=1, max=Inf, step=1, width=150)), '',
      selectInput(inputId = ns("gen.con"), label = "Display by", choices = c("Gene"="gene", "Condition"="con", "None"="none"), selected = cfg$lis.par$shm.img['display.by', 'default'], multiple = FALSE, width = NULL), '',
      numericInput(inputId=ns('scale.shm'), label='Scale plots', value=as.numeric(cfg$lis.par$shm.img['scale.plots', 'default']), min=0.1, max=Inf, step=0.1, width=150), '',
      numericInput(inputId=ns('title.size'), label='Title size', value=as.numeric(cfg$lis.par$shm.img['title.size', 'default']), min=0.1, max=Inf, step=0.5, width=150), '',
      div(style = "margin-top: 23px;",
      dropdownButton(inputId=ns('dropdown'), label='Color key', circle=FALSE, icon= icon("sliders"), status='primary', inline=FALSE, width=250,
      fluidRow(splitLayout(cellWidths=c('1%', '60%', '35%'), '', textInput(ns("color"), "Color scheme", cfg$lis.par$shm.img['color', 'default'], placeholder=paste0('Eg: ', cfg$lis.par$shm.img['color', 'default']), width=200),
      actionButton(ns("col.but"), "Confirm", icon=NULL, style = "margin-top: 24px;"))), 
      radioButtons(inputId=ns('cs.v'), label='Color key based on', choices=c("Selected rows", "All rows"), selected=cfg$lis.par$shm.img['color.scale', 'default'], inline=TRUE)
      )) 
      )), # fluidRow
      textOutput(ns('h.w.c')), textOutput(ns('msg.col')),
      fluidRow(splitLayout(cellWidths=c('100%'),  checkboxGroupInput(inputId=ns("tis"), label="Select features to be transparent", choices='', selected='', inline=TRUE)))
      ), # tabPanel
      tabPanel("Value legend",
      column(9, offset=0, style='padding-left:50px; padding-right:80px; padding-top:0px; padding-bottom:5px',
      fluidRow(splitLayout(cellWidths=c('1%', '32%', '1%', '13%', '1%', '17%', '1%', '16%', '1%', '28%'), '', 
      actionButton(ns("val.lgd"), "Add/Remove", icon=icon("refresh"), style = "margin-top: 24px;"), '',  
      numericInput(inputId=ns('val.lgd.row'), label='Rows', value=as.numeric(cfg$lis.par$shm.img['value.legend.rows', 'default']), min=1, max=Inf, step=1, width=150), '',
      numericInput(inputId=ns('val.lgd.key'), label='Key size', value=as.numeric(cfg$lis.par$shm.img['value.legend.key', 'default']), min=0.0001, max=1, step=0.01, width=150), '',
      numericInput(inputId=ns('val.lgd.text'), label='Text size', value=as.numeric(cfg$lis.par$shm.img['value.legend.text', 'default']), min=0.0001, max=Inf, step=1, width=140), '',
      radioButtons(inputId=ns('val.lgd.feat'), label='Include features', choices=c('No', 'Yes'), selected=cfg$lis.par$shm.img['include.feature', 'default'], inline=TRUE)
      ))) # column
      ), # tabPanel
      tabPanel("Shape outline",
      splitLayout(cellWidths=c('1%', '15%', '1%', '13%'), '', 
      selectInput(ns('line.color'), label='Line color', choices=c('grey70', 'black', 'red', 'green', 'blue'), selected=cfg$lis.par$shm.img['line.color', 'default']), '', 
      numericInput(inputId=ns('line.size'), label='Line size', value=as.numeric(cfg$lis.par$shm.img['line.size', 'default']), min=0.05, max=Inf, step=0.05, width=150) 
      )), # tabPanel
     tabPanel("Download",

#     tags$div(title="Download the spatial heatmaps and legend plot.",
      # h1(strong("Download paramters:"), style = "font-size:20px;"),
      fluidRow(splitLayout(cellWidths=c('30%', '1%', '18%', '1%', '13%', '1%', '19%', '1%', '16%'),
      tags$div(class='tp', span(class='tpt', strong('Select a file type to download.')),
      radioButtons(inputId=ns('ext'), label='File type', choices=c('NA', "png", "jpg", "pdf"), selected=cfg$lis.par$shm.img['file.type', 'default'], inline=TRUE)), '', 
      numericInput(inputId=ns('res'), label='Resolution (dpi)', value=as.numeric(cfg$lis.par$shm.img['dpi', 'default']), min=10, max=Inf, step=10, width=150), ''
      )), # fluidRow

      fluidRow(splitLayout(cellWidths=c('18%', '1%', '25%', '1%', '18%'), 
      numericInput(inputId=ns('lgd.w'), label='Legend plot width', value=as.numeric(cfg$lis.par$shm.img['legend.width', 'default']), min=0, max=1, step=0.1, width=150), '',
      tags$div(title="Adjust the aspect ratio of legend plot.",
      numericInput(inputId=ns('lgd.ratio'), label='Legend plot aspect ratio', value=as.numeric(cfg$lis.par$shm.img['legend.aspect.ratio', 'default']), min=0.0001, max=Inf, step=0.1, width=140)
      ), '', downloadButton(ns("dld.shm"), "Download", style = "margin-top: 24px;")
      )) # fluidRow 
      ), # tabPanel 
     
     # navbarMenu("More", # Create dropdown menu.
     tabPanel("Relative size",
       tags$div(title="Adjust the relative size in multiple aSVGs.",
       numericInput(inputId=ns('relaSize'), label='Relative sizes', value=as.numeric(cfg$lis.par$shm.img['relative.size', 'default']), min=0.01, max=Inf, step=0.1, width=140))
      ), # tabPanel

      tabPanel(title="Re-match features", value='rematch',
        column(12, fluidRow(splitLayout(cellWidths=c('0.2%', "40%", '20%', "30%"), '',
          uiOutput(ns('svg'), style = 'margin-left:-5px'), '', 
          actionButton(ns("match"), "Confirm re-matching", icon=icon("refresh"), style="color: #fff; background-color:#3498DB;border-color: #2e6da4;margin-top: 24px;")
        ))), verbatimTextOutput(ns('msg.match')),
        column(12, uiOutput(ns('ft.match')))
        # verbatimTextOutput('ft.re'),
      # fluidRow(splitLayout(cellWidths=c("1%", "49%", "49%", "1%"), "", uiOutput('sam'), uiOutput('sf'), "")),
      )
      #) # navbarMenu
      ), # navbarPage 
                          
      verbatimTextOutput(ns('msg.shm')),
      fluidRow(splitLayout(id = 'barSHM', cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput(ns("bar1")),
      div(style = 'overflow-y:scroll;height:400px;',
      plotOutput(ns("shm"), height='100%', width = '100%')
      ), ""))
      ), # tabPanel 

      tabPanel(title='Animation', value='shm2', 
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '16%', '1%', '12%', '1%', '12%', '1%', '13%'), '', 
      radioButtons(inputId=ns("ggly.but"), label="Show animation", choices=c("Yes", "No"), selected=as.numeric(cfg$lis.par$shm.img['show', 'default']), inline=TRUE), '',
      uiOutput(ns('tran.t')), '',
      uiOutput(ns('anm.h')), '', uiOutput(ns('anm.w')), '', uiOutput(ns('dld.anm.but'))
      )), textOutput(ns('tran')), uiOutput(ns('sld.fm')),
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput(ns("bar2")), htmlOutput(ns("ggly")), ""))
      ),

      tabPanel(title='Video', value='shm3',
      navbarPage('Parameters:',
      tabPanel("Dimension",
      fluidRow(splitLayout(cellWidths=c('1%', '8%', '1%', '8%', '1%', '14%'), '', 
      numericInput(inputId=ns('vdo.height'), label='Height', value=as.numeric(cfg$lis.par$shm.video['height', 'default']), min=0.1, max=0.99, step=0.1, width=270), '',
      numericInput(inputId=ns('vdo.width'), label='Width', value=as.numeric(cfg$lis.par$shm.video['width', 'default']), min=0.1, max=0.92, step=0.1, width=270), '', uiOutput(ns('video.dim'))
      )) # fluidRow
      ),
      tabPanel("Key",
      fluidRow(splitLayout(cellWidths=c('1%', '8%', '1%', '10%', '1%', '13%'), '', 
      numericInput(inputId=ns('vdo.key.row'), label='Key rows', value=as.numeric(cfg$lis.par$shm.video['key.rows', 'default']), min=1, max=Inf, step=1, width=270), '',
      numericInput(inputId=ns('vdo.key.size'), label='Key size', value=as.numeric(cfg$lis.par$shm.video['key.size', 'default']), min=0.01, max=Inf, step=0.1, width=270), '',
      radioButtons(inputId=ns("vdo.val.lgd"), label="Key value", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['value.legend', 'default'], inline=TRUE) 
      ))),
      tabPanel("Label",
      splitLayout(cellWidths=c("15%", "1%", '10%'),
      radioButtons(inputId=ns("vdo.label"), label="Feature label", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['feature.label', 'default'], inline=TRUE), '',
      numericInput(inputId=ns('vdo.lab.size'), label='Label size', value=as.numeric(cfg$lis.par$shm.video['label.size', 'default']), min=0, max=Inf, step=0.5, width=150)
      )), # tabPanel
      tabPanel("More",
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '13%'), '', 
      numericInput(inputId=ns('vdo.itvl'), label='Transition time (s)', value=as.numeric(cfg$lis.par$shm.video['transition', 'default']), min=0.1, max=Inf, step=1, width=270), '',
      numericInput(inputId=ns('vdo.res'), label='Resolution (dpi)', value=as.numeric(cfg$lis.par$shm.video['dpi', 'default']), min=1, max=1000, step=5, width=270)
      )) # fluidRow
      )
      ), # navbarPage
      textOutput(ns('tran.vdo')), htmlOutput(ns('ffm')), 
      radioButtons(inputId=ns("vdo.but"), label="Show/update video", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['show', 'default'], inline=TRUE),
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput(ns("bar3")), uiOutput(ns('video')), ""))
      ) # tabPanel
      
    )
      ) # column 
  
  })

# addPopover(session=session, id="height", title="", content="Check 'Yes' to preserve the aspect ratio defined in the aSVG file.", placement = "bottom", trigger = "hover", options = NULL)  
  
 output$shms.o <- renderUI({
    ns <- session$ns 
    if (is.null(svg.path1())) return(NULL)
    if (length(svg.path1()$svg.na)==1) return(NULL)
    svg.pa <- svg.path1()[['svg.na']]
    selectInput(ns('shms.in'), label='Select plots', choices=svg.pa, selected=svg.pa[1])
  })

  observe({
    ipt$fileIn; ipt$geneInpath; lis.par <- cfg$lis.par
    updateRadioButtons(session, inputId='cs.v', label='Color key based on', choices=c("Selected rows", "All rows"), selected=cfg$lis.par$shm.img['color.scale', 'default'], inline=TRUE)
    # updateNumericInput(session, inputId="height", label="Overall height", value=as.numeric(cfg$lis.par$shm.img['height', 'default']), min=0.1, max=Inf, step=NA)
    # updateNumericInput(session, inputId="width", label="Overall width", value=as.numeric(cfg$lis.par$shm.img['width', 'default']), min=0.1, max=Inf, step=NA)
    updateNumericInput(session, inputId="col.n", label="Columns", value=as.numeric(cfg$lis.par$shm.img['columns', 'default']), min=1, max=Inf, step=1)
  })
  return(list(gID = gID))
})} # shm_server

shm.mod.lis <- shm_server('shmAll', ipt, cfg, dat.mod.lis)

network_server <- function(id, ipt, cfg, dat.mod.lis, shm.mod.lis) {
  
  ipt.dat$dat <- dat.mod.lis$ipt.dat; A <- dat.mod.lis$A
  P <- dat.mod.lis$P; CV1 <- dat.mod.lis$CV1
  CV2 <- dat.mod.lis$CV2
  gID <- shm.mod.lis$gID; geneIn <- dat.mod.lis$geneIn
  moduleServer(id, function(input, output, session) {

  observe({
    geneIn(); ipt$adj.modInpath; input$A; input$P; input$CV1
    input$CV2; input$min.size; input$net.type
    input$measure; input$cor.abs; input$thr; input$mhm.v
    updateRadioButtons(session, "mat.scale", "Scale by: ", c("No", "Column", "Row"), "Row", inline=TRUE)
  })

  observe({  
    ipt$fileIn; geneIn(); input$adj.modInpath; input$A; input$P; input$CV1; input$CV2; ipt.dat$dat$dt_rows_selected
    updateActionButton(session, inputId='mhm.but', label='Update', icon=icon("refresh"))
    #updateRadioButtons(session, inputId="mhm.but", label="Show plot:", choices=c("Yes", "No"), selected=cfg$lis.par$mhm['show', 'default'], inline=TRUE)
  })
# Avoid unnecessay correlation/distance computation if geneIn() updates due to column reordering. For example, correlation/distance, extracting coordinates, which only depend on the df.aggr.
  df.net <- reactive({
    if (is.null(geneIn())) return()
    return(list(df.aggr = geneIn()[['df.aggr']], df.met = geneIn()[['df.met']]))
  })
  # Calculate whole correlation or distance matrix.
  cor.dis <- reactive({
    cat('Correlation/distance matrix ... \n')
    if (ipt$fileIn %in% cfg$na.cus & is.null(input$svgInpath1) & is.null(input$svgInpath2)) return()
    if (is.null(input$mhm.but)) return() 
    if (is.null(df.net()$df.aggr)|input$mhm.but=='No') return()
    if ((ipt$fileIn %in% cfg$na.cus & is.null(df.net()))|ipt$fileIn=="none") return(NULL)
    withProgress(message="Compute similarity/distance matrix: ", value = 0, {
      incProgress(0.5, detail="please wait ...")
      gene <- df.net()$df.aggr
      # Too many genes may crash the app.
      if (nrow(gene)>15000 & input$mhm.but==0) return()
      if (input$measure=='correlation') {
        m <- cor(x=t(gene)); cat('Done! \n')
        if (input$cor.abs==TRUE) { m <- abs(m) }; return(m)
      } else if (input$measure=='distance') { return(-as.matrix(dist(x=gene))) }
    })
  })

  # Subset nearest neighbours for target genes based on correlation or distance matrix.
  submat <- reactive({ 
    cat('Subsetting nearest neighbors ... \n')
    if (ipt$fileIn=="none") return()
    if (is.null(input$mhm.but)) return()
    if (is.null(cor.dis())|input$mhm.but=='No') return()
    gene <- df.net()$df.aggr; rna <- rownames(gene)
    gen.tar<- gID$geneSel; mat <- cor.dis()
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
      if (input$measure=='distance' & input$thr=='v') arg['v'] <- -arg[['v']]
      if (!all(gen.tar %in% rownames(mat))) return()    
      validate(need(try(ncol(gene)>4), 'The "sample__condition" variables in the Data Matrix are less than 5, so no coexpression analysis is applied!'))
      gen.na <- do.call(sub_na, c(mat=list(mat), ID=list(gen.tar), arg)); cat('Done! \n')
      if (any(is.na(gen.na))) return() 
      validate(need(try(length(gen.na)>=2), paste0('Only ', gen.na, ' selected!'))); return(gene[gen.na, ])
    })
  })
  mhm <- reactiveValues(hm=NULL)
  # Plot matrix heatmap.
  observe({
    cat('Initial matrix heatmap ... \n')
    if (is.null(input$mhm.but)) return() # Matrix heatmap sections is removed.
    if (input$mhm.but!=0) return(); if (is.null(submat())) return()
    gene <- df.net()$df.aggr; rna <- rownames(gene)
    gen.tar <- gID$geneSel; if (length(gen.tar)>1) return()
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="Plotting ...")
      if (input$mat.scale=="Column") scale.hm <- 'column' else if (input$mat.scale=="Row") scale.hm <- 'row' else scale.hm <- 'no'
      mhm$hm <- matrix_hm(ID=gen.tar, data=submat(), scale=scale.hm, main='Target Genes and Their Nearest Neighbours', title.size=10, static=FALSE)
      cat('Done! \n')
    })
  })
  hmly <- eventReactive(input$mhm.but, {
    cat('Matrix heatmap ... \n')
    #if (is.null(submat())|input$mhm.but=='No') return()
    if (is.null(submat())) return()
    gene <- df.net()$df.aggr; rna <- rownames(gene)
    gen.tar<- gID$geneSel
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="plotting ...")
      if (input$mat.scale=="Column") scale.hm <- 'column' else if (input$mat.scale=="Row") scale.hm <- 'row' else scale.hm <- 'no'  
      cat('Done!')
      matrix_hm(ID=gen.tar, data=submat(), scale=scale.hm, main='Target Genes and Their Nearest Neighbours', title.size=10, static=FALSE)
    })
  })

  output$HMly <- renderPlotly({ 
    if (is.null(ipt.dat$dat$dt_rows_selected)) return()
    if (is.null(gID$geneSel)|is.null(submat())) return()
    if (gID$geneSel[1]=='none'|is.na(gID$geneSel[1])) return()
    if (input$mhm.but!=0) hmly() else if (input$mhm.but==0) mhm$hm else return() 
  })

  adj.mod <- reactive({ 
    if (ipt$fileIn=="customComputedData" & !is.null(input$adj.modInpath)) {

      name <- input$adj.modInpath$name; path <- input$adj.modInpath$datapath
      path1 <- path[name=="adj.txt"]; path2 <- path[name=="mod.txt"]

      withProgress(message="Loading: ", value = 0, {
        incProgress(0.5, detail="adjacency matrix and module definition ...")
        adj <- fread(path1, sep="\t", header=TRUE, fill=TRUE); c.na <- colnames(adj)[-ncol(adj)]
        r.na <- as.data.frame(adj[, 1])[, 1];  adj <- as.data.frame(adj)[, -1] 
        rownames(adj) <- r.na; colnames(adj) <- c.na

        mcol <- fread(path2, sep="\t", header=TRUE, fill=TRUE); c.na <- colnames(mcol)[-ncol(mcol)]
        r.na <- as.data.frame(mcol[, 1])[, 1]; mcol <- as.data.frame(mcol)[, -1] 
        rownames(mcol) <- r.na; colnames(mcol) <- c.na

      }); return(list(adj=adj, mcol=mcol))

    }

  })

    #gene <- geneIn()[["df.aggr.tran"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    adj.mods <- reactiveValues(lis=NULL)
    observe({
      cat('Initial adjacency matrix and modules ...\n')
      if (ipt$fileIn=="none") return()
      if (is.null(input$cpt.nw)) return() # Network section is removed.
      if (is.null(submat())|input$cpt.nw!=0|length(gID$geneSel)>1) return()

    if (ipt$fileIn=="customData"|any(ipt$fileIn %in% cfg$na.def)) {

      gene <- df.net()$df.aggr; if (is.null(gene)) return()
      type <- input$net.type; sft <- if (type=='distance') 1 else 6
      withProgress(message="Computing: ", value = 0, {
        incProgress(0.3, detail="adjacency matrix ...")
        incProgress(0.5, detail="topological overlap matrix ...")
        incProgress(0.1, detail="dynamic tree cutting ...")
        adj.mods$lis <- adj_mod(data=submat(), type=type, minSize=input$min.size, dir=NULL)
        cat('Done! \n')
      })
    }
    })

  # er <- eventReactive(exp, {}). If its reactive value "er()" is called before eventReactive is triggered, the code execution stops where "er()" is called.
  observeEvent(input$cpt.nw, {
    cat('Adjacency and modules ... \n')
    if (ipt$fileIn %in% cfg$na.cus & is.null(input$svgInpath1) & is.null(input$svgInpath2)) return()
    #gene <- geneIn()[["df.aggr.tran"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    if (is.null(submat())|input$cpt.nw==0) return()
    if (ipt$fileIn=="customData"|any(ipt$fileIn %in% cfg$na.def)) {

      gene <- df.net()$df.aggr; if (is.null(gene)) return()
      type <- input$net.type; sft <- if (type=='distance') 1 else 6
      withProgress(message="Computing: ", value = 0, {
        incProgress(0.3, detail="adjacency matrix ...")
        incProgress(0.5, detail="topological overlap matrix ...")
        incProgress(0.1, detail="dynamic tree cutting ...")
        adj.mods$lis <- adj_mod(data=submat(), type=type, minSize=input$min.size)
        cat('Done! \n')
      })
    }
  })

  observe({
    if (is.null(geneIn())) return(NULL)
    r.na <- rownames(geneIn()[["df.aggr.tran"]]); gens.sel <- r.na[ipt.dat$dat$dt_rows_selected]
    if (length(gens.sel)==0) return()
    updateSelectInput(session, inputId="gen.sel", label="Select a target gene:", choices=c("None", gens.sel), selected=gens.sel[1])
  })
  observe({ 
    input$gen.sel; input$measure; input$cor.abs; input$thr; input$mhm.v
    updateSelectInput(session, 'ds', "Module splitting sensitivity level:", 3:2, selected=cfg$lis.par$network['ds', 'default'])
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

  visNet <- reactive({
    input$cpt.nw; if (ipt$fileIn=="none") return()
    if (ipt$fileIn=='customComputedData' & is.null(df.net())) return()
    # if (input$adj.in=="None") return(NULL)
    if (ipt$fileIn=="customComputedData") { adj <- adj.mod()[['adj']]; mods <- adj.mod()[['mcol']] } else if (ipt$fileIn=="customData"|ipt$fileIn %in% cfg$na.def) { 
      adj <- adj.mods$lis[['adj']]; mods <- adj.mods$lis[['mod']]
    }
    if (ipt$fileIn=='customComputedData') gene <- df.net()$df.aggr else gene <- submat()
    if (is.null(input$gen.sel)) return() # Matrix heatmap section is removed.
    if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    lab <- mods[, input$ds][rownames(gene)==input$gen.sel]
    validate(need(try(length(lab)==1 & !is.na(lab) & nrow(mods)==nrow(gene)), 'Click "Update" to display new network!'))
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
        
      }; meta <- df.net()$df.met
      if (ncol(meta) > 0) node <- cbind(node, title=meta[node$id, ], borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      if (ncol(meta) == 0) node <- cbind(node, borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      net.lis <- list(node=node, link=link1, adjs=adjs, lins=lins)

    }); net.lis

  })
  # The order of reactive expression matters so "updateSelectInput" of "adj.in" should be after visNet().
  observe({
 
    if (ipt$fileIn=="none") return()
    geneIn(); gID$geneSel; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$P; input$CV1; input$CV2; input$min.size; input$net.type
    input$gen.sel; input$measure; input$cor.abs; input$thr; input$mhm.v; input$cpt.nw
     #if ((input$adj.in==1 & is.null(visNet()[["adjs1"]]))|(input$cpt.nw!=cfg$lis.par$network['max.edges', 'default'] & is.null(visNet()[["adjs1"]]))) { updateSelectInput(session, "adj.in", "Adjacency threshold:", sort(seq(0, 1, 0.002), decreasing=TRUE), visNet()[["adjs"]]) } else if (!is.null(visNet()[["adjs1"]])) updateSelectInput(session, "adj.in", "Adjacency threshold:", sort(seq(0, 1, 0.002), decreasing=TRUE), visNet()[["adjs1"]])
    lins <- visNet()[["lins"]]
    if (is.null(input$adj.in)) return() # Network section is removed. 
    if (input$adj.in==1|is.null(lins)|is.numeric(lins)) updateSelectInput(session, "adj.in", "Adjacency threshold (the  smaller, the more edges):", sort(seq(0, 1, 0.002), decreasing=TRUE), as.numeric(visNet()[["adjs"]])) 
  
  })
  output$bar.net <- renderPlot({ 
    cat('Network bar ... \n')
    #if (input$adj.in=="None"|input$cpt.nw=="No") return(NULL)
    if (input$adj.in=="None") return(NULL)
    if (length(color.net$col.net=="none")==0) return(NULL)
    gene <- df.net()$df.aggr; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    if(input$col.but.net==0) color.net$col.net <- colorRampPalette(col_sep(cfg$lis.par$network['color', 'default']))(len.cs.net) # color.net$col.net is changed alse outside renderPlot, since it is a reactive value.

      withProgress(message="Color scale: ", value = 0, {
      incProgress(0.25, detail="preparing data, please wait ...")
      incProgress(0.75, detail="plotting, please wait ...")
      node <- visNet()[["node"]]; if (is.null(node)) return()
      node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
      cs.net <- col_bar(geneV=v.net, cols=color.net$col.net, width=1); 
      cat('Done! \n'); return(cs.net) # '((max(v.net)-min(v.net))/len.cs.net)*0.7' avoids bar overlap.

      })

  })

  observeEvent(visNet(), {

    output$edge <- renderUI({ 
      cat('Remaining edges ... \n')
      if (input$adj.in=="None"|is.null(visNet())) return(NULL)
      if (ipt$fileIn=="none"|(ipt$fileIn=="customData" & is.null(geneIn()))|
      input$gen.sel=="None") return(NULL)
      span(style = "color:black;font-weight:NULL;", HTML(paste0("Remaining edges: ", dim((visNet()[["link"]]))[1])))
      cat('Done! \n')
    })

  })
  vis.net <- reactive({ 
    
    cat('Network ... \n')
    #if (input$adj.in=="None"|input$cpt.nw=="No") return(NULL)
    if (input$adj.in=="None") return(NULL)
    gene <- df.net()$df.aggr; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    withProgress(message="Network:", value=0.5, {
    incProgress(0.3, detail="prepare for plotting ...")
    # Match colours with gene connectivity by approximation.
    node <- visNet()[["node"]]; if (is.null(node)) return() 
    node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
    col.nod <- NULL; for (i in node$value) {

      ab <- abs(i-v.net); col.nod <- c(col.nod, color.net$col.net[which(ab==min(ab))[1]])

    }; node$color <- col.nod
    cat('Done! \n')
    visNetwork(node, visNet()[["link"]], height="300px", width="100%", background="", main=paste0("Network Module Containing ", input$gen.sel), submain="", footer= "") %>% visIgraphLayout(physics=FALSE, smooth=TRUE) %>% visOptions(highlightNearest=list(enabled=TRUE, hover=TRUE), nodesIdSelection=TRUE)

    })
    
  })

  output$vis <- renderVisNetwork({

    if (is.null(ipt.dat$dat$dt_rows_selected)) return()
    if (ipt$fileIn=="none"|is.null(vis.net())) return(NULL)
    # if (input$cpt.nw=="No") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="plotting ...")
      cat('Rendering network...\n'); vis.net()

    })

  })

})

}
network_server('net', ipt, cfg, dat.mod.lis, shm.mod.lis)


deg_server <- function(id, ipt, cfg, dat.mod.lis, shm.mod.lis) {

  ipt.dat$dat <- dat.mod.lis$ipt.dat; A <- dat.mod.lis$A
  P <- dat.mod.lis$P; CV1 <- dat.mod.lis$CV1
  CV2 <- dat.mod.lis$CV2; gID <- shm.mod.lis$gID
  geneIn <- dat.mod.lis$geneIn; geneIn1 <- dat.mod.lis$geneIn1
  moduleServer(id, function(input, output, session) {
                   
    output$dt.rep <- renderDataTable({
    cat('Presenting data matrix (DEG) ... \n')
    if (is.null(geneIn())) return()
    if ((ipt$fileIn %in% cfg$na.cus & is.null(geneIn()))|ipt$fileIn=="none") return(NULL)

    withProgress(message="Data table: ", value = 0, {

      incProgress(0.5, detail="displaying, please wait ...")
      if (ipt$fileIn!="none") {
      gene <- geneIn0(); df.met <- gene[["df.met"]][, , drop=FALSE]
      gene.rep <- gene[["df.rep"]][, , drop=FALSE]
      if (nrow(df.met) > 0) gene.dt <- cbind.data.frame(df.met, gene.rep, stringsAsFactors=FALSE) else gene.dt <- gene.rep
   }
   #if (is.null(sear$id)) sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) else sel <- sear$id
   #if (length(sel)==1 & sel[1]==as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) & nrow(gene.dt)>1) sel <- sel else if (nrow(gene.dt)==1)  sel <- 1 else if (length(sel)>1) sel <- seq_along(sel)

   d.tab <- datatable(gene.dt, selection='none',
   filter="top", extensions=c('Scroller'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE, searchHighlight=FALSE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs = list(list(targets = c(1), render = JS("$.fn.dataTable.render.ellipsis(5, false)")))), 
   class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
   formatRound(colnames(geneIn()[["df.aggr.tran"]]), 2)
   cat('Done! \n'); d.tab
    })

  })

  ssg.dat <- reactive({
    # if (input$fileIn!="customData") return(NULL)
    cat('DEG: SummarizedExperiment, features, factors ... \n')
    if (is.null(geneIn0())) return(NULL)
    gen.lis <- geneIn0()
    df.rep <- as.matrix(gen.lis[['df.rep']])
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
  cna.sel <- NULL
  for (i in cons) { cna.sel <- c(cna.sel, paste0(sams, '__', i)) }
  se <- se[, colnames(se) %in% cna.sel]; return(se)
}

  se <- reactive({
    cat('Subsetting SE with input features/factors ... \n')
    # if (input$fileIn!="customData") return(NULL)
    dat <- ssg.dat(); sam.con <- input$sam.con 
    if (is.null(geneIn0())|is.null(dat)|is.null(sam.con)) return(NULL)
    sam <- input$ssg.sam; con <- input$ssg.con; se <- dat$se
    if (is.null(sam)|is.null(con)) return() 
    if ('all' %in% sam) sam <- unique(dat$sams)
    if ('all' %in% con) con <- unique(dat$cons)
    se <- sub_se(se, sam, con)
    if (sam.con=='feature') fct <- 'sample' else if (sam.con=='factor') fct <- 'condition' else if (sam.con=='feature__factor') fct <- 'samCon'
    cat('Done! \n')
    return(list(se=se, fct=fct))
  })

  se.nor.rok.dis <- reactive({
    cat('Normalizing data for ROKU/distinct ... \n')
    # if (input$fileIn!="customData") return(NULL)
    dat <- ssg.dat() 
    if (is.null(geneIn0())|is.null(dat)) return(NULL)
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
    # Reference, query.
    if (sam.con %in% c('feature', 'factor')) {
      if (sam.con == 'feature') {  
        under <- con; ref <- setdiff(sam, tar)
      } else if (sam.con == 'factor') {
        under <- sam; ref <- setdiff(con, tar)
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
    d.tab <- datatable(vs, selection = 'none', extensions='Scroller', plugins = "ellipsis", class='cell-border strip hover', options = list(dom = 't', scrollX = TRUE)); cat('Done! \n')
    d.tab
  })
  edg0 <- reactive({
    cat('edgeR all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (is.null(geneIn0())|is.null(se())) return(NULL)
    withProgress(message="edgeR: ", value=0, {
    incProgress(0.5, detail="in progress ...")
    edg <- edgeR(se=se()$se, method.norm=sub('.*-', '', input$edg.lim.nor), com.factor=se()$fct, method.adjust='BH', return.all=TRUE); cat('Done! \n')
    edg
    })
  })
  edg <- eventReactive(input$ssg.update, {
    cat('edgeR log2/fc ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (!'edgeR' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn0())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(edg0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct]) 
    up.dn <- up_dn(sam.all=lvl, df.all=edg0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='logFC', fdr.na='FDR'); cat('Done! \n')
    up.dn
  })

  dsq0 <- reactive({
    cat('DESeq2 all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (is.null(geneIn0())|is.null(se())) return(NULL)
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
    if (is.null(geneIn0())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(dsq0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct])
    up.dn <- up_dn(sam.all=lvl, df.all=dsq0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='log2FoldChange', fdr.na='padj'); cat('Done! \n'); up.dn
  })

  lim0 <- reactive({
    cat('limma all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    if (is.null(geneIn0())|is.null(se())) return(NULL)
    withProgress(message="limma: ", value=0, {
    incProgress(0.5, detail="in progress ...") 
    lim <- limma(se=se()$se, method.norm=sub('.*-', '', input$edg.lim.nor), m.array=FALSE, com.factor=se()$fct, method.adjust='BH', return.all=TRUE); cat('Done! \n'); lim
    })
  })
  lim <- eventReactive(input$ssg.update, {
    cat('limma log2/fc ... \n')
    # if (input$fileIn!="customData") return(NULL)
    # if (!'limma' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn0())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(lim0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct])
    up.dn <- up_dn(sam.all=lvl, df.all=lim0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='logFC', fdr.na='adj.P.Val'); cat('Done! \n'); up.dn
  })

  rok0 <- reactive({
    cat('ROKU all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    lis <- se.nor.rok.dis()
    if (is.null(geneIn0())|is.null(se())|is.null(lis)) return(NULL)
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
    if (is.null(geneIn0())|is.null(se())|is.null(rok0())) return(NULL)
    cat('ROKU up/down ... \n')
    up_dn_roku(rok0())
  })

  dis0 <- reactive({
    cat('distinct all ... \n')
    # if (input$fileIn!="customData") return(NULL)
    lis <- se.nor.rok.dis()
    if (is.null(geneIn0())|is.null(se())|is.null(lis)) return(NULL)
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
    if (is.null(geneIn0())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(dis0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct])
    up.dn <- up_dn(sam.all=lvl, df.all=dis0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='log2FC', fdr.na='FDR'); cat('Done! \n'); up.dn
  })
  
  output$upset1 <- renderPlot({
    cat('Upset up ... \n')
    if (is.null(geneIn0())|is.null(input$ssg.tar)) return(NULL)
    meth.sel <- c(!is.null(edg()), !is.null(lim()), !is.null(dsq()), !is.null(rok()), !is.null(dis()))
    validate(need(try(sum(meth.sel)>1), 'At least 2 methods should be selected!'))
    lis <- list(edgeR=edg(), limma=lim(), DESeq2=dsq(), ROKU=rok(), distinct=dis())[meth.sel]
    deg.lis <- deg_lis(lis, sam=input$ssg.tar, 'up')
    validate(need(try(length(deg.lis)>1), 'Over-expressed genes detected in less than two methods, so upset plot is unavailable for the selected target!'))
    # if (length(deg.lis)<2) return()
    up <- upset(fromList(deg.lis), order.by="degree", nintersects=40, point.size=3, line.size=1, mb.ratio=c(0.6, 0.4), text.scale=1.5); cat('Done! \n'); up
  })
  output$upset2 <- renderPlot({
    cat('Upset down ... \n')
    if (is.null(geneIn0())|is.null(input$ssg.tar)) return(NULL)
    meth.sel <- c(!is.null(edg()), !is.null(lim()), !is.null(dsq()), !is.null(rok()), !is.null(dis()))
    validate(need(try(sum(meth.sel)>1), 'At least 2 methods should be selected!'))
    lis <- list(edgeR=edg(), limma=lim(), DESeq2=dsq(), ROKU=rok(), distinct=dis())[meth.sel]
    deg.lis <- deg_lis(lis, sam=input$ssg.tar, 'down')
    validate(need(try(length(deg.lis)>1), 'Under-expressed genes detected in less than two methods, so upset plot is unavailable for the selected target!'))
    up <- upset(fromList(deg.lis), order.by="degree", nintersects=40, point.size=3, line.size=1, mb.ratio=c(0.6, 0.4), text.scale=1.5); cat('Done! \n'); up
  })

  output$ovl1 <- renderPlot({
    cat('Overlap up DEG ... \n')
    if (is.null(geneIn0())|is.null(input$ssg.tar)) return(NULL)
    meth.sel <- c(!is.null(edg()), !is.null(lim()), !is.null(dsq()), !is.null(rok()), !is.null(dis()))
    validate(need(try(sum(meth.sel)>1), 'At least 2 methods should be selected!'))
    lis <- list(edgeR=edg(), limma=lim(), DESeq2=dsq(), ROKU=rok(), distinct=dis())[meth.sel]
    deg.lis <- deg_lis(lis, sam=input$ssg.tar, 'up')
    validate(need(try(length(deg.lis)>1), 'Under-expressed genes detected in less than two methods, so upset plot is unavailable for the selected target!')) 
    names(deg.lis) <- sub('\\.up$|\\.down$', '', names(deg.lis))
    mat <- vapply(names(deg.lis), function(x) vapply(names(deg.lis), function(y) length(intersect(deg.lis[[x]], deg.lis[[y]])), numeric(1)), numeric(length(deg.lis)))
    mel <- reshape2::melt(mat)
    g <- ggplot(data=mel, aes(x=Var1, y=Var2, fill=value))+geom_tile(colour="white")+scale_fill_gradient(low="lightcyan3", high="darkorange")+theme_minimal()+theme(axis.text=element_text(angle=45, vjust=1, size=10, hjust=1))+coord_fixed()+geom_text(aes(Var2, Var1, label=value), color="black", size=4)+theme(axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.ticks=element_blank()); cat('Done! \n'); g
  })

  output$ovl2 <- renderPlot({
    cat('Overlap down DEG ... \n')
    if (is.null(geneIn0())|is.null(input$ssg.tar)) return(NULL)
    meth.sel <- c(!is.null(edg()), !is.null(lim()), !is.null(dsq()), !is.null(rok()), !is.null(dis()))
    validate(need(try(sum(meth.sel)>1), 'At least 2 methods should be selected!'))
    lis <- list(edgeR=edg(), limma=lim(), DESeq2=dsq(), ROKU=rok(), distinct=dis())[meth.sel]
    deg.lis <- deg_lis(lis, sam=input$ssg.tar, 'down')
    validate(need(try(length(deg.lis)>1), 'Under-expressed genes detected in less than two methods, so upset plot is unavailable for the selected target!'))
    names(deg.lis) <- sub('\\.up$|\\.down$', '', names(deg.lis)) 
    mat <- vapply(names(deg.lis), function(x) vapply(names(deg.lis), function(y) length(intersect(deg.lis[[x]], deg.lis[[y]])), numeric(1)), numeric(length(deg.lis)))
    mel <- reshape2::melt(mat)
    g <- ggplot(data=mel, aes(x=Var1, y=Var2, fill=value))+geom_tile(colour="white")+scale_fill_gradient(low="lightcyan3", high="darkorange")+theme_minimal()+theme(axis.text=element_text(angle=45, vjust=1, size=10, hjust=1))+coord_fixed()+geom_text(aes(Var2, Var1, label=value), color="black", size=4)+theme(axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.ticks=element_blank()); cat('Done! \n'); g
  })
venn_inter <- function(lis.all, suffix) {
  library(gplots)
  gen.all <- unique(unlist(lis.all))
  zero <- rep(0, length(gen.all))
  df.all <- as.data.frame(matrix(rep(0, length(gen.all)*length(lis.all)), ncol=length(lis.all)))
  cna <- names(lis.all)
  suf <- unique(gsub('(.*)\\.(.*)', '\\2', cna))
  meth <- unique(gsub('(.*)\\.(.*)', '\\1', cna))
  names(lis.all) <- colnames(df.all) <- meth
  rownames(df.all) <- gen.all
  # Retrieve all overlaps.
  lis <- venn(lis.all, show.plot=FALSE) 
  inter <- attributes(lis)$intersections
  # Translate all overlaps into a data frame.
  for (i in seq_along(inter)) {
    lis0 <- inter[[i]]; na0 <- strsplit(names(inter[i]), ':')[[1]]
    w1 <- which(gen.all %in% lis0); w2 <- which(meth %in% na0)  
    df.all[w1, w2] <- 1
  }

  df.all <- cbind(type=suf, total=rowSums(df.all), df.all)
  df.all <- df.all[order(df.all$total, decreasing=TRUE), ]
  return(df.all)
}


  dt.deg <- reactive({
    cat('DEG data frame ... \n')
    if (is.null(geneIn0())|is.null(geneIn1())|is.null(input$ssg.tar)) return(NULL)
    meth.sel <- c(!is.null(edg()), !is.null(lim()), !is.null(dsq()), !is.null(rok()), !is.null(dis()))
    validate(need(try(sum(meth.sel)>1), 'At least 2 methods should be selected!'))
    cat('Preparing DEG table ... \n')
    lis <- list(edgeR=edg(), limma=lim(), DESeq2=dsq(), ROKU=rok(), distinct=dis())[meth.sel]
    deg.lis1 <- deg_lis(lis, sam=input$ssg.tar, 'up')
    deg.lis2 <- deg_lis(lis, sam=input$ssg.tar, 'down')
    df.deg <- rbind(venn_inter(deg.lis1), venn_inter(deg.lis2))
    df.met <- geneIn1()[["df.met"]][, , drop=FALSE]
    na.deg <- rownames(df.deg); df.aggr.tran <- round(geneIn1()$df.aggr.tran, 2)
    if (nrow(df.met) > 0) df.deg <- cbind.data.frame(df.met[na.deg, , drop = FALSE], df.deg, stringsAsFactors=FALSE)
    df.deg <- cbind.data.frame(df.deg, df.aggr.tran[na.deg, , drop = FALSE], stringsAsFactors=FALSE)
    cat('Done! \n'); df.deg 
  })
    output$dld.ssg.tab <- downloadHandler(
      filename=function(){ "tissue_specific_genes.txt" },  content=function(file=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/tissue_specific_genes.txt')){
      write.table(dt.deg(), file, sep='\t', col.names=TRUE, row.names=TRUE) }
    )
    
  output$dt.deg <- renderDataTable({
      cat('DEG summary table ... \n')
      if (is.null(dt.deg())) return()
      d.tab <- datatable(rbind(dt.deg()), selection=list(mode="multiple", target="row", selected=NULL), filter="top", extensions='Scroller', plugins = "ellipsis",
      options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE, searchHighlight = TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching = TRUE, columnDefs = list(list(targets = c(1), render = JS("$.fn.dataTable.render.ellipsis(5, false)")))), 
      class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer'); cat('Done! \n')
      d.tab
  }) 

  # gID is a reactiveValue. It is shared between modules, so no need to pass the updated gID to shm_server.
  observeEvent(list(input$dt.deg_rows_selected), {
    cat('Selected genes in DEG summary table ... \n')
    if (is.null(input$dt.deg_rows_selected)) return()
    r.na <- rownames(dt.deg()); gID$geneSel <- r.na[input$dt.deg_rows_selected]
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

})}

deg.mod.lis <- deg_server('deg', ipt, cfg, dat.mod.lis, shm.mod.lis)
    #id.r <- rownames(se)[1]

    #plot_gen(se=se.nor, id=id.r)
  onStop(function() { ggly_rm(); vdo_rm(); cat("Session stopped! \n") })

})



