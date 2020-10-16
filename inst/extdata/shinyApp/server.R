
options(stringsAsFactors=FALSE) 

# source('~/tissue_specific_gene/function/fun.R')

# Right before submit the package the following functions will be deleted, and they will be imported as above. They are listed here now for the convenience of functionality development.

# Import internal functions.
sort_gen_con <- get('sort_gen_con', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

read_hdf5 <- get('read_hdf5', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

matrix_hm <- get('matrix_hm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Function to extract nearest genes.
sub_na <- get('sub_na', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

adj_mod <- get('adj_mod', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

filter_data <- get('filter_data', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

col_bar <- get('col_bar', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

lay_shm <- get('lay_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

nod_lin <- get('nod_lin', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Break combined path to a group (g=TRUE) or siblings (g=FALSE).
path_br <- get('path_br', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# The outline or tissue nodes are checked for combines paths. If combined paths are detected, those outside a group are broken to a group while those inside a group are broken as siblings.  
path_br_all <- get('path_br_all', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# 'a' nodes are not removed.
svg_attr <- get('svg_attr', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

svg_df <- get('svg_df', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

grob_list <- get('grob_list', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Separate SHMs of grobs and ggplot. Different SHMs of same 'gene_condition' are indexed with suffixed of '_1', '_2', ...
grob_gg <- get('grob_gg', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Subset data matrix by correlation or distance measure.
submatrix <- get('submatrix', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Adjust legend key size and rows in ggplot.
gg_lgd <- get('gg_lgd', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Add value keys SHMs.
gg_2lgd <- get('gg_2lgd', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Prepare interactive SHMs in html.
html_ly <- get('html_ly', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Make videos.
video <- get('video', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Import input matrix.
fread.df <- function(input, isRowGene=TRUE, header=TRUE, sep='auto', fill=TRUE, rep.aggr='mean', check.names=FALSE) {
  
  if (!is(input, 'data.frame') & !is(input, 'matrix')) { 
    
    df0 <- fread(input=input, header=header, sep=sep, fill=fill, check.names=check.names)
    cna <- make.names(colnames(df0))
    if (cna[1]=='V1') cna <- cna[-1] else cna <- cna[-ncol(df0)] 
    df1 <- as.data.frame(df0); rownames(df1) <- df1[, 1]
    # Subsetting identical column names in a matrix will not trigger appending numbers.
    df1 <- as.matrix(df1[, -1]); colnames(df1) <- cna
    if(isRowGene==FALSE) df1 <- t(df1)
    cna <- colnames(df1); rna <- rownames(df1) 
  
  } else { df1 <- input; rna <- rownames(df1); cna <- colnames(df1) }
  
  # Isolate data and annotation.
  na <- vapply(seq_len(ncol(df1)), function(i) { tryCatch({ as.numeric(df1[, i]) }, warning=function(w) { return(rep(NA, nrow(df1))) }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(df1)) )
  na <- as.data.frame(na); rownames(na) <- rna
  idx <- colSums(apply(na, 2, is.na))!=0
  gene2 <- na[!idx]; colnames(gene2) <- cna <- cna[!idx]
  gene3 <- as.data.frame(df1)[idx]
  form <- grepl("__", cna); if (sum(form)==0) { colnames(gene2) <- paste0(cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE
  if (ncol(gene3)>0) { colnames(gene3)[1] <- 'ann'; gene3 <- gene3[1] }
  if(sum(is.na(as.numeric(as.matrix(gene2))))>=1) return('Make sure all values in data matrix are numeric.')
  
  gen.rep <- gene2; rna <- rownames(gen.rep); gen.rep <-apply(gen.rep, 2, as.numeric); rownames(gen.rep) <- rna
  # Aggregate replicates.
  if (any(duplicated(cna)) & !is.null(rep.aggr)) {

    # To keep colnames, "X" should be a character, not a factor.
    if (rep.aggr=='mean') gene2 <- sapply(X=unique(cna), function(x) rowMeans(gene2[, cna==x, drop=FALSE]))
    if (rep.aggr=='median') {

      gene2 <- sapply(X=unique(cna), function(x) Biobase::rowMedians(gene2[, cna==x, drop=FALSE]))
      rownames(gene2) <- rna

    }

  }; gene2 <-apply(gene2, 2, as.numeric); rownames(gene2) <- rna 
  return(list(gene2=as.data.frame(gene2), gene3=as.data.frame(gene3), gen.rep=as.data.frame(gen.rep), con.na=con.na))

}

# Separate colour ingredients.
col_sep <- function(color) {

  color <- gsub(' |\\.|-|;|,|/', '_', color)
  color <- strsplit(color, '_')[[1]]
  color <- color[color!='']; return(color)

}

# Check/process "sample__condition" in the se extracted from tar.
se_from_db <- function(se) {
  dat <- assay(se); cold <- colData(se)
  form <- grepl("__", colnames(dat))
  if (sum(form)==0) {
    if (all(c('sample', 'condition') %in% colnames(cold))) { 
      cna <- colnames(dat) <- paste0(cold$sample, '__', cold$condition)
      if (any(duplicated(cna))) return('Duplicated "sample__condition" replicates are detected in the selected dataset!')
    } else if ('sample' %in% colnames(cold)) {
      if (any(duplicated(cold$sample))) return('The "sample" should not be duplicated in the absence of "condition"!')
      colnames(dat) <- cold$sample
    }
  }; return(dat)
}

# Extract target svgs in tar into tmp folder, and return the paths. 
extr_svg <- function(file, name) {

  dir <- paste0(tempdir(check=TRUE), '/svg_shm')
  if (!dir.exists(dir)) dir.create(dir)
  sys <- system(paste0('tar -xf ', file, ' -C ', dir, ' ', name))
  if (sys==0) return(paste0(dir, '/', name)) else return()

}

# Extract svg path/na from uploaded or internal tar files.
svg_pa_na <- function(svg.path, pa.svg.upl) {
  svg.na <- NULL; for (i in seq_along(svg.path)) {
    # Extract svg names. 
    str <- strsplit(svg.path[[i]], '/')[[1]]
    na0 <- str[length(str)]
    if (!grepl('\\.svg$', na0)) return('No aSVG file is detected! Solution: 1) select another aSVG and rematch it to data; 2) add an aSVG file for the selected data in the backend aSVG tar file or uploaded aSVG tar file.')
    svg.na <- c(svg.na, na0)
    # Complete uploaded svg paths.
    if (!grepl('example/', svg.path[[i]])) {
      # The data/svg precedence: uploaded tar > internal tar > default examples. The duplicated data/svgs are removed according to this precedence when processing data/svg upstream.
      pa0 <- NULL; if (!is.null(pa.svg.upl)) pa0 <- extr_svg(file=pa.svg.upl, name=na0)
      if (is.null(pa.svg.upl)|is.null(pa0)) {
        tar.all <- list.files('example', pattern='\\.tar$', full.names=TRUE)
        tar.svg <- tar.all[!grepl('data_shm.tar$', tar.all)][1]
        pa0 <- extr_svg(file=tar.svg, name=na0)
      }; if (is.null(pa0)) return(paste0("This aSVG file is not detected: ", na0, "!")) else svg.path[[i]] <- pa0
    }
  }; return(list(svg.path=svg.path, svg.na=svg.na))
}

# Check suffixes if multiple svgs.
svg_suffix <- function(svg.path, svg.na) {
  if (length(svg.na)>1) {
    if (!all(grepl('_shm\\d+\\.svg$', svg.na, perl=TRUE))) return("Suffixes of aSVGs should be indexed as '_shm1.svg', '_shm2.svg', '_shm3.svg', ...")
    ord <- order(gsub('.*_(shm.*)$', '\\1', svg.na))
    svg.path <- svg.path[ord]; svg.na <- svg.na[ord]  
  }; return(list(svg.path=svg.path, svg.na=svg.na))
}

# enableWGCNAThreads()
shinyServer(function(input, output, session) {

  cfg <- reactiveValues(lis.dat=NULL, lis.dld=NULL, lis.par=NULL, na.def=NULL, dat.def=NULL, svg.def=NULL, pa.upl=NULL, pa.dat.upl=NULL, pa.svg.upl=NULL, na.cus=NULL)

  observe({

    withProgress(message="Loading dependencies: ", value=0, {
    incProgress(0.3, detail="in progress...")
    library(SummarizedExperiment); library(shiny); library(shinydashboard); library(grImport); library(rsvg); library(ggplot2); library(DT) 
    incProgress(0.6, detail="in progress...")
    library(gridExtra); library(ggdendro); library(WGCNA); library(grid); library(xml2); library(plotly); library(data.table); library(genefilter); library(flashClust); library(visNetwork); 
    incProgress(0.9, detail="in progress...")
    library(reshape2); library(igraph); library(animation); library(av); library(shinyWidgets); library(yaml); library(HDF5Array); library(sortable)
  })

    lis.cfg <- yaml.load_file('config/config.yaml')
    lis.dat <- lis.cfg[grepl('^dataset\\d+', names(lis.cfg))]
    lis.dld <- lis.cfg[grepl('download_single|download_multiple', names(lis.cfg))]
    if (is.null(input$config)) lis.par <- lis.cfg[!grepl('^dataset\\d+|download_single|download_multiple', names(lis.cfg))] else lis.par <- yaml.load_file(input$config$datapath[1])
    upl.size <- toupper(lis.par$max.upload.size)
    num <- as.numeric(gsub('(\\d+)(G|M)', '\\1', upl.size))
    if (grepl('\\d+G$', upl.size)) max.size <- num*1024^3
    if (grepl('\\d+M$', upl.size)) max.size <- num*1024^2 
    options(shiny.maxRequestSize=max.size) 
    if (!any(lis.par$hide.legend %in% c('Yes', 'No'))) lis.par$hide.legend <- ifelse(lis.par$hide.legend==TRUE, 'Yes', 'No')
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
    if (!is.null(df.tar)) validate(need(try(sum(tar.num)==2), 'Two tar files of data and aSVGs respectively are expected!'))
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
    na.cus <- c('customData', 'customComputedData')
    dat.def <- c(dat.upl, dat.ipt[na.def]); svg.def <- c(svg.upl, svg.ipt[na.def])
    # If data/svg are duplicated between the server and upload, the data/svg on server side is removed.
    dat.def <- dat.def[unique(names(dat.def))]; svg.def <- svg.def[unique(names(svg.def))]
    cfg$lis.dat <- lis.dat; cfg$lis.dld <- lis.dld; cfg$lis.par <- lis.par; cfg$na.def <- names(dat.def); cfg$svg.def <- svg.def; cfg$dat.def <- dat.def; cfg$na.cus <- na.cus
    output$spatialHeatmap <- renderText({ lis.par$title['title', 'default'] })
    output$title.w <- renderText({ lis.par$title['width', 'default'] })
    dat.nas <- c('none', 'customData', 'customComputedData', names(dat.def))
    updateSelectInput(session, 'fileIn', 'Step 1: data sets', dat.nas, lis.par$default.dataset)
    updateRadioButtons(session, inputId='dimName', label='Step 4: is column or row gene?', choices=c("None", "Row", "Column"), selected=lis.par$col.row.gene, inline=TRUE)
    updateNumericInput(session, inputId="A", label="Value (A) to exceed:", value=as.numeric(lis.par$data.matrix['A', 'default']))
    updateNumericInput(session, inputId="P", label="Proportion (P) of samples with values >= A:", value=as.numeric(lis.par$data.matrix['P', 'default']))
    updateNumericInput(session, inputId="CV1", label="Min coefficient of variation (CV1):", value=as.numeric(lis.par$data.matrix['CV1', 'default']))
    updateNumericInput(session, inputId="CV2", label="Max coefficient of variation (CV2):", value=as.numeric(lis.par$data.matrix['CV2', 'default']))
    updateRadioButtons(session, inputId='log', label='Log/exp:', choices=c("No", "log2", "exp2"), selected=lis.par$log.exp, inline=TRUE)
    updateRadioButtons(session, inputId='scale', label='Scale by:', choices=c('No', 'Row', 'Column'), selected=lis.par$data.matrix.scale, inline=TRUE)
    updateRadioButtons(session, inputId='hide.lgd', label="Hide legend:", choices=c('Yes', 'No'), selected=lis.par$hide.legend, inline=TRUE)
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

    dld.exp <- reactiveValues(sgl=NULL, mul=NULL)
    dld.exp$sgl <- cfg$lis.dld$download_single
    dld.exp$mul <- cfg$lis.dld$download_multiple

    output$dld.cfg <- downloadHandler(
      filename=function(){ "config_par.yaml" },
      content=function(file=paste0(tempdir(), '/config_par.yaml')){ 
 
        lis.cfg <- yaml.load_file('config/config.yaml')
        lis.par <- lis.cfg[c("max.upload.size", "default.dataset", "col.row.gene", "separator", "hide.legend", "data.matrix", "shm.img", "shm.anm", "shm.video", "legend", "mhm", "network")]
        write_yaml(lis.par, file)

      }
    )
    output$dld.sgl <- downloadHandler(
      filename=function(){ "single_aSVG_data.zip" },
      content=function(file=paste0(tempdir(), '/single_aSVG_data.zip')){ zip(file, c(dld.exp$sgl$data, dld.exp$sgl$svg)) }
    )
    output$dld.mul <- downloadHandler(   
      filename=function(){ "multiple_aSVG_data.zip" },
      content=function(file=paste0(tempdir(), '/multiple_aSVG_data.zip')){ zip(file, c(dld.exp$mul$data, dld.exp$mul$svg)) }
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

  # Filter parameters.
  fil <- reactiveValues(P=0, A=0, CV1=-Inf, CV2=Inf)

  observe({

    input$fileIn; input$geneInpath
    updateRadioButtons(session, inputId="dimName", label="Step 4: is column or row gene?", 
    inline=TRUE, choices=c("None", "Row", "Column"), selected="None")
    updateRadioButtons(session, inputId='log', label='Log/exp transform:', choices=c("No", "log2", "exp2"), selected=cfg$lis.par$data.matrix['log.exp', 'default'], inline=TRUE)
    updateRadioButtons(session, 'scale', label='Scale by:', choices=c('No', 'Row', 'Column'), selected=cfg$lis.par$data.matrix['scale', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId='cs.v', label='Color scale based on:', choices=c("Selected rows", "All rows"), selected=cfg$lis.par$shm.img['color.scale', 'default'], inline=TRUE)
    updateNumericInput(session, inputId="height", label="Overall height:", value=as.numeric(cfg$lis.par$shm.img['height', 'default']), min=0.1, max=Inf, step=NA)
    updateNumericInput(session, inputId="width", label="Overall width:", value=as.numeric(cfg$lis.par$shm.img['width', 'default']), min=0.1, max=Inf, step=NA)
    updateNumericInput(session, inputId="col.n", label="Columns:", value=as.numeric(cfg$lis.par$shm.img['columns', 'default']), min=1, max=Inf, step=1)

  })

  observe({

    input$fileIn; input$geneInpath; input$log
    updateNumericInput(session, inputId="A", label="Value (A) to exceed:", value=as.numeric(cfg$lis.par$data.matrix['A', 'default'])) 
    updateNumericInput(session, inputId="P", label="Proportion (P) of samples with values >= A:", value=as.numeric(cfg$lis.par$data.matrix['P', 'default']), min=0, max=1)
    updateNumericInput(session, inputId="CV1", label="Min coefficient of variation (CV1):", value=as.numeric(cfg$lis.par$data.matrix['CV1', 'default']))
    updateNumericInput(session, inputId="CV2", label="Max coefficient of variation (CV2):", value=as.numeric(cfg$lis.par$data.matrix['CV2', 'default'])) 
    fil$P <- 0; fil$A <- 0; fil$CV1 <- -Inf; fil$CV2 <- Inf

  })

  observeEvent(input$fil.but, {

    if (input$fileIn=="none") return(NULL)  
    fil$P <- input$P; fil$A <- input$A; fil$CV1 <- input$CV1; fil$CV2 <- input$CV2
  
  })

  output$fil.par <- renderText({
    
    if (input$fileIn=="none") return(NULL)  
    P <- input$P
    validate(need(try(P<=1 & P>=0), 'P should be between 0 to 1 !'))

  })

  geneIn0 <- reactive({

    if (input$fileIn=="none") return(NULL)
    withProgress(message="Loading data: ", value = 0, {
    if (any(input$fileIn %in% cfg$na.def)) {

      incProgress(0.5, detail="Loading matrix. Please wait.")
      dat.na <- cfg$dat.def[input$fileIn]
      if ('example' %in% strsplit(dat.na, '/')[[1]]) df.te <- fread.df(input=dat.na, isRowGene=TRUE) else { 
        # Exrtact data from uploaded tar.
        dat <- NULL; if (!is.null(cfg$pa.dat.upl)) if (file.exists(cfg$pa.dat.upl)) {
          cat('Extracting uploaded data... \n')
          # The prefix is not input$fileIn. The returned value of read_hdf5 is either data.frame or SE.
          dat <- read_hdf5(cfg$pa.dat.upl, prefix=dat.na)[[1]]
        }
        if (is.null(dat)|is.character(dat)|is.null(input$tar)) {
          cat('Extracting internal data... \n')
          dat <- read_hdf5('example/data_shm.tar', dat.na)[[1]]
        }
        validate(need(try(is(dat, 'data.frame')|is(dat, 'SummarizedExperiment')), 'The selected data is empty! Solution: 1) select another data, then rematch its samples to the selected aSVG; 2) include a data for the selected aSVG file in the backend database or the uploaded database.'))
        if (is(dat, 'SummarizedExperiment')) {
          dat <- se_from_db(dat)
          validate(need(try(!is.character(dat)), dat))
        }; df.te <- fread.df(input=dat)
      }; return(df.te)
    }
    if (any(input$fileIn %in% cfg$na.cus) & 
    !is.null(input$geneInpath) & input$dimName!="None") {

      incProgress(0.25, detail="Importing matrix. Please wait.")
      geneInpath <- input$geneInpath
      df.upl <- fread.df(input=geneInpath$datapath, isRowGene=(input$dimName=='Row')); return(df.upl)
   
    }

    })

  })

  # Transform data.
  geneIn1 <- reactive({

    if (is.null(geneIn0())|is.null(input$scale)) return(NULL)
    gene1 <- gene2 <- geneIn0()[['gene2']] 
    if (input$log=='log2') {
   
      g.min <- min(gene2)
      if (g.min<0) gene2 <- gene2-g.min+1; if (g.min==0) gene2 <- gene2+1; gene2 <- log2(gene2)

    }; if (input$log=='exp2') gene2 <- 2^gene2
    # Scale by row/column
    if (input$scale=='Row') { gene2 <- t(scale(t(gene2))) } else if (input$scale=='Column') { gene2 <- scale(gene2) }

    gene3 <- geneIn0()[['gene3']]; gen.rep <- geneIn0()[['gen.rep']]
    return(list(gene1=gene1, gene2=gene2, gene3=gene3, gen.rep=gen.rep))

  })

  output$col.order <- renderUI({

    if (is.null(geneIn1())) return()
    col.nas <- colnames(geneIn1()[['gene2']])
    dropdownButton(inputId='dropdown', label='Re-order columns', circle=FALSE, icon=NULL, status='primary',
    actionButton("col.cfm", "Confirm", icon=icon("refresh")), 
    selectizeInput(inputId="col.na", label='', choices=col.nas, selected=col.nas, multiple=TRUE, options= list(plugins=list('remove_button', 'drag_drop')))
    )
  
  })
  sear <- reactiveValues(id=NULL)
  observeEvent(input$fileIn, { sear$id <- NULL })
  observeEvent(input$search.but, {
    
    if (is.null(geneIn1())) return()
    if (input$search=='') sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) else {

      gens <- strsplit(gsub(' |,', '_', input$search), '_')[[1]]
      pat <- paste0('^', gens, '$', collapse='|')
      sel <- which(grepl(pat, x=rownames(geneIn1()[['gene2']]), ignore.case=TRUE, perl=TRUE))
      if (length(sel)==0) sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default'])

     }; sear$id <- sel

  })
  geneIn <- reactive({
    
    if (is.null(geneIn1())) return(NULL)    
    gene1 <- geneIn1()[['gene1']]; gene2 <- geneIn1()[['gene2']]; gene3 <- geneIn1()[['gene3']]; input$fil.but
    if (!identical(sort(input$col.na), sort(colnames(gene2)))) return()
    # Input variables in "isolate" will not triger re-excution, but if the whole reactive object is trigered by "input$fil.but" then code inside "isolate" will re-excute.
    isolate({
  
      se <- SummarizedExperiment(assays=list(expr=as.matrix(gene2)), rowData=gene3)
      if (ncol(gene3)>0) ann.col <- colnames(gene3)[1] else ann.col <- NULL
      # If scaled by row, sd is 1, mean is 0, cv is Inf.
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=c(ifelse(input$scale=='Row', -Inf, fil$CV1), ifelse(input$scale=='Row', Inf, fil$CV2)), dir=NULL)
      if (nrow(se)==0) { validate(need(try(nrow(se)>0), 'All rows are filtered out!')); return() }
      # In case of all rows are filtered, the app continues to work without refreshing after the filter parameters are reduced.
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=c(ifelse(input$scale=='Row', -Inf, fil$CV1), ifelse(input$scale=='Row', Inf, fil$CV2)), dir=NULL)

      gene2 <- as.data.frame(assay(se), stringsAsfactors=FALSE); colnames(gene2) <- make.names(colnames(gene2))
      gene1 <- gene1[rownames(gene2), ]
      gene3 <- as.data.frame(rowData(se))[, , drop=FALSE]
      
    })
    cat('Preparing data matrix... \n')
    if (is.null(sear$id)) rows <- seq_len(nrow(gene2)) else rows <- sear$id
    if (length(rows)==1 & rows[1]==as.numeric(cfg$lis.par$data.matrix['row.selected', 'default'])) rows <- seq_len(nrow(gene2))
    return(list(gene1=gene1[rows, input$col.na], gene2=gene2[rows, input$col.na], gene3=gene3[rows, , drop=FALSE]))

  })
  output$dt <- renderDataTable({

    if (is.null(geneIn())) return()
    if ((any(input$fileIn %in% cfg$na.cus) & is.null(geneIn()))|input$fileIn=="none") return(NULL)

    withProgress(message="Data table: ", value = 0, {

      incProgress(0.5, detail="Displaying. Please wait.")
      if (input$fileIn!="none") {

      gene <- geneIn(); gene.dt <- cbind.data.frame(gene[["gene2"]][, , drop=FALSE], gene[["gene3"]][, , drop=FALSE], stringsAsFactors=FALSE) 

   }; cat('Presenting data matrix... \n')
   if (is.null(sear$id)) sel <- as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) else sel <- sear$id
   if (length(sel)==1 & sel[1]==as.numeric(cfg$lis.par$data.matrix['row.selected', 'default']) & nrow(gene.dt)>1) sel <- sel else if (nrow(gene.dt)==1)  sel <- 1 else if (length(sel)>1) sel <- seq_along(sel)

   datatable(gene.dt, selection=list(mode="multiple", target="row", selected=sel),
   filter="top", extensions=c('Scroller'), options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE, searchHighlight=FALSE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
    formatRound(colnames(geneIn()[["gene2"]]), 2)

    })

  })

  gID <- reactiveValues(geneID="none", new=NULL, all=NULL)
  observe({ input$geneInpath; input$fileIn; gID$geneID <- "none" })
  observe({ if (is.null(geneIn())) gID$geneID <- "none" })
  # To make the "gID$new" and "gID$all" updated with the new "input$fileIn", since the selected row is fixed (3rd row), the "gID$new" is not updated when "input$fileIn" is changed, and the downstream is not updated either. The shoot/root examples use the same data matrix, so the "gID$all" is the same (pre-selected 3rd row) when change from the default "shoot" to others like "organ". As a result, the "gene$new" is null and downstream is not updated. Also the "gene$new" is the same when change from shoot to organ, and downstream is not updated, thus "gene$new" and "gene$all" are both set NULL above upon new "input$fileIn".  
  observeEvent(input$fileIn, {

    if (is.null(input$dt_rows_selected)) return()
    gID$all <- gID$new <- NULL
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    # Avoid multiple selected rows from last input$fileIn. Must be behind gID$geneID. 
    if (length(input$dt_rows_selected)>1) return()
    gID$new <- setdiff(gID$geneID, gID$all); gID$all <- c(gID$all, gID$new)
    if (is.null(r.na)) gID$geneID <- "none"
    })
  observeEvent(list(input$dt_rows_selected, geneIn()), {
 
    if (is.null(input$dt_rows_selected)) return()
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    if (any(is.na(gID$geneID))) gID$geneID <- "none"
    gID$new <- setdiff(gID$geneID, gID$all); gID$all <- c(gID$all, gID$new)
    
  })
  observeEvent(list(input$search.but), {
  
    if (is.null(input$search.but)|is.null(sear$id)|is.null(geneIn())) return()
    gID$geneID <- rownames(geneIn()[["gene2"]])[sear$id]
    if (any(is.na(gID$geneID))) gID$geneID <- "none"
 
  })

  geneV <- reactive({

    if (any(is.na(gID$geneID))) return()
    if (is.null(geneIn())|sum(gID$geneID[1]!='none')==0) return(NULL)
    if (input$cs.v=="Selected rows" & is.null(input$dt_rows_selected)) return(NULL)
    if (input$fileIn!="none") { if (input$cs.v=="Selected rows") gene <- geneIn()[["gene2"]][gID$geneID, ]
    if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]] }
    seq(min(gene), max(gene), len=1000) # len must be same with that from the function "spatial_hm()". Otherwise the mapping of a gene value to the colour bar is not accurate. 

  })

  col.sch <- reactive({ 

    if(input$color=="") return(NULL)
    col <- gsub(' |\\.|-|;|,|/', '_', input$color)
    col <- strsplit(col, '_')[[1]]
    col <- col[col!='']; col1 <- col[!col %in% colors()]
    if (length(col1>0)) validate(need(try(col1 %in% colors()), paste0('Colors not valid: ', col1, ' !'))); col

  })
  
  color <- reactiveValues(col="none")
  observe({
    
    col0 <- cfg$lis.par$shm.img['color', 'default']
    if (is.null(input$col.but)|is.null(col0)) return()
    if(input$col.but==0) color$col <- colorRampPalette(col_sep(col0))(length(geneV()))

  })
  # As long as a button is used, observeEvent should be used. All variables inside 'observeEvent' trigger code evaluation, not only 'eventExpr'.  
  observeEvent(input$col.but, {

    if (is.null(col.sch())) return (NULL)
    if (input$fileIn!="none") { color$col <- colorRampPalette(col.sch())(length(geneV())) }

  })

  shm.bar <- reactive({

    if (is.null(gID$all)) return(NULL)
    if ((any(input$fileIn %in% cfg$na.def) & !is.null(geneIn()))|(any(input$fileIn %in% cfg$na.cus) & (!is.null(input$svgInpath)|!is.null(input$svgInpath1)) & !is.null(geneIn()))) {

      if (length(color$col=="none")==0|input$color==""|is.null(geneV())) return(NULL)

      withProgress(message="Color scale: ", value = 0, {

        incProgress(0.75, detail="Plotting. Please wait.")
        cat('Colour key... \n')
        cs.g <- col_bar(geneV=geneV(), cols=color$col, width=1); return(cs.g)

      })

    }

  })
  # One output can only be used once in ui.R.
  output$bar1 <- renderPlot({ if (!is.null(shm.bar)) shm.bar() })
  output$bar2 <- renderPlot({ if (!is.null(shm.bar)) shm.bar() })

  observe({

    if (is.null(geneIn())) return(NULL)
    r.na <- rownames(geneIn()[["gene2"]]); gens.sel <- r.na[input$dt_rows_selected]
    if (length(gens.sel)==0) return()
    updateSelectInput(session, inputId="gen.sel", label="Select a target gene:", choices=c("None", gens.sel), selected=gens.sel[1])

  })

  svg.na.mat <- reactiveValues(svg.path=NULL, svg.na=NULL)
  svg.path <- reactive({
    if (input$fileIn=='none') return()
    if (any(input$fileIn %in% cfg$na.cus)) {
      if (is.null(input$svgInpath1)) svgIn.df <- input$svgInpath else svgIn.df <- input$svgInpath1
      svg.path <- svgIn.df$datapath; svg.na <- svgIn.df$name
    } else {
      # Extract svg path and name: single or multiple svg paths are treated same way.
      lis <- svg_pa_na(cfg$svg.def[[input$fileIn]], cfg$pa.svg.upl)
      validate(need(try(!is.character(lis)), lis))
      svg.path <- lis$svg.path; svg.na <- lis$svg.na
    }; cat('Access aSVG path... \n')
    # If multiple svgs, check suffixes.
    lis <- svg_suffix(svg.path, svg.na)
    validate(need(try(!is.character(lis)), lis)); return(lis)
  })

  sam <- reactive({ 
    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\1", c.na) else return(NULL) 
  })

  output$svg <- renderUI({
    nas <- names(cfg$svg.def)
    selectInput('svg', label='aSVGs to re-match:', choices=nas, selected=input$fileIn)
  })

  # The object scope in reactive expression is restricted to the reactive evironment, thus sf.sep and sam.sep are moved outside 'observe'.
  sf.sep <- 'only_above_features_are_re-matched'
  sam.sep <- 'only_above_samples_are_re-matched'
  observeEvent(input$svg, {
    if (input$fileIn=='none'|is.null(cfg$svg.def)|is.null(input$svg)) return()
    if (!(any(input$fileIn %in% cfg$na.def) & is.null(input$svgInpath))) return()
    # Single or multiple svg paths are treated same way.
    lis <- svg_pa_na(cfg$svg.def[[input$svg]], cfg$pa.svg.upl)
    validate(need(try(!is.character(lis)), lis))
    svg.path <- lis$svg.path; svg.na <- lis$svg.na
    cat('Access aSVG path for re-matching... \n')
    # If multiple svgs, check suffixes.
    lis <- svg_suffix(svg.path, svg.na)
    validate(need(try(!is.character(lis)), lis))
    svg.path <- lis$svg.path; svg.na <- lis$svg.na
 
  withProgress(message="Tissue heatmap: ", value=0, {  
    incProgress(0.5, detail="Extracting coordinates. Please wait.") 
    # Whether a single or multiple SVGs, all are returned in a list.
    sf.all <- NULL; for (i in seq_along(svg.na)) { 
      cat('Extract all spatial features for re-matching:', svg.na[i], '\n')
      df_tis <- svg_df(svg.path=svg.path[i], feature=sam())
      validate(need(!is.character(df_tis), paste0(svg.na[i], ': ', df_tis)))
      sf.all <- c(sf.all, df_tis$tis.path)
    }
  })
  # paths and gs are dropped to bottom. 
  sf.all <- unique(sf.all); pas.idx <- grepl('^path\\d+|^g\\d+', sf.all)
  sf.all <- c(sf.all[!pas.idx], sf.all[pas.idx])
  # Matching samples are raised to top.
  sam.all <- unique(sam()); inter <- intersect(sam.all, sf.all)
  cat('Adding separators to samples and features for re-matching... \n')
  if (length(inter)!=0) { 
    int.idx1 <- sam.all %in% inter; inter1 <- sam.all[int.idx1]
    sam.all <- c(inter1, sam.sep, sam.all[!int.idx1])
    int.idx2 <- sf.all %in% inter; sf.all <- c(inter1, sf.sep, sf.all[!int.idx2])
  } else { sf.all <- c(sf.sep, sf.all); sam.all <- c(sam.sep, sam.all) }
  output$sam <- renderUI({
    rank_list(text="Drag samples in any desired order to match right", labels=sam.all, input_id="sam")
  })

  output$sf <- renderUI({
    rank_list(text="Drag spatial features in any desired order to match left", labels=sf.all, input_id="sf"
    )
  }); svg.na.mat$svg.path <- svg.path; svg.na.mat$svg.na <- svg.na

})
    
  observeEvent(input$fileIn, {
    cna.match$cna <- NULL
    svg.na.mat$svg.path <- svg.na.mat$svg.na <- NULL
  })
  # Reactive object in "observeEvent" is not accessible outside observeEvent. The solution is eventReactive. 
  svg.path1 <- reactive({
    if (!is.null(svg.na.mat$svg.path) & !is.null(svg.na.mat$svg.na)) { svg.path <- svg.na.mat$svg.path; svg.na <- svg.na.mat$svg.na } else { svg.path <- svg.path()$svg.path; svg.na <- svg.path()$svg.na }
    return(list(svg.path=svg.path, svg.na=svg.na))
  })

  cna.match <- reactiveValues(cna=NULL)
  observeEvent(input$match, {
  
    sam <- input$sam; sf <- input$sf
    cna <- colnames(geneIn()[["gene2"]])
    w.sf <- which(sf %in% sf.sep); w.sam <- which(sam %in% sam.sep)
    # If no "renderText", the validate message is not assigned to any reactive values and thus not seen on the user interface.
    output$msg.match <- renderText({
    validate(need(try(w.sf==w.sam), 'Samples and features should be one-to-one re-matched!')); NULL })
    sf <- sf[-w.sf]; w.sf <- w.sf-1; sam <- sam[-w.sam]; w.sam <- w.sam-1
    output$msg.match <- renderText({
    validate(need(try(w.sf!=0|w.sam!=0), 'No samples or features are re-matched!')); NULL })
    if (w.sf>0) {
      cat('Re-match samples and features... \n')
      sf.mat <- sf[seq_len(w.sf)]
      cna <- colnames(geneIn()[["gene2"]]);
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

    if ((any(input$fileIn %in% cfg$na.cus) & 
    (!is.null(input$svgInpath)|!is.null(input$svgInpath1)))|(any(input$fileIn %in% cfg$na.def) & is.null(input$svgInpath))) {

      withProgress(message="Tissue heatmap: ", value=0, {
    
        incProgress(0.5, detail="Extracting coordinates. Please wait.")
          svg.path <- svg.path1()$svg.path; svg.na <- svg.path1()$svg.na
          # Whether a single or multiple SVGs, all are returned in a list.
         svg.df.lis <- NULL; for (i in seq_along(svg.na)) {
         
            cat('Coordinate:', svg.na[i], '\n')
            df_tis <- svg_df(svg.path=svg.path[i], feature=sam())
            validate(need(!is.character(df_tis), paste0(svg.na[i], ': ', df_tis)))
            svg.df.lis <- c(svg.df.lis, list(df_tis))
   
          }; names(svg.df.lis) <- svg.na; return(svg.df.lis)

      })

    }

  })


  observe({

    input$fileIn; geneIn(); input$adj.modInpath; svg.df(); input$hide.lgd 
    tis.tran <- NULL; for (i in seq_along(svg.df())) { tis.tran <- c(tis.tran, svg.df()[[i]][['tis.path']]) }
    updateCheckboxGroupInput(session, inputId="tis", label='Select tissues to be transparent:', choices=intersect(unique(sam()), unique(tis.tran)), selected='', inline=TRUE)

  })

  con <- reactive({ 

    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\3", c.na) else return(NULL) 

  })

  # General selected gene/condition pattern.
  pat.con <- reactive({ con.uni <- unique(con()); if (is.null(con.uni)) return(NULL); paste0(con.uni, collapse='|') })
  pat.gen <- reactive({ if (is.null(gID$geneID)) return(); if (gID$geneID[1]=='none') return(NULL);  paste0(gID$geneID, collapse='|') })
  pat.all <- reactive({ if (is.null(pat.con())|is.null(pat.gen())) return(NULL); paste0('(', pat.gen(), ')_(', pat.con(), ')') })

  grob <- reactiveValues(all=NULL, all1=NULL, gg.all=NULL, gg.all1=NULL, lgd.all=NULL)
  observeEvent(input$fileIn, { grob$all <- grob$gg.all1 <- grob$gg.all1 <- grob$gg.all <- grob$lgd.all <- NULL })
  # Avoid repetitive computation under input$cs.v=='w.mat'.
  gs.new <- reactive({ 
    if.con <- is.null(svg.df())|is.null(geneIn())|is.null(gID$new)|length(gID$new)==0|is.null(gID$all)|is.null(input$dt_rows_selected)|color$col[1]=='none'
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)

    if (input$cs.v=="Selected rows") ID <- gID$geneID
    if (input$cs.v=="w.mat") ID <- gID$new
    if (is.null(ID)|length(gID$new)>1|length(ID)>1|ID[1]=='none') return()
    # Avoid repetitive computation.  
    pat.new <- paste0('^', gID$new, '_(', pat.con(), ')_\\d+$')
    if (any(grepl(pat.new, names(grob$all)))) return()
    withProgress(message="Tissue heatmap: ", value=0, {
 
      incProgress(0.25, detail="preparing data.")
      gene <- geneIn()[["gene2"]]
      svg.df.lis <- svg.df(); grob.lis.all <- w.h.all <- NULL
      # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
      for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      svg.na <- names(svg.df.lis)
      for (i in seq_along(svg.df.lis)) {

        svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]; w.h <- svg.df[['w.h']]
        tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']]
        if (input$pre.scale=='Yes') mar <- (1-w.h/w.h.max*0.99)/2 else mar <- NULL
        cat('New grob/ggplot:', ID, ' \n')
        if (!is.null(cna.match$cna)) { 
		  if (ncol(gene)==length(cna.match$cna)) colnames(gene) <- cna.match$cna 
        }
        grob.lis <- grob_list(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=ID, legend.col=fil.cols, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=18, mar.lb=mar, legend.nrow=2, legend.key.size=0.04) # Only gID$new is used.
        msg <- paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')
        if (is.null(grob.lis)) cat(msg, '\n')
        output$msg.shm <- ({ validate(need(!is.null(grob.lis), msg)) })
        grob.lis.all <- c(grob.lis.all, list(grob.lis))

      }; names(grob.lis.all) <- svg.na; return(grob.lis.all)

    })

  })

  # Extension of 'observeEvent': any of 'input$log; input$tis; input$col.but; input$cs.v' causes evaluation of all code. 
  # input$tis as an argument in "grob_list" will not cause evaluation of all code, thus it is listed here.
  # Use "observeEvent" to replace "observe" and list events (input$log, input$tis, ...), since if the events are in "observe", every time a new gene is clicked, "input$dt_rows_selected" causes the evaluation of all code in "observe", and the evaluation is duplicated with "gs.new".
  col.reorder <- reactiveValues(col.re='Y')
  observeEvent(input$col.na, { if (input$col.cfm>0) col.reorder$col.re <- 'N' })
  observeEvent(list(input$log, input$tis, input$col.but, input$cs.v, input$pre.scale, input$col.cfm, input$scale, input$match), {
    
    grob$all <- grob$gg.all <- grob$lgd.all <- NULL; gs.all <- reactive({ 

      if.con <- is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none'|is.null(input$pre.scale)|gID$geneID[1]=='none'

      if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
      withProgress(message="Spatial heatmap: ", value=0, {
      incProgress(0.25, detail="preparing data.")
      #if (input$cs.v=="Selected rows") gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
      #if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]]
      gene <- geneIn()[["gene2"]][gID$geneID, ]
      svg.df.lis <- svg.df(); grob.lis.all <- w.h.all <- NULL
      # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
      for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      for (i in seq_along(svg.df.lis)) {
        
        svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]
        tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']]; w.h <- svg.df[['w.h']]
        if (input$pre.scale=='Yes') mar <- (1-w.h/w.h.max*0.99)/2 else mar <- NULL
        cat('All grob/ggplot:', gID$geneID, ' \n')
        svg.na <- names(svg.df.lis)
        incProgress(0.75, detail=paste0('preparing ', paste0(gID$geneID, collapse=';')))
        if (!is.null(cna.match$cna)) { 
		  if (ncol(gene)!=length(cna.match$cna)) return()
          colnames(gene) <- cna.match$cna 
        }
        grob.lis <- grob_list(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=gID$geneID, legend.col=fil.cols, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=18, mar.lb=mar, legend.nrow=2, legend.key.size=0.04) # All gene IDs are used.
        msg <- paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')
        if (is.null(grob.lis)) cat(msg, '\n')
        output$msg.shm <- ({ validate(need(!is.null(grob.lis), msg)) })
        grob.lis.all <- c(grob.lis.all, list(grob.lis))

      }; names(grob.lis.all) <- svg.na; return(grob.lis.all)

      })

    }); grob.gg.lis <- grob_gg(gs=gs.all())
    grob$all <- grob.gg.lis[['grob']]; grob$gg.all <- grob.gg.lis[['gg']]; grob$lgd.all <- grob.gg.lis[['lgd.all']]

  })
  # Avoid repetitive computation under input$cs.v=='gen.sel'.
  observeEvent(list(gID$geneID), {
    
    if.con <-  is.null(input$cs.v)|gID$geneID[1]=='none'|input$cs.v=='w.mat'

    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    ID <- gID$geneID
    grob$all <- grob$gg.all <- grob$lgd.all <- NULL; gs.all <- reactive({ 

      if.con <- is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none'|is.null(input$pre.scale)
      if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
      withProgress(message="Spatial heatmap: ", value=0, {
      incProgress(0.25, detail="preparing data.")
      gene <- geneIn()[["gene2"]][gID$geneID, ]
      svg.df.lis <- svg.df(); grob.lis.all <- w.h.all <- NULL
      # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
      for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      for (i in seq_along(svg.df.lis)) {

        svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]
        tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']]; w.h <- svg.df[['w.h']]
        if (input$pre.scale=='Yes') mar <- (1-w.h/w.h.max*0.99)/2 else mar <- NULL
        cat('All grob/ggplot of row selection:', ID, ' \n')
        svg.na <- names(svg.df.lis)
        incProgress(0.75, detail=paste0('preparing ', paste0(ID, collapse=';')))
        if (!is.null(cna.match$cna)) { 
		  if (ncol(gene)!=length(cna.match$cna)) return()
          colnames(gene) <- cna.match$cna 
        }
        grob.lis <- grob_list(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=ID, legend.col=fil.cols, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=18, mar.lb=mar, legend.nrow=2, legend.key.size=0.04) # All gene IDs are used.
        msg <- paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')
        if (is.null(grob.lis)) cat(msg, '\n')
        output$msg.shm <- ({ validate(need(!is.null(grob.lis), msg)) })
        grob.lis.all <- c(grob.lis.all, list(grob.lis))

      }; names(grob.lis.all) <- svg.na; return(grob.lis.all)

      })

    }); grob.gg.lis <- grob_gg(gs=gs.all())
    grob$all <- grob.gg.lis[['grob']]; grob$gg.all <- grob.gg.lis[['gg']]; grob$lgd.all <- grob.gg.lis[['lgd.all']]

  })
  # when 'color <- reactiveValues(col="none")', upon the app is launched, 'gs.new' is evaluated for 3 time. In the 1st time, 'gID$new'/'gID$all' are NULL, so 'gs.new' is NULL. In the 2nd time, 'color$col[1]=='none'' is TRUE, so NULL is returned to 'gs.new', but 'gID$new'/'gID$all' are 'HRE2'. In the third time, 'color$col[1]=='none'' is FALSE, so 'gs.new' is not NULL, but 'gID$new' is still 'HRE2', so it does not triger evaluation of 'observeEvent' and hence SHMs and legend plot are not returned upon being launched. The solution is to assign colors to 'color$col' in 'observe' upon being launched so that in the 2nd time 'gs.new' is not NULL, and no 3rd time.
  observeEvent(gs.new(), { 

    if (is.null(svg.df())|is.null(gID$new)|length(gID$new)==0|is.null(gID$all)|is.null(gs.new())) return(NULL)
    cat('Updating grobs/ggplots/lgds... \n')
    grob.gg.lis <- grob_gg(gs=gs.new())
    grobs <- grob.gg.lis[['grob']]
    grob.rm <- !names(grob$all) %in% names(grobs)
    grob$all <- c(grob$all[grob.rm], grobs)
    ggs <- grob.gg.lis[['gg']]
    gg.rm <- !names(grob$gg.all) %in% names(ggs)
    grob$gg.all <- c(grob$gg.all[gg.rm], ggs) 
    lgd0 <- grob.gg.lis[['lgd.all']]
    grob$lgd.all <- c(grob$lgd.all, lgd0[!names(lgd0) %in% names(grob$lgd.all)])
  
  })
 
  output$h.w.c <- renderText({
    
    if (is.null(geneIn())|is.null(input$dt_rows_selected)|is.null(svg.df())|is.null(grob$all)) return(NULL)

    height <- input$height; width <- input$width
    col.n <- input$col.n;
    validate(need(height>=0.1 & !is.na(height), 'Height should be a positive numeric !'))
    validate(need(width>=0.1 & !is.na(width), 'Width should be a positive numeric !'))
    validate(need(col.n>=1 & as.integer(col.n)==col.n & !is.na(col.n), 'No. of columns should be a positive integer !'))

  })
  observeEvent(list(input$lgd.key.size, input$lgd.row, input$tis, input$lgd.label, input$lgd.lab.size), {

    lgd.key.size <- input$lgd.key.size; lgd.row <- input$lgd.row
    lgd.label <- input$lgd.label; label.size <- input$lgd.lab.size
    if (is.null(grob$lgd.all)|is.null(lgd.key.size)|is.null(lgd.row)|is.null(lgd.label)) return()
    cat('Adjust legend size/rows... \n')
    grob$lgd.all <- gg_lgd(gg.all=grob$lgd.all, size.key=lgd.key.size, size.text.key=NULL, row=lgd.row, sam.dat=sam(), tis.trans=input$tis, position.text.key='right', label=(lgd.label=='Yes'), label.size=label.size)
 
  })
  observeEvent(list(grob.all=grob$all, gen.con=input$gen.con), {
  
    if (is.null(gID$all)|is.null(grob$all)|is.null(grob$gg.all)) return()
    cat('Reordering grobs/ggplots... \n')
    na.all <- names(grob$all); pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', na.all)); if (length(con)==0) return()
    na.all <- sort_gen_con(ID.sel=gID$all, na.all=na.all, con.all=con, by=input$gen.con)
    grob$all1 <- grob$all[na.all]; grob$gg.all1 <- grob$gg.all[na.all]

  })
  # Add value legend to SHMs.
  # 'observeEvent' is able to avoid infinite cycles while 'observe' may cause such cycles. E.g. in the latter, 'is.null(grob$gg.all)' and 'grob$gg.all1 <- gg.all <- gg_2lgd()' would induce each other and form infinit circles.
  observeEvent(list(val.lgd=input$val.lgd, row=input$val.lgd.row, key=input$val.lgd.key, text=input$val.lgd.text, feat=input$val.lgd.feat), {
    
    validate(need(try(as.integer(input$val.lgd.row)==input$val.lgd.row & input$val.lgd.row>0), 'Legend key rows should be a positive integer!'))
    validate(need(try(input$val.lgd.key>0), 'Legend key size should be a positive numeric!'))
    validate(need(try(input$val.lgd.text>0), 'Legend text size should be a positive numeric!'))
    
    cat('Adding value legend... \n')
    if.con <- is.null(grob$gg.all)|is.null(sam())|is.null(input$val.lgd)|is.null(input$val.lgd.feat)|input$val.lgd==0
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    gg.all <- grob$gg.all1
    if ((input$val.lgd %% 2)==1) {

      gg.all <- gg_2lgd(gg.all=gg.all, sam.dat=sam(), tis.trans=input$tis, position.2nd='bottom', legend.nrow.2nd=input$val.lgd.row, legend.key.size.2nd=input$val.lgd.key, legend.text.size.2nd=input$val.lgd.text, add.feature.2nd=(input$val.lgd.feat=='Yes'))
      grob$all1 <- gg.all <- lapply(gg.all, function(x) { x+theme(legend.position="bottom") } )
      tmp <- tempfile()
      png(tmp); grob$all1 <- lapply(gg.all, ggplotGrob)
      dev.off(); if (file.exists(tmp)) do.call(file.remove, list(tmp))
    
    } else if ((input$val.lgd %% 2)==0) { 

      cat('Remove value legend... \n')
      grob$gg.all1 <- gg.all <- lapply(gg.all, function(x) { x+theme(legend.position="none") })
      tmp <- tempfile(); png(tmp); grob$all1 <- lapply(gg.all, ggplotGrob) 
      dev.off(); if (file.exists(tmp)) do.call(file.remove, list(tmp))

    }

  })

  observeEvent(input$col.cfm, { col.reorder$col.re <- 'Y' })
  # In "observe" and "observeEvent", if one code return (NULL), then all the following code stops. If one code changes, all the code renews.
  observe({
    if.con <- is.null(geneIn())|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all1)
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    
    output$shm <- renderPlot(width=as.numeric(input$width)/2*as.numeric(input$col.n), height=as.numeric(input$height)*length(input$dt_rows_selected), {
    
    if (col.reorder$col.re=='N') return()
    if.con <-  is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all1)
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    if (is.na(color$col[1])|length(color$col=="none")==0|input$color=="") return(NULL)
    r.na <- rownames(geneIn()[["gene2"]])
    grob.na <- names(grob$all1)
    # Select target grobs.
    # Use definite patterns and avoid using '.*' as much as possible. Try to as specific as possible.
    pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    grob.lis.p <- grob$all1[grepl(pat.all, grob.na)] # grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', names(grob.lis.p))); if (length(con)==0) return()
    cat('Plotting spatial heatmaps... \n')
    lay <- input$gen.con; ID <- gID$geneID; ncol <- input$col.n
    shm.lay <- lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, width=input$width, height=input$height, shiny=TRUE); shm <- shm.lay$shm
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
      lgd.lis <- grob$lgd.all; lgd.lis <- gg_lgd(gg.all=lgd.lis, sam.dat=sam(), tis.trans=input$tis, label=FALSE)
      lgd.lis <- gg_lgd(gg.all=lgd.lis, size.key=input$lgd.key.size*0.5, size.text.key=NULL, label.size=input$lgd.lab.size, row=input$lgd.row, sam.dat=sam(), tis.trans=input$tis, position.text.key='right', label=(input$lgd.label=='Yes'))
      if (input$lgd.w>0) {
  
        grob.lgd.lis <- lapply(lgd.lis, ggplotGrob)
        lgd.tr <- lapply(grob.lgd.lis, grobTree)
    # In 'arrangeGrob', if numbers in 'layout_matrix' are more than items in 'grobs', there is no difference. The width/height of each subplot is decided by 'widths' and 'heights'.
    lgd.arr <- arrangeGrob(grobs=lgd.tr, layout_matrix=matrix(seq_along(lgd.lis), ncol=1), widths=unit(1, "npc"), heights=unit(rep(1/length(lgd.lis)/input$lgd.ratio, length(lgd.lis)), "npc"))
    w.lgd <- (1-0.08)/(ncol+1)*input$lgd.w # Legend is reduced.
    png(paste0(tempdir(check=TRUE), '/tmp.png')); shm1 <- grid.arrange(cs.arr, shm, lgd.arr, ncol=3, widths=unit(c(0.08-0.005, 1-0.08-w.lgd, w.lgd), 'npc')); dev.off() } else { png(paste0(tempdir(check=TRUE), '/tmp.png')); shm1 <- grid.arrange(cs.arr, shm, ncol=2, widths=unit(c(0.08-0.005, 1-0.08), 'npc')); dev.off() }
    ggsave(paste0(tempdir(check=TRUE), '/shm.', input$ext), plot=shm1, device=input$ext, width=input$width/72, height=input$height/72, dpi=input$res, unit='in') }
    
    })

  })

  output$dld.shm <- downloadHandler(
    filename=function() { paste0('shm.', input$ext) },
    content=function(file) { file0 <- paste0(tempdir(check=TRUE), '/shm.', input$ext); 
    cat("Downloading 'shm' from", tempdir(check=TRUE), '...\n')
    file.copy(file0, file, overwrite=TRUE) }
  )

  observe({
  
    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$dt_rows_selected; input$tis; input$gen.con  
    updateRadioButtons(session, inputId='ext', label='File type:', choices=c('NA', "png", "jpg", "pdf"), selected=cfg$lis.par$shm.img['file.type', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId="ggly.but", label="Show animation:", choices=c("Yes", "No"), selected=cfg$lis.par$shm.anm['show', 'default'], inline=TRUE)
    updateRadioButtons(session, inputId="vdo.but", label="Show video:", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['show', 'default'], inline=TRUE)

  })

  observe({
   input$lgd.key.size; input$lgd.row; input$tis; input$lgd.label; input$lgd.lab.size
   updateRadioButtons(session, inputId="vdo.but", label="Show video:", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['show', 'default'], inline=TRUE)
  })


 output$shm.ui <- renderUI({

      box(title="Spatial Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, width=ifelse(input$hide.lgd=='No', 9, 12), height=NULL, 
      tabBox(title="", width=12, id='shm_all', selected='shm1', side='right', 
      tabPanel(title='Video', value='shm3', 
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '15%', '1%', '8%', '1%', '8%', '1%', '13%', '1%', '16%'), '', 
      radioButtons(inputId="vdo.but", label="Show video:", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['show', 'default'], inline=TRUE), '',
      numericInput(inputId='vdo.itvl', label='Transition time (s):', value=as.numeric(cfg$lis.par$shm.video['transition', 'default']), min=0.1, max=Inf, step=1, width=270), '',
      numericInput(inputId='vdo.height', label='Height:', value=as.numeric(cfg$lis.par$shm.video['height', 'default']), min=0.1, max=0.99, step=0.1, width=270), '',
      numericInput(inputId='vdo.width', label='Width:', value=as.numeric(cfg$lis.par$shm.video['width', 'default']), min=0.1, max=0.92, step=0.1, width=270), '',
      numericInput(inputId='vdo.res', label='Resolution (dpi):', value=as.numeric(cfg$lis.par$shm.video['dpi', 'default']), min=1, max=1000, step=5, width=270), '',
      radioButtons(inputId="vdo.val.lgd", label="Add values to legend:", choices=c("Yes", "No"), selected=cfg$lis.par$shm.video['value.legend', 'default'], inline=TRUE)

      )), uiOutput('video.dim'), textOutput('tran.vdo'), htmlOutput('ffm'), 
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar3"), uiOutput('video'), ""))),
      
      tabPanel(title='Animation', value='shm2', 
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '15%', '1%', '10%', '1%', '10%', '1%', '13%'), '', 
      radioButtons(inputId="ggly.but", label="Show animation:", choices=c("Yes", "No"), selected=as.numeric(cfg$lis.par$shm.img['show', 'default']), inline=TRUE), '',
      numericInput(inputId='t', label='Transition time (s):', value=as.numeric(cfg$lis.par$shm.anm['transition', 'default']), min=0.1, max=Inf, step=NA, width=270), '',
      uiOutput('anm.h'), '', uiOutput('anm.w'), '', uiOutput('dld.anm.but')
      )), textOutput('tran'), uiOutput('sld.fm'),
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar2"), htmlOutput("ggly"), ""))
      ),

      tabPanel(title="Image", value='shm1',  
      fluidRow(column(10, splitLayout(cellWidths=c('14%', '1%', '14%', '1%', '9%', '1%', '32%', '1%', '15%'),
      numericInput(inputId='height', label='Overall height:', value=as.numeric(cfg$lis.par$shm.img['height', 'default']), min=1, max=Inf, step=NA, width=170), '',
      numericInput(inputId='width', label='Overall width:', value=as.numeric(cfg$lis.par$shm.img['width', 'default']), min=1, max=Inf, step=NA, width=170), '',
      numericInput(inputId='col.n', label='Columns:', value=as.numeric(cfg$lis.par$shm.img['columns', 'default']), min=1, max=Inf, step=1, width=150), '',
      radioButtons(inputId="gen.con", label="Display by:", choices=c("Gene"="gene", "Condition"="con", "None"="none"), selected=cfg$lis.par$shm.img['display.by', 'default'], inline=TRUE), '', 
     radioButtons(inputId="pre.scale", label="Preserve.scale:", choices=c("Yes", "No"), selected=cfg$lis.par$shm.img['preserve.scale', 'default'], inline=TRUE)
      )),
      column(1,
      dropdownButton(inputId='dropdown', label='Color key', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      fluidRow(splitLayout(cellWidths=c('1%', '70%', '25%'), '', textInput("color", "Color scheme:", cfg$lis.par$shm.img['color', 'default'], placeholder=paste0('Eg: ', cfg$lis.par$shm.img['color', 'default']), width=200),
      actionButton("col.but", "Go", icon=icon("refresh")))), 
      radioButtons(inputId='cs.v', label='Color scale based on:', choices=c("Selected rows", "All rows"), selected=cfg$lis.par$shm.img['color.scale', 'default'], inline=TRUE)
      ))
      ), textOutput('h.w.c'), textOutput('msg.col'),

      fluidRow(splitLayout(cellWidths=c('100%'),  checkboxGroupInput(inputId="tis", label="Select tissues to be transparent:", choices='', selected='', inline=TRUE))),

      fluidRow(column(1, offset=0, style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      dropdownButton(inputId='dropdown', label='Download', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=800,
      fluidRow(splitLayout(cellWidths=c('27%', '1%', '14%', '1%', '11%', '1%', '17%', '1%', '14%'), 
      radioButtons(inputId='ext', label='File type:', choices=c('NA', "png", "jpg", "pdf"), selected=cfg$lis.par$shm.img['file.type', 'default'], inline=TRUE), '',
      numericInput(inputId='res', label='Resolustion (dpi):', value=as.numeric(cfg$lis.par$shm.img['dpi', 'default']), min=10, max=Inf, step=10, width=150), '',
      numericInput(inputId='lgd.w', label='Legend width:', value=as.numeric(cfg$lis.par$shm.img['legend.width', 'default']), min=0, max=1, step=0.1, width=150), '',
      numericInput(inputId='lgd.ratio', label='Legend aspect ratio:', value=as.numeric(cfg$lis.par$shm.img['legend.aspect.ratio', 'default']), min=0.0001, max=Inf, step=0.1, width=140), '', downloadButton("dld.shm", "Download")
      )))
      ),  
      column(1, offset=0, style='padding-left:40px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      dropdownButton(inputId='value.lgd', label='Value legend', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=500,
      fluidRow(splitLayout(cellWidths=c('1%', '25%', '1%', '13%', '1%', '17%', '1%', '14%', '1%', '28%'), '', 
      actionButton("val.lgd", "Add/Remove", icon=icon("refresh")), '',  
      numericInput(inputId='val.lgd.row', label='Rows:', value=as.numeric(cfg$lis.par$shm.img['value.legend.rows', 'default']), min=1, max=Inf, step=1, width=150), '',
      numericInput(inputId='val.lgd.key', label='Key size:', value=as.numeric(cfg$lis.par$shm.img['value.legend.key', 'default']), min=0.0001, max=1, step=0.01, width=150), '',
      numericInput(inputId='val.lgd.text', label='Text size:', value=as.numeric(cfg$lis.par$shm.img['value.legend.text', 'default']), min=0.0001, max=Inf, step=1, width=140), '',
      radioButtons(inputId='val.lgd.feat', label='Include feature:', choices=c('No', 'Yes'), selected=cfg$lis.par$shm.img['include.feature', 'default'], inline=TRUE)

      ))
      ))
      ), verbatimTextOutput('msg.shm'),
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar1"), plotOutput("shm", height='auto'), "")))

      ))
  
  })

  output$shms.o <- renderUI({
 
    if (is.null(svg.path1())) return(NULL)
    if (length(svg.path1()$svg.na)==1) return(NULL)
    svg.pa <- svg.path1()[['svg.na']]
    selectInput('shms.in', label='aSVG for legend:', choices=svg.pa, selected=svg.pa[1])

  })

  output$lgd <- renderPlot(width='auto', height="auto", {
    validate(need(try(as.integer(input$lgd.row)==input$lgd.row & input$lgd.row>0), 'Legend key rows should be a positive integer!'))
    validate(need(try(input$lgd.key.size>0&input$lgd.key.size<1), 'Legend key size should be between 0 and 1!'))
    svg.path <- svg.path1()
    if (is.null(svg.path1())|is.null(grob$lgd.all)|(length(svg.path$svg.na)>1 & is.null(input$shms.in))) return(ggplot())
      # Width and height in original SVG.
      if (length(svg.path$svg.na)>1) svg.na <- input$shms.in else svg.na <- 1

      w.h <- svg.df()[[svg.na]][['w.h']]
      w.h <- as.numeric(gsub("^(\\d+\\.\\d+|\\d+).*", "\\1", w.h)); r <- w.h[1]/w.h[2]
      if (is.na(r)) return(); cat('Plotting legend plot... \n')
      g.lgd <- grob$lgd.all[[svg.na]]; g.lgd <- g.lgd+coord_fixed(ratio=r); return(g.lgd)

  })

  output$lgd.ui <- renderUI({ 
    
    if (is.null(input$hide.lgd)) return(NULL) 
    if (input$hide.lgd=='Yes') return(NULL)
    box(title="Legend Plot", status="primary", solidHeader=TRUE, collapsible=TRUE, uiOutput('shms.o'),
    splitLayout(cellWidths=c("43%", "1%", '43%'),
    numericInput(inputId='lgd.row', label='Legend key rows:', value=as.numeric(cfg$lis.par$legend['key.row', 'default']), min=1, max=Inf, step=1, width=150), '',
    numericInput(inputId='lgd.key.size', label='Legend key size:', value=as.numeric(cfg$lis.par$legend['key.size', 'default']), min=0, max=1, step=0.02, width=150)
    ),
    splitLayout(cellWidths=c("43%", "1%", '43%'),
    radioButtons(inputId="lgd.label", label="Label feature:", choices=c("Yes", "No"), selected=cfg$lis.par$legend['label', 'default'], inline=TRUE), '',
    numericInput(inputId='lgd.lab.size', label='Label size:', value=as.numeric(cfg$lis.par$legend['label.size', 'default']), min=0, max=Inf, step=0.5, width=150)
    ),
    splitLayout(cellWidths=c("99%", "1%"), plotOutput("lgd"), ""), width=3) 

  })


  output$tran <- renderText({
    
    if (is.null(geneIn())|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)
    validate(need(try(input$t>=0.1), 'Transition time should be at least 0.1 second!'))

  })

  observeEvent(list(fineIn=input$fileIn, log=input$log, tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, pre.scale=input$pre.scale), {
    if (dir.exists('www/ggly/')) { system('rm -fr www/ggly/lib'); system('rm -f www/ggly/*html') } else dir.create('www/ggly/')
    if (dir.exists('html_shm/')) { system('rm -fr html_shm/lib'); system('rm -f html_shm/*html') } else dir.create('html_shm/')
    if (dir.exists('www/video/')) system('rm -fr www/video/*.mp4') else dir.create('www/video/')
  })

  observeEvent(list(width.ly=input$width.ly, height.ly=input$height.ly), {
    if (dir.exists('html_shm/')) { system('rm -fr html_shm/lib'); system('rm -f html_shm/*html') } else dir.create('html_shm/')
  })
  observeEvent(list(log=input$log, tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, pre.scale=input$pre.scale, ggly.but=input$ggly.but, gID.new=gID$new), {

    if (is.null(input$ggly.but)) return() 
    if (input$ggly.but=='No') return()
    if (is.null(geneIn())|is.null(gID$new)|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$gg.all1)|input$ggly.but=='No') return(NULL)
    if (length(color$col=="none")==0|input$color=="") return(NULL)

    withProgress(message="Animation: ", value=0, {
    incProgress(0.25, detail="preparing frames...") 
    gg.all <- grob$gg.all1; na <- names(gg.all)
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]; gg.all <- gg.all[na]
    for (i in seq_along(gg.all)) {

      na0 <- paste0(na[i], ".html")
      if (length(list.files('www/ggly/', na0))>0) next
      gly <- ggplotly(gg.all[[i]], tooltip='text') %>% layout(showlegend=FALSE)
      gly$sizingPolicy$padding <- 0
      cat('Animation: saving', na0, '\n')
      saveWidget(gly, na0, selfcontained=FALSE, libdir="lib")
      system(paste0('mv ', na0, ' www/ggly/'))

    }; if (!dir.exists('www/ggly/lib')) system('mv lib/ www/ggly/') else if (dir.exists('lib/')) system('rm -rf lib')

    })

  })

  output$sld.fm <- renderUI({
 
    if (is.null(grob$gg.all)|is.null(pat.all())|is.null(gID$geneID)) return(NULL) 
    gen.con.pat <- paste0('^', pat.all(), '_\\d+$') 
    sliderInput(inputId='fm', 'Frames', min=1, max=sum(grepl(gen.con.pat, names(grob$gg.all1))), step=1, value=1, animate=animationOptions(interval=input$t*10^3, loop=FALSE, playButton=icon('play'), pauseButton=icon('pause')))
  
  })

  # As long as the variable of 'reactive' is used in the 'ui.R', changes of elements in 'reactive' would cause chain change all the way to 'ui.R'. E.g. the change in "input$ggly.but=='No'" leads to changes in 'output$ggly' and 'ui.R', not necessarily changes in 'output$ggly' call changes in 'gly.url'.
  gly.url <- reactive({ 
    
    if (is.null(input$ggly.but)) return() 
    if (is.null(grob$gg.all1)|input$ggly.but=='No'|gID$geneID[1]=='none'|is.null(pat.all())) return(NULL)
    gg.all <- grob$gg.all1; na <- names(gg.all)
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
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
    
    gg.all <- grob$gg.all1; na <- names(gg.all)
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
    output$anm.w <- renderUI({

      numericInput(inputId='width.ly', label='Width:', value=w, min=1, max=Inf, step=NA, width=170)

    })
    output$anm.h <- renderUI({

      numericInput(inputId='height.ly', label='Height:', value=h, min=1, max=Inf, step=NA, width=170)

    })
  output$dld.anm.but <- renderUI({ downloadButton("dld.anm", "Download") })

  })

  
  observeEvent(list(log=input$log, tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, pre.scale=input$pre.scale, ggly.but=input$ggly.but, fm=input$fm), {
  
  output$ggly <- renderUI({ 

    if (input$ggly.but=='No'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
    withProgress(message="Animation: ", value=0, {
    incProgress(0.75, detail="plotting...")
    gly.url <- gly.url(); cat('Animation: plotting', gly.url, '\n')
    tags$iframe(src=gly.url, height=input$height.ly, width=input$width.ly, scrolling='auto')
  
    })

  })

  })

  anm.dld <- reactive({
    
    if (input$ggly.but=='No'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL) 
    withProgress(message="Downloading animation: ", value=0, {
    incProgress(0.1, detail="in progress...")
    gg.all <- grob$gg.all1; na <- names(gg.all)
    gg.na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]
    gg <- gg.all[gg.na]
    pro <- 0.1; for (i in seq_along(gg.na)) {
    incProgress(pro+0.2, detail=paste0('preparing ', gg.na[i], '.html...'))
    html_ly(gg=gg[i], cs.g=shm.bar(), tis.trans=input$tis, sam.uni=sam(), anm.width=input$width.ly, anm.height=input$height.ly, out.dir='.') }

   })

  })
  
  # This step leaves 'fil.na' in 'output$dld.anm' being a global variable.
  output$dld.anm <- downloadHandler( 
    # The rest code will run only after 'anm.dld()' is done.
    filename=function(){ anm.dld(); "html_shm.zip" },
    fil.na <- paste0(tempdir(), '/html_shm.zip'),
    content=function(fil.na){ cat('Downloading animation... \n'); zip(fil.na, 'html_shm/') }
  )

  output$video.dim <- renderUI({

    selectInput("vdo.dim", label="Fixed dimension:", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected=cfg$lis.par$shm.video['dimension', 'default'], width=110)

  })

  output$ffm <- renderText({
    
    ffm <- tryCatch({ system('which ffmpeg', intern=TRUE) }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
    if (!grepl('ffmpeg', ffm)) paste("<span style=\"color:red\">Error: \"ffmpeg\" is not detected!\"</span>")

  })

  observeEvent(list(log=input$log, tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, pre.scale=input$pre.scale, vdo.but=input$vdo.but, vdo.dim=input$vdo.dim, vdo.itvl=input$vdo.itvl, vdo.height=input$vdo.height, vdo.width=input$vdo.width, vdo.res=input$vdo.res, vdo.val.lgd=input$'vdo.val.lgd'), {

    if (is.null(input$vdo.but)) return(NULL) 
    if (input$vdo.but=='No'|is.null(pat.all())) return(NULL)
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
    validate(need(try(!is.na(input$vdo.itvl)&input$vdo.itvl>0), 'Transition time should be a positive numeric!'))
    validate(need(try(!is.na(input$vdo.height)&input$vdo.height>0&input$vdo.height<=0.99), 'Height should be between 0.1 and 0.99!'))
    validate(need(try(!is.na(input$vdo.width)&input$vdo.width>0&input$vdo.width<=0.92), 'Width should be between 0.1 and 0.92!'))
    validate(need(try(!is.na(input$vdo.res)&input$vdo.res>=1&input$vdo.res<=700), 'Resolution should be between 1 and 700!'))
    
    withProgress(message="Video: ", value=0, {
    incProgress(0.75, detail="in progress...")
    gg.all <- grob$gg.all1; na <- names(gg.all)
    pat <- paste0('^', pat.all(), '_\\d+$'); na <- na[grepl(pat, na)]
    gg.all1 <- gg.all[na]
    cat('Making video... \n')
    res <- input$vdo.res; dim <- input$vdo.dim
    if (dim %in% c('1280x800', '1280x1024', '1280x720')&res>450) res <- 450
    if (dim=='1920x1080'&res>300) res <- 300
    # selectInput("vdo.dim", label="Fixed dimension:", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected='640x480', width=110)
    vdo <- video(gg=gg.all1, cs.g=shm.bar(), sam.uni=sam(), tis.trans=input$tis, lgd.key.size=input$lgd.key.size, lgd.text.size=NULL, position.text.key='right', legend.value.vdo=(input$'vdo.val.lgd'=='Yes'), label=(input$lgd.label=='Yes'), label.size=input$lgd.lab.size, sub.title.size=8, bar.value.size=6, lgd.row=input$lgd.row, width=input$vdo.width, height=input$vdo.height, video.dim=dim, interval=input$vdo.itvl, res=res, out.dir='./www/video'); if (is.null(vdo)) return()
    cat('Presenting video... \n')
    incProgress(0.95, detail="Presenting video...")
    w.h <- as.numeric(strsplit(input$vdo.dim, 'x')[[1]])
    output$video <-renderUI({ tags$video(id="video", type="video/mp4", src="video/shm.mp4", width=w.h[1], height=w.h[2], controls="controls") })

    })

  })

  observe({

    geneIn(); input$adj.modInpath; input$A; input$p; input$cv1
    input$cv2; input$min.size; input$net.type
    input$measure; input$cor.abs; input$thr; input$mhm.v
    updateRadioButtons(session, "mat.scale", "Scale by: ", c("No", "Column", "Row"), "No", inline=TRUE)

  })


  observe({
  
    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$dt_rows_selected
    updateActionButton(session, inputId='mhm.but', label='Update', icon=icon("refresh"))
    #updateRadioButtons(session, inputId="mhm.but", label="Show plot:", choices=c("Yes", "No"), selected=cfg$lis.par$mhm['show', 'default'], inline=TRUE)

  })

  # Calculate whole correlation or distance matrix.
  cor.dis <- reactive({

    if (is.null(geneIn())|input$mhm.but=='No') return()
    if ((any(input$fileIn %in% cfg$na.cus) & is.null(geneIn()))|input$fileIn=="none") return(NULL)
    
    withProgress(message="Compute similarity/distance matrix: ", value = 0, {

      incProgress(0.5, detail="Please wait...")
      gene <- geneIn()[['gene1']]
      cat('Correlation/distance matrix...\n')
      if (input$measure=='correlation') {
      
        m <- cor(x=t(gene))
        if (input$cor.abs==TRUE) { m <- abs(m) }; return(m)

      } else if (input$measure=='distance') { return(-as.matrix(dist(x=gene))) }

    })

  })

  # Subset nearest neighbours for target genes based on correlation or distance matrix.
  submat <- reactive({
  
    if (input$fileIn=="None") return()
    if (is.null(cor.dis())|input$mhm.but=='No') return()
    gene <- geneIn()[["gene1"]]; rna <- rownames(gene)
    gen.tar<- gID$geneID; mat <- cor.dis()
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

      incProgress(0.5, detail="Please wait...")
      arg <- list(p=NULL, n=NULL, v=NULL)
      arg[names(arg) %in% input$thr] <- input$mhm.v
      if (input$measure=='distance' & input$thr=='v') arg['v'] <- -arg[['v']]
      if (!all(gen.tar %in% rownames(mat))) return()    
      cat('Subsetting nearest neighbors...\n')
      validate(need(try(ncol(gene)>4), 'The "sample__condition" variables in the Data Matrix are less than 5, so no coexpression analysis is applied!'))
      gen.na <- do.call(sub_na, c(mat=list(mat), ID=list(gen.tar), arg))
      if (any(is.na(gen.na))) return() 
      validate(need(try(length(gen.na)>=2), paste0('Only ', gen.na, ' selected!'))); return(gene[gen.na, ])

    })

  })
  mhm <- reactiveValues(hm=NULL)
  # Plot matrix heatmap.
  observe({

    if (input$mhm.but!=0) return() 
    if (is.null(submat())) return()
    gene <- geneIn()[["gene1"]]; rna <- rownames(gene)
    gen.tar <- gID$geneID; if (length(gen.tar)>1) return()
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="Plotting...")
      if (input$mat.scale=="Column") scale.hm <- 'column' else if (input$mat.scale=="Row") scale.hm <- 'row' else scale.hm <- 'no'
      cat('Initial matrix heatmap...\n')
      mhm$hm <- matrix_hm(ID=gen.tar, data=submat(), scale=scale.hm, main='Target Genes and Their Nearest Neighbours', title.size=10, static=FALSE)

    })

  })
  hmly <- eventReactive(input$mhm.but, {
    
    #if (is.null(submat())|input$mhm.but=='No') return()
    if (is.null(submat())) return()
    gene <- geneIn()[["gene1"]]; rna <- rownames(gene)
    gen.tar<- gID$geneID
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="Plotting...")
      if (input$mat.scale=="Column") scale.hm <- 'column' else if (input$mat.scale=="Row") scale.hm <- 'row' else scale.hm <- 'no'
      cat('Matrix heatmap...\n')
      matrix_hm(ID=gen.tar, data=submat(), scale=scale.hm, main='Target Genes and Their Nearest Neighbours', title.size=10, static=FALSE)

    })

  })

  output$HMly <- renderPlotly({ 
   
    if (is.null(input$dt_rows_selected)) return()
    if (is.null(gID$geneID)|is.null(submat())) return()
    if (gID$geneID[1]=='none'|is.na(gID$geneID[1])) return()
    if (input$mhm.but!=0) hmly() else if (input$mhm.but==0) mhm$hm else return() 
  
  })

  adj.mod <- reactive({ 

    if (input$fileIn=="customComputedData" & !is.null(input$adj.modInpath)) {

      name <- input$adj.modInpath$name; path <- input$adj.modInpath$datapath
      path1 <- path[name=="adj.txt"]; path2 <- path[name=="mod.txt"]

      withProgress(message="Loading: ", value = 0, {
        incProgress(0.5, detail="adjacency matrix and module definition.")
        adj <- fread(path1, sep="\t", header=TRUE, fill=TRUE); c.na <- colnames(adj)[-ncol(adj)]
        r.na <- as.data.frame(adj[, 1])[, 1];  adj <- as.data.frame(adj)[, -1] 
        rownames(adj) <- r.na; colnames(adj) <- c.na

        mcol <- fread(path2, sep="\t", header=TRUE, fill=TRUE); c.na <- colnames(mcol)[-ncol(mcol)]
        r.na <- as.data.frame(mcol[, 1])[, 1]; mcol <- as.data.frame(mcol)[, -1] 
        rownames(mcol) <- r.na; colnames(mcol) <- c.na

      }); return(list(adj=adj, mcol=mcol))

    }

  })

    #gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    adj.mods <- reactiveValues(lis=NULL)
    observe({
    
      if (input$fileIn=="None") return()
      if (is.null(submat())|input$cpt.nw!=0|length(gID$geneID)>1) return()

    if (input$fileIn=="customData"|any(input$fileIn %in% cfg$na.def)) {

      gene <- geneIn()[["gene1"]]; if (is.null(gene)) return()
      type <- input$net.type; sft <- if (type=='distance') 1 else 6
      withProgress(message="Computing: ", value = 0, {
        incProgress(0.3, detail="adjacency matrix.")
        incProgress(0.5, detail="topological overlap matrix.")
        incProgress(0.1, detail="dynamic tree cutting.")
        cat('Initial adjacency matrix and modules...\n')
        adj.mods$lis <- adj_mod(data=submat(), type=type, minSize=input$min.size, dir=NULL)

      })

    }
    
    })

  # er <- eventReactive(exp, {}). If its reactive value "er()" is called before eventReactive is triggered, the code execution stops where "er()" is called.
  observeEvent(input$cpt.nw, { 
    #gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    if (is.null(submat())|input$cpt.nw==0) return()
    if (input$fileIn=="customData"|any(input$fileIn %in% cfg$na.def)) {

      gene <- geneIn()[["gene1"]]; if (is.null(gene)) return()
      type <- input$net.type; sft <- if (type=='distance') 1 else 6
      withProgress(message="Computing: ", value = 0, {
        incProgress(0.3, detail="adjacency matrix.")
        incProgress(0.5, detail="topological overlap matrix.")
        incProgress(0.1, detail="dynamic tree cutting.")
        cat('Adjacency and modules... \n')
        adj.mods$lis <- adj_mod(data=submat(), type=type, minSize=input$min.size)

      })

    }

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
  
   # geneIn(); gID$geneID; input$adj.in; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
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
 
    input$cpt.nw; if (input$fileIn=="None") return()
    if (input$fileIn=='customComputedData' & is.null(geneIn())) return()
    # if (input$adj.in=="None") return(NULL)
    if (input$fileIn=="customComputedData") { adj <- adj.mod()[['adj']]; mods <- adj.mod()[['mcol']] } else if (input$fileIn=="customData"|input$fileIn %in% cfg$na.def) { 
      adj <- adj.mods$lis[['adj']]; mods <- adj.mods$lis[['mod']]
    }
    if (input$fileIn=='customComputedData') gene <- geneIn()$gene2 else gene <- submat()
    if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    lab <- mods[, input$ds][rownames(gene)==input$gen.sel]
    validate(need(try(length(lab)==1 & !is.na(lab) & nrow(mods)==nrow(gene)), 'Click "Update" to display new network!'))
    if (length(lab)>1|is.na(lab)) return() # When input$fileIn is changed, gene is changed also, but mods is not since it is controled by observeEvent.
    validate(need(try(lab!='0'), 'Warning: the selected gene is not assigned to any module. Please select a different one or adjust the "Minmum module size"!'))
    idx.m <- mods[, input$ds]==lab; adj.m <- adj[idx.m, idx.m]; gen.na <- colnames(adj.m) 
    idx.sel <- grep(paste0("^", input$gen.sel, "$"), gen.na); gen.na[idx.sel] <- paste0(input$gen.sel, "_target")
    colnames(adj.m) <- rownames(adj.m) <- gen.na
    withProgress(message="Computing network:", value=0, { 
      incProgress(0.8, detail="making network data frame")
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
        
      }; ann <- geneIn()[['gene3']]
      if (!is.null(ann)) node <- cbind(node, title=ann[node$id, ], borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      if (is.null(ann)) node <- cbind(node, borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      net.lis <- list(node=node, link=link1, adjs=adjs, lins=lins)

    }); net.lis

  })
  # The order of reactive expression matters so "updateSelectInput" of "adj.in" should be after visNet().
  observe({
 
    if (input$fileIn=="None") return()
    geneIn(); gID$geneID; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    input$gen.sel; input$measure; input$cor.abs; input$thr; input$mhm.v; input$cpt.nw
     #if ((input$adj.in==1 & is.null(visNet()[["adjs1"]]))|(input$cpt.nw!=cfg$lis.par$network['max.edges', 'default'] & is.null(visNet()[["adjs1"]]))) { updateSelectInput(session, "adj.in", "Adjacency threshold:", sort(seq(0, 1, 0.002), decreasing=TRUE), visNet()[["adjs"]]) } else if (!is.null(visNet()[["adjs1"]])) updateSelectInput(session, "adj.in", "Adjacency threshold:", sort(seq(0, 1, 0.002), decreasing=TRUE), visNet()[["adjs1"]])
     lins <- visNet()[["lins"]]
     if (input$adj.in==1|is.null(lins)|is.numeric(lins)) updateSelectInput(session, "adj.in", "Adjacency threshold (the  smaller, the more edges):", sort(seq(0, 1, 0.002), decreasing=TRUE), as.numeric(visNet()[["adjs"]])) 
  
  })
  output$bar.net <- renderPlot({  

    #if (input$adj.in=="None"|input$cpt.nw=="No") return(NULL)
    if (input$adj.in=="None") return(NULL)
    if (length(color.net$col.net=="none")==0) return(NULL)
    gene <- geneIn()[["gene1"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    if(input$col.but.net==0) color.net$col.net <- colorRampPalette(col_sep(cfg$lis.par$network['color', 'default']))(len.cs.net) # color.net$col.net is changed alse outside renderPlot, since it is a reactive value.

      withProgress(message="Color scale: ", value = 0, {
      incProgress(0.25, detail="Preparing data. Please wait.")
      incProgress(0.75, detail="Plotting. Please wait.")
      node <- visNet()[["node"]]; if (is.null(node)) return()
      node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
      cat('Network bar... \n')
      cs.net <- col_bar(geneV=v.net, cols=color.net$col.net, width=1); return(cs.net) # '((max(v.net)-min(v.net))/len.cs.net)*0.7' avoids bar overlap.

      })

  })

  observeEvent(visNet(), {

    output$edge <- renderUI({ 

      if (input$adj.in=="None"|is.null(visNet())) return(NULL)
      if (input$fileIn=="none"|(input$fileIn=="Your own" & is.null(geneIn()))|
      input$gen.sel=="None") return(NULL)
      cat('Remaining edges... \n')
      span(style = "color:black;font-weight:NULL;", HTML(paste0("Remaining edges: ", dim((visNet()[["link"]]))[1])))

    })

  })
  vis.net <- reactive({ 
    
    #if (input$adj.in=="None"|input$cpt.nw=="No") return(NULL)
    if (input$adj.in=="None") return(NULL)
    gene <- geneIn()[["gene1"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    withProgress(message="Network:", value=0.5, {
    incProgress(0.3, detail="prepare for plotting.")
    # Match colours with gene connectivity by approximation.
    node <- visNet()[["node"]]; if (is.null(node)) return() 
    node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
    col.nod <- NULL; for (i in node$value) {

      ab <- abs(i-v.net); col.nod <- c(col.nod, color.net$col.net[which(ab==min(ab))[1]])

    }; node$color <- col.nod
    cat('Network... \n')
    visNetwork(node, visNet()[["link"]], height="300px", width="100%", background="", main=paste0("Network Module Containing ", input$gen.sel), submain="", footer= "") %>% visIgraphLayout(physics=FALSE, smooth=TRUE) %>% visOptions(highlightNearest=list(enabled=TRUE, hover=TRUE), nodesIdSelection=TRUE)

    })
    
  })

  output$vis <- renderVisNetwork({

    if (is.null(input$dt_rows_selected)) return()
    if (input$fileIn=="none"|is.null(vis.net())) return(NULL)
    # if (input$cpt.nw=="No") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="plotting.")
      cat('Rendering network...\n'); vis.net()

    })

  })

  onStop(function() { 

    if (dir.exists('www/ggly/')) {  cat("Removing animation files in 'www/ggly/' ... \n"); system('rm -fr www/ggly/lib/*'); system('rm -f www/ggly/*html') }
    if (dir.exists('html_shm/lib')) {  cat("Removing animation files in 'html_shm/lib/' ... \n"); system('rm -fr html_shm/lib/*'); system('rm -f html_shm/*html') }
    if (dir.exists('www/video/')) {  cat("Removing video file in 'www/video/' ... \n"); system('rm -fr www/video/.mp4*') }
    cat("Session stopped\n")

  })

})



