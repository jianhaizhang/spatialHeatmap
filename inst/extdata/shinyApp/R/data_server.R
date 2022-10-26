# Module for processing data.
data_server <- function(id, sch, lis.url, ids, deg = FALSE, upl.mod.lis, scell.mod.lis=NULL, session) {
  moduleServer(id, function(input, output, session) {
    ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg
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
    cat('Import data, row metadata, targets files ... \n')
    fileIn <- ipt$fileIn; geneInpath <- ipt$geneInpath; dimName <- ipt$dimName
    svgInpath1 <- ipt$svgInpath1; svgInpath2 <- ipt$svgInpath2
    # validate/need: Suspend the process and no return, can cause errors, so "if" is choosen.
    if (fileIn=="none") return(NULL)
    withProgress(message="Loading data: ", value = 0, {
    if (any(fileIn %in% cfg$na.def)) {
      incProgress(0.5, detail="loading matrix, please wait ...")
      dat.na <- cfg$dat.def[fileIn]
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
    }; dimNa <- dimName
    if (fileIn %in% cfg$na.cus) {
      if (is.null(dimNa)) return()
      if (is.null(geneInpath) | dimNa=="None") return()
      incProgress(0.25, detail="importing matrix, please wait ...")
      geneInpath <- geneInpath$datapath; targetInpath <- ipt$target$datapath; metInpath <- ipt$met$datapath
      # Keep replicates unchaged, and compared with targets/metadata files.
      df.upl <- fread_df(read_fr(geneInpath), isRowGene=(dimNa=='Row'), rep.aggr=NULL)
      df.rep <- df.upl$df.rep; df.met <- df.upl$df.met
      if (!is.null(targetInpath)) {
        df.tar <- read_fr(targetInpath)
        # If errors detected, modalDialog terminates the app, while validate only stops the current step. 
        if (nrow(df.tar) != ncol(df.rep)) showModal(modal(msg = 'Ensure "columns" in the data matrix corresponds with "rows" in the targets file respectively!'))
        validate(need(try(nrow(df.tar) == ncol(df.rep)), 'Ensure "columns" in the data matrix corresponds with "rows" in the targets file respectively!'))
        # Check feature/variable columns in targets file.
        cna.tar <- colnames(df.tar) <- tolower(colnames(df.tar))
        if (all(c('feature', 'variable') %in% cna.tar)) cna <- paste0(df.tar$feature, '__', df.tar$variable) else if ('feature' %in% cna.tar) cna <- paste0(df.tar$feature, '__', 'con')
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
      if ((!is.null(svgInpath1) | !is.null(svgInpath2))) 
      cat('Done! \n'); return(df.upl) 
    }

    })

  })
  thr <- reactiveValues()
  output$msg.sig.thr <- renderText({
    sig.max <- input$sig.max; sig.min <- input$sig.min
    er.wa1 <- er.wa2 <- 1
    if (sig.min!='') er.wa1 <- check(sig.min, as.numeric)
    if (sig.max!='') er.wa2 <- check(sig.max, as.numeric)
    if (!is.numeric(c(er.wa1, er.wa2))) return('Signal thresholds must be numeric!') 
    if (!is.null(thr$max)) return(thr$max)
    if (!is.null(thr$min)) return(thr$min)
    if (!is.null(thr$vs)) return(thr$vs)
  })
  # Transform data.
  geneIn1 <- reactive({
    cat('Scale/tranform data ... \n')
    norm.meth <- input$normDat; scale.dat <- input$scaleDat
    df.aggr.thr <- NULL
    if (deg == FALSE) {
      # A default single-cell example should be considered in the "if" statement.
      if (ipt$fileIn!='customSingleCellData') { # Bulk data. 
        if (is.null(geneIn0())|is.null(norm.meth)|is.null(scale.dat)|is.null(input$log)) return() 
        df.aggr <- df.aggr.tran <- geneIn0()[['df.aggr']]
        # Normalization
        if (grepl('^CNF-', norm.meth)) {
          norm.fun <- 'CNF'
          par.lis <- list(method=sub('.*-', '', norm.meth))
        } else { norm.fun <- norm.meth; par.lis <- NULL }
        if (all(df.aggr.tran==round(df.aggr.tran))) df.aggr.tran <- norm_data(data=df.aggr.tran, norm.fun=norm.fun, parameter.list=par.lis, log2.trans=FALSE)
      } else if (ipt$fileIn=='customSingleCellData') { # Single-cell data.
        if (is.null(scell.mod.lis)) return()
        df.aggr <- df.aggr.tran <- scell.mod.lis$df.lis()$df.aggr
        if (is.null(df.aggr)) return()
    } else return()

    if (input$log=='log2') { 
      g.min <- min(df.aggr.tran)
      if (g.min<0) df.aggr.tran <- df.aggr.tran-g.min+1; if (g.min==0) df.aggr.tran <- df.aggr.tran+1; df.aggr.tran <- log2(df.aggr.tran)
    } else if (input$log=='exp2') df.aggr.tran <- 2^df.aggr.tran

    input$sig.but; isolate({
      cat('Threshold ... \n')
      thr$max <- thr$min <- thr$vs <- NULL
      max.v <- max(df.aggr.tran); min.v <- min(df.aggr.tran)
      max.v1 <- min.v1 <- NULL
      sig.max <- input$sig.max; sig.min <- input$sig.min
      if (sig.min!='') { er.wa <- check(sig.min, as.numeric)
        if (is.numeric(er.wa)) { min.v1 <- er.wa
          if (min.v1 >= max.v) thr$min <- paste0('The min threshold should be < the max un-scaled value: ', max.v, '!')
        }
      }
      if (sig.max!='') { er.wa <- check(sig.max, as.numeric)
        if (is.numeric(er.wa)) { max.v1 <- er.wa
          if (min.v >= max.v1) thr$max <- paste0('The max threshold should be > the min un-scaled value: ', min.v, '!')
        }
      }
      # validate/need can give rise to errors to downstream processes while if condition not.
      # validate(need(min.v < max.v, ''))
      if (!is.null(max.v1)) if (max.v1 < max.v) max.v <- max.v1 
      if (!is.null(min.v1)) if (min.v1 > min.v) min.v <- min.v1
      if (!is.null(max.v1) & !is.null(min.v1)) {
        if (max.v1 <= min.v1) thr$vs <- 'The max threshold must be > min threshold!'
      }
      if (max.v <= min.v) thr$vs <- 'The max threshold must be > min threshold!'
      if (max.v > min.v) {
      df.aggr.tran[df.aggr.tran >= max.v] <- max.v
      df.aggr.tran[df.aggr.tran <= min.v] <- min.v
      }
    })
    # The only difference between "df.aggr.thr" and "df.aggr.tran" is scaling. The former is used to indicate un-scaled range in the colour bar title.
    df.aggr.thr <- df.aggr.tran
    # Scale by row/column
    if (scale.dat=='Row') { df.aggr.tran <- t(scale(t(df.aggr.tran))) } else if (scale.dat=='Selected') {
      # validate(need(!is.null(ids$sel), ''))
      idx.sel <- rownames(df.aggr.tran) %in% ids$sel
      df.sel <- df.aggr.tran[idx.sel, ]
      df.aggr.tran[idx.sel, ] <- scale_all(df.sel) 
    } else if (scale.dat=='All') { df.aggr.tran <- scale_all(df.aggr.tran) }
    } else if (deg == TRUE) {   
      if (is.null(geneIn0())) return()
      df.aggr <- df.aggr.tran <- geneIn0()[['df.aggr']]
    }
    if (ipt$fileIn!='customSingleCellData') { 
      df.met <- geneIn0()[['df.met']]; df.rep <- geneIn0()[['df.rep']]
      con.na <- geneIn0()$con.na
    } else {
      if (is.null(scell.mod.lis)) return()
      df.lis <- scell.mod.lis$df.lis()
      df.met <- df.lis$df.met; df.rep <- df.lis$df.rep
      con.na <- df.lis$con.na
    }
    cat('Done! \n'); return(list(df.aggr=df.aggr, df.aggr.tran=df.aggr.tran, df.aggr.thr=df.aggr.thr, df.met=df.met, df.rep=df.rep, con.na=con.na))

  })

  output$col.order <- renderUI({
    if (is.null(geneIn1())) return()
    col.nas <- colnames(geneIn1()[['df.aggr.tran']])
    column(12, actionButton("col.cfm", "Confirm", icon=icon("sync")), 
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
    df.aggr <- gene.lis[['df.aggr']]
    if (is.null(df.aggr)) return()
    df.aggr.tran <- gene.lis[['df.aggr.tran']]
    df.aggr.thr <- gene.lis$df.aggr.thr
    df.met <- gene.lis[['df.met']]; df.rep <- gene.lis[['df.rep']]
    input$fil.but; lis.url
    if (deg==FALSE) {
      if (is.null(df.aggr.thr)) return()
      # Filter un-scaled values.
      df2fil <- df.aggr.thr
      # When ipt$fileIn changes, col.reorder$col.na.re of last session still exists and causes errors.
      if (!identical(sort(col.reorder$col.na.re), sort(colnames(df2fil)))) return()
    } else df2fil <- df.rep
    # Input variables in "isolate" will not triger re-excution, but if the whole reactive object is trigered by "input$fil.but" then code inside "isolate" will re-excute.
    scale.dat <- input$scaleDat
    isolate({ 
      se <- SummarizedExperiment(assays=list(expr=as.matrix(df2fil)), rowData=df.met)
      ann.col <- NULL
      if (!is.null(df.met)) if (ncol(df.met)>0 & 'metadata' %in% colnames(df.met)) ann.col <- 'metadata'
      # If scaled by row, sd is 1, mean is 0, cv is Inf.
      if (deg == FALSE) CVs <- c(ifelse(scale.dat %in% c('Row', 'Selected'), -Inf, fil$CV1), ifelse(scale.dat %in% c('Row', 'Selected'), Inf, fil$CV2)) else CVs <- c(fil$CV1, fil$CV2)
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV = CVs, dir=NULL, verbose=FALSE)
      if (fil$P!=0 | fil$A!=0 | fil$CV1!=-10000 | fil$CV2!=10000) showNotification("Filtering is applied on un-scaled values.", duration=2)
      # In case all rows are filtered, the app continues to work without refreshing after the filter parameters are reduced.
      if (nrow(se)==0) { validate(need(try(nrow(se)>0), 'All rows are filtered out!')); return() }
      # se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=CVs, dir=NULL, verbose=FALSE)
      df2fil <- as.data.frame(assay(se), stringsAsfactors=FALSE)
      colnames(df2fil) <- make.names(colnames(df2fil))
      rna.fil <- rownames(df2fil)
      df.aggr <- df.aggr[rna.fil, , drop=FALSE]
      df.aggr.tran <- df.aggr.tran[rna.fil, , drop=FALSE]
      df.met <- as.data.frame(rowData(se))[, , drop=FALSE]
    })
    cat('Preparing data matrix ... \n')
    if (deg == FALSE) {
    if (is.null(sear$id)) rows <- seq_len(nrow(df.aggr.tran)) else rows <- sear$id
    if (length(rows)==1 & rows[1]==as.numeric(cfg$lis.par$data.matrix['row.selected', 'default'])) rows <- seq_len(nrow(df.aggr.tran))
    df.aggr.tran.order <- df.aggr.tran[rows, col.reorder$col.na.re, drop=FALSE]
    df.aggr.thr <- df2fil
    df.aggr.tran <- df.aggr.tran[rows, , drop=FALSE]
    df.aggr.thr <- df.aggr.thr[rows, , drop=FALSE]
    df.aggr <- df.aggr[rows, , drop=FALSE]
    df.met <- df.met[rows, , drop = FALSE]
    } else {
      # Only "df.rep" is used.
      df.aggr.tran.order <- df.aggr.tran
      df.rep <- df2fil  
    }
    if (!is.null(df.aggr.thr)) df.aggr.thr <- as.data.frame(as.matrix(df.aggr.thr))
    # df.aggr.tran is used for SHMs.
    cat('Done! \n'); return(list(df.aggr = as.data.frame(as.matrix(df.aggr)), df.aggr.tran = as.data.frame(as.matrix(df.aggr.tran)), df.aggr.tran.order = as.data.frame(as.matrix(df.aggr.tran.order)), df.aggr.thr = df.aggr.thr, df.met=df.met, df.rep = as.data.frame(as.matrix(df.rep)), con.na=gene.lis$con.na))

  })
  dt.shm <- reactive({
    cat('Preparing data matrix ... \n')
    if (is.null(geneIn())|length(ids$sel)==0) return()
    if ((ipt$fileIn %in% cfg$na.cus & is.null(geneIn()) & is.null(scell.mod.lis))|ipt$fileIn=="none") return()
    withProgress(message="Data table: ", value = 0, {
      incProgress(0.5, detail="Preparing data matrix, please wait ...")
      if (ipt$fileIn!="none") {
      df.all <- geneIn()
      if (deg == FALSE) gene.dt <- cbind.data.frame(df.all[["df.met"]][, , drop=FALSE], df.all[["df.aggr.tran"]][, , drop=FALSE], stringsAsFactors=FALSE) else gene.dt <- cbind.data.frame(df.all[["df.met"]][, , drop=FALSE], df.all[["df.rep"]][, , drop=FALSE], stringsAsFactors=FALSE)
      # Remove '__con' only in the data table, not in the downstream (shm, network).
      if (df.all$con.na==FALSE) colnames(gene.dt) <- sub('__con$', '', colnames(gene.dt))
   }
   cat('Done!\n'); gene.dt
    })
  })

  if (deg==FALSE) output$selProf <- renderPlot(height=800, {
    cat('Profile of selected genes ... \n')
    gene.dt <- dt.shm(); df.all <- geneIn()
    df.aggr.thr <- df.all$df.aggr.thr
    if (is.null(gene.dt)|length(ids$sel)==0|is.null(df.aggr.thr)) return()
    validate(need(length(ids$sel)<=100, 'Due to space limitation, profiles of 100+ genes are not plotted!'))
    dt.sel <- gene.dt[ids$sel, , drop=FALSE]
    dt.sel <- dt.sel[, !colnames(dt.sel) %in% c('metadata', 'link'), drop=FALSE]
    df.aggr.thr <- df.aggr.thr[ids$sel, , drop=FALSE]
    scale.dat <- input$scaleDat; if (is.null(scale.dat)) return()
    if (!is.null(scale.dat)) if (scale.dat=='No') title <- 'No scaling' else if (scale.dat=='Row') title <- 'Scaled by row' else if (scale.dat=='Selected') title <- 'Scaled across selected genes' else if (scale.dat=='All') title <- 'Scaled across all genes' else title <- '' 
    # print(list('sel.prof', dim(df.all$df.aggr), dim(gene.dt)))
    g1 <- graph_line(dt.sel, y.title=paste0(title, ' (', round(min(df.aggr.thr), 2), '-', round(max(df.aggr.thr), 2), ')'), text.size=12)
    dt.sel.ori <- df.all$df.aggr[ids$sel, , drop=FALSE]
    if (df.all$con.na==FALSE) colnames(dt.sel.ori) <- sub('__con$', '', colnames(dt.sel.ori))
    g2 <- graph_line(dt.sel.ori, y.title=paste0('Un-scaled values (', round(min(dt.sel.ori), 2), '-', round(max(dt.sel.ori), 2), ')'))
    cat('Done! \n'); grid.arrange(g1, g2, nrow=2)
  })
  if (deg==FALSE) output$dtSel <- renderDataTable({
    cat('Preparing selected data matrix ... \n')
    gene.dt <- dt.shm()
    if (is.null(gene.dt)|length(ids$sel)==0) return()
    # Decimals.
    df.num <- gene.dt[, !colnames(gene.dt) %in% c('metadata', 'link'), drop=FALSE]
    geneIn <- geneIn(); cna <- colnames(geneIn$df.aggr.tran)
    if (geneIn$con.na==FALSE) cna <- sub('__con$', '', cna)
    if (all(df.num==round(df.num))) deci <- 0 else deci <- 2
    dt.sel <- gene.dt[ids$sel, , drop=FALSE]
    # Tooltip on metadata.
    col1 <- list(list(targets = c(1), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    # In case no metadata column.
    if (colnames(dt.sel)[1]!='metadata') col1 <- NULL
    dtab <- datatable(dt.sel, selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs=col1, dom='t', fixedColumns = list(leftColumns=2)), 
   class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
   formatRound(cna, deci)
    cat('Done! \n'); dtab
  })
  dt.sel <- reactiveValues(val='none')
  if (deg==FALSE) output$dtAll <- renderDataTable({
    cat('Preparing complete data matrix ... \n')
    gene.dt <- dt.shm(); page.h <- input$page
    if (is.null(gene.dt)|!is.numeric(page.h)) return()
    # Decimals.
    df.num <- gene.dt[, !colnames(gene.dt) %in% c('metadata', 'link'), drop=FALSE]
    geneIn <- geneIn(); cna <- colnames(geneIn$df.aggr.tran)
    if (geneIn$con.na==FALSE) cna <- sub('__con$', '', cna)
    if (all(df.num==round(df.num))) deci <- 0 else deci <- 2
    # Tooltip on metadata.
    col1 <- list(list(targets = c(1), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    if (colnames(gene.dt)[1]!='metadata') col1 <- NULL
    dtab <- datatable(gene.dt, selection=list(mode="multiple", target="row", selected=dt.sel$val), escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=page.h, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs=col1), 
   class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
   formatRound(cna, deci)
    cat('Done! \n'); dtab
  })

  if (deg==TRUE) output$dtRep <- renderDataTable({
    cat('Preparing complete data matrix ... \n')
    gene.dt <- dt.shm(); if (is.null(gene.dt)) return()
    col1 <- list(list(targets = c(1), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    # In case no metadata column.
    if (colnames(gene.dt)[1]!='metadata') col1 <- NULL
    dtab <- datatable(gene.dt, selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs=col1, dom = 't', fixedColumns = list(leftColumns=2)), 
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
  observeEvent(list(input$tran.scale.but.sel, input$tran.scale.but.prof), ignoreInit=TRUE, {
    updateTabsetPanel(session, "dtab.shm", selected='dTabAll')
  })

  observeEvent(scell.mod.lis$covis.man$match.mod.lis$but.match$val, ignoreInit=TRUE, {
    updateTabsetPanel(session, "dtab.shm", selected='dTabScell')
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
  sig.but <- reactive({ input$sig.but })
  onBookmark(function(state) { state })
  return(list(geneIn0 = geneIn0, geneIn1 = geneIn1, geneIn = geneIn, sear = sear, ipt.dat = ipt.dat, col.reorder = col.reorder, col.cfm = col.cfm, col.na = col.na, scaleDat=scaleDat, log = log, A = A, P = P, CV1 = CV1, CV2 = CV2, search.but = search.but, sig.but=sig.but))
  })

}
