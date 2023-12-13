# Module for processing data.
data_server <- function(id, sch, lis.url, ids, upl.mod.lis, deg.mod.lis=NULL, scell.mod.lis=NULL, shm.mod=NULL, session, parent=NULL) {
  moduleServer(id, function(input, output, session) {
  ipt <- upl.mod.lis$ipt; covis.pa <- upl.mod.lis$covis.pa
  ns <- session$ns; cfg <- upl.mod.lis$cfg
  con.na <- reactiveValues(v=FALSE)
  con.na.cell <- reactiveValues(v=FALSE)
  # test <- reactive({ req(''); 1})
  # Error: observe({ print(test()) })
  # "test2" will be suspended: observe({ print('test1'); test(); print('test2') })
  observeEvent(list(deg.mod.lis$input$eSHMBut), { # Should not be merged with below.
   fileIn <- ipt$fileIn; req(!dat.no %in% fileIn)
   query.res <- check_exp(check_obj(deg.mod.lis$query.res()))
   if (!is(query.res, 'character') & !is.null(query.res) & !grepl(na.sgl, fileIn)) {
     cho <- c("Complete"='all', 'Spatial enrichment'='enr')
     sel <- 'enr' 
     updateSelectInput(session, inputId='datIn', choices=cho, selected=sel)
   }
  })
  observeEvent(list(scell.mod.lis$covis.auto$but.covis, scell.mod.lis$covis.man$match.mod.lis$but.match$val), { # Should not be merged with above.
   fileIn <- ipt$fileIn; req(!dat.no %in% fileIn)
   if (grepl(na.sgl, fileIn)) {
     cho <- c("Co-visualization"='covis'); sel <- 'covis'
     updateSelectInput(session, inputId='datIn', choices=cho, selected=sel)
   }
  })
  observeEvent(list(ipt$fileIn), {
   updateSelectInput(session, inputId='datIn', choices=c('Complete'='all'), selected='all') 
  })
  observeEvent(list(ipt$fileIn, dat(), input$datIn), {
    dat <- dat(); lis.par <- cfg$lis.par
    if (!check_obj(list(dat, lis.par))) return()
    def <- lis.par$data.matrix['norm', 'default']
    assay <- assay(dat$se.rep)
    if (all(round(assay)==assay)) sel <- def else sel <- 'None' 
    updateSelectInput(session, inputId='normDat', selected=sel)
  })
  # Import data, row metadata, targets file, aggregate replicates.
  dat <- reactive({
    cat('Import data, row metadata, targets files ... \n')
    fileIn <- ipt$fileIn; geneInpath <- ipt$geneInpath; dimName <- ipt$dimName
    svgInpath1 <- ipt$svgInpath1; svgInpath2 <- ipt$svgInpath2
    req(!dat.no %in% fileIn)
    if (na.cus.covis %in% fileIn) {
      if (!check_obj(list(covis.pa$dat, !is.null(covis.pa$svg1)|!is.null(covis.pa$svg2)))) return()
    }
    # validate/need: Suspend the process and no return, can cause errors, so "if" is choosen.
    if (grepl(na.sgl, fileIn)) {
     sce.res <- scell.mod.lis$sce.res
     if (is.null(sce.res)) req('')
     sce.lis <- sce.res()$sce.lis
     message('Done!'); return(sce.lis)
    }
    withProgress(message="Loading data: ", value = 0, {
    if (fileIn %in% cfg$na.def & !grepl(na.sgl, fileIn)) {
      incProgress(0.5, detail="loading matrix, please wait ...")
      dat.na <- cfg$dat.def[fileIn] 
      # All default example data have genes in rows.
      if (!'data_shm.tar' %in% basename(dat.na)) {
        if (grepl('\\.rds$', dat.na)) {
          # SummarizedExperiment
          dat.ex <- readRDS(dat.na)
          if (is(dat.ex, 'SummarizedExperiment')) {
            dat.ex <- check_se(dat.ex)
            lgc.se <- is.character(dat.ex)
            if (lgc.se) showModal(modal(msg = dat.ex)); req(!lgc.se)
            se <- dat.ex$se; rdat <- rowData(se); rowData(se) <- link_dat(rdat, link.only=FALSE)
            # Assay metadata. 
            df.meta <- metadata(se)$df.meta
            if (!is.null(df.meta)) {
              if (ncol(data.frame(df.meta))<2) { 
                lgc.mt <- FALSE; msg <- 'The "df.meta" in the "metadata" slot should be a "data.frame" with at least two columns!';  
               if (!lgc.mt) showModal(modal(msg = msg)); req(!lgc.mt)
              } else {
                 metadata(se)$df.meta <- link_dat(df.meta, link.only=FALSE)
              }
            }
            dat.ex <- list(se.rep=se, se.aggr=NULL, con.na=dat.ex$con.na)
          }
        } else {
          dat.ex <- fread_df(input=dat.na, isRowGene=TRUE, rep.aggr=NULL) 
        }  
      } else { 
        # Exrtact data from uploaded tar.
        dat <- NULL; if (!is.null(cfg$pa.dat.upl)) if (file.exists(cfg$pa.dat.upl)) {
          cat('Extracting uploaded data... \n')
          # The prefix is not input$fileIn. The returned value of read_hdf5 is either data.frame or SE.
          dat <- read_hdf5(cfg$pa.dat.upl, dat.na)[[1]]
        }
        if (is.null(dat)|is.character(dat)|is.null(ipt$tar)) {
          cat('Extracting data from internal database ... \n')
          # "spFeature" and "variable" are are already checked when creating "data_shm.tar".
          dat <- readRDS(read_hdf5(dat.na, fileIn)[[1]]$data)
          lgc.se <- is.character(dat)
          if (lgc.se) showModal(modal(msg = dat)); req(!lgc.se)
        }
        dat <- check_se(dat); lgc.se <- is.character(dat)
        if (lgc.se) showModal(modal(msg = dat)); req(!lgc.se)
        se <- dat$se; rdat <- rowData(se); rowData(se) <- link_dat(rdat, link.only=FALSE)
        dat.ex <- list(se.rep=se, se.aggr=NULL, con.na=dat$con.na)
      }
      incProgress(0.3, detail="loading assay matrix ...")
      message('Done!'); return(dat.ex)
    }; dimNa <- dimName
    if (fileIn %in% cfg$na.cus & !grepl(na.sgl, fileIn)) {
      if (is.null(dimNa)) req('')
      if ((is.null(svgInpath1) & is.null(svgInpath2))) req('')
      if (is.null(geneInpath) | dimNa=="None") req('')
      incProgress(0.25, detail="importing data matrix ...")
      geneInpath <- geneInpath$datapath; targetInpath <- ipt$target$datapath
      metInpath <- ipt$met$datapath; asymetp <- ipt$asymet$datapath
      # Keep replicates unchaged, and compared with targets/metadata files.
      if (grepl('\\.rds$', geneInpath)) {
        dat <- readRDS(geneInpath)
        dat <- check_se(dat); lgc.se <- is.character(dat)
        if (lgc.se) showModal(modal(msg = dat)); req(!lgc.se)
        se <- dat$se; rdat <- rowData(se); rowData(se) <- link_dat(rdat, link.only=FALSE)
        dat <- list(se.rep=se, se.aggr=NULL, con.na=dat$con.na)
        return(dat)
      } else dat.cus <- fread_df(read_fr(geneInpath), isRowGene=(dimNa=='Row'), rep.aggr=NULL)
      se.rep <- dat.cus$se.rep; # df.met <- dat.cus$df.met
      if (!is.null(targetInpath)) {
        df.tar <- read_fr(targetInpath)
        # If errors detected, modalDialog terminates the app, while validate only stops the current step.
        lgc.tar <- nrow(df.tar) == ncol(se.rep) 
        if (!lgc.tar) showModal(modal(msg = 'Ensure "columns" in the assay matrix corresponds with "rows" in the targets file respectively!')); req(lgc.tar)
        # Check feature/variable columns in targets file.
        cna.tar <- colnames(df.tar)
        if (all(c('spFeature', 'variable') %in% cna.tar)) { 
          cna <- paste0(df.tar$spFeature, '__', df.tar$variable)
          ft <- df.tar$spFeature; vari <- df.tar$variable
        } else if ('spFeature' %in% cna.tar) { 
          cna <- paste0(df.tar$spFeature, '__', 'con')
          ft <- df.tar$spFeature; vari <- 'con'
        }
          se.rep$spFeature <- ft; se.rep$variable <- vari
          colnames(se.rep) <- cna
      } 
      if (!is.null(metInpath)) { 
        df.met <- read_fr(metInpath); lgc.met <- nrow(df.met) == nrow(se.rep)
        if (!lgc.met) showModal(modal(msg = 'Ensure "rows" in the assay matrix corresponds with "rows" in the row metadata file respectively!')); req(lgc.met)
        rownames(df.met) <- rownames(se.rep)
        rdat <- rowData(se.rep)
        rdat <- cbind(DataFrame(df.met[, !colnames(df.met) %in% colnames(rdat), drop=FALSE]), DataFrame(rdat))
        rowData(se.rep) <- rdat; dat.cus$se.rep <- se.rep
        if (!is.null(dat.cus$se.aggr)) rowData(dat.cus$se.aggr) <- rdat
      }
      if (!is.null(asymetp)) { 
        df.meta <- read_fr(asymetp)
        if (!is.null(df.meta)) {
          if (ncol(data.frame(df.meta))<2) { 
            lgc.mt <- FALSE; msg <- 'The "df.meta" in the "metadata" slot should be a "data.frame" with at least two columns!';  
            if (!lgc.mt) showModal(modal(msg = msg)); req(!lgc.mt)
          } else {
            df.meta <- link_dat(df.meta, link.only=FALSE)
            metadata(dat.cus$se.rep)$df.meta <- df.meta
            if (!is.null(dat.cus$se.aggr)) metadata(dat.cus$se.aggr)$df.meta <- df.meta
          }
        }
      }
      # Aggregate replicates after targets file is processed.
      # dat.cus <- fread_df(se.rep, isRowGene=(dimNa=='Row'), rep.aggr = NULL)
      incProgress(0.3, detail="loading assay matrix ...")
      cat('Done! \n'); return(dat.cus) 
    }
    })
  })
 
  nor.par <- reactiveValues()
  observeEvent(list(input$run, dat(), input$datIn, deg.mod.lis$input$eSHMBut, scell.mod.lis$covis.auto$but.covis, scell.mod.lis$covis.man$match.mod.lis$but.match$val, shm.mod$ipt$profile), {
  run <- input$run; dat <- dat(); datIn <- input$datIn
  fileIn <- ipt$fileIn; profile <- shm.mod$ipt$profile
  if (!check_obj(list(input$normDat, run, dat, datIn, !dat.no %in% fileIn))) return()
  pars <- list(input$normDat, dat, datIn)
  query.res <- check_exp(deg.mod.lis$query.res())
  if (!is(query.res, 'character') & !is.null(query.res) & !grepl(na.sgl, fileIn)) {
     pars <- list(normDat=input$normDat, dat=dat, query.res=query.res, datIn=datIn)
   } else if (grepl(na.sgl, fileIn)) {
     if (!check_obj(profile)) return()
     if (!'idp' %in% profile) { 
       sce.res <- scell.mod.lis$sce.res
       if (is.null(sce.res)) req('')
       se.aggr <- sce.res()$sce.lis$se.aggr
       pars <- list(normDat=input$normDat, se.aggr=se.aggr, datIn=datIn, profile=profile)
     } else {
       covis.type <- scell.mod.lis$sce.upl$covis.type 
       if (!check_obj(covis.type)) return() 
       if (covis.type %in% c('toBulkAuto', 'toCellAuto')) { 
         blk.sc <- scell.mod.lis$covis.auto$res
         # By default, data in co-clustering have no variables.
         colnames(blk.sc) <- paste0(colnames(blk.sc), '__con')
         blk.sc$variable <- 'con'
         con.na$v <- scell.mod.lis$covis.auto$con.na
       } else if (covis.type %in% c('toBulk', 'toCell')) {
         # Format features and variables in bulk.
         label <- scell.mod.lis$covis.man$covisGrp; bulk <- scell.mod.lis$covis.man$bulk
         if (!check_obj(list(label, bulk))) return()
         cdat.blk <- colData(bulk) 
         if (!'variable' %in% colnames(cdat.blk)) {
           bulk$variable <- 'con'; cdat.blk <- colData(bulk)
         } else con.na$v <- TRUE
         colnames(bulk) <- paste0(cdat.blk[, label], '__', cdat.blk[, 'variable'])
         vars.blk <- unique(colData(bulk)[, 'variable'])
         # Format features and variables in cell.
         cell <- scell.mod.lis$covis.man$dimred; var.cell <- NULL
         if ('variable' %in% colnames(colData(cell))) var.cell <- 'variable'
         sce.lis <- check_sce(sce=cell, label, var.cell)
         lgc.sc <- is(sce.lis, 'list')
         if (!lgc.sc) show_mod(lgc.sc, msg=sce.lis); req(lgc.sc)
         cell <- sce.lis$sce; var.cell <- sce.lis$var.cell
         con.na.cell$v <- sce.lis$con.na.cell
         vars.cell <- unique(colData(cell)[, var.cell]) 
         
         lgc.var <- identical(sort(vars.blk), sort(vars.cell))
         if (!lgc.var) show_mod(lgc.var, msg='Variables between bulk and single cell data should be the same!')
         req(lgc.var)
         # To combine bulk and cell. Bulk do not have reduced dimensions.
         reducedDim(cell, 'PCA') <- reducedDim(cell, 'UMAP') <- reducedDim(cell, 'TSNE') <- NULL
         reducedDim(bulk, 'PCA') <- reducedDim(bulk, 'UMAP') <- reducedDim(bulk, 'TSNE') <- NULL
         # By default, bulk and cell data have the same column names in colData, since they are stored in the same SCE. 
         # req(identical(sort(colnames(colData(bulk))), sort(colnames(colData(cell)))))
         blk.sc <- cbind_se(bulk, cell)
       }
       pars <- list(normDat=input$normDat, blk.sc=blk.sc, datIn=datIn, profile=profile) 
     }
   }; nor.par$pars <- pars
  })
  dat.nor <- eventReactive(nor.par$pars, {
    message('SHM: normalizing ... ')
    fileIn <- ipt$fileIn; dat <- dat(); normDat <- input$normDat
    datIn <- input$datIn; lis.par <- cfg$lis.par
    req(check_obj(list(fileIn, dat, normDat, datIn, lis.par, !dat.no %in% fileIn)))
    se <- dat$se.rep
    if ('enr' %in% datIn & !grepl(na.sgl, fileIn)) {
      se <- nor.par$pars$query.res
    } else if ('covis' %in% datIn & grepl(na.sgl, fileIn)) {
      prof <- nor.par$pars$profile
      if (!check_obj(prof)) return()
      if (! 'idp' %in% prof) se <- nor.par$pars$se.aggr else se <- nor.par$pars$blk.sc
      return(se)
    }; req(check_obj(list(se)))
    # Organize rowData.
    rdat <- rowData(se); idx.met <- grep(met.pat, colnames(rdat))
    rdat.sel <- rdat[, idx.met, drop=FALSE]
    cna.sel <- colnames(rdat.sel)
    idx.md <- grep('^metadata$', cna.sel)
    if (length(idx.md)>0) rdat.sel <- cbind(rdat.sel[, idx.md, drop=FALSE], rdat.sel[, -idx.md, drop=FALSE])
    rdat <- cbind(rdat.sel, rdat[, -idx.met, drop=FALSE])
    rowData(se) <- rdat 

    assay <- assay(se); lgc.as <- all(round(assay)==assay)
    # Must be before req(lgc.as).
    if ('None' %in% normDat) {
      if (!lgc.as) { message('Done!'); return(se) } else normDat <- lis.par$data.matrix['norm', 'default'] # Force to normalize count data. 
    }
    if (!lgc.as) {
      showNotification(HTML('Spatial Heatmap -> Data Table -> Settings: <br> normalization is skipped, since the input data are not count matrix.'), duration=2, closeButton = TRUE) 
      updateSelectInput(session, inputId='normDat', selected='None')
      se.aggr <- aggr_rep(data=se, assay.na=NULL, sam.factor='spFeature', con.factor='variable', aggr='mean')
      message('Done!'); return(se.aggr)
    }; req(lgc.as)
    withProgress(message="Normalizing: ", value = 0, {
      incProgress(0.5, detail="please wait ...")
      if (grepl('^CNF-', normDat)) {
        norm.fun <- 'CNF'
        par.lis <- list(method=sub('.*-', '', normDat))
      } else { norm.fun <- normDat; par.lis <- NULL }
        se.nor <- norm_data(data=se, norm.fun=norm.fun, par.list=par.lis, log2.trans=TRUE)
        se.aggr <- aggr_rep(data=se.nor, assay.na=NULL, sam.factor='spFeature', con.factor='variable', aggr='mean')
        incProgress(0.4, detail="...")
        message('Done!'); return(se.aggr)
    })
  })

  tran.par <- reactiveValues()
  observeEvent(list(input$run, dat.nor()), {
    run <- input$run; dat.nor <- dat.nor()
    if (!check_obj(list(input$log, run, dat.nor))) req('')
    pars <- list(input$log, dat.nor)
    tran.par$pars <- pars
  })
  dat.tran <- eventReactive(tran.par$pars, {
    message('SHM: log2/exp2 ... ')
    dat.nor <- dat.nor(); log <- input$log
    if (!check_obj(list(dat.nor, log))) req('')
    assay <- assay(dat.nor)
    lgc.as <- all(round(assay)==assay)
    if (lgc.as & 'exp2' %in% log) { 
      showNotification(HTML('Spatial Heatmap -> Data Table -> Settings: <br> exponent transformation is skipped, since the input assay data in this step are integers.'), duration=2, closeButton = TRUE) 
      message('Done!'); return(dat.nor)
    }
    if (!lgc.as & 'log2' %in% log) { 
      showNotification(HTML('Spatial Heatmap -> Data Table -> Settings: <br> exponent transformation is skipped, since the input assay data in this step are integers.'), duration=2, closeButton = TRUE) 
     message('Done!'); return(dat.nor)
    }
    withProgress(message="Log/exponent transformation: ", value = 0, {
      incProgress(0.2, detail="please wait ...")
      incProgress(0.4, detail="...")
      if ('log2' %in% log) {
        lgc.pos <- (min(assay) >= 0)
        if (!lgc.pos) { 
          showModal(modal(msg = 'Only non-negative intergers are accepted in log2 transformation!')); req('')
        }
        if (min(assay)==0) assay <- assay +1
        assay(dat.nor) <- log2(assay)
      } else if ('exp2' %in% log) { assay(dat.nor) <- 2^assay }
      incProgress(0.2, detail="...")
      message('Done!'); return(dat.nor)
    })
  })
 
  fil.par <- reactiveValues()
  observeEvent(list(input$run, dat.tran()), {
    A <- input$A; P <- input$P; CV1 <- input$CV1
    CV2 <- input$CV2; dat.tran <- dat.tran()
    if (!check_obj(list(A, P, CV1, CV2, dat.tran))) req('')
    pars <- list(A, P, CV1, CV2, dat.tran)
    fil.par$pars <- pars
  })
   se.fil <- eventReactive(list(fil.par$pars), {
     message('SHM: filtering ... ')
     A <- input$A; P <- input$P; CV1 <- input$CV1
     CV2 <- input$CV2; dat.tran <- dat.tran()
     if (!check_obj(list(A, P, CV1, CV2, dat.tran))) req('')
     p.lgc <- (P >= 0 & P <=1)
     show_mod(p.lgc, 'P should be between 0-1!'); req(p.lgc)
     cv.lgc <- (CV1 < CV2)
     show_mod(cv.lgc, 'CV1 should be less than CV2!'); req(cv.lgc)
     se.fil <- filter_data(data=dat.tran, sam.factor=NULL, con.factor=NULL, pOA=c(P, A), CV = c(CV1, CV2), verbose=FALSE)
     # se.lgc <- (nrow(se.fil) >= 5)
     # show_mod(se.lgc, 'Less than 5 rows remain!'); req(se.lgc)
     message('Done!'); se.fil
   })
  thr.par <- reactiveValues()
  observeEvent(list(input$run, se.fil()), {
    run <- input$run; se.fil <- se.fil()
    sig.max <- input$sig.max; sig.min <- input$sig.min
    if (!check_obj(list(sig.max, sig.min, run, se.fil))) req('')
    pars <- list(sig.min, sig.max, se.fil)
    thr.par$pars <- pars
  })
   se.thr <- eventReactive(list(thr.par$pars), {
     message('SHM: thresholding ... ')
     sig.max <- input$sig.max; sig.min <- input$sig.min
     se.fil <- se.fil()
     if (!check_obj(list(sig.max, sig.min))) req('')
     assay <- assay(se.fil)
     assay <- thrsd(thr.min=sig.min, thr.max=sig.max, data=assay)
     lgc.as <- !is(assay, 'character')
     if (!lgc.as) { msg <- assay; show_mod(lgc.as, msg) }
     req(lgc.as)
     assay(se.fil) <- round(assay, 2); message('Done!'); se.fil
   })
   observeEvent(list(ipt$fileIn, se.thr(), input$datIn), {
    fileIn <- ipt$fileIn; se.thr <- se.thr(); datIn <- input$datIn
    if (!check_obj(list(fileIn, se.thr, datIn))) return()
    if (grepl(na.sgl, fileIn) | !'all' %in% datIn) {
      shinyjs::hide(id='uplRefD') 
    } else shinyjs::show(id='uplRefD')
   })

   observeEvent(list(ipt$fileIn, se.thr(), input$datIn), {
    fileIn <- ipt$fileIn; se.thr <- se.thr(); datIn <- input$datIn
    if (!check_obj(list(fileIn, se.thr, datIn))) return()
    updateSelectInput(session, 'ref', selected='No')
    cna <- colnames(colData(se.thr))
    if (grepl(na.sgl, fileIn) | !'reference' %in% cna | !'all' %in% datIn) {
      shinyjs::hide(id='refD'); shinyjs::hide(id='ref')
    } else if ('reference' %in% cna & 'all' %in% datIn) {
      shinyjs::show(id='refD'); shinyjs::show(id='ref')
    }
   })
   uplref <- eventReactive(list(input$uplRef, input$refSel, ipt$fileIn, se.thr(), input$datIn), {
     message('Importing uploaded references ...')
     pa <- input$uplRef$datapath; yes <- input$refSel
     fileIn <- ipt$fileIn; se.thr <- se.thr()
     datIn <- input$datIn
     if (!check_obj(list(pa, yes, fileIn, se.thr, datIn))) return()
     if (grepl(na.sgl, fileIn) | !'all' %in% datIn) return()
     if (TRUE %in% yes) {
       ref <- tryCatch({ read_fr(pa) }, warning=function(w) { 'w' }, error=function(e) { 'e' }
       ); lgc.ref <- is(ref, 'character')
       if (lgc.ref) {
         msg <- 'The uploaded table cannot be imported!'
         show_mod(!lgc.ref, msg)  
       }; req(!lgc.ref)
       lgc.idt <- identical(colnames(se.thr), rownames(ref))
       if (!lgc.idt) {
         msg <- 'Ensure the rownames are correct!'
         show_mod(lgc.idt, msg)  
       }; req(lgc.idt); message('Done!'); ref
     } else return()
  })

   ref.par <- reactiveValues()
   observeEvent(list(input$ref, se.thr()), {
     pars <- list(input$ref, se.thr())
     if (!check_obj(pars)) return(); ref.par$pars <- pars
   })
   se.ref <- eventReactive(list(ref.par$pars, uplref()), {
     message('SHM: relative expressions ... ')
     ref <- input$ref; se.thr <- se.thr(); datIn <- input$datIn
     if (!check_obj(list(ref, se.thr, datIn))) return()
     uplref <- uplref()
     if (check_obj(uplref)) colData(se.thr)[, 'reference'] <- uplref[, 1] 
     if (!'Yes' %in% ref | !'reference' %in% colnames(colData(se.thr)) | !'all' %in% datIn) return(se.thr)
     se <- data_ref(se.thr)
     lgc.ref <- !is(se, 'character') & 'Yes' %in% ref
     if (!lgc.ref) { msg <- se; show_mod(lgc.ref, msg)
     return() }  
     message('Done!'); se
   })

  output$dldRef <- downloadHandler(# Download example references 
    filename=function(){ "mouse_organ_reference.txt" },
content=function(file=paste0(tmp.dir, '/mouse_organ_reference.txt')){
    ref <- read_fr('data/mouse_organ_reference.txt')
    write.table(ref, file, col.names=TRUE, row.names=TRUE, sep='\t')
  }
  ) 

  scl.par <- reactiveValues()
  observeEvent(list(input$run, input$ref, se.thr(), se.ref()), {
    scl <- input$scl; ref <- input$ref
    run <- input$run; se.thr <- se.thr(); se.ref <- se.ref()
    if (!check_obj(list(scl, ref, run, se.thr, se.ref))) req('')
    pars <- list(scl, ref, se.thr, se.ref)
    scl.par$pars <- pars
  })
   se.scl <- eventReactive(list(scl.par$pars), {
     message('SHM: scaling ... ')
     scl <- input$scl; ref <- input$ref
     run <- input$run; se.thr <- se.thr(); se.ref <- se.ref()
     if (!check_obj(list(scl, ref, run, se.thr, se.ref))) req('') 
     # if ('Yes' %in% ref) se <- se.ref else se <- se.thr
     se <- se.ref; assay <- assay(se)
     # Scale by row/column
     if (scl=='Row') { assay <- t(scale(t(assay))) 
     } else if (scl=='All') { assay <- scale_all(assay) }
     assay(se) <- assay; cna <- colnames(se)
     # Co-clustering unlabeled cells.
     idx <- grepl('^none$|^none__', colnames(se))
     se <- cbind(se[, !idx], se[, idx])
     message('Done!'); se
   })
  scl.sel.par <- reactiveValues()
  observeEvent(list(input$run, input$ref, se.scl(), ids$sel), {
    scl <- input$scl; run <- input$run; ref <- input$ref
    se.scl <- se.scl()
    if (!check_obj(list(scl, run, ref, se.scl, ids$sel))) return()
    pars <- list(scl, se.scl, ref, ids$sel)
    scl.sel.par$pars <- pars
  })
   se.scl.sel <- eventReactive(list(scl.sel.par$pars), {
     message('SHM: selected data ... ')
     scl <- input$scl; run <- input$run; se.scl <- se.scl()
     ref <- input$ref
     if (!check_obj(list(scl, run, se.scl, ref, ids$sel))) return()
     if (!all(ids$sel %in% rownames(se.scl))) return()
     se.scl.sel <- se.scl[ids$sel, ]
     assay.sel <- assay(se.scl.sel)
     if (scl=='Selected' & !'Yes' %in% ref) { assay.sel <- scale_all(assay.sel) }
     assay(se.scl.sel) <- assay.sel
     message('Done!'); se.scl.sel
   })

  sear <- reactiveValues(id=NULL)
  observeEvent(ipt$fileIn, { sear$id <- NULL })
  observeEvent(sch$but, {
    se.scl <- se.scl(); lis.par <- cfg$lis.par
    if (!check_obj(list(se.scl, lis.par))) req('')
    if (sch$sch=='') sel <- as.numeric(lis.par$data.matrix['row.selected', 'default']) else {
      gens <- strsplit(gsub(' |,', '_', sch$sch), '_')[[1]]
      pat <- paste0('^', gens, '$', collapse='|')
      sel <- which(grepl(pat, x=rownames(se.scl), ignore.case=TRUE, perl=TRUE))
      if (length(sel)==0) sel <- as.numeric(lis.par$data.matrix['row.selected', 'default'])
     }; sear$id <- sel
  })
  dt.shm <- eventReactive(list(se.scl(), input$spk), {
    cat('Prepaing data table ... \n')
    se.scl <- se.scl(); dat <- dat(); spk <- input$spk
    if (!check_obj(list(se.scl, dat, spk))) req('')
    assay <- assay(se.scl); rdat <- rowData(se.scl)
    withProgress(message="Data table: ", value = 0, {
      incProgress(0.5, detail="please wait ...")
      rdat <- rdat[, grep(met.pat, colnames(rdat)), drop=FALSE]
      if (!all(assay==round(assay))) assay <- round(assay, 2)
      if ('Yes' %in% spk) {
      spk.code <- NULL
      for (i in seq_len(nrow(assay))) {
        spk.code <- c(spk.code, spk_chr(setNames(unlist(assay[i, ]), NULL), lineColor = 'black', fillColor = '#ccc', chartRangeMin = 0, chartRangeMax = 8, width = 80, height = 30, highlightLineColor = 'orange', highlightSpotColor = 'orange'))
      }; df.spk <- cbind.data.frame(Sparklines=spk.code, assay) 
      }
      if (ncol(rdat) > 0) {
        if ('Yes' %in% spk) df.tab <- cbind.data.frame(rdat, df.spk, stringsAsFactors=FALSE) else df.tab <- cbind.data.frame(rdat, assay, stringsAsFactors=FALSE)
      } else {   
        if ('Yes' %in% spk) df.tab <- df.spk else df.tab <- assay
      }
      # Remove '__con' only in the data table, not in the downstream (shm, network).
      if (dat$con.na==FALSE) colnames(df.tab) <- sub('__con$', '', colnames(df.tab))
      incProgress(0.4, detail="...")
      cat('Done!\n'); df.tab
    })
  }); observe({ dt.shm() })
  output$selProf <- renderPlot({
    cat('Profile of selected rows ... \n') 
    scl <- input$scl; se.scl.sel <- se.scl.sel(); dat <- dat()
    if (!check_obj(list(se.scl.sel, dat, scl))) req('')
    # validate(need(length(ids$sel)<=100, 'Due to space limitation, profiles of 100+ genes are not plotted!')) 
    dt.sel <- assay(se.scl.sel)
    if (dat$con.na==FALSE) colnames(dt.sel) <- sub('__con$', '', colnames(dt.sel)) 
    if (scl=='No') title <- 'No scaling' else if (scl=='Row') title <- 'Scaled by row' else if (scl=='Selected') title <- 'Scaled across selected rows' else if (scl=='All') title <- 'Scaled across all rows' else title <- '' 
    lgd.guide <- guides(color=guide_legend(nrow=2, byrow=FALSE, title=NULL))
    # grid.arrange(g1, g2, nrow=2)
    if ('No' %in% scl) { 
      y.title=paste0('Un-scaled values (', round(min(dt.sel), 2), '-', round(max(dt.sel), 2), ')')
    } else {
      y.title=paste0(title, ' (', round(min(dt.sel), 2), '-', round(max(dt.sel), 2), ')')
    }
    g <- graph_line(dt.sel, y.title=y.title, text.size=12, lgd.guide=lgd.guide); message('Done!'); g
  })

  observeEvent(list(input$run), ignoreInit=FALSE, {
    updateNavbarPage(session, inputId='settNav', selected = 'dat')
  })

  observeEvent(dt.shm(), {
    gene.dt <- dt.shm(); if (!check_obj(list(gene.dt))) return() 
    updateNumericInput(session, inputId="r2", value=nrow(gene.dt)) 
    updateNumericInput(session, inputId="c2", value=ncol(gene.dt)) 
  })
  subdat <- reactiveValues(r1=1, r2=500, c1=1, c2=20)
  observeEvent(list(input$run+1, dt.shm()), ignoreInit=FALSE, {
    r1 <- input$r1; r2 <- input$r2
    c1 <- input$c1; c2 <- input$c2
    gene.dt <- dt.shm()
    if (!check_obj(list(r1, r2, c1, c2, gene.dt))) return() 
    if (0 %in% input$run) {
      subdat$r2 <- nrow(gene.dt)
      if (ncol(gene.dt) >= 30) subdat$c2 <- 30 else subdat$c2 <- ncol(gene.dt)
      return()
    }
    if (r1 < 1) r1 <- 1; if (c1 < 1) c1 <- 1
    lgc.r <- r2 > r1; if (!lgc.r) {
      show_mod(lgc.r, msg='Row End should > Row Start!')
    }; req(lgc.r)
    lgc.c <- c2 > c1; if (!lgc.c) {
      show_mod(lgc.c, msg='Column End should > Column Start!')
    }; req(lgc.c)
    if (nrow(gene.dt) < r2) r2 <- nrow(gene.dt)
    if (ncol(gene.dt) < c2) c2 <- ncol(gene.dt)
    subdat$r1 <- r1; subdat$r2 <- r2
    subdat$c1 <- c1; subdat$c2 <- c2
  })
  output$dtSel <- renderDataTable({
    cat('Preparing selected data matrix ... \n')
    gene.dt <- dt.shm()
    if (!check_obj(list(gene.dt, ids$sel))) req('')
    withProgress(message="Data table (selected): ", value = 0, {
      incProgress(0.5, detail="please wait ...") 
      # Tooltip on metadata.
      col1 <- list(list(targets = c(1), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
      # In case no metadata column.
      if (colnames(gene.dt)[1]!='metadata') col1 <- NULL
      dtab <- datatable(gene.dt[ids$sel, seq(subdat$c1, subdat$c2, 1), drop=FALSE], selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=FALSE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=col1, dom='t', fixedColumns = list(leftColumns=3), class='cell-border strip hover', 
   fnDrawCallback = htmlwidgets::JS('function(){HTMLWidgets.staticRender()}')
   )) %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% spk_add_deps()
   # formatRound(colnames(assay.sel), deci)
   incProgress(0.4, detail="please wait ...")
    cat('Done! \n'); dtab
  })
  })
  dt.sel <- reactiveValues(val='none')
  output$dtAll <- renderDataTable({
    cat('Preparing complete data matrix ... \n')
    gene.dt <- dt.shm(); page.h <- input$page
    if (is.null(gene.dt)|!is.numeric(page.h)) return()
    # Decimals.
    # Tooltip on metadata.
    col1 <- NULL; cna <- colnames(gene.dt)
    idx.met <- which(cna %in% 'metadata')
    if (length(idx.met)>0) col1 <- idx.met
    col1 <- list(list(targets = col1, render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    # if (colnames(gene.dt)[1]!='metadata') col1 <- NULL
    # dat <- gene.dt[seq(subdat$r1, subdat$r2, 1), seq(subdat$c1, subdat$c2, 1), drop=FALSE]
    colnames(gene.dt) <- sub('__', '_', colnames(gene.dt))
    dtab <- datatable(gene.dt[seq(subdat$r1, subdat$r2, 1), seq(subdat$c1, subdat$c2, 1), drop=FALSE], selection=list(mode="multiple", target="row", selected=dt.sel$val), escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=page.h, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, class='cell-border strip hover', columnDefs=col1, fixedColumns = list(leftColumns=3), 
   fnDrawCallback = htmlwidgets::JS('function(){HTMLWidgets.staticRender()}')
   )) %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% spk_add_deps()
   cat('Done! \n'); dtab
   # formatRound(idx.num, deci); cat('Done! \n'); dtab
  })

  output$expDsg <- renderDataTable({
    cat('Preparing Experiment design ... \n')
    se.scl <- se.scl(); dat <- dat()
    if (!check_obj(list(se.scl, dat))) req('')
    cdat <- colData(se.scl)
    withProgress(message="Experiment design: ", value = 0, {
      incProgress(0.5, detail="please wait ...") 
      # Tooltip on metadata.
      col1 <- list(list(targets = seq_len(ncol(cdat)), render = DT::JS("$.fn.dataTable.render.ellipsis(50, false)")))
      rownames(cdat) <- sub('__', '_', rownames(cdat))
      if ('reference' %in% colnames(cdat)) cdat$reference <- sub('__', '_', cdat$reference)
      dtab <- datatable(data.frame(cdat), selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=FALSE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, class='cell-border strip hover', columnDefs=col1, fixedColumns =NULL, 
   fnDrawCallback = htmlwidgets::JS('function(){HTMLWidgets.staticRender()}')
   )) %>% formatStyle(0, backgroundColor="orange", cursor='pointer')
   # formatRound(colnames(assay.sel), deci)
   incProgress(0.4, detail="please wait ...")
    cat('Done! \n'); dtab
  })
  })
  output$over <- renderDT({
    cat('Preparing assay metadata ... \n')
    se.scl <- se.scl(); dat <- dat()
    if (!check_obj(list(se.scl, dat))) req('')
    meta <- metadata(se.scl)$df.meta
    if (is.null(meta)) return() 
    withProgress(message="Assay/image overview: ", value = 0, {
      incProgress(0.5, detail="please wait ...")
      dat <- datatable(data.frame(meta), escape=FALSE, selection='none',
      options = list(
      deferRender = TRUE, scrollY = TRUE, scrollX = TRUE, scroller = TRUE, autoWidth = TRUE, columnDefs = list(list(width = '50%', targets = grep('^description$', colnames(meta)))),
      fnDrawCallback = htmlwidgets::JS('function(){HTMLWidgets.staticRender()}')
      )) %>% formatStyle(0, backgroundColor="orange", cursor='pointer') # %>% formatStyle(columns = c(1), width='10%')
      incProgress(0.4, detail="please wait ...")
      cat('Done! \n'); dat
  })
  })
  # ids from URL are only used once. 
  url.id.sel <- reactiveValues(init=TRUE)
  observeEvent(list(lis.url$par$ids, dt.shm()), {
    gene.dt <- dt.shm(); # selRow <- input$selRow
    id.url <- lis.url$par$ids
    req(check_obj(list(gene.dt, id.url, url.id.sel$init)))
    rna <- rownames(gene.dt)
    id.url.no <- setdiff(id.url, rna); lgc.no <- length(id.url.no)==0
    if (!lgc.no) showModal(modal(msg = paste0('Invalid IDs in the URL: ', paste0(id.url.no, collapse=',')))); req(lgc.no)
    dt.sel$val <- which(rna %in% id.url); # ids$sel <- id.url
    url.id.sel$init <- FALSE
  })
  observeEvent(list(input$selRow), { # Select genes in table.
    gene.dt <- dt.shm(); if (is.null(gene.dt)) return()
    ids$sel <- rownames(gene.dt)[input$dtAll_rows_selected]
    dt.sel$val <- input$dtAll_rows_selected
    lgc.ids <- check_obj(ids$sel)
    if (!lgc.ids) showModal(modal(msg = 'No genes are selected!')); req(lgc.ids)
    tabTop <- parent$input$tabTop
    if ('shmPanelAll' %in% tabTop & check_obj(ids$sel)) {
      runjs('document.getElementById("tabTop").scrollIntoView()') 
    }
  })
  observeEvent(list(input$deSel), { # Select genes in table.
    gene.dt <- dt.shm(); if (is.null(gene.dt)) return()
    ids$sel <- NULL; dt.sel$val <- 'none'
  })
  observeEvent(list(ipt$fileIn, deg.mod.lis$input$eSHMBut, scell.mod.lis$covis.man$match.mod.lis$but.match$val, scell.mod.lis$covis.auto$but.covis, input$datIn), { # Select genes in table.
    id.url <- lis.url$par$ids; but.sgl <- ids$but.sgl
    but.mul <- ids$but.mul; selRow <- input$selRow
    fileIn <- gsub('"', '', lis.url$par$'upl-fileIn')
    # If bookmarked url, ids$sel is not set NULL by "covis.auto$but.covis".  
    if (check_obj(id.url) & (is.null(but.sgl)|0 %in% but.sgl) & (is.null(but.mul) | 0 %in% but.mul) & (is.null(selRow) |0 %in% selRow) & ipt$fileIn %in% fileIn) return()
    ids$sel <- NULL; dt.sel$val <- 'none'
    # print(list('clear', id.url, but.sgl, but.mul, selRow, ids$sel, dt.sel$val))
  })
  observeEvent(list(parent$input$btnInf), {
    btnInf <- parent$input$btnInf
    if (!check_obj(btnInf)) return()
    if (btnInf > 0) updateTabsetPanel(session, "settNav", selected='over')
  })

  observe({
    ipt$fileIn; ipt$geneInpath; lis.par <- cfg$lis.par; lis.url
    req(check_obj(list(lis.par, lis.url)))
    updateSelectInput(session, 'log', selected=url_val('dat-log', lis.url, def=lis.par$data.matrix['log.exp', 'default']))
    updateSelectInput(session, 'scl', selected=url_val('dat-scale', lis.url, def=lis.par$data.matrix['scale', 'default']))
  })

  observe({
    ipt$fileIn; ipt$geneInpath; input$log; lis.url
    lis.par <- cfg$lis.par; req(check_obj(list(lis.par, lis.url)))
    updateSelectInput(session, 'normDat', selected=url_val('dat-normDat', lis.url, def=lis.par$data.matrix['norm', 'default']))
    updateNumericInput(session, "A", value=url_val('dat-A', lis.url, def=as.numeric(lis.par$data.matrix['A', 'default'])))  
    updateNumericInput(session, inputId="P", value=url_val('dat-P', lis.url, def=as.numeric(lis.par$data.matrix['P', 'default'])))
    updateNumericInput(session, inputId="CV1", value=url_val('dat-CV1', lis.url, def=as.numeric(lis.par$data.matrix['CV1', 'default'])))
    updateNumericInput(session, inputId="CV2", value=url_val('dat-CV2', lis.url, def=as.numeric(lis.par$data.matrix['CV2', 'default']))) 
  })

  observeEvent(scell.mod.lis$covis.man$match.mod.lis$but.match$val, ignoreInit=TRUE, {
    updateTabsetPanel(session, "dtab.shm", selected='dTabScell')
  })
  ipt.dat <- reactiveValues()
  observe({ ipt.dat$dt_rows_selected <- input$dt_rows_selected })
  col.cfm <- reactive({ input$col.cfm })
  scl <- reactive({ input$scl })
  log <- reactive({ input$log }); A <- reactive({ input$A })
  CV1 <- reactive({ input$CV1 }); CV2 <- reactive({ input$CV2 })
  P <- reactive({ input$P })
  search.but <- reactive({ sch$but })
  sig.but <- reactive({ input$sig.but })
  observe({
    dat <- dat(); profile <- shm.mod$ipt$profile; fileIn <- ipt$fileIn
    if (!check_obj(list(dat, profile, fileIn, !dat.no %in% fileIn))) return()
    if (! 'idp' %in% profile | !grepl(na.sgl, fileIn)) con.na$v <- dat$con.na     
  })
  onBookmark(function(state) { state })
  return(list(sear = sear, ipt.dat = ipt.dat, col.cfm = col.cfm, scaleDat=scl, log = log, A = A, P = P, CV1 = CV1, CV2 = CV2, search.but = search.but, sig.but=sig.but, dat=dat, se.scl=se.scl, se.scl.sel=se.scl.sel, con.na=con.na, con.na.cell=con.na.cell, se.thr=se.thr, sn=session))
  })

}



