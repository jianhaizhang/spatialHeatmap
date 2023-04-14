# Module for spatial enrichment.
deg_server <- function(id, sch, lis.url, url.id, ids, upl.mod.lis, dat.mod.lis, shm.mod.lis, session) {
  moduleServer(id, function(input, output, session) {
  ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg
  gID <- shm.mod.lis$gID; datIn <- dat.mod.lis$dat

  cnt.q.res <- reactiveValues(v=0)
  observeEvent(input$degAll, {
    if (!'dt.deg' %in% input$degAll|cnt.q.res$v > 2) return()
    show_mod(FALSE, 'Please click "Enrichment SHMs".')
    cnt.q.res$v <- cnt.q.res$v + 1
  })
  # dat.deg.mod.lis <- data_server('datDEG', sch, lis.url, ids, deg=TRUE, upl.mod.lis)
  dat <- reactive({
    cat('DEG: SummarizedExperiment, features, variables ... \n')
    datIn <- datIn(); if (is.null(datIn)|grepl(na.sgl, ipt$fileIn)) return()
    se.rep <- datIn$se.rep; df.rep <- as.matrix(assay(se.rep))
    int <- all(round(df.rep) == df.rep)
    show_mod(int, 'Only count matrix is accepted!', title='Spatial Enrichment'); req(int)
    rows <- nrow(se.rep) >= 50
    show_mod(rows, 'Make sure count matrix includes at least 10 biomolecules!'); req(rows)
    cna <- colnames(se.rep)
    sams <- gsub('^(.*)__(.*)$', '\\1', cna)
    cons <- gsub('^(.*)__(.*)$', '\\2', cna)
    sf.var.lgc <- length(unique(sams))>1 & length(unique(cons))>1
    show_mod(sf.var.lgc, 'At least 2 spatial features and 2 variables are required!!'); req(sf.var.lgc)
    cat('Done! \n')
    return(list(se=se.rep, sams=sams, cons=cons))
  })
  fil.par <- reactiveValues() 
  observeEvent(list(input$run, dat()), {  
    A <- input$A; P <- input$P; CV1 <- input$CV1
    CV2 <- input$CV2; dat <- dat(); run <- input$run
    if (!check_obj(list(A, P, CV1, CV2, dat, run))) req('')
    # if (run==0) req('') 
    pars <- list(A, P, CV1, CV2, dat)  
    fil.par$pars <- pars 
  })   
  se.fil <- eventReactive(list(fil.par$pars), {
    message('DEG: filtering ... ')
    A <- input$A; P <- input$P; CV1 <- input$CV1
    CV2 <- input$CV2; dat <- dat()
    if (!check_obj(list(A, P, CV1, CV2, dat))) req('')
    se <- dat$se
    p.lgc <- (P >= 0 & P <=1)
    show_mod(p.lgc, 'P should be between 0-1!'); req(p.lgc)
    cv.lgc <- (CV1 < CV2)
    show_mod(cv.lgc, 'CV1 should be less than CV2!'); req(cv.lgc)
    se.fil <- filter_data(data=se, sam.factor=NULL, con.factor=NULL, pOA=c(P, A), CV = c(CV1, CV2), verbose=FALSE)
    lgc.se <- (nrow(se.fil)>=50)
    show_mod(lgc.se, 'Less than 50 rows remain!'); req(lgc.se)
    message('Done!')
    return(list(se=se.fil, sams=dat$sams, cons=dat$cons))
   })
  # Normalizing once on the complete data: avoid repetitive normalizations when features/variables change. 
  nor.par <- reactiveValues() 
  observeEvent(list(input$run, se.fil()), {
    norMeth <- input$norMeth; se.fil <- se.fil()
    if (!check_obj(list(norMeth, se.fil))) req('') 
    pars <- list(norMeth, se.fil)  
    nor.par$pars <- pars 
  })   
  se.nor <- eventReactive(nor.par$pars, {
    cat('DEG: normalizing data ... \n')
    se.fil <- se.fil(); norMeth <- input$norMeth
    if (!check_obj(list(se.fil, norMeth))) req('')
    if ('none' %in% norMeth) return(se.fil$se)
    se.nor <- norm_data(data=se.fil$se, norm.fun='CNF', par.list=list(method=norMeth), log2.trans=TRUE)
    cat('Done! \n'); se.nor
  })

  output$ssg.sam <- renderUI({
    ns <- session$ns
    dat <- dat(); if (is.null(dat)) return()
    cho <- c('all', unique(dat$sams))
    selectInput(ns('ssg.sam'), label='Select spatial features', choices=cho, selected=cho[2:3], multiple=TRUE)
  })
  output$ssg.con <- renderUI({
    ns <- session$ns
    dat <- dat(); if (is.null(dat)) return()
    cho <- c('all', unique(dat$cons))
    selectInput(ns('ssg.con'), label='Select variables', choices=cho, selected=cho[2:3], multiple=TRUE)
  })

  se.sub <- reactive({
    cat('Subsetting SE with input features/variables ... \n')
    se.fil <- se.fil(); comBy <- input$comBy
    sam <- input$ssg.sam; con <- input$ssg.con
    if (!check_obj(list(se.fil, comBy, sam, con))) return()
    if ('all' %in% sam) sam <- unique(se.fil$sams)
    if ('all' %in% con) con <- unique(se.fil$cons)
    if (comBy=='feature') fct <- 'ft' else if (comBy=='variable') fct <- 'var' else if (comBy=='feature__variable') fct <- 'ft.var'
    se <- se.fil$se
# save(se, sam, con, fct, file='sscff')
    se.sub <- sf_var(data=se.fil$se, feature='spFeature', ft.sel=sam, variable='variable', var.sel=con, com.by=fct)
    if (is(se.sub, 'character')) return()
    # Replicates >= 2.
    fct.tab <- table(colData(se.sub)$com.by)
    fct.na <- names(fct.tab)[fct.tab==1]
    rep.lgc <- length(fct.na)==0
    msg <- paste0('At least 2 replicates are required: ', paste0(fct.na, collapse=', '))
    if (!rep.lgc) showModal(modal(msg = msg))
    validate(need(rep.lgc, ''))
    cat('Done! \n'); se.sub
  })
  output$query <- renderUI({
    ns <- session$ns; se.sub <- se.sub()
    if (is.null(se.sub)) return()
    selectInput(ns('query'), label='Select a query', choices=sort(unique(se.sub$com.by)), selected=NULL)
  }) 

  # Pairwise comparison coefficients.
  output$dt.vs1 <- output$dt.vs2 <- renderDataTable({
    cat('Pairwise comparison coefficients ... \n')
    dat <- dat(); sam <- input$ssg.sam; con <- input$ssg.con
    comBy <- input$comBy; tar <- input$query
    if (is.null(dat)|!is.character(sam)|!is.character(con)|!is.character(comBy)|!is.character(tar)) return()
    if ('all' %in% sam) sam <- unique(dat$sams)
    if ('all' %in% con) con <- unique(dat$cons)
    vs <- data.frame()
    # save(dat, sam, con, comBy, tar, file='dscstf')
    # Reference, query.
    if (comBy %in% c('feature', 'variable')) {
      if (comBy == 'feature') {  
        under <- con; ref <- setdiff(sam, tar)
        num.sf <- length(ref) > 0
        msg <- 'If compare by "spatial feature", at least 2 features should be selected!'
        if (!num.sf) showModal(modal(msg = msg))
        validate(need(num.sf, ''))
      } else if (comBy == 'variable') {
        under <- sam; ref <- setdiff(con, tar)
        num.var <- length(ref) > 0
        msg <- 'If compare by "variable", at least 2 variables should be selected!'
        if (!num.var) showModal(modal(msg = msg))
        validate(need(num.var, ''))
      }
      tar.all <- paste0(tar, '__', under)
      ref.all <- paste0(unlist(lapply(ref, function(i) paste0(i, '__', under))))
      vs <- rbind(vs, c(tar.all, 'VS', ref.all))
      colnames(vs) <- c(paste0('query', seq_along(tar.all)), 'VS', paste0('reference', seq_along(ref.all))) 
    } else if (comBy == 'feature__variable') {
      coms <- unlist(lapply(sam, function(i) {paste0(i, '__', con)} ))
      ref <- setdiff(coms, tar)
      vs <- rbind(vs, c(tar, 'VS', ref))
      colnames(vs) <- c('query', 'VS', paste0('reference', seq_along(ref))) 
    }
    d.tab <- datatable(vs, selection='none', extensions='Scroller', plugins = "ellipsis", class='cell-border strip hover', options = list(dom = 't', scrollX = TRUE))
    cat('Done! \n'); d.tab
  })
  edg0 <- reactive({
    cat('edgeR all ... \n'); se.sub <- se.sub()
    norMeth <- input$norMeth
    if (!check_obj(list(se.sub, norMeth))) return()
    if ('none' %in% norMeth) req('')
    withProgress(message="edgeR: ", value=0, {
      incProgress(0.5, detail="in progress ...")
      edg <- edgeR(se=se.sub, method.norm=norMeth, com.factor='com.by', method.adjust='BH', return.all=TRUE)
      incProgress(0.4, detail="in progress ...")
      cat('Done! \n'); edg
    })
  })
  edg <- eventReactive(input$run, {
    cat('edgeR log2/fc ... \n')
    if (!'edgeR' %in% input$meth) req('')
    se.sub <- se.sub(); fc <- input$ssg.fc; fdr <- input$ssg.fdr
    edg0 <- edg0(); outlier <- input$outlier
    if (!check_obj(list(se.sub, fc, fdr, edg0, outlier))) return()
    sam.sub <- sort(unique(se.sub$com.by))
    up.dn <- up_dn(sam.all=sam.sub, df.all=edg0, log.fc=abs(fc), fdr=fdr, log.na='logFC', fdr.na='FDR', method='edgeR', outliers=outlier)
    message('Done!'); up.dn
  }); observe({ edg() })

  dsq0 <- reactive({
    cat('DESeq2 all ... \n'); se.sub <- se.sub()
    norMeth <- input$norMeth
    if (!check_obj(list(se.sub, norMeth))) return() 
    withProgress(message="DESeq2: ", value=0, {
    incProgress(0.5, detail="in progress ...")
    dsq <- deseq2(se=se.sub, com.factor='com.by', method.adjust='BH', return.all=TRUE); cat('Done! \n')
    incProgress(0.4, detail="in progress ...")
    cat('Done! \n'); dsq
    })
  }) 
  dsq <- eventReactive(input$run, {
    cat('DESeq2 log2/fc ... \n')
    if (!'DESeq2' %in% input$meth) return()
    se.sub <- se.sub(); fc <- input$ssg.fc; fdr <- input$ssg.fdr
    dsq0 <- dsq0(); outlier <- input$outlier
    if (!check_obj(list(se.sub, fc, fdr, dsq0, outlier))) return()
    sam.sub <- sort(unique(se.sub$com.by))
    up.dn <- up_dn(sam.all=sam.sub, df.all=dsq0, log.fc=abs(fc), fdr=fdr, log.na='log2FoldChange', fdr.na='padj', method='DESeq2', outliers=outlier)
    cat('Done! \n'); up.dn
  }); observe({ dsq() })

  lim0 <- reactive({
    cat('limma all ... \n'); se.sub <- se.sub()
    norMeth <- input$norMeth
    if (!check_obj(list(se.sub, norMeth))) return() 
    if ('none' %in% norMeth) req('')
    withProgress(message="limma: ", value=0, {
    incProgress(0.5, detail="in progress ...") 
    lim <- limma(se=se.sub, method.norm=norMeth, m.array=FALSE, com.factor='com.by', method.adjust='BH', return.all=TRUE)
    incProgress(0.4, detail="in progress ...")
    cat('Done! \n'); lim
    })
  })
  lim <- eventReactive(input$run, {
    cat('limma log2/fc ... \n')
    if (!'limma' %in% input$meth) return()
    se.sub <- se.sub(); fc <- input$ssg.fc; fdr <- input$ssg.fdr
    lim0 <- lim0(); outlier <- input$outlier
    if (!check_obj(list(se.sub, fc, fdr, lim0, outlier))) return()
    sam.sub <- sort(unique(se.sub$com.by))
    up.dn <- up_dn(sam.all=sam.sub, df.all=lim0, log.fc=abs(fc), fdr=fdr, log.na='logFC', fdr.na='adj.P.Val', method='limma', outliers=outlier)
    cat('Done! \n'); up.dn
  }); observe({ lim() })
  dis0 <- reactive({
    cat('distinct all ... \n'); se.sub <- se.sub()
    norMeth <- input$norMeth
    if (!check_obj(list(se.sub, norMeth))) return()
    if ('none' %in% norMeth) req('')
    withProgress(message="distinct: ", value=0, {
      incProgress(0.5, detail="this method takes longer time ...")
      dis <- distt(se.sub, norm.fun='CNF', par.list=list(method=norMeth), com.factor='com.by', return.all=TRUE)
      incProgress(0.4, detail="this method takes longer time ...")
      cat('Done! \n'); dis
    })
  })
  dis <- eventReactive(input$run, {
    cat('distinct log2/fc ... \n')
    if (!'distinct' %in% input$meth) return()
    se.sub <- se.sub(); fc <- input$ssg.fc; fdr <- input$ssg.fdr
    dis0 <- dis0(); outlier <- input$outlier
    if (!check_obj(list(se.sub, fc, fdr, dis0, outlier))) return()
    sam.sub <- sort(unique(se.sub$com.by))
    up.dn <- up_dn(sam.all=sam.sub, df.all=dis0, log.fc=abs(fc), fdr=fdr, log.na='log2FC', fdr.na='FDR', method='distinct', outliers=outlier)
    cat('Done! \n'); up.dn
  }); observe({ dis() })
  # All enrichment results: up/down results, normalized data.
  res <- eventReactive(list(edg(), dsq(), lim(), dis(), se.nor()), {
    se.nor <- se.nor(); se.sub <- se.sub(); meth <- input$meth
    if (!check_obj(list(meth))) return()
    if ('edgeR' %in% meth) up.dn <- edg()
    if ('DESeq2' %in% meth) up.dn <- dsq()
    if ('limma' %in% meth) up.dn <- lim()
    if ('distinct' %in% meth) up.dn <- dis()
    if (!check_obj(list(se.nor, up.dn, se.sub))) return()
    up.dn.all <- sum(unlist(lapply(lapply(up.dn, function(i) do.call('rbind', i)), function(i) nrow(i))))
    if (up.dn.all==0) {
      msg <- 'No enriched/depleted biomolecules detected!'
      showModal(modal(msg = msg)); return()
    }
    se.nor.sub <- se.nor[, colnames(se.sub)]
    colData(se.nor.sub) <- colData(se.sub)
    se.nor.sub <- aggr_rep(se.nor.sub, sam.factor=NULL, con.factor=NULL, aggr='mean')
    list(result=up.dn, data=se.nor.sub)
  })
  output$upset <- renderPlot({
    meth <- input$meth; enrType <- input$enrType; res <- res()
    if (!check_obj(list(res, meth, enrType, res))) return()
    ovl_enrich(res, type=enrType, plot='upset')
  })
  output$matrix <- renderPlot({
    meth <- input$meth; enrType <- input$enrType; res <- res()
    if (!check_obj(list(res, meth, enrType, res))) return()
    ovl_enrich(res, type=enrType, plot='matrix')
  })
  output$venn <- renderPlot({
    meth <- input$meth; enrType <- input$enrType; res <- res()
    if (!check_obj(list(res, meth, enrType, res))) return()
    sams <- unique(res$data$com.by)
    sams.lgc <- length(sams) <= 5
    if (!sams.lgc) {
      msg <- 'Venn diagrams are used for max 5-way overlaps.'
      showModal(modal(msg = msg)); return()
    }
    ovl_enrich(res, type=enrType, plot='venn')
  })
  
  query.res <- eventReactive(list(input$query, res()), {
    query <- input$query; meth <- input$meth; res <- res()
    if (!check_obj(list(res, query, meth, res))) return()
    q.res <- query_enrich(res=res, query=query)
    if (is.character(q.res)) {
      showModal(modal(msg = msg)) 
      validate(need(try(!is.character(q.res))), '')
    }; q.res
  }) 

  observeEvent(res(), ignoreInit=FALSE, ignoreNULL=TRUE, {
    if (is.null(res())) return()
    updateNavbarPage (session, "degAll", selected='ovl')
  })
  dt.deg <- reactive({
    cat('DEG data frame ... \n')
    dat <- dat(); query.res <- query.res()
    se.sub <- se.sub(); se.nor <- se.nor()
    if (!check_obj(list(dat, query.res, se.sub, se.nor))) req('')
    df.met <- rowData(dat$se)
    df.deg <- rowData(query.res)[, c('type', 'total', 'method')]
    # The data (geneIn1) before filtering is used, so even though all genes are filtered out, the DEG SHMs can still be displayed.
    na.deg <- rownames(df.deg)
    df.fil <- round(assay(se.nor), 2)
    # Switch data sets.
    if (!all(na.deg %in% rownames(df.fil))) req('')
    if (ncol(df.met) > 0) df.deg <- cbind.data.frame(df.met[na.deg, , drop = FALSE], df.deg, stringsAsFactors=FALSE)
    df.deg <- cbind.data.frame(df.deg, df.fil[na.deg, colnames(se.sub), drop = FALSE], stringsAsFactors=FALSE)
    cat('Done! \n'); df.deg 
  })
    output$dld.ssg.tab <- downloadHandler(
      filename=function(){ "spatial_enrichment.txt" },  content=function(file=paste0(tmp.dir, '/spatial_enrichment.txt')){
      write.table(dt.deg(), file, sep='\t', col.names=TRUE, row.names=TRUE) }
    )
    
  output$dt.deg <- renderDataTable({
    cat('DEG summary table ... \n'); dt.deg <- dt.deg()
    if (is.null(dt.deg)) return()
    col1 <- list(list(targets=c(1), render=DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    # In case no metadata column.
    if (colnames(dt.deg)[1]!='metadata') col1 <- NULL
    d.tab <- datatable(dt.deg, selection=list(mode="multiple", target="row"), escape=FALSE, filter="top", extensions='Scroller', plugins = "ellipsis",
    options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=col1), 
    class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer'); cat('Done! \n')
    d.tab
  })
  dat.all <- eventReactive(list(se.nor(), input$run, dat()), {
    dat <- dat(); run <- input$run
    if (!check_obj(list(run, dat))) req('')
    if (run==0) dat.all <- dat$se else {
      se.nor <- se.nor(); norMeth <- input$norMeth
      if (!check_obj(list(se.nor, norMeth))) req('')
      dat.all <- se.nor
    }; df.met <- rowData(dat.all)
    if ('metadata' %in% colnames(df.met)) dat.all <- cbind.data.frame(df.met[, 'metadata', drop=FALSE], round(assay(dat.all), 2), stringsAsFactors=FALSE) else dat.all <- round(assay(dat.all), 2)
    dat.all
  })
 output$datAll <- renderDataTable({
    message('DEG: complete data matrix ...')
    dat.all <- dat.all(); run <- input$run
    if (!check_obj(list(run))) return()
    if (!check_obj(list(dat.all)) & run > 0) return()
    col1 <- list(list(targets = c(1), render = DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    # In case no metadata column.
    if (colnames(dat.all)[1]!='metadata') col1 <- NULL
    dtab <- datatable(dat.all, selection='none', escape=FALSE, filter="top", extensions=c('Scroller', 'FixedColumns'), plugins = "ellipsis",
   options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=TRUE, columnDefs=col1, fixedColumns = list(leftColumns=2)),
   class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer')
    message('Done!'); dtab
  })

  dat.mod.lis.deg <- reactiveValues()
  #sch.mod.lis <- reactiveValues()
  #observe({
  #  cat('Preparing search box in DEG section ... \n')
  #  dat <- dat(); dat.deg <- dt.deg() 
  #  if (!check_obj(list(dat, dat.deg))) req('')
  #  deg.rna <- rownames(dat.deg)
  #  dat.mod.lis.deg$se.scl <- reactive({ dat$se[deg.rna, ] })
  #  sch.mod.lis$val <- search_server('deg', ids, cfg, lis.url, url.id, dat.mod.lis.deg)
  #  cat('Done! \n')
  #})
  #but.sgl <- reactive({ sch.mod.lis$val$ids$but.sgl })
  #but.mul <- reactive({ sch.mod.lis$val$ids$but.mul })
  onBookmark(function(state) { state })
  # return(list(but.sgl=but.sgl, but.mul=but.mul, query.res=query.res, input=input))
  return(list(query.res=query.res, input=input))
})}

