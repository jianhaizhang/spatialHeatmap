# Module for spatial enrichment.
deg_server <- function(id, sch, lis.url, url.id, ids, upl.mod.lis, dat.mod.lis, shm.mod.lis, session) {

  moduleServer(id, function(input, output, session) {
  ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg

  gID <- shm.mod.lis$gID
  cat('Presenting data matrix (DEG) ... \n')
  dat.deg.mod.lis <- data_server('datDEG', sch, lis.url, ids, deg=TRUE, upl.mod.lis)
  cat('Done! \n')
  geneIn1 <- dat.mod.lis$geneIn1 # Take transformed values in DEG table.
  geneIn <- dat.deg.mod.lis$geneIn # Take filted matrix with replicates. 
  ssg.dat <- reactive({
    # if (input$fileIn!="customBulkData") return(NULL)
    cat('DEG: SummarizedExperiment, features, variables ... \n')
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
    selectInput(ns('ssg.con'), label='Select variables', choices=cho, selected=cho[2:3], multiple=TRUE)
  })

  # Subset SE.
  sub_se <- function(se, sams, cons) {
    cna.sel <- NULL; for (i in cons) { 
      cna.sel <- c(cna.sel, paste0(sams, '__', i)) }
    se <- se[, colnames(se) %in% cna.sel]; return(se)
  }

  se <- reactive({
    cat('Subsetting SE with input features/variables ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
    dat <- ssg.dat(); sam.con <- input$sam.con 
    if (is.null(geneIn())|is.null(dat)|is.null(sam.con)) return(NULL)
    sam <- input$ssg.sam; con <- input$ssg.con; se <- dat$se
    if (is.null(sam)|is.null(con)) return() 
    if ('all' %in% sam) sam <- unique(dat$sams)
    if ('all' %in% con) con <- unique(dat$cons)
    se <- sub_se(se, sam, con)
    if (sam.con=='feature') fct <- 'sample' else if (sam.con=='variable') fct <- 'condition' else if (sam.con=='feature__variable') fct <- 'samCon'
    fct.tab <- table(colData(se)[, fct])
    fct.na <- names(fct.tab)[fct.tab==1] 
    validate(need(length(fct.na)==0, paste0('At least 2 replicates are required: ', paste0(fct.na, collapse=', '))))
    cat('Done! \n'); return(list(se=se, fct=fct))
  })

  se.nor.rok.dis <- reactive({
    cat('Normalizing data for ROKU/distinct ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
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
    # if (input$fileIn!="customBulkData") return(NULL)
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
    if (sam.con %in% c('feature', 'variable')) {
      if (sam.con == 'feature') {  
        under <- con; ref <- setdiff(sam, tar)
        validate(need(length(ref)>0, 'If compare by "feature", at least 2 features should be selected!'))
      } else if (sam.con == 'variable') {
        under <- sam; ref <- setdiff(con, tar)
        validate(need(length(ref)>0, 'If compare by "variable", at least 2 variables should be selected!'))
      }
      tar.all <- paste0(tar, '__', under)
      ref.all <- paste0(unlist(lapply(ref, function(i) paste0(i, '__', under))))
      vs <- rbind(vs, c(tar.all, 'VS', ref.all))
      colnames(vs) <- c(paste0('target', seq_along(tar.all)), 'VS', paste0('reference', seq_along(ref.all))) 
    } else if (sam.con == 'feature__variable') {
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
    # if (input$fileIn!="customBulkData") return(NULL)
    if (is.null(geneIn())|is.null(se())) return(NULL)
    withProgress(message="edgeR: ", value=0, {
    incProgress(0.5, detail="in progress ...")
    edg <- edgeR(se=se()$se, method.norm=sub('.*-', '', input$edg.lim.nor), com.factor=se()$fct, method.adjust='BH', return.all=TRUE); cat('Done! \n'); edg
    })
  })
  edg <- eventReactive(input$ssg.update, {
    cat('edgeR log2/fc ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
    if (!'edgeR' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(edg0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct]) 
    up.dn <- up_dn(sam.all=lvl, df.all=edg0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='logFC', fdr.na='FDR'); cat('Done! \n')
    up.dn
  })

  dsq0 <- reactive({
    cat('DESeq2 all ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
    if (is.null(geneIn())|is.null(se())) return(NULL)
    withProgress(message="DESeq2: ", value=0, {
    incProgress(0.5, detail="in progress ...")
    dsq <- deseq2(se=se()$se, com.factor=se()$fct, method.adjust='BH', return.all=TRUE); cat('Done! \n')
    dsq
    })
  })
  dsq <- eventReactive(input$ssg.update, {
    cat('DESeq2 log2/fc ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
    if (!'DESeq2' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(dsq0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct])
    up.dn <- up_dn(sam.all=lvl, df.all=dsq0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='log2FoldChange', fdr.na='padj'); cat('Done! \n'); up.dn
  })

  lim0 <- reactive({
    cat('limma all ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
    if (is.null(geneIn())|is.null(se())) return(NULL)
    withProgress(message="limma: ", value=0, {
    incProgress(0.5, detail="in progress ...") 
    lim <- limma(se=se()$se, method.norm=sub('.*-', '', input$edg.lim.nor), m.array=FALSE, com.factor=se()$fct, method.adjust='BH', return.all=TRUE); cat('Done! \n'); lim
    })
  })
  lim <- eventReactive(input$ssg.update, {
    cat('limma log2/fc ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
    if (!'limma' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(input$ssg.fc)|is.null(input$ssg.fdr)|is.null(lim0())) return(NULL)
    lvl <- unique(colData(se()$se)[, se()$fct])
    up.dn <- up_dn(sam.all=lvl, df.all=lim0(), log.fc=abs(input$ssg.fc), fdr=input$ssg.fdr, log.na='logFC', fdr.na='adj.P.Val'); cat('Done! \n'); up.dn
  })

  rok0 <- reactive({
    cat('ROKU all ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
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
    # if (input$fileIn!="customBulkData") return(NULL)
    if (!'ROKU' %in% input$ssg.meth) return(NULL)
    if (is.null(geneIn())|is.null(se())|is.null(rok0())) return(NULL)
    cat('ROKU up/down ... \n')
    up_dn_roku(rok0())
  })

  dis0 <- reactive({
    cat('distinct all ... \n')
    # if (input$fileIn!="customBulkData") return(NULL)
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
    # if (input$fileIn!="customBulkData") return(NULL)
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
      na0 <- names(deg.lis1)
      deg.lis1 <- c(deg.lis1[1], deg.lis1[1])
      names(deg.lis1) <- c(na0, sub('\\.', '1\\.', na0))
    }
    if (length(deg.lis2) == 1) { # Only one method is selected.
      na0 <- names(deg.lis2)
      deg.lis2 <- c(deg.lis2[1], deg.lis2[1])
      names(deg.lis2) <- c(na0, sub('\\.', '1\\.', na0))
    }
    cat('Done! \n'); return(list(up.lis = deg.lis1, down.lis = deg.lis2))
  })

  observeEvent(deg.lis(), ignoreInit=FALSE, ignoreNULL=TRUE, {
    if (is.null(deg.lis())) return()
    updateNavbarPage (session, "degAll", selected='ovlMeth')
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
      filename=function(){ "tissue_specific_genes.txt" },  content=function(file=paste0(tmp.dir, '/tissue_specific_genes.txt')){
      write.table(dt.deg(), file, sep='\t', col.names=TRUE, row.names=TRUE) }
    )
    
  output$dt.deg <- renderDataTable({
    cat('DEG summary table ... \n'); gene.dt <- dt.deg()
    if (is.null(gene.dt)) return()
    col1 <- list(list(targets=c(1), render=DT::JS("$.fn.dataTable.render.ellipsis(40, false)")))
    # In case no metadata column.
    if (colnames(gene.dt)[1]!='metadata') col1 <- NULL
    d.tab <- datatable(gene.dt, selection='none', escape=FALSE, filter="top", extensions='Scroller', plugins = "ellipsis",
      options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=300, scroller=TRUE, searchHighlight=TRUE, search=list(regex=TRUE, smart=FALSE, caseInsensitive=TRUE), searching=FALSE, columnDefs=col1), 
      class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer'); cat('Done! \n')
      d.tab
  }) 
  dat.mod.lis.deg <- reactiveValues()
  sch.mod.lis <- reactiveValues()
  observe({
    cat('Preparing search box in DEG section ... \n')
    dat.deg <- dt.deg(); if (is.null(dat.deg)) return()
    dat.lis <- dat.mod.lis; if (is.null(dat.lis)) return()
    deg.rna <- rownames(dat.deg); gen.lis <- dat.lis$geneIn()
    dat.mod.lis.deg$geneIn <- reactive({  
    return(list(df.aggr=gen.lis$df.aggr[deg.rna, , drop=FALSE], df.aggr.tran=gen.lis$df.aggr.tran[deg.rna, , drop=FALSE], df.aggr.tran.order=gen.lis$df.aggr.tran.order[deg.rna, , drop=FALSE], df.met=gen.lis$df.met[deg.rna, , drop=FALSE], df.rep=gen.lis$df.rep[deg.rna, , drop=FALSE]))
    })
   sch.mod.lis$val <- search_server('deg', ids, cfg, lis.url, url.id, dat.mod.lis.deg)
   cat('Done! \n')
  })
  onBookmark(function(state) { state })
  but.sgl <- reactive({ sch.mod.lis$val$ids$but.sgl })
  but.mul <- reactive({ sch.mod.lis$val$ids$but.mul })
  return(list(but.sgl=but.sgl, but.mul=but.mul))
})}

