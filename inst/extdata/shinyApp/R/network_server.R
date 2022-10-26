# Module for matrix heatmap and network analysis.
network_server <- function(id, upl.mod.lis, dat.mod.lis, shm.mod.lis, sch.mod.lis, session) {
  
  moduleServer(id, function(input, output, session) {
    ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg

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
    updateActionButton(session, inputId='mhm.but', icon=icon("sync"))
    #updateRadioButtons(session, inputId="mhm.but", label="Show plot:", choices=c("Yes", "No"), selected=cfg$lis.par$mhm['show', 'default'], inline=TRUE)
  })
# Avoid unnecessay correlation/distance computation if geneIn() updates due to column reordering. For example, correlation/distance, extracting coordinates, which only depend on the df.aggr.
  df.net <- reactive({
    if (is.null(geneIn())) return()
    df.aggr <- geneIn()[['df.aggr']]
    if (all(df.aggr==round(df.aggr))) df.aggr <- log2(df.aggr+1) 
    return(list(df.aggr=df.aggr, df.met=geneIn()[['df.met']]))
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
    # ggsave(paste0(tmp.dir, '/shm_all.', ext), plot=shm.all, device=ext, width=width/72, height=height/72, dpi=res, unit='in', limitsize = FALSE) 
      cat('Initial adjacency matrix and modules ...\n')
      if (ipt$fileIn=="none"|tab.net$val!='yes') return()
      if (is.null(input$cpt.nw)) return() # Network section is removed.
      if (is.null(submat())|input$cpt.nw!=0|length(gID$geneSel)>1) return()
    if (ipt$fileIn=="customBulkData"|any(ipt$fileIn %in% cfg$na.def)) {
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
       library(WGCNA); lis <- adj_mod(data=sub.mat, type=input$net.type, minSize=input$min.size); cat('Done! \n'); lis
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

   # if ((input$cpt.nw!=0|!is.null(adj.mods$lis)) & (input$fileIn=="customBulkData"|any(input$fileIn %in% cfg$na.def))) { 

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
    # if (ipt$fileIn=="customComputedData") { adj <- adj.mod()[['adj']]; mod <- adj.mod()[['mcol']] } else 
    # if (ipt$fileIn=="customBulkData"|ipt$fileIn %in% cfg$na.def) { 
      adj <- adj.mods.lis[['adj']]; mod <- adj.mods.lis[['mod']]
    # }
    # if (ipt$fileIn=='customComputedData') gene <- df.net()$df.aggr else {
    gene <- submat(); if (is.null(gene)) return()
    if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    ds <- input$ds; # save(mod, ds, file='md')
    lab <- mod[, input$ds][rownames(gene)==input$gen.sel]
    # validate(need(try(length(lab)==1 & !is.na(lab) & nrow(mod)==nrow(gene)), 'Click "Update" to display new network!'))
    if (length(lab)==0) return()
    if (length(lab)>1|is.na(lab)) return() # When input$fileIn is changed, gene is changed also, but mod is not since it is controled by observeEvent.
    validate(need(try(lab!='0'), 'Warning: the selected gene is not assigned to any module. Please select a different one or adjust the "Minmum module size"!'))
    idx.m <- mod[, input$ds]==lab; adj.m <- adj[idx.m, idx.m]; gen.na <- colnames(adj.m) 
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
          nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mod, adj=adj, geneID=ID, adj.min=adjs)
          lin <- nrow(nod.lin[['link']])
          vec0 <- adjs; names(vec0) <- lin
          adj.lin.vec <- c(adj.lin.vec, vec0)
          if (adjs==0) break

      }; cat('Adjacency-edge pairs done! \n')
      # The first version of links computed from the min adj or the input adj, which is subjected to the following check.
      nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mod, adj=adj, geneID=ID, adj.min=ifelse(input$adj.in==1, adjs, input$adj.in))
      link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'
      # If the links are 0 due to the input adj, change the "adjs" to the value bringing 1 or 2 links.
      lins <- NULL; if (nrow(link1)==0) {

        adjs <- adj.lin.vec[names(adj.lin.vec)>=1][1]
        nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mod, adj=adj, geneID=ID, adj.min=adjs)
        link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)

      } else if (nrow(link1)>input$max.edg) {
       
        # If the links are larger than max links due to the input adj, change the "adjs" to the value producing max links.
        adjs <- adjs+0.002
        nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mod, adj=adj, geneID=ID, adj.min=adjs)
        link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)

      } else if (nrow(link1)<=input$max.edg & input$adj.in!=1) {
        
        # If 0<link total<max links, use the input adj.
        adjs <- input$adj.in
        nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mod, adj=adj, geneID=ID, adj.min=adjs)
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
      if (ipt$fileIn=="none"|(ipt$fileIn=="customBulkData" & is.null(geneIn()))|
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
