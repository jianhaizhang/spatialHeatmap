

net_server <- function(id, submat, ids, gID, ipt.up, ipt.ana, cfg, data, df.net, clus.h=NULL, clus.k=NULL, session) {  
moduleServer(id, function(input, output, session) {
  ns <- session$ns
  observe({
    library(flashClust); library(visNetwork)
  })
  # er <- eventReactive(exp, {}). If its reactive value "er()" is called before eventReactive is triggered, the code execution stops where "er()" is called.
  data <- eventReactive(list(ipt.ana$showBut, ipt.ana$showButNet), {
    sub.mat <- submat(); method <- ipt.ana$method; cl <- NULL
    if (!check_obj(list(sub.mat, method))) return() 
    if ('net' %in% method & id=='netS3') return()
    if ('net' %in% method) return(sub.mat)
    lgc.op <- !'net' %in% ipt.ana$method & 'shmAll-net-netOp' %in% ipt.ana$clusNav
    if (lgc.op) {
    hclus <- clus.h$hclus
    if ('hcl' %in% method & !is.null(hclus)) { cl <- hclus() }
    kclus <- clus.k$lis
    if ('kmean' %in% method & !is.null(kclus)) cl <- kclus()$cluster
    if (!check_obj(list(cl))) return()
    sub.mat[cl, , drop=FALSE]
    }
  }); observe({ data() })
  adj.mod.par <- reactiveValues()
  observeEvent(list(ipt.ana$showBut, ipt.ana$showButNet), ignoreInit=FALSE, {
    if (!'net' %in% ipt.ana$method & !'shmAll-net-netOp' %in% ipt.ana$clusNav) return()
    if ('net' %in% ipt.ana$method & 'shmAll-net-netOp' %in% ipt.ana$clusNav) return()
    pars <- list(query=ipt.ana$query, ds=input$ds, type=input$net.type, min.size=input$min.size, data())
    adj.mod.par$pars <- pars
  })
  adj.mods <- eventReactive(adj.mod.par$pars, {
    cat('Adjacency and modules ... \n')
    fileIn <- ipt.up$fileIn; req(!dat.no %in% fileIn)
    if (fileIn %in% cfg$na.cus & is.null(ipt.up$svgInpath1) & is.null(ipt.up$svgInpath2)) return()
    #gene <- geneIn()[["df.aggr.tran"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    sub.mat <- data(); if (is.null(sub.mat)) return()
    withProgress(message="Network modules: ", value = 0, {
      incProgress(0.3, detail="adjacency matrix ...")
      incProgress(0.5, detail="topological overlap matrix ...")
       library(WGCNA); lis <- adj_mod(data=sub.mat, type=input$net.type, minSize=input$min.size, ds=2:3)
      incProgress(0.1, detail="dynamic tree cutting ...")
      cat('Done!\n'); lis
      })
    }); observe({ adj.mods() })

  #observe({
    # if (is.null(geneIn())) return(NULL)
  #  if (length(ids$sel)==0) return()
  #  gens.sel <- ids$sel
  #  updateSelectInput(session, inputId="gen.sel", choices=c("None", gens.sel), selected=gens.sel[1])
  #})
  observe({ 
    ipt.ana$query; # input$measure; input$thr; input$mhm.v
    lis.par <- cfg$lis.par; req(check_obj(lis.par))
    updateSelectInput(session, 'ds', choices=3:2, selected=lis.par$network['ds', 'default'])
  })

  col.sch.net <- reactive({ 
    if(input$color.net=="") return() 
    col <- gsub(' |\\.|-|;|,|/', '_', input$color.net)
    col <- strsplit(col, '_')[[1]]
    col <- col[col!='']; col1 <- col[!col %in% colors()]
    if (length(col1>0)) { 
      msg <- paste0('Invalid colors: ', paste0(col1, collapse=', '), '!') 
      showModal(modal(msg = msg))
      validate(need(try(col1 %in% colors()), '')); 
    }; col
  }); color.net <- reactiveValues(col.net="none")
  len.cs.net <- 350
  observe({
    lis.par <- cfg$lis.par; but <- input$col.but.net
    req(check_obj(list(lis.par, but)))
    if(but==0) color.net$col.net <- colorRampPalette(col_sep(lis.par$network['color', 'default']))(len.cs.net)
  })
  observeEvent(input$col.but.net, {
    col.sch <- col.sch.net()
    if (is.null(col.sch)) return ()
    color.net$col.net <- colorRampPalette(col.sch)(len.cs.net)
  })
  nod.edg.par <- reactiveValues()
  observeEvent(list(ipt.ana$showBut, ipt.ana$showButNet, adj.mods()), ignoreInit=FALSE, {
    if (!'net' %in% ipt.ana$method & !'shmAll-net-netOp' %in% ipt.ana$clusNav) return()
    if ('net' %in% ipt.ana$method & 'shmAll-net-netOp' %in% ipt.ana$clusNav) return()
    pars <- list(ipt.ana$query, input$ds, input$net.type, input$min.size, input$max.edg, input$adj.in, adj.mods())
    nod.edg.par$pars <- pars
  })
  nod.edg <- eventReactive(nod.edg.par$pars, {
    cat('Network nodes and edges ... \n')
    adj.mods.lis <- adj.mods(); gene <- data()
    ds <- input$ds; query <- ipt.ana$query
    max.edg <- input$max.edg; adj.in <- input$adj.in
    if (!check_obj(list(adj.mods.lis, ds, query, gene, max.edg, adj.in))) return()
    if (query=='None') return()
    adj <- adj.mods.lis[['adj']]; mod <- adj.mods.lis[['mod']]
    if (!(query %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    lab <- mod[, ds][rownames(gene)==query]
    if (length(lab)==0) return()
    # When input$fileIn is changed, gene is changed also, but mod is not since it is controled by observeEvent.
    if (length(lab)>1|is.na(lab)) return()
    if (lab=='0') {
      msg <- 'The query biomolecule is not assigned to any module!' 
      showModal(modal(msg = msg))
      validate(need(try(lab!='0'), ''))
    }
    idx.m <- mod[, ds]==lab; adj.m <- adj[idx.m, idx.m]
    gen.na <- colnames(adj.m) 
    idx.sel <- grep(paste0("^", query, "$"), gen.na)
    gen.na[idx.sel] <- paste0(query, "_target")
    colnames(adj.m) <- rownames(adj.m) <- gen.na
    withProgress(message="Module properties:", value=0, { 
      incProgress(0.8, detail="Nodes and edges ...")
      cat('Extracting nodes and edges... \n')
      # Identify adjcency threshold with edges < max number (e.g. 300) 
      ID <- query; adjs <- 1; lin <- 0; adj.lin.vec <- NULL
      if (round(max.edg)!=max.edg) {
        showModal(modal(msg = 'The number of edges should be an integer!'))
        validate(need(try(round(max.edg)==max.edg), ''))
      }
      # Compute the min adj.
      while (lin < max.edg) {    
          adjs <- adjs-0.002; if (adjs<=10^-15) adjs <- 0
          nod.lin <- nod_lin(ds=ds, lab=lab, mods=mod, adj=adj, geneID=ID, adj.min=adjs)
          lin <- nrow(nod.lin[['link']])
          vec0 <- adjs; names(vec0) <- lin
          adj.lin.vec <- c(adj.lin.vec, vec0)
          if (adjs==0) break
      }; message('Adjacency-edge pairs done!')
      # The first version of links computed from the min adj or the input adj, which is subjected to the following check.
      nod.lin <- nod_lin(ds=ds, lab=lab, mods=mod, adj=adj, geneID=ID, adj.min=ifelse(adj.in==1, adjs, adj.in))
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
      } else if (nrow(link1)<=max.edg & adj.in!=1) {
        # If 0<link total<max links, use the input adj.
        adjs <- adj.in
        nod.lin <- nod_lin(ds=ds, lab=lab, mods=mod, adj=adj, geneID=ID, adj.min=adjs)
        link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)
      }
      node <- nod.lin[['node']]; colnames(node) <- c('id', 'value') 
      if (nrow(link1)!=0) {    
        # 'length' is not well indicative of adjacency value, so replaced by 'value'.
        link1$title <- link1$value; link1$color <- 'lightblue' 
      }; df.met <- rowData(df.net()$se.thr)
      if ('metadata' %in% colnames(df.met)) meta <- df.met[, 'metadata', drop=FALSE] else meta <- data.frame()
      # if (ncol(meta)>0) node <- cbind(node, title=meta[sub('_target$', '', node$id), ], borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      df.other <- cbind(borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      if (ncol(meta)>0) node <- cbind(node, title=meta[sub('_target$', '', node$id), ], df.other) else node <- cbind(node, df.other) 
      # if (ncol(meta)==0) node <- cbind(node, borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      net.lis <- list(node=node, link=link1, adjs=adjs, lins=lins)
    }); cat('Done! \n'); net.lis
  })
  mod <- reactive({ 
    nod.edg <- nod.edg(); if (!check_obj(list(nod.edg))) return()    
    mod <- sub('_target$', '', nod.edg$node$id); mod
  }); observe({ mod() }) # reactive: only after called, will it be updated.
  # The order of reactive expression matters so "updateSelectInput" of "adj.in" should be after visNet().
  observe({ 
    data(); gID$geneSel; ipt.ana$query; input$ds; input$adj.modInpath; input$A; input$P; input$CV1; input$CV2; input$min.size; input$net.type
    ipt.ana$showBut; ipt.ana$showButNet
    # input$measure; input$thr; input$mhm.v
     #if ((input$adj.in==1 & is.null(visNet()[["adjs1"]]))|(input$cpt.nw!=cfg$lis.par$network['max.edges', 'default'] & is.null(visNet()[["adjs1"]]))) { updateSelectInput(session, "adj.in", "Adjacency threshold:", sort(seq(0, 1, 0.002), decreasing=TRUE), visNet()[["adjs"]]) } else if (!is.null(visNet()[["adjs1"]])) updateSelectInput(session, "adj.in", "Adjacency threshold:", sort(seq(0, 1, 0.002), decreasing=TRUE), visNet()[["adjs1"]])
    # if (is.null(visNet$val)) return(); if (is.null(visNet$val())) return()
    nd.eg <- nod.edg(); adj.in <- input$adj.in
    if (!check_obj(list(nd.eg, adj.in))) return()
    lins <- nd.eg[["lins"]]
    if (input$adj.in==1|is.null(lins)|is.numeric(lins)) updateSelectInput(session, "adj.in", choices=sort(seq(0, 1, 0.002), decreasing=TRUE), selected=as.numeric(nd.eg[["adjs"]])) 
  })
  observeEvent(list(nod.edg(), color.net$col.net), {
  output$bar.net <- renderPlot({ 
    message('Network bar ...')
    nd.eg <- nod.edg(); query <- ipt.ana$query
    col.net <- color.net$col.net
    se.thr <- df.net()$se.thr
    if (!check_obj(list(nd.eg, query, col.net, se.thr))) req('')
    if (query=="None"|input$adj.in=="None"|col.net[1]=='none') return()
    gene <- assay(se.thr)
    if (!(query %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    withProgress(message="Network color key: ", value = 0, {
    incProgress(0.25, detail="preparing data ...")
      node <- nd.eg[["node"]]; if (is.null(node)) return()
      node.v <- node$value
      v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
      cs.net <- col_bar(geneV=v.net, cols=col.net, width=1)
      incProgress(0.75, detail="plotting ...")
      message('Done!'); return(cs.net) # '((max(v.net)-min(v.net))/len.cs.net)*0.7' avoids bar overlap.
    })
  })
  })
  observeEvent(nod.edg(), ignoreInit=FALSE, {
    nd.eg <- nod.edg(); if (is.null(nd.eg)) return()
    output$edge <- renderUI({ 
      if (input$adj.in=="None") return()
      if ((ipt.up$fileIn=="customBulkData" & is.null(data()))|
      ipt.ana$query=="None") return()
      span(style = "color:black;font-weight:NULL;", HTML(paste0("Remaining edges: ", dim((nd.eg[["link"]]))[1])))
    })
  })
  vis.net <- eventReactive(list(nod.edg(), color.net$col.net), {
    message('Network ...')
    if (input$adj.in=="None") return(NULL)
    col.net <- color.net$col.net; nd.eg <- nod.edg()
    query <- ipt.ana$query; se.thr <- df.net()$se.thr
    if (!check_obj(list(nd.eg, query, se.thr, col.net))) return()
    gene <- assay(se.thr)
    if (query=='None'|col.net[1]=="none") return()
    if (!(query %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    withProgress(message="Plotting module:", value=0.5, {
    incProgress(0.3, detail="preparing ...")
    # Match colours with gene connectivity by approximation.
    node <- nd.eg[["node"]]; if (is.null(node)) return() 
    node.v <- node$value
    v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
    col.nod <- NULL; for (i in node$value) {
      ab <- abs(i-v.net); col.nod <- c(col.nod, col.net[which(ab==min(ab))[1]])
    }; node$color <- col.nod
    message('Done!'); # save(node, nd.eg, file='nn')
    visNetwork(node, nd.eg[["link"]], height="300px", width="100%", background="", main=paste0("Module Containing ", query), submain="", footer= "") %>% visIgraphLayout(physics=FALSE, smooth=TRUE) %>% visOptions(highlightNearest=list(enabled=TRUE, hover=TRUE), nodesIdSelection=TRUE)
    }) 
  })
  output$vis <- renderVisNetwork({
    if (!check_obj(list(ids$sel, vis.net()))) return()
    withProgress(message="Plotting module:", value=0.5, {
      incProgress(0.3, detail="plotting ...")
      message('Rendering network ... '); vis.net()
    })
  }) 
  
  observeEvent(list(mod()), {
  output$dld <- renderUI({  
    query <- ipt.ana$query; ds <- input$ds; mod <- mod()
    adj.mods <- adj.mods(); nod.edg <- nod.edg()
    if (!check_obj(list(query, ds, mod, adj.mods, nod.edg))) return()
    if (query=='None') return()

    pa <- 'tmp/module.txt'
    write(paste0('# Module containing the query: ', paste0(mod, collapse=','), '\n# Connectivities of the module containing the query'), file.path('www', pa))
    # Connectivity.
    if ('title' %in% colnames(nod.edg$node)) {
      df.con <- nod.edg$node[, c('id', 'value', 'title')]
      colnames(df.con) <- c('node', 'connectivity', 'description')
    } else {
      df.con <- nod.edg$node[, c('id', 'value')]
      colnames(df.con) <- c('node', 'connectivity')
    }
    df.con <- rbind(colnames(df.con), df.con)
    write.table(df.con, file.path('www', pa), col.names=FALSE, row.names=FALSE, sep='\t', append=TRUE)
    # Adjacency.
    write('# Adjacencies of the module containing the query', file.path('www', pa), append=TRUE)
    df.adj <- nod.edg$link[, c('from', 'to', 'value')]
    colnames(df.adj) <- c('node1', 'node2', 'adjacency')
    df.adj <- rbind(colnames(df.adj), df.adj)
    write.table(df.adj, file.path('www', pa), col.names=FALSE, row.names=FALSE, sep='\t', append=TRUE)
    # All modules.
    df.mod <- as.data.frame(adj.mods$mod)
    df.mod <- cbind(node=rownames(df.mod), df.mod)
    df.m <- df.mod[, c('node', as.character(ds))]
    colnames(df.m)[2] <- 'module'
    df.m <- df.m[order(-df.m$module), ]
    # All modules.
    write('# All modules', file.path('www', pa), append=TRUE)
    df.m <- rbind(colnames(df.m), df.m)
    write.table(df.m, file.path('www', pa), col.names=FALSE, row.names=FALSE, sep='\t', append=TRUE)
    a(href=pa, target='blank', 'Download') 
  })
  })

  onBookmark(function(state) { state })
  return(list(mod=mod))
})}


kmean_server <- function(id, ipt.ana, submat, session) {  
moduleServer(id, function(input, output, session) {
  dat.scl <- eventReactive(ipt.ana$showBut, ignoreInit=FALSE, {
    if (! 'kmean' %in% ipt.ana$method) return()
    query <- ipt.ana$query; submat <- submat()
    if (!check_obj(list(submat, query))) return()
    if (query=='None') return()
    scale(submat)
  }); observe({ dat.scl() })

  op.par <- reactiveValues()
  observeEvent(list(ipt.ana$showBut, dat.scl()), ignoreInit=FALSE, {
    if (! 'kmean' %in% ipt.ana$method) return()
    pars <- list(input$maxK, dat.scl())
    op.par$pars <- pars
  })

  gg.op <- eventReactive(list(op.par$pars), {
    set.seed(input$seed)
    message('Optimal k ... ')
    maxK <- input$maxK; dat.scl <- dat.scl()
    if (!check_obj(list(op.par$pars, maxK, dat.scl))) return()
    # Reactive expression of sub.mat() is NULL, but sub.mat is not NULL.
    withProgress(message="K-means clustering across a range of k:", value=0, {
      incProgress(0.5, detail=" ...")
      gg <- optimal_k(dat.scl, maxK)
      incProgress(0.4, detail=" ...")
      message('Done!'); gg
    })
  })
  output$elbow <- renderPlot({
    gg.op <- gg.op()
    if (!check_obj(list(gg.op))) return()
    withProgress(message="Plotting results of the elbow method:", value=0, {
      incProgress(0.7, detail="...")
      gg.op
    })
  })
  
  clu.par <- reactiveValues()
  observeEvent(list(ipt.ana$showBut, dat.scl()), ignoreInit=FALSE, {
    if (! 'kmean' %in% ipt.ana$method) return()
    pars <- list(input$k, dat.scl())
    clu.par$pars <- pars
  })
  clus <- eventReactive(list(clu.par$pars), {
    message('Kmeans clustering ... ')
    k <- input$k; dat.scl <- dat.scl()
    query <- ipt.ana$query; 
    if (!check_obj(list(clu.par$pars, k, dat.scl, query))) return()
    if (query=='None') return()
    # Reactive expression of sub.mat() is NULL, but sub.mat is not NULL.
    withProgress(message="K-means clustering with chosen k:", value=0, {
      incProgress(0.5, detail=" ...")
      clus <- kmeans(dat.scl, k)
      clus.all <- clus$cluster 
      cl.q <- clus.all[names(clus.all)==query]
      cl <- unique(names(clus.all[clus.all==cl.q]))
      incProgress(0.4, detail=" ...")
      message('Done!'); list(res=clus, cluster=cl)
    })
  })
  output$cluster <- renderPlot({
    clus <- clus(); dat.scl <- dat.scl()   
    query <- ipt.ana$query; dimred <- input$dimred
    if (!check_obj(list(clus, dat.scl, query, dimred))) return()
    if (query=='None') return()
    withProgress(message="Plotting clusters:", value=0, {
      incProgress(0.7, detail=" ...")
      gg <- plot_kmeans(data=dat.scl, res=clus$res, query=query, dimred=dimred)
      incProgress(0.2, detail=" ...")
      message('Done!'); gg
    })
  }) 
  observeEvent(clus()$cluster, {
  output$dld <- renderUI({  
    query <- ipt.ana$query; lis <- clus()
    if (!check_obj(list(lis, query))) return()
    if (query=='None') return()
    # All clusters in a data frame.
    clus.all <- lis$res$cluster; bm <- names(clus.all)
    clus.all <- paste0('clus', lis$res$cluster)
    names(clus.all) <- bm
    df.clus <- data.frame(biomolecule=names(clus.all), cluster=clus.all)
    w <- which(names(clus.all)==query)
    cl.query <- clus.all[w]
    cl.qu.w <- which(clus.all==cl.query)
    df.clus$cluster[cl.qu.w] <- paste0(cl.query, '.query')
    df.clus <- df.clus[order(df.clus$cluster), ]
    
    pa <- 'tmp/kcluster.txt'
    writeLines(paste0('# Cluster containing the query: ', paste0( lis$cluster, collapse=','), '\n# All clusters'), file.path('www', pa))
    # Append all clusters in a data frame.
    df.clus <- rbind(colnames(df.clus), df.clus)
    write.table(df.clus, file.path('www', pa), col.names=FALSE, row.names=FALSE, sep='\t', append=TRUE)
    a(href=pa, target='blank', 'Download')
  })
  })

  onBookmark(function(state) { state })
  return(list(lis=clus))
})}


mhm_server <- function(id, ipt.ana, submat, session) {  
moduleServer(id, function(input, output, session) {
  mhm.par <- reactiveValues()
  observeEvent(ipt.ana$showBut, ignoreInit=FALSE, {
    if (! 'hcl' %in% ipt.ana$method) return()
    # Include all upstream parameter changes. Since correlation/distance paremeter changes are included in submat.par$pars, only submat.par$pars represents upstream parameter change.
    pars <- list(mat.scale=input$mat.scale, submat(), input$cutH)
    mhm.par$pars <- pars
  })
  # mhm <- reactiveValues(hm=NULL)
  # Plot matrix heatmap.
  mhm <- eventReactive(list(mhm.par$pars), {
    cat('Initial matrix heatmap ... \n')
    mat.scale <- input$mat.scale; sub.mat <- submat()
    query <- ipt.ana$query; h <- input$cutH
    if (!check_obj(list(mhm.par$pars, mat.scale, sub.mat, query, h))) return()
    if (query=='None') return()
    # In submat, there is a validate expression. Only if submat is included in reactive({}) not observe({}), can the message in validate be propagated to the user interface.
    # Reactive expression of sub.mat() is NULL, but sub.mat is not NULL.
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="Plotting ...")
      if (mat.scale=="Column") scale.hm <- 'column' else if (mat.scale=="Row") scale.hm <- 'row' else scale.hm <- 'no'
      lis <- matrix_hm(ID=query, data=sub.mat, col=c('yellow', 'red'), scale=scale.hm, cut.h=h, main='', title.size=10, static=FALSE)
      cat('Done! \n'); lis
    })
  }); observe({ mhm() })

  output$HMly <- renderPlotly({
    # if (length(ids$sel)==0) return()
    if (is.null(mhm())) return()
    #if (is.null(gID$geneSel)|is.null(submat())) return()
    #if (gID$geneSel[1]=='none'|is.na(gID$geneSel[1])) return()
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="plotting ...")
      mhm()$plot
    })
  })

  #cut.par <- reactiveValues()
  #observeEvent(list(ipt.ana$showBut, mhm()), ignoreInit=FALSE, {
  #  if (!'hcl' %in% ipt.ana$method) return()
  #  query <- ipt.ana$query; lis <- mhm(); h <- input$cutH
  #  if (!check_obj(list(query, lis, h))) return()
  # if (query=='None') return()
  #  cut.par$pars <- list(lis, h, query)
  #})
  # An eventReactive is not evaluated unless used somewhere even though triggered by an action button.
  hclus <- eventReactive(mhm()$cluster, {
    #message('Cutting row dendrograms ... ')
    #query <- ipt.ana$query; lis <- mhm(); h <- input$cutH
    #if (!check_obj(list(query, lis, h))) return()
    #if (query=='None') return()
    #dendro <- lis$dendro.row; cl <- cut_dendro(dendro, h, query)
    cl <- mhm()$cluster
    if (!check_obj(list(cl))) return()
    message('Cluster size ', length(cl)); message('Done!')
    msg <- paste0('Size of the cluster containing the query at current cutting height: ', length(cl), '.')
    showNotification(msg, closeButton = TRUE, duration = 60)
    cl
  }); observe({ hclus() })
  
  observeEvent(mhm()$cluster, {
  output$dld <- renderUI({  
    lis <- mhm(); if (!check_obj(lis)) return()
    pa <- 'tmp/hcluster.txt'
    writeLines(paste0('Cluster containing the query: ', paste0( lis$cluster, collapse=',')), file.path('www', pa))
    a(href=pa, target='blank', 'Download')
  })
  })

  onBookmark(function(state) { state })
  return(list(hclus=hclus))
})}

fun_enrich_server <- function(id, ipt.ana, clus.h, clus.k, mod, mod.op=NULL, session) {  
moduleServer(id, function(input, output, session) {
  enr.par <- reactiveValues()
  observeEvent(input$FunBut, ignoreInit=FALSE, {
    stpClus <- input$stpClus
    method <- ipt.ana$method; if (!check_obj(list(method, stpClus))) return()
    pars <- NULL; if (stpClus=='Step2') {
    if ('hcl' %in% method) pars <- list(sp=input$species, id=input$funID, cl=clus.h$hclus())
    if ('net' %in% method) pars <- list(sp=input$species, id=input$funID, cl=mod$mod())
    if ('kmean' %in% method) pars <- list(sp=input$species, id=input$funID, cl=clus.k$lis()$cluster)
    } else if (stpClus=='Step3') {
    if (!'net' %in% method) pars <- list(sp=input$species, id=input$funID, cl=mod.op$mod())
    }
    enr.par$pars <- pars
  })
  entr.id <- eventReactive(enr.par$pars, {
    cat('Converting ids to "ENTREZID" ... \n')
    # library(AnnotationDbi)
    species <- input$species; funID <- input$funID
    cl <- enr.par$pars$cl
    if (!check_obj(list(species, funID, cl))) return()
    if (species=='hum') { 
      library(org.Hs.eg.db); db <- org.Hs.eg.db 
    } else if (species=='mus') { 
      library(org.Mm.eg.db); db <- org.Mm.eg.db
    }  else if (species=='arab') {
      library(org.At.tair.db); db <- org.At.tair.db
    } else if (species=='fish') {
      library(org.Dr.eg.db); db <- org.Dr.eg.db
    } else if (species=='fly') {
      library(org.Dm.eg.db); db <- org.Dm.eg.db
    }

    if (funID=='ense') keytype <- 'ENSEMBL'
    if (funID=='unpr') keytype <- 'SYMBOL'
    if (funID=='tair') keytype <- 'TAIR'
    if (funID=='entr') { cat('Done! \n'); return(cl) }
    df.id <- check_exp(select(db, keys=cl, keytype=keytype, columns=c(keytype, 'ENTREZID')))
    if ('e' %in% df.id) { 
      msg <- 'No ids are mapped. Ensure the "species" and "Source ID" are correct!'
      showModal(modal(msg=msg, easyClose=TRUE)); return()
    }; cat('Done! \n'); unique(df.id$ENTREZID)
  })

  go.keg.par <- reactiveValues()
  observeEvent(list(input$FunBut, entr.id()), ignoreInit=FALSE, {
    pars <- list(input$ont, input$Padj, input$Pcut, entr.id())
    go.keg.par$pars <- pars
  })
  library(clusterProfiler)
  go <- eventReactive(list(go.keg.par$pars), {
    message('GO enrichment ... ')
    entriz <- entr.id(); ont <- input$ont; Padj <- input$Padj
    Pcut <- input$Pcut; species <- input$species
    if (!check_obj(list(entriz, ont, Padj, Pcut, species))) return()
    if (species=='hum') db <- 'org.Hs.eg.db'
    if (species=='mus') db <- 'org.Mm.eg.db'
    if (species=='arab') db <- 'org.At.tair.db'
    if (species=='fish') db <- 'org.Dr.eg.db'
    if (species=='fly') db <- 'org.Dm.eg.db'
    withProgress(message="GO enrichment:", value=0.2, {
    incProgress(0.3, detail="in progress ...")
    # save(entriz, db, ont, Padj, Pcut, file='go')
    go <- enrichGO(gene = entriz, keyType = "ENTREZID", OrgDb = db, ont = ont, pAdjustMethod = Padj, pvalueCutoff  = Pcut, qvalueCutoff  = 0.05, readable = TRUE)
    incProgress(0.3, detail=" done!")
    }); message('Done!') 
    if (nrow(as.data.frame(go))==0) { 
      msg <- 'No GO terms enriched for the current settings!'
      showModal(modal(msg = msg))
    }; go 
  })

  observeEvent(list(go(), input$terms), {
    go <- go(); terms <- input$terms
    if (!check_obj(list(go, terms))) return()
    output$go <- renderPlot({
      if (nrow(as.data.frame(go))==0) return()
      barplot(go, showCategory=terms)
    })
  })

  keg <- eventReactive(list(go.keg.par$pars), {
    message('KEGG enrichment ... ')
    entriz <- entr.id(); Padj <- input$PadjKeg
    Pcut <- input$PcutKeg; species <- input$species
    if (!check_obj(list(entriz, Padj, Pcut, species))) return()
    if (species=='hum') sp <- 'hsa'
    if (species=='mus') sp <- 'mmu'
    if (species=='arab') sp <- 'ath'
    if (species=='fish') sp <- 'dre'
    if (species=='fly') sp <- 'dme'
    withProgress(message="KEGG enrichment:", value=0.2, {
    incProgress(0.3, detail="in progress ...")
    # save(entriz, db, ont, Padj, Pcut, file='go')
    kk <- enrichKEGG(gene = entriz, keyType='kegg', organism = sp, pvalueCutoff = Pcut, pAdjustMethod=Padj, qvalueCutoff=0.05)
    incProgress(0.3, detail=" done!")
    }); message('Done!')
    if (nrow(as.data.frame(kk))==0) { 
      msg <- 'No KEGG terms enriched for the current settings!'
      showModal(modal(msg = msg))
    }; kk 
  })
  observeEvent(list(keg(), input$termsKeg), {
    keg <- keg(); terms <- input$termsKeg
    if (!check_obj(list(keg, terms))) return()
    output$keg <- renderPlot({
      if (nrow(as.data.frame(keg))==0) return()
      barplot(keg, showCategory=terms)
    })
  })

  observeEvent(list(go(), keg()), {
  output$dld <- renderUI({  
    go <- go(); keg <- keg()
    if (is.null(go) & is.null(keg)) return()
    
    pa <- 'tmp/function_enrich.txt'
    write('# Go enrichment', file.path('www', pa))
    # Connectivity.
    go <- as.data.frame(go); go <- rbind(colnames(go), go)
    write.table(go, file.path('www', pa), col.names=FALSE, row.names=FALSE, sep='\t', append=TRUE)
    write('# KEGG enrichment', file.path('www', pa), append=TRUE)
    keg <- as.data.frame(keg); keg <- rbind(colnames(keg), keg)
    write.table(keg, file.path('www', pa), col.names=FALSE, row.names=FALSE, sep='\t', append=TRUE)
    a(href=pa, target='blank', 'Download') 
  })
  })
  onBookmark(function(state) { state })
  return(list(sn=session))
})}

# Module for large-scale analysis.
analysis_server <- function(id, upl.mod.lis, dat.mod.lis, shm.mod.lis, ids, session) { 
  moduleServer(id, function(input, output, session) {
  ns <- session$ns; ipt <- upl.mod.lis$ipt; cfg <- upl.mod.lis$cfg
  ipt.dat <- reactiveValues()
  ipt.dat$dat <- dat.mod.lis$ipt.dat; sear <- dat.mod.lis$sear
  col.reorder <- reactiveValues(); col.reorder <- dat.mod.lis$col.reorder
  se.thr <- dat.mod.lis$se.thr
  data <- dat.mod.lis$dat
  ipt.dat$dat <- dat.mod.lis$ipt.dat; A <- dat.mod.lis$A
  P <- dat.mod.lis$P; CV1 <- dat.mod.lis$CV1
  CV2 <- dat.mod.lis$CV2
  gID <- shm.mod.lis$gID; 
  observe({
    se.thr(); ipt$adj.modInpath; input$A; input$P; input$CV1
    input$CV2; input$min.size; input$net.type
    input$measure; input$thr; input$mhm.v
    updateSelectInput(session, "mat.scale", choices=c("No", "Row", "Column"), selected="Row")
  })

  observe({  
    ipt$fileIn; se.thr(); input$adj.modInpath; input$A; input$P; input$CV1; input$CV2; ids$sel
    updateActionButton(session, inputId='showBut', icon=icon("sync"))
  })
  observe({
    if (!check_obj(list(se.thr(), ids$sel))) return()
    gens.sel <- ids$sel
    updateSelectInput(session, inputId="query", choices=c("None", gens.sel), selected=gens.sel[1])
  })
  observe({
    method <- input$method; if (!check_obj(method)) return()
    if (method=='hcl') { showElement('mhm'); hideElement('net'); hideElement('kmean') }
    if (method=='net') { showElement('net'); hideElement('mhm'); hideElement('kmean') }
    if (method=='kmean') { showElement('kmean'); hideElement('mhm'); hideElement('net') }
  })
# Avoid unnecessay correlation/distance computation if geneIn() updates due to column reordering. For example, correlation/distance, extracting coordinates, which only depend on the df.aggr.
  df.net <- reactive({
    if (is.null(se.thr())) return()
    se.thr <- se.thr(); assay <- assay(se.thr)
    if (all(assay==round(assay))) {
      if (min(assay)==0) assay <- assay + 1 
      assay(se.thr) <- log2(assay) 
    }; return(list(se.thr=se.thr))
  })
  # Combine all relevant parameters to pars. After all parameters are adjusted, they only are controlled by buttons, which avoids execution after each parameter is adjusted.
  # A series of steps: the first "ignoreInit" needs to be TRUE. Otherwise, the combined pars are not able to trigger the following steps.
  cor.dis.par <- reactiveValues()
  observeEvent(list(input$showBut), ignoreInit=TRUE, {
    measure <- input$measure; se.thr <- df.net()$se.thr
    req(check_obj(list(measure, se.thr)))
    pars <- list(measure=measure, se.thr=se.thr)
    cor.dis.par$pars <- pars
  })
  # cor.dis <- reactiveValues(val=NULL)
  # Calculate whole correlation or distance matrix.
  cor.dis <- eventReactive(list(cor.dis.par$pars), ignoreInit=FALSE, valueExpr={
    cat('Correlation/distance matrix ... \n')
    measure <- input$measure
    if (!check_obj(list(ipt$fileIn, measure))) return()
    if (ipt$fileIn %in% cfg$na.cus & is.null(ipt$svgInpath1) & is.null(ipt$svgInpath2)) return()
    se.thr <- df.net()$se.thr; if (is.null(se.thr)) return()
    gene <- assay(se.thr)
    if ((ipt$fileIn %in% cfg$na.cus & is.null(gene))|ipt$fileIn=="none") return(NULL)
    withProgress(message="Compute similarity/distance matrix: ", value = 0, {
      incProgress(0.5, detail="please wait ...")
      # Too many genes may crash the app.
      if (nrow(gene)>15000) showModal(modal(msg=strong('Too many genes (e.g. 15,000+) may crash the app!'), easyClose=TRUE))
      # if (nrow(gene)>15000 & input$showBut==0) return()
      if (input$measure=='cor') {
        m <- cor(x=t(gene)); cat('Done! \n')
        if (measure=='absCor') { m <- abs(m) }; return(m)
      } else if (input$measure=='dis') { cat('Done! \n'); return(-as.matrix(dist(x=gene))) }
    })
  })
  # Subset nearest neighbours for target genes based on correlation or distance matrix.
  # observe({
    # The submat reactive expression is only accessible inside the same obeserve environment.
    # submat <- reactive({})
   #  if (tab.mhm$val!='yes') { return() }
  submat.par <- reactiveValues()
  observeEvent(list(input$showBut, cor.dis()), ignoreInit=FALSE, {
    query <- input$query; req(query)
    if ('None' %in% query) return()
    pars <- list(input$thr, input$mhm.v, cor.dis.par$pars, input$query, cor.dis())
    submat.par$pars <- pars
  })
  # Avoid calling eventReactive with ignoreInit=TRUE in another eventReactive also with ignoreInit=TRUE, since the latter will be slient even though the former is triggered.
  submat <- eventReactive(list(submat.par$pars), { 
    cat('Subsetting nearest neighbors ... \n')
    corDis <- cor.dis(); query <- input$query
    se.thr <- df.net()$se.thr
    if (!check_obj(list(corDis, query, se.thr))) return()
    if (query=='None') return()
    gene <- assay(se.thr)
    msg <- 'Subsetting nearest neighbors: columns in the assay data are less than 5!'
    if (ncol(gene)<5) showNotification(msg, duration=2, closeButton = TRUE)
    # Validate filtering parameters in matrix heatmap. 
    measure <- input$measure; mhm.v <- input$mhm.v; thr <- input$thr
    if (thr=='p') {
      validate(need(try(mhm.v>0 & mhm.v<=1), 'Proportion should be between 0 to 1 !'))
    } else if (thr=='n') {
      validate(need(try(mhm.v>=1 & as.integer(mhm.v)==mhm.v & !is.na(mhm.v)), 'Number should be a positive integer !'))
      g.n <- nrow(corDis); if (mhm.v >= g.n) mhm.v <- g.n
    } else if (thr=='v' & measure=='cor') {
      validate(need(try(mhm.v>-1 & mhm.v <1), 'Correlation value should be between -1 to 1 !'))
    } else if (thr=='v' & measure=='dis') {
      validate(need(try(mhm.v>=0), 'Distance value should be non-negative !'))
    }
    withProgress(message="Selecting nearest neighbours: ", value = 0, {
      incProgress(0.5, detail="please wait ...")
      arg <- list(p=NULL, n=NULL, v=NULL)
      arg[names(arg) %in% thr] <- mhm.v
      if (input$measure=='dis' & thr=='v') arg['v'] <- -arg[['v']]
      if (!all(query %in% rownames(corDis))) return()    
      gen.na <- do.call(sub_na, c(mat=list(corDis), ID=list(query), arg))
      if (any(is.na(gen.na))) return()
      if (length(gen.na)==1) { 
        msg <- paste0('Only ', gen.na, ' is remaining!')
        showModal(modal(msg = msg))
        validate(need(try(length(gen.na)>=2), '')) 
      }; cat('Done! \n'); return(gene[gen.na, , drop=FALSE])
    })
  }); observe({ submat() })
 # })

  #mhm.mod <- reactiveValues()
  #observe({
  #  mhm.mod$v <- mhm_server('mhm', ipt=input, submat=submat, df.net=df.net)
  #})
  #hclu <- reactive({ mhm.mod$v$hclus(); })
  h.clus <- mhm_server('mhm', ipt.ana=input, submat=submat)
  k.lis <- kmean_server('kmean', ipt.ana=input, submat=submat)
  net.mod <- net_server('net', submat=submat, gID=gID, ids=ids, ipt.up=ipt, ipt.ana=input, cfg=cfg, data=data, df.net=df.net)
  net.op.mod <- net_server('netS3', submat=submat, gID=gID, ids=ids, ipt.up=ipt, ipt.ana=input, cfg=cfg, data=data, df.net=df.net, clus.h=h.clus, clus.k=k.lis)
  fun.mod <- fun_enrich_server('goKeg', ipt.ana=input, clus.h=h.clus, clus.k=k.lis, mod=net.mod, mod.op=net.op.mod)
  observe({
    method <- input$method; cl <- NULL
    if ('hcl' %in% method) cl <- h.clus$hclus()
    if ('net' %in% method) cl <- net.mod$mod()
    if ('kmean' %in% method) cl <- k.lis$lis()$cluster
    id <- fun.mod$sn$ns('funEnr')
    # print(list('id', id, method, cl, fun.mod, fun.mod$sn, fun.mod$sn$ns))
    if (!check_obj(list(method, cl))) shinyjs::hide(id = id) else shinyjs::showElement(id = id) 
  })

  output$help <- renderUI({ 
    tags$iframe(seamless="seamless", src= "html/shm_shiny_manual.html#24_Data_mining", width='100%', height='100%')   
  })
  observe({
    disable(selector='a[data-value="shmAll-net-netOp"]')
    hclus <- h.clus$hclus; kclus <- k.lis$lis
    method <- input$method
    if (!check_obj(method)) return()
    if ('net' %in% method) return()
    hc <- check_exp(check_obj(list(hclus())))
    kc <- check_exp(check_obj(list(kclus()$cluster)))
    if ('hcl' %in% method & !('e' %in% hc[1])) { 
      enable(selector='a[data-value="shmAll-net-netOp"]')
    } else if ('kmean' %in% method & !('e' %in% kc[1])) { 
      enable(selector='a[data-value="shmAll-net-netOp"]')
    }
  })
  onBookmark(function(state) { state })

})

}
