options(shiny.maxRequestSize=7*1024^3, stringsAsFactors=FALSE) 

# Import internal functions.
#filter_data <- get('filter_data', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#svg_df <- get('svg_df', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#nod_lin <- get('nod_lin', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#grob_list <- get('grob_list', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#col_bar <- get('col_bar', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#lay_shm <- get('lay_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#adj_mod <- get('adj_mod', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#matrix_hm <- get('matrix_hm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#network <- get('network', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# source('~/tissue_specific_gene/function/fun.R')

# Right before submit the package the following functions will be deleted, and they will be imported as above. They are listed here now for the convenience of functionality development.

sort_gen_con <- function(ID.sel, na.all, con.all, by='gene') {

  # Sort vector of letter and number mixture.
  sort_mix <- function(vec) {

    w <- is.na(as.numeric(gsub('\\D', '', vec))); let <- vec[w]; num <- vec[!w]
    let.num <- as.numeric(gsub('\\D', '', vec)); vec[!w] <- num[order(let.num[!w])]
    vec[w] <- sort(let); return(vec)
        
  }  

  pat.all <- paste0(paste0('^(', paste0(ID.sel, collapse='|'), ')'), '_', paste0('(', paste0(con.all, collapse='|'), ')$'))
  svg.idx <- unique(gsub('(.*)(_\\d+)$', '\\2', con.all))
  
  if (by=='gene') {

    # Sort conditions under each gene.
    con.pat1 <- paste0('_(', paste0(con.all, collapse='|'), ')$') # Avoid using '.*' as much as possible.
    na.sort <- NULL; for (i in sort_mix(ID.sel)) {
        
      na0 <- na.all[grepl(paste0('^', i, con.pat1), na.all)]
      if (length(na0)==0) next; con1 <- gsub(pat.all, '\\2', na0)
      # Sort conditions after isolating svg.idx.
      for (j in svg.idx) {

        con2 <- con1[grepl(paste0(j, '$'), con1)]
        con0 <- gsub('(.*)(_\\d+)$', '\\1', con2)
        na.sort <- c(na.sort, paste0(i, '_', sort_mix(con0), j))

      }

    }

  } else if (by=='con') {

    # Sort conditions and genes.
    gen.pat1 <- paste0('^(', paste0(ID.sel, collapse='|'), ')_') # Avoid using '.*' as much as possible.
    # Sort conditions after isolating svg.idx.
    con.all1 <- NULL; for (j in svg.idx) {

      con1 <- con.all[grepl(paste0(j, '$'), con.all)]
      con0 <- gsub('(.*)(_\\d+)$', '\\1', con1)
      con.all1 <- c(con.all1, paste0(sort_mix(con0), j))

    }
    na.sort <- NULL; for (i in sort_mix(con.all1)) {
      
      na0 <- na.all[grepl(paste0(gen.pat1, i, '$'), na.all)]
      if (length(na0)==0) next; gen1 <- gsub(pat.all, '\\1', na0)
      na.sort <- c(na.sort, paste0(sort_mix(gen1), '_', i))

    }

  } else if (by=='none') {
    
    # Sort conditions under each gene.
    con.pat1 <- paste0('_(', paste0(con.all, collapse='|'), ')$') # Avoid using '.*' as much as possible.
    na.sort <- NULL; for (i in ID.sel) {
        
      na0 <- na.all[grepl(paste0('^', i, con.pat1), na.all)]
      if (length(na0)==0) next; na.sort <- c(na.sort, na0)

    }
        
  }; return(na.sort)

}

matrix_hm <- function(ID, data, scale='no', col=c('purple', 'yellow', 'blue'), main=NULL, title.size=10, cexCol=1, cexRow=1, angleCol=45, angleRow=45, sep.color="black", sep.width=0.02, static=TRUE, margin=c(10, 10), arg.lis1=list(), arg.lis2=list()) {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data)) 
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    gene <- na[!idx]; colnames(gene) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) { gene <- assay(data) }
  mod <- as.matrix(gene)
 
  if (static==TRUE) {

    tmp <- tempdir(check=TRUE); pa <- paste0(tmp, '/delete_hm.png')
    png(pa); hm <- heatmap.2(x=mod, scale=scale, main=main, trace="none"); dev.off()
    do.call(file.remove, list(pa))
    # Select the row of target gene.  
    idx <- which(rev(colnames(hm$carpet) %in% ID))
    # If colour codes are more than 500, the colour key is blank.
    lis1 <- c(arg.lis1, list(x=mod, scale=scale, main=main, margin=margin, col=colorRampPalette(col)(500), rowsep=c(idx-1, idx), cexCol=cexCol, cexRow=cexRow, srtRow=angleRow, srtCol=angleCol, dendrogram='both', sepcolor=sep.color, sepwidth=c(sep.width, sep.width), key=TRUE, trace="none", density.info="none", Rowv=TRUE, Colv=TRUE))
    do.call(heatmap.2, lis1)

  } else if (static==FALSE) {

     x <- x1 <- x2 <- y <- y1 <- y2 <- xend <- yend <- value <- NULL 
     dd.gen <- as.dendrogram(hclust(dist(mod))); dd.sam <- as.dendrogram(hclust(dist(t(mod))))
     d.sam <- dendro_data(dd.sam); d.gen <- dendro_data(dd.gen)

     g.dengra <- function(df) {

       ggplot()+geom_segment(data=df, aes(x=x, y=y, xend=xend, yend=yend))+labs(x="", y="")+theme_minimal()+theme(axis.text= element_blank(), axis.ticks=element_blank(), panel.grid=element_blank())

     }

     p.gen <- g.dengra(d.gen$segments)+coord_flip(); p.sam <- g.dengra(d.sam$segments)
     gen.ord <- order.dendrogram(dd.gen); sam.ord <- order.dendrogram(dd.sam); mod.cl <- mod[gen.ord, sam.ord]
     if (scale=="column") mod.cl <- scale(mod.cl); if (scale=="row") mod.cl <- t(scale(t(mod.cl)))
     mod.cl <- data.frame(mod.cl); mod.cl$gene <- rownames(mod.cl)
     mod.m <- reshape2::melt(mod.cl, id.vars='gene', measure.vars=colnames(mod)); colnames(mod.m) <- c('gene', 'sample', 'value')
     # Use "factor" to re-order rows and columns as specified in dendrograms. 
     mod.m$gene <- factor(mod.m$gene, levels=rownames(mod.cl)); mod.m$sample <- factor(mod.m$sample, levels=colnames(mod.cl))
     # Plot the re-ordered heatmap.
     lis2 <- c(arg.lis2, list(data=mod.m, mapping=aes(x=sample, y=gene))) 
     g <- do.call(ggplot, lis2)+geom_tile(aes(fill=value), colour="white")+scale_fill_gradient(low=col[1], high=col[2])+theme(axis.text.x=element_text(size=cexRow*10, angle=angleCol), axis.text.y=element_text(size=cexCol*10, angle=angleRow))
     # Label target row/gene.
     g.idx <- which(rownames(mod.cl) %in% ID)
     g <- g+geom_hline(yintercept=c(g.idx-0.5, g.idx+0.5), linetype="solid", color=sep.color, size=sep.width*25)
     ft <- list(family = "sans serif", size=title.size, color='black')
     subplot(p.sam, ggplot(), g, p.gen, nrows=2, shareX=TRUE, shareY=TRUE, margin=0, heights=c(0.2, 0.8), widths=c(0.8, 0.2)) %>% plotly::layout(title=main, font=ft)

   }

}



sub_na <- function(mat, ID, p=0.3, n=NULL, v=NULL) {

  len <- nrow(mat)
  na <- NULL; for (i in ID) {

    if (!is.null(p)) {
  
      vec <- sort(mat[, i]); thr <- vec[len-floor(len*p)+1]
      na0 <- names(vec[vec >= thr]); na <- c(na, na0)

    } else if (!is.null(n)) {

      vec <- sort(mat[, i]); thr <- vec[len-n+1]
      na0 <- names(vec[vec >= thr]); na <- c(na, na0)
    } else if (!is.null(v)) {
  
      vec <- mat[, i]; na0 <- names(vec[vec >= v]); na <- c(na, na0)

    }

  }; na <- unique(na); return(na)

}

submatrix <- function(data, ann=NULL, ID, p=0.3, n=NULL, v=NULL, fun='cor', cor.absolute=FALSE, arg.cor=list(method="pearson"), arg.dist=list(method="euclidean")) {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data)) 
    if (any(duplicated(cna))) stop('Please use function \'aggr_rep\' to aggregate replicates!')
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0; ann <- data[idx]
    data <- na[!idx]; colnames(data) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) { 

    ann <- rowData(data)[ann]; data <- as.data.frame(assay(data)) 
    if (any(duplicated(rownames(data)))) stop('Please use function \'aggr_rep\' to aggregate replicates!')

  }; if (nrow(data)<5) cat('Warning: variables of sample/condition are less than 5! \n')

  na <- NULL; len <- nrow(data)
  if (len>=50000) cat('More than 50,000 rows are detected in data. Computation may take a long time! \n')
  if (sum(c(!is.null(p), !is.null(n), !is.null(v)))>1) return('Please only use one of \'p\', \'n\', \'v\' as the threshold!')

  if (fun=='cor') {

    m <- do.call(cor, c(x=list(t(data)), arg.cor))
    if (cor.absolute==TRUE) { m1 <- m; m <- abs(m) }
    na <- sub_na(mat=m, ID=ID, p=p, n=n, v=v)
    if (cor.absolute==TRUE) m <- m1

  } else if (fun=='dist') { 
    
    m <- -as.matrix(do.call(dist, c(x=list(data), arg.dist)))
    if (!is.null(v)) v <- -v
    na <- sub_na(mat=m, ID=ID, p=p, n=n, v=v); m <- -m

  }; sub.m <- m[na, na] 
  return(list(sub_matrix=cbind(data[na, ], ann[na, , drop=FALSE]), cor_dist=m))

}

adj_mod <- function(data, type='signed', power=if (type=='distance') 1 else 6, arg.adj=list(), TOMType='unsigned', arg.tom=list(), method='complete', minSize=15, arg.cut=list(), dir=NULL) {

  options(stringsAsFactors=FALSE)
  # Get data matrix.
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data)) 
    if (any(duplicated(cna))) stop('Please use function \'aggr_rep\' to aggregate replicates!')
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    data <- na[!idx]; colnames(data) <- cna[!idx]; data <- t(data)

  } else if (is(data, 'SummarizedExperiment')) { data <- t(assay(data)) 

    if (any(duplicated(rownames(data)))) stop('Please use function \'aggr_rep\' to aggregate replicates!')

  }; if (nrow(data)<5) cat('Warning: variables of sample/condition are less than 5! \n')

  if (ncol(data)>10000) cat('More than 10,000 rows are detected in data. Computation may take a long time! \n')
  
  # Compute adjacency matrix.
  arg.adj <- c(list(datExpr=data, power=power, type=type), arg.adj)
  adj <- do.call(adjacency, arg.adj)
  # Compute TOM and hierarchical clustering.
  arg.tom <- c(list(adjMat=adj, TOMType=TOMType), arg.tom)
  tom <- do.call(TOMsimilarity, arg.tom)
  dissTOM <- 1-tom; tree.hclust <- flashClust(d=as.dist(dissTOM), method=method)
  # Cut the tree to get modules.
  cutHeight <- quantile(tree.hclust[['height']], probs=seq(0, 1, 0.05))[19]
  arg.cut1 <- list(dendro=tree.hclust, pamStage=FALSE, cutHeight=cutHeight, distM=dissTOM)
  na.cut1 <- names(arg.cut1); na.cut <- names(arg.cut)
  if ('deepSplit' %in% na.cut) arg.cut <- arg.cut[!(na.cut %in% 'deepSplit')]
  w <- na.cut1 %in% na.cut
  arg.cut1 <- c(arg.cut1[!w], arg.cut)
  mcol <- NULL; for (ds in 2:3) {
         
    min <- as.numeric(minSize)-3*ds; if (min < 5) min <- 5
    arg.cut <- c(list(minClusterSize=min, deepSplit=ds), arg.cut1)
    tree <- do.call(cutreeHybrid, arg.cut)
    mcol <- cbind(mcol, tree$labels); arg.cut <- list()

  }; colnames(mcol) <- as.character(2:3); rownames(mcol) <- colnames(adj) <- rownames(adj) <- colnames(data)
  if (!is.null(dir)) { 

    path <- paste0(dir, "/local_mode_result/"); if (!dir.exists(path)) dir.create(path)
    write.table(adj, paste0(path, "adj.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    write.table(mcol, paste0(path, "mod.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    
  }; return(list(adj=adj, mod=mcol))

}



filter_data <- function(data, pOA=c(0, 0), CV=c(-Inf, Inf), ann=NULL, sam.factor, con.factor,  dir=NULL) {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data))
    if (!identical(cna, colnames(data))) cat('Syntactically valid column names are made! \n')
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    row.meta <- data[idx]; expr <- na[!idx]; colnames(expr) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) {

    expr <- assay(data); col.meta <- as.data.frame(colData(data))
    row.meta <- as.data.frame(rowData(data), stringsAsFactors=FALSE)[, , drop=FALSE]
    # Factors teated by paste0/make.names are vectors.
    if (!is.null(sam.factor) & !is.null(con.factor)) { cna <- paste0(col.meta[, sam.factor], '__', col.meta[, con.factor]) } else if (!is.null(sam.factor) & is.null(con.factor)) { cna <- as.vector(col.meta[, sam.factor]) } else if (is.null(sam.factor) & !is.null(con.factor)) { cna <- as.vector(col.meta[, con.factor]) } else cna <- colnames(expr)
    colnames(expr) <- make.names(cna); if (!identical(cna, make.names(cna))) cat('Syntactically valid column names are made! \n')

  }

  ffun <- filterfun(pOverA(p=pOA[1], A=pOA[2]), cv(CV[1], CV[2]))
  filtered <- genefilter(expr, ffun); expr <- expr[filtered, , drop=FALSE] # Subset one row in a matrix, the result is a numeric vector not a matrix, so drop=FALSE.
  row.meta <- row.meta[filtered, , drop=FALSE]

  if (!is.null(dir)) { 

    path <- paste0(dir, "/local_mode_result/"); if (!dir.exists(path)) dir.create(path)

    if (is(data, 'data.frame')|is(data, 'matrix')) {

      expr1 <- cbind.data.frame(expr, row.meta, stringsAsFactors=FALSE)

    }  else if (is(data, 'SummarizedExperiment')) {
      
      if (ncol(row.meta)>0 & !is.null(ann)) {
        
        expr1 <- cbind.data.frame(expr, row.meta[, ann], stringsAsFactors=FALSE)
        colnames(expr1)[ncol(expr1)] <- ann
      
      } else expr1 <- expr

    }
    write.table(expr1, paste0(path, "processed_data.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
      
  }

  if (is(data, 'data.frame')|is(data, 'matrix')) { return(cbind(expr, row.meta)) } else if (is(data, 'SummarizedExperiment')) {
  
    rownames(col.meta) <- NULL # If row names present in colData(data), if will become column names of assay(data).
    expr <- SummarizedExperiment(assays=list(expr=expr), rowData=row.meta, colData=col.meta); return(expr)

  }

}


norm_data <- function(data, norm.fun='CNF', parameter.list=NULL, data.trans='none') { 
  
  if (is(data, 'data.frame')|is(data, 'matrix')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- colnames(data) 
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    ann <- data[idx]; expr <- na[!idx]; colnames(expr) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) { expr <- SummarizedExperiment::assay(data) }
  
  if (is.null(norm.fun)) norm.fun <- 'no'
  if (!(norm.fun %in% c("CNF", "ESF", "VST", "rlog"))) { 

    if (data.trans=='log2') { 

      v.min <- min(expr); if (v.min<0) expr <- expr-v.min
      expr <- log2(expr+1) 

  } else if (data.trans=='exp2') expr <- 2^expr }
  
  if (min(expr)>=0 & all(round(expr)==expr)) {

    if (norm.fun=='CNF') {

      na <- names(parameter.list); if (!('method' %in% na)|is.null(parameter.list)) {

        parameter.list <- c(list(method='TMM'), parameter.list)

      }
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      y <- DGEList(counts=expr); 
      y <- do.call(calcNormFactors, c(list(object=y), parameter.list))
      expr <- cpm(y, normalized.lib.sizes=TRUE, log=(data.trans=='log2'))

    } else dds <- DESeqDataSetFromMatrix(countData=expr, colData=data.frame(col.dat=colnames(expr)), design=~1) # "design" does not affect "rlog" and "varianceStabilizingTransformation".

    if (norm.fun=='ESF') {

      na <- names(parameter.list); if (!('type' %in% na)|is.null(parameter.list)) { parameter.list <- c(list(type='ratio'), parameter.list) }
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      dds <- do.call(estimateSizeFactors, c(list(object=dds), parameter.list))
      expr <- counts(dds, normalized=TRUE)
      if (data.trans=='log2') expr <- log2(expr+1)

    } else if (norm.fun=='VST') { 
    
      # Apply A Variance Stabilizing Transformation (VST) To The Count Data.
      if (is.null(parameter.list)) parameter.list <- list(fitType='parametric', blind=TRUE)
      na <- names(parameter.list); if (!('fitType' %in% na)) { parameter.list <- c(list(fitType='parametric'), parameter.list) }
      if (!('blind' %in% na)) { parameter.list <- c(list(blind=TRUE), parameter.list) }    
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      # Returns log2-scale data. 
      vsd <- do.call(varianceStabilizingTransformation, c(list(object=dds), parameter.list))
      expr <- SummarizedExperiment::assay(vsd)
      if (data.trans=='exp2') expr <- 2^expr

    } else if (norm.fun=='rlog') { 
  
      if (is.null(parameter.list)) parameter.list <- list(fitType='parametric', blind=TRUE)
      na <- names(parameter.list); if (!('fitType' %in% na)) { parameter.list <- c(list(fitType='parametric'), parameter.list) }
      if (!('blind' %in% na)) { parameter.list <- c(list(blind=TRUE), parameter.list) }    
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      # Apply A 'Regularized Log' Transformation. 
      rld <- do.call(rlog, c(list(object=dds), parameter.list))
      expr <- SummarizedExperiment::assay(rld) 
      if (data.trans=='exp2') expr <- 2^expr  

    }

  } else cat('Nornalisation only applies to data matrix with all non-negative values! \n')

  if (is(data, 'data.frame')|is(data, 'matrix')) { return(cbind(expr, ann)) } else if (is(data, 'SummarizedExperiment')) { SummarizedExperiment::assay(data) <- expr; return(data) }

} 



aggr_rep <- function(data, sam.factor, con.factor, aggr='mean') {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data))
    if (!identical(cna, colnames(data))) cat('Syntactically valid column names are made! \n')
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    ann <- data[idx]; mat <- na[!idx]; fct <- colnames(mat) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) {

    mat <- assay(data); col.meta <- as.data.frame(colData(data))
    # Factors teated by paste0/make.names are vectrs.
    if (!is.null(sam.factor) & !is.null(con.factor)) { fct <- paste0(col.meta[, sam.factor], '__', col.meta[, con.factor]) } else if (!is.null(sam.factor) & is.null(con.factor)) {  fct <- as.vector(col.meta[, sam.factor]) } else if (is.null(sam.factor) & !is.null(con.factor)) { fct <- as.vector(col.meta[, con.factor]) } else fct <- colnames(mat)
    if (!identical(fct, make.names(fct))) cat('Syntactically valid column names are made! \n')
    fct <- colnames(mat) <- make.names(fct) 

  }
  # To keep colnames, "X" should be a character, not a factor.
  if (aggr=='mean') mat <- sapply(X=unique(fct), function(x) rowMeans(mat[, fct==x, drop=FALSE]))
  if (aggr=='median') {
  
    mat <- sapply(X=unique(fct), function(x) Biobase::rowMedians(mat[, fct==x, drop=FALSE]))
    rownames(mat) <- rownames(data)

  }
  
  if (is(data, 'data.frame')|is(data, 'matrix')) { return(cbind(mat, ann)) } else if (is(data, 'SummarizedExperiment')) { 
  
    col.meta <- col.meta[!duplicated(fct), ]; rownames(col.meta) <- NULL
    data <- SummarizedExperiment(assays=list(expr=mat), rowData=rowData(data), colData=col.meta); return(data)

  } 

}


col_bar <- function(geneV, cols, width, bar.title.size=0, bar.value.size=10, mar=c(0.07, 0.01, 0.05, 0.1)) {        

  color_scale <- y <- NULL
  cs.df <- data.frame(color_scale=geneV, y=1)
  cs.g <- ggplot(data=cs.df)+geom_bar(aes(x=color_scale, y=y), fill=cols, orientation='x', stat="identity", width=((max(geneV)-min(geneV))/length(geneV))*width)+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=bar.value.size), plot.margin=margin(t=mar[1], r=mar[2], b=mar[3], l=mar[4], "npc"), panel.grid=element_blank(), panel.background=element_blank(), plot.title= element_text(hjust=0.1, size=bar.title.size))+coord_flip()+labs(title="value", x=NULL, y=NULL)
  if (max(abs(geneV))<=10000 & max(abs(geneV))>=0.0001) cs.g <- cs.g+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0)) else cs.g <- cs.g+scale_x_continuous(expand=c(0,0), labels=function(x) format(x, scientific=TRUE))+scale_y_continuous(expand=c(0,0))+theme(axis.text.y=element_text(angle=45))
  return(cs.g)

}


lay_shm <- function(lay.shm, con, ncol, ID.sel, grob.list, width, height, shiny) {
  
  width <- as.numeric(width); height <- as.numeric(height); ncol <- as.numeric(ncol); con <- unique(con)
  grob.all.na <- names(grob.list)
  if (lay.shm=="gene"|lay.shm=="none") {

    all.cell <- ceiling(length(con)/ncol)*ncol
    cell.idx <- c(seq_len(length(con)), rep(NA, all.cell-length(con)))
    m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
    lay <- NULL; for (i in seq_len(length(ID.sel))) { lay <- rbind(lay, m+(i-1)*length(con)) }
    # Sort conditions under each gene.
    na.sort <- sort_gen_con(ID.sel=ID.sel, na.all=grob.all.na, con.all=con, by=lay.shm)
    grob.list <- grob.list[na.sort]
    if (shiny==TRUE & length(grob.list)>=1) return(grid.arrange(grobs=grob.list, layout_matrix=lay, newpage=TRUE))
       
    g.tr <- lapply(grob.list[seq_len(length(grob.list))], grobTree)
    n.col <- ncol(lay); n.row <- nrow(lay)
    g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc"))

  } else if (lay.shm=="con") {

    all.cell <- ceiling(length(ID.sel)/ncol)*ncol
    cell.idx <- c(seq_len(length(ID.sel)), rep(NA, all.cell-length(ID.sel)))
    m <- matrix(cell.idx, ncol=ncol, byrow=TRUE)
    lay <- NULL; for (i in seq_len(length(con))) { lay <- rbind(lay, m+(i-1)*length(ID.sel)) }
    na.sort <- sort_gen_con(ID.sel=ID.sel, na.all=grob.all.na, con.all=con, by='con')
    grob.list <- grob.list[na.sort]
    if (shiny==TRUE & length(grob.list)>=1) return(grid.arrange(grobs=grob.list, layout_matrix=lay, newpage=TRUE))
    g.tr <- lapply(grob.list, grobTree); g.tr <- g.tr[names(grob.list)]
    n.col <- ncol(lay); n.row <- nrow(lay)
    g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc")) 

  }; return(g.arr)

}

nod_lin <- function(ds, lab, mods, adj, geneID, adj.min) {

  from <- to <- NULL
  idx.m <- mods[, ds]==lab; adj.m <- adj[idx.m, idx.m]; gen.na <- colnames(adj.m) 
  idx.sel <- grep(paste0("^", geneID, "$"), gen.na); gen.na[idx.sel] <- paste0(geneID, "_target")
  colnames(adj.m) <- rownames(adj.m) <- gen.na; idx = adj.m > as.numeric(adj.min)
  link <- data.frame(from=rownames(adj.m)[row(adj.m)[idx]], to=colnames(adj.m)[col(adj.m)[idx]], width=adj.m[idx], stringsAsFactors=FALSE)
  # Should not exclude duplicate rows by "length".
  node.pas <- NULL; for (i in seq_len(nrow(link))) { node.pas <- c(node.pas, paste0(sort(c(link[i, 'from'], link[i, 'to'])), collapse='')) }
  w <- which(duplicated(node.pas)); link <- link[-w, ]
  link1 <- subset(link, from!=to, stringsAsFactors=FALSE); link1 <- link1[order(-link1$width), ]
  node <- data.frame(label=colnames(adj.m), size=colSums(adj.m), stringsAsFactors=FALSE)
  node <- node[order(-node$size), ]; return(list(node=node, link=link1))

}


# Break combined path to a group (g=TRUE) or siblings (g=FALSE).

path_br <- function(node, g=TRUE) {

    na <- xml2::xml_name(node); if (na!='g') {

      d <- xml2::xml_attr(node, 'd') 
      if (grepl('m ', d)) return('Please use absolute coordinates for all paths!')
      if (grepl('Z M', d)) {

        z <- paste0(strsplit(d, 'Z')[[1]], 'Z')
        ids <- paste0(xml2::xml_attr(node, 'id'), '_', seq_along(z))
        # Make node empty.
        xml2::xml_attr(node, 'd') <- NA
        
        # Break the combined path to a group.
        if (g==TRUE) {
        
          # Isolate 'title' node.
          na.chil <- xml_name(xml_children(node))
          w <- which(na.chil=='title')
          if (length(w)>0) { tit <- xml_children(node)[[w]]; xml_remove(xml_children(node)[w], free=FALSE) }

          # Add the empty node to itself as the first child.
          xml_add_child(node, node)
          # Copy the first child for length(z)-1 times.
          node1 <- xml_children(node)[[1]]
          for (j in seq_len(length(z)-1)) { xml_add_child(node, node1) }
          node.chl <- xml_children(node) # Function applies to 'nodeset' recusively. 
          # Set d and id for all childrend of node.
          xml_set_attr(node.chl, 'd', z)
          xml_set_attr(node.chl, 'id', ids)  
          # Name node 'g'.
          xml2::xml_name(node) <- 'g'; xml2::xml_attr(node, 'd') <- NULL
          if (length(w)>0) xml_add_child(node, tit, .where=0)

        } else {

          for (j in seq_along(z)[-1]) {

            # Copy node as its own siblings.
            xml_set_attr(node, 'd', z[j]); xml_set_attr(node, 'id', ids[j])
            xml_add_sibling(node, node, 'after')
            # Change 'd' in node at last.  
            xml_set_attr(node, 'd', z[1]); xml_set_attr(node, 'id', ids[1])

          }

        }

      }
    
    }

  }



# The outline or tissue nodes are checked for combines paths. If combined paths are detected, those outside a group are broken to a group while those inside a group are broken as siblings.  
path_br_all <- function(node.parent) {

  len <- xml_length(node.parent); chdn <- xml_children(node.parent)
  for (i in seq_len(len)) {

    nod0 <- chdn[[i]]; na <- xml_name(nod0)
      if (na=='a') next
      if (na!='g') path_br(nod0, g=TRUE) else {

        nod0.chl <- xml_children(nod0); nas <- xml_name(nod0.chl)
        if ('g' %in% nas) return(paste0('Nested group detected in ', xml_attr(nod0, 'id'), '!'))
        if ('use' %in% nas) return(paste0('use node detected in ', xml_attr(nod0, 'id'), '!'))
        for (j in seq_along(nod0.chl)) {

          nod1 <- nod0.chl[[j]]; if (xml_name(nod1)=='a') next
          path_br(nod1, g=FALSE)

        }

      }

  }

}


# 'a' nodes are not removed.
svg_attr <- function(doc, feature) {

  options(stringsAsFactors=FALSE); name <- NULL
  len <- xml_length(doc); out <- xml_children(doc)[[len-1]]; ply <- xml_children(doc)[[len]]
  # Break combined path to a group or siblings.
  path_br_all(out); path_br_all(ply)


  # If out is not a group, it is assigned an empty node.
  if (xml_name(out)!='g') { xml_add_child(out, 'empty', .where=0); out1 <- xml_children(out)[[1]]; xml_remove(xml_children(out)[[1]], free=FALSE); out <- out1 }
  # If ply is not a group, it is assigned an empty node.
  if (xml_name(ply)!='g') { xml_add_child(ply, 'empty', .where=0); ply1 <- xml_children(ply)[[1]]; xml_remove(xml_children(ply)[[1]], free=FALSE); ply <- ply1 }
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  
  ## Exrtact basic attributes into a data frame.
  idx <- seq_len(length(chdn.out)+length(chdn.ply))
  idx1 <- c(seq_len(length(chdn.out)), seq_len(length(chdn.ply)))
  parent <- c(rep(xml_attr(out, 'id'), length(chdn.out)), rep(xml_attr(ply, 'id'), length(chdn.ply)))
  nas <- c(xml_name(chdn.out), xml_name(chdn.ply))
  ids <- make.names(c(xml_attr(chdn.out, 'id'), xml_attr(chdn.ply, 'id')))
  if (any(duplicated(ids))) return(paste0('Duplicated node ids detected: ', paste0(ids[duplicated(ids)], collapse=' '), '!'))
  title <- c(xml_text(chdn.out), xml_text(chdn.ply)) # Use original names, no 'make.names'. '' after applied to 'make.names' becomes 'X'.
  w <- which(title=='X'|title==''); title[w] <- ids[w]
  dup <- duplicated(title); if (any(dup)) {
   
    tit.dup <- unique(title[dup]) 
    if (length(intersect(tit.dup, unique(feature)))>0) return(paste0('Duplicated title text detected: ', paste0(tit.dup, collapse=' '), '!')) else {

      w <- title %in% tit.dup
      title[w] <- paste0(title[w], seq_len(sum(w)))

    }
  
  }
  # Style inside groups are ignored. 
  sty <- c(xml_attr(chdn.out, 'style'), xml_attr(chdn.ply, 'style'))
  sty[!grepl('fill:', sty)] <- 'none'
  w1 <- grepl(';', sty); st <- sty[w1]; st <- strsplit(st, ';')
  st1 <- NULL; for (i in st) { st1 <- c(st1, i[grepl('fill:', i)]) }; sty[w1] <- st1
  # If only keep part of the string, the pattern should cover everything in the string, e.g. the '.*' on both ends.
  fil.cols <- gsub('.*(fill:)(.*).*', '\\2', sty)
  df.attr <- data.frame(index=idx, index1=idx1, parent=parent, name=nas, id=ids, title=title, fil.cols=fil.cols)
  df.attr <- subset(df.attr, name!='a')    
  return(list(df.attr=df.attr, out=out, ply=ply))

}


svg_df <- function(svg.path, feature) {

  # Make sure the style is correct. If the stroke width is not the same across polygons such as '0.0002px', '0.216px', some stroke outlines cannot be recognised by 'PostScriptTrace'. Then some polygons are missing. Since the ggplot is based on 'stroke' not 'fill'.
  options(stringsAsFactors=FALSE)
  doc <- read_xml(svg.path); spa <- xml_attr(doc, 'space')
  if (!is.na(spa)) if (spa=='preserve') xml_set_attr(doc, 'xml:space', 'default')
  # w.h <- c(xml_attr(doc, 'width'), xml_attr(doc, 'height'))
  # w.h <- as.numeric(gsub("^(\\d+\\.\\d+|\\d+).*", "\\1", w.h))
  # names(w.h) <- c('width', 'height')
  svg.attr <- svg_attr(doc, feature=feature); if (is(svg.attr, 'character')) return(svg.attr)
  df.attr <- svg.attr[['df.attr']]; df.attr$title <- make.names(df.attr$title)
  out <- svg.attr[['out']]; ply <- svg.attr[['ply']]
  # Paths in 'a' node are recognised in .ps.xml file, so all 'a' nodes in or out groups are removed. 
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  chdn.all <- c(chdn.out, chdn.ply)
  for (i in chdn.all) {

    na <- xml_name(i)
    if (na=='a') xml_remove(i, free=FALSE) else if (na=='g') {

      chil <- xml_children(i); for (j in chil) {

        na1 <- xml_name(j); if (na1=='a') xml_remove(j, free=FALSE) else if (na1=='g') return(paste0('Nested group detected in ', xml_attr(i, 'id'), '!'))

      }

    }
 
  }
  # Renew the children after deletion of 'a' nodes.
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  chdn.all <- c(chdn.out, chdn.ply)

  # Get ids and titles for every path, including paths inside groups, except for 'a' nodes.
  tit <- id.all <- NULL; for (i in seq_along(chdn.all)) {

    if (df.attr[i, 'name']=='g') {

     na <- xml_name(xml_children(chdn.all[[i]]))
     tit0 <- rep(df.attr[i, 'title'], xml_length(chdn.all[[i]])-sum(na=='title')); tit0 <- paste0(tit0, '_', seq_along(tit0)); tit <- c(tit, tit0)
     id0 <- rep(df.attr[i, 'id'], xml_length(chdn.all[[i]])-sum(na=='title')); id.all <- c(id.all, id0)
     # If the styles in paths of a group are different with group style, they can lead to messy 'fill' and 'stroke' in '.ps.xml', so they are set NULL. This step is super important.
     xml_set_attr(xml_children(chdn.all[[i]]), 'style', NULL)

    } else if (df.attr[i, 'name']=='use') {

      ref <- paste0('#', df.attr[, 'id'])
      w <- which(ref %in% xml_attr(chdn.all[[i]], 'href'))
      # If reference is inside a group. A group contains no groups, so the use node has 1 shape.
      if (length(w)==0) { tit <- c(tit, df.attr[i, 'title']); id.all <- c(id.all, df.attr[i, 'id']) }
      # If reference is outside a group.
      if (length(w)>0) if (df.attr[w, 'name']=='g') {

        na <- xml_name(xml_children(chdn.all[[w]]))
        tit0 <- rep(df.attr[i, 'title'], xml_length(chdn.all[[w]])-sum(na=='title')); tit0 <- paste0(tit0, '_', seq_along(tit0)); tit <- c(tit, tit0)
        id0 <- rep(df.attr[i, 'id'], xml_length(chdn.all[[w]])-sum(na=='title')); id.all <- c(id.all, id0)

      } else { tit <- c(tit, df.attr[i, 'title']); id.all <- c(id.all, df.attr[i, 'id']) }

    } else { tit <- c(tit, df.attr[i, 'title']); id.all <- c(id.all, df.attr[i, 'id']) }

  }; tis.path <- gsub("_\\d+$", "", tit)  

  style <- 'fill:#46e8e8;fill-opacity:1;stroke:#000000;stroke-width:3;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1' # 'fill' is not necessary. In Inkscape, "group" or move an object adds transforms (relative positions), and this can lead to related polygons uncolored in the spatial heatmaps. Solution: ungroup and regroup to get rid of transforms and get absolute positions.
  # Change 'style' of all polygons.
  xml_set_attr(chdn.out, 'style', style); xml_set_attr(chdn.ply, 'style', style)  
  # xml_set_attr(out, 'style', style); xml_set_attr(ply, 'style', style)  
  # Export internal SVG.
  tmp <- tempdir(check=TRUE); svg.inter <- paste0(tmp, '/internal.svg')
  if (grepl("~", svg.inter)) svg.inter <- normalizePath(svg.inter)
  write_xml(doc, file=svg.inter)
  
  # SVG file conversion. 
  rsvg_ps(svg.inter, file=sub("svg$", "ps", svg.inter)) # Only the paths inside canvas of SVG are valid.
  p1 <- sub("svg$", "ps", svg.inter); p2 <- paste0(sub("svg$", "ps", svg.inter), ".xml"); PostScriptTrace(p1, p2) 
  chdn1 <- xml_children(read_xml(p2)) # Use internal svg to get coordinates.
     
  # Detect groups that use relative coordinates ("transform", "matrix" in Inkscape.), which leads to some plygons missing in ".ps.xml" file.
  # EBI SVG, if the outline shapes and tissue shapes are separate, they must be in two layers NOT two groups. Otherwise, 'fill' and 'stroke' in '.ps.xml' can be  messy.
  fil.stk <- xml_attr(chdn1[-length(chdn1)], 'type'); tab <- table(fil.stk)
  w <- which(fil.stk=='fill')%%2==0
  if (any(w) & tab['fill'] > tab['stroke']) { 
 
    # Index of wrong path.
    w1 <- which(w)[1]
    # Wrong path and related group.
    tis.wrg <- paste0(id.all[c(w1-1, w1)], collapse='; ')
    return(paste0("Error detected in '", tis.wrg, "' in SVG image. Please ungroup and regroup the respective group they belong to.")) 

  }
  
  # Get coordinates from '.ps.xml'.
  stroke <- chdn1[which(xml_attr(chdn1, 'type')=='stroke')]
  df <- NULL; for (i in seq_along(stroke)) {

    xy <- xml_children(stroke[[i]])[-1]
    x <- as.numeric(xml_attr(xy, 'x'))
    y <- as.numeric(xml_attr(xy, 'y'))
    df0 <- cbind(tissue=tit[i], data.frame(x=x, y=y), stringsAsFactors=TRUE) # The coordinates should not be factor.
    df <- rbind(df, df0)

  }; fil.cols <- df.attr$fil.cols; names(fil.cols) <- df.attr$title
  w.h <- c(max(abs(df$x)), max(abs(df$y))) 
  names(w.h) <- c('width', 'height')
  lis <- list(df=df, tis.path=sub('_\\d+$', '', tit), fil.cols=fil.cols, w.h=w.h); return(lis)

}


grob_list <- function(gene, con.na=TRUE, geneV, coord, ID, cols, tis.path, tis.trans=NULL, sub.title.size, sam.legend='identical', legend.col, legend.title=NULL, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, legend.title.size=8, line.size=0.2, line.color='grey70', mar.lb=NULL, ...) {

#con.na=TRUE; sam.legend='identical'; legend.title=NULL; legend.ncol=NULL; legend.nrow=NULL; legend.position='bottom'; legend.direction=NULL; legend.key.size=0.5; legend.text.size=8; legend.title.size=8; line.size=0.2; line.color='grey70'; mar.lb=NULL

# save(con.na, tis.trans, sam.legend, legend.col, legend.title, legend.ncol, legend.nrow, legend.position, legend.direction, legend.key.size, legend.text.size, legend.title.size, line.size, line.color, mar.lb, gene, geneV, coord, ID, cols, tis.path, sub.title.size, file='all')



  g_list <- function(con, lgd=FALSE, ...) {

    x <- y <- tissue <- NULL; tis.df <- unique(coord[, 'tissue'])
    if (lgd==FALSE) { 

      g.col <- NULL; con.idx <- grep(paste0("^", con, "$"), cons)
      tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]

      for (i in tis.path) {

        tis.idx <- which(tis.col1 %in% i); if (length(tis.idx)==1) { g.col <- c(g.col, scol1[tis.idx])
        } else if (length(tis.idx)==0) { g.col <- c(g.col, NA) }

      }; names(g.col) <- tis.df
    # Make selected tissues transparent by setting their colours as NA.
    if (!is.null(tis.trans)) for (i in tis.df) { if (sub('_\\d+$', '', i) %in% tis.trans) g.col[i] <- NA }
    
    } else g.col <- rep(NA, length(tis.path))
    names(g.col) <- tis.df # The colors might be internally re-ordered alphabetically during mapping, so give them names to fix the match with tissues. E.g. c('yellow', 'blue') can be re-ordered to c('blue', 'yellow'), which makes tissue mapping wrong. Correct: colours are not re-ordered. The 'tissue' in 'data=coord' are internally re-ordered according to a factor. Therfore, 'tissue' should be a factor with the right order. Otherwise, disordered mapping can happen.

    # Show selected or all samples in legend.
    if (length(sam.legend)==1) if (sam.legend=='identical') sam.legend <- intersect(sam.uni, unique(tis.path)) else if (sam.legend=='all') sam.legend <- unique(tis.path)
    
    if (lgd==FALSE) {
    
      sam.legend <- setdiff(sam.legend, tis.trans) 
      leg.idx <- !duplicated(tis.path) & (tis.path %in% sam.legend)
      # Legends are set for each SHM and then removed in 'ggplotGrob', but a copy with legend is saved separately for later used in video.
      scl.fil <- scale_fill_manual(values=g.col, breaks=as.character(tis.df)[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=legend.title, ncol=legend.ncol, nrow=legend.nrow))
   
    } else { 

      # Select legend key colours if identical samples between SVG and matrix have colors of "none".
      legend.col <- legend.col[sam.legend] 
      if (any(legend.col=='none')) {
       
         n <- sum(legend.col=='none'); col.all <- grDevices::colors()[grep('honeydew|aliceblue|white|gr(a|e)y', grDevices::colors(), invert=TRUE)]
         col.none <- col.all[seq(from=1, to=length(col.all), by=floor(length(col.all)/n))]
         legend.col[legend.col=='none'] <- col.none[seq_len(n)]

       }
       # Map legend colours to tissues.
       sam.legend <- setdiff(sam.legend, tis.trans) 
       leg.idx <- !duplicated(tis.path) & (tis.path %in% sam.legend)
       legend.col <- legend.col[sam.legend] # Exclude transparent tissues. 
       for (i in seq_along(g.col)) {

         g.col0 <- legend.col[sub('_\\d+', '', names(g.col)[i])]
         if (!is.na(g.col0)) g.col[i] <- g.col0

       }; scl.fil <- scale_fill_manual(values=g.col, breaks=as.character(tis.df)[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=legend.title, ncol=legend.ncol, nrow=legend.nrow)) 

    }
    lgd.par <- theme(legend.position=legend.position, legend.direction=legend.direction, legend.background = element_rect(fill=alpha(NA, 0)), legend.key.size=unit(legend.key.size, "npc"), legend.text=element_text(size=legend.text.size), legend.title=element_text(size=legend.title.size), legend.margin=margin(l=0.1, r=0.1, unit='npc'))
    ## Add 'feature' and 'value' to coordinate data frame, since the resulting ggplot object is used in 'ggplotly'. Otherwise, the coordinate data frame is applied to 'ggplot' directly by skipping the following code.
    coord$gene <- k; coord$condition <- con; coord$value <- NA
    ft.pat <- paste0('(', paste0(unique(tis.path), collapse='|'), ')(_\\d+)$')
    coord$feature <- gsub(ft.pat, '\\1', coord$tissue)    
    # Assign values to each tissue.
    col.na <- paste0(coord$feature, '__', coord$condition)
    idx1 <- col.na %in% colnames(gene); df0 <- coord[idx1, ]
    df0$value <- unlist(gene[df0$gene[1], col.na[idx1]])
    coord[idx1, ] <- df0

    # If "data" is not in ggplot(), g$data slot is empty.
    g <- ggplot(data=coord, aes(x=x, y=y, value=value, group=tissue, text=paste0('feature: ', feature, '\n', 'value: ', value)), ...)+geom_polygon(aes(fill=tissue), color=line.color, size=line.size, linetype='solid')+scl.fil+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), axis.title.x=element_text(size=16, face="bold"), plot.title=element_text(hjust=0.5, size=sub.title.size))+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))+lgd.par
    if (is.null(mar.lb)) g <- g+theme(plot.margin=margin(0.005, 0.005, 0.005, 0.005, "npc")) else g <- g+theme(plot.margin=margin(mar.lb[2], mar.lb[1], mar.lb[2], mar.lb[1], "npc"))
    if (con.na==FALSE) g.tit <- ggtitle(k) else g.tit <- ggtitle(paste0(k, "_", con)); g <- g+g.tit

    if (lgd==TRUE) {

      g <- g+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.005, 0.005, 0.2, 0, "npc"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=15, face="bold"))+lgd.par+ggtitle('Legend')

    }; return(g)

  }
  # Map colours to samples according to expression level.
  cname <- colnames(gene); form <- grep('__', cname) # Only take the column names with "__".
  cons <- gsub("(.*)(__)(.*)", "\\3", cname[form]); con.uni <- unique(cons)
  sam.uni <- unique(gsub("(.*)(__)(.*)", "\\1", cname)); tis.trans <- make.names(tis.trans)
  grob.na <- grob.lis <- g.lis.all <- NULL; for (k in ID) {

    scol <- NULL; for (i in gene[k, ]) { 
      ab <- abs(i-geneV); col.ind <- which(ab==min(ab))[1]; scol <- c(scol, cols[col.ind])
    }

    idx <- grep("__", cname); c.na <- cname[idx]
    tis.col <- gsub("(.*)(__)(.*)", "\\1", c.na) 
    tis.col.uni <- unique(tis.col); tis.path.uni <- unique(tis.path)
    tis.tar <- tis.col.uni[tis.col.uni %in% tis.path.uni]
    if (length(tis.tar)==0) return(NULL)
    c.na1 <- c.na[grepl(paste0('^(', paste0(tis.tar,collapse='|'), ')__'), c.na)]
    # Only conditions paired with valid tissues (have matching samples in data) are used. 
    con.vld <- gsub("(.*)(__)(.*)", "\\3", c.na1); con.vld.uni <- unique(con.vld)
    g.lis <- NULL; grob.na0 <- paste0(k, "_", con.vld.uni); g.lis <- lapply(con.vld.uni, g_list)
    # Repress popups by saving it to a png file, then delete it.
    tmp <- tempfile()
    png(tmp); grob <- lapply(g.lis, function(x) { x <- x+theme(legend.position="none"); ggplotGrob(x) })
    dev.off(); if (file.exists(tmp)) do.call(file.remove, list(tmp))
    names(g.lis) <- names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob); g.lis.all <- c(g.lis.all, g.lis)

  }; g.lgd <- g_list(con=NULL, lgd=TRUE)
  return(list(grob.lis=grob.lis, g.lgd=g.lgd, g.lis.all=g.lis.all))

}

# Separate SHMs of grobs and ggplot. Different SHMs of same 'gene_condition' are indexed with suffixed of '_1', '_2', ...
grob_gg <- function(gs) {
 
  # One SVG and multiple SVGs are treated same way. The list of gs slots are named with SVG file names.
  grob.all <- gg.all <- lgd.all <- NULL; idx <- seq_along(gs)
  for (i in idx) {
    
    grob0 <- gs[[i]][['grob.lis']]; names(grob0) <- paste0(names(grob0), '_', i)
    grob.all <- c(grob.all, grob0); gg0 <- gs[[i]][['g.lis.all']]
    names(gg0) <- names(grob0); gg.all <- c(gg.all, gg0)
    lgd.all <- c(lgd.all, list(gs[[i]][['g.lgd']]))

  }; names(lgd.all) <- names(gs)
  return(list(grob=grob.all, gg=gg.all, lgd.all=lgd.all))

}

# Subset data matrix by correlation or distance measure.
submatrix <- function(data, ann=NULL, ID, p=0.3, n=NULL, v=NULL, fun='cor', cor.absolute=FALSE, arg.cor=list(method="pearson"), arg.dist=list(method="euclidean")) {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data)) 
    if (any(duplicated(cna))) stop('Please use function \'aggr_rep\' to aggregate replicates!')
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0; ann <- data[idx]
    data <- na[!idx]; colnames(data) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) { 

    ann <- rowData(data)[ann]; data <- as.data.frame(assay(data)) 
    if (any(duplicated(rownames(data)))) stop('Please use function \'aggr_rep\' to aggregate replicates!')

  }; if (nrow(data)<5) cat('Warning: variables of sample/condition are less than 5! \n')

  na <- NULL; len <- nrow(data)
  if (len>=50000) cat('More than 50,000 rows are detected in data. Computation may take a long time! \n')
  if (sum(c(!is.null(p), !is.null(n), !is.null(v)))>1) return('Please only use one of \'p\', \'n\', \'v\' as the threshold!')
   
  # Function to extract nearest genes.
  sub_na <- function() {

    na <- NULL; for (i in ID) {

      if (!is.null(p)) {
  
        vec <- sort(m[, i]); thr <- vec[len-floor(len*p)+1]
        na0 <- names(vec[vec >= thr]); na <- c(na, na0)

      } else if (!is.null(n)) {

        vec <- sort(m[, i]); thr <- vec[len-n+1]
        na0 <- names(vec[vec >= thr]); na <- c(na, na0)

      } else if (!is.null(v)) {
  
        vec <- m[, i]; na0 <- names(vec[vec >= v]); na <- c(na, na0)

      }

    }; na <- unique(na); return(na)

  }

  if (fun=='cor') {

    m <- do.call(cor, c(x=list(t(data)), arg.cor))
    if (cor.absolute==TRUE) { m1 <- m; m <- abs(m) }
    na <- sub_na()
    if (cor.absolute==TRUE) m <- m1

  } else if (fun=='dist') { 
    
    m <- -as.matrix(do.call(dist, c(x=list(data), arg.dist)))
    if (!is.null(v)) v <- -v
    na <- sub_na(); m <- -m

  }; sub.m <- m[na, na] 
  return(list(sub_matrix=cbind(data[na, ], ann[na, , drop=FALSE]), cor_dist=m))

}

# Adjust legend key size and rows in ggplot.
gg_lgd <- function(gg.all, size.key=NULL, size.text.key=8, angle.text.key=NULL, position.text.key=NULL, sub.title.size=NULL, row=NULL, label=FALSE, label.size=3, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, sam.dat, tis.trans=NULL) {

  # Function to remove feature labels. 
  rm_label <- function(g) {
        
    g.layer <- g$layer; if (length(g.layer)==1) return(g) 
    for (k in rev(seq_along(g.layer))) {

      na.lay <- unique(names(as.list(g.layer[[k]])$geom_params))
      if (all(c('check_overlap', 'angle', 'size') %in% na.lay)) g$layers[[k]] <- NULL

    }; return(g)

  }
  for (i in seq_along(gg.all)) {
  
    g <- gg.all[[i]] 
    if (!is.null(size.key)) g <- g+theme(legend.key.size=unit(size.key, "npc"), legend.text=element_text(size=ifelse(is.null(size.text.key), 8*size.key*33, size.text.key)))
    if (!is.null(sub.title.size)) g <- g+theme(plot.title=element_text(hjust=0.5, size=sub.title.size))
    if (!is.null(row)|label==TRUE|opacity!=1|!is.null(angle.text.key)|!is.null(position.text.key)) {

      lay.dat <- layer_data(g)
      dat <- g$data; g.col <- lay.dat$fill
      names(g.col) <- dat$tissue; df.tis <- as.vector(dat$tissue)
      g.col <- g.col[!duplicated(names(g.col))]; tis.path <- dat$feature
      sam.legend <- intersect(unique(sam.dat), unique(tis.path))
      sam.legend <- setdiff(sam.legend, tis.trans) 
      leg.idx <- !duplicated(tis.path) & (tis.path %in% sam.legend)
      df.tar <- df.tis[leg.idx]; path.tar <- tis.path[leg.idx]
      if (opacity!=1) g.col <- alpha(g.col, opacity)
      if (key==TRUE) gde <- guide_legend(title=NULL, nrow=row, label.theme=element_text(angle=angle.text.key, size=g$theme$legend.text$size), label.position=position.text.key)
      if (key==FALSE) gde <- FALSE
      if (!is.null(row)|opacity!=1|key==FALSE|!is.null(angle.text.key)|!is.null(position.text.key)) g <- g+scale_fill_manual(values=g.col, breaks=df.tar, labels=path.tar, guide=gde)
      if (label==TRUE) {

        dat$x0 <- dat$y0 <- dat$label <- NA
        lab.idx <- dat$feature %in% path.tar
        dat1 <- dat[lab.idx, ]; dat1$label <- dat1$feature
        df.lab <- data.frame() 
        for (j in unique(dat1$tissue)) {

          df0 <- subset(dat1, tissue==j)
          x <- mean(df0$x); y <- mean(df0$y)
          df0$x0 <- x; df0$y0 <- y
          df.lab <- rbind(df.lab, df0)
       
         } 
         g <- rm_label(g)+geom_text(data=df.lab, aes(label=label, x=x0, y=y0), check_overlap=TRUE, size=label.size, angle=label.angle, hjust=hjust, vjust=vjust)

      }; gg.all[[i]] <- rm_label(g)

    }; if (label==FALSE) { gg.all[[i]] <- rm_label(g) }

  }; return(gg.all)

}


# Prepare interactive SHMs in html.
html_ly <- function(gg, cs.g, tis.trans, sam.uni, anm.width, anm.height, selfcontained=FALSE, out.dir) {

  gg.na <- names(gg); cs.lis <- gg2list(cs.g, tooltip='color_scale')
  cs.lis$layout$title$text <- NULL 
  csly <- as_widget(cs.lis, tooltip='color_scale') 
  dir <- paste0(normalizePath(out.dir), '/html_shm')
  if (!dir.exists(dir)) dir.create(dir)
  rd1 <- '1. Double click the "html" files to display the interactive spatial heatmaps in a web browser.'
  rd2 <- '2. All files in the "lib" folder are required to display the spatial heatmaps, so the "html" files cannot work independently.'
  writeLines(text=c(rd1, rd2), con=paste0(dir, '/README.txt'))
  for (i in gg.na) {
    
    na.hl <- paste0(i, '.html')
    cat('Preparing', paste0("'", na.hl, "'"), '... \n')
    g <- gg[[i]]; lay.dat <- layer_data(g)
    dat <- g$data; g.col <- lay.dat$fill; names(g.col) <- dat$tissue
    g.col <- g.col[!duplicated(names(g.col))]; tis.path <- dat$feature
    sam.legend <- intersect(sam.uni, unique(tis.path))
    sam.legend <- setdiff(sam.legend, tis.trans) 
    leg.idx <- !duplicated(tis.path) & (tis.path %in% sam.legend)
    tis.show <- as.vector(dat$tissue)[leg.idx]
    tis.show1 <- tis.path[leg.idx]
    g2l <- gg2list(g, tooltip="text")
    cat('Preparing legend for', paste0("'", na.hl, "'"), '... \n')
    for (i in seq_along(g2l$data)) {
     
      idx.tis <- tis.show %in% g2l$data[[i]]$name
      if (any(idx.tis)) g2l$data[[i]]$name <- tis.show1[idx.tis] else g2l$data[[i]]$showlegend <- FALSE

    }; ggly <- as_widget(g2l)
    subly <- subplot(csly, ggly, nrows=1, shareX=F, shareY=F, margin=0, widths=c(0.05, 0.95))
    subly$width <- anm.width; subly$height <- anm.height
    saveWidget(subly, na.hl, selfcontained=selfcontained, libdir="lib") 
    system(paste0('mv ', na.hl, ' ', dir))

  }; system(paste0('rm -fr ', dir, '/lib')) # Use '-fr' to avoid errors if no target directory.
  if (dir.exists('lib')) system(paste0('mv lib ', dir))

}

# Make videos.
video <- function(gg, cs.g, sam.uni, tis.trans, sub.title.size=NULL, bar.value.size=NULL, lgd.key.size=0.02, lgd.text.size=8, angle.text.key=NULL, position.text.key=NULL, lgd.row=2, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, width=0.92, height=0.99, video.dim='640x480', res=500, interval=1, framerate=1, out.dir) {

  # Test if "av" works.
  test <- function() {
    av_capture_graphics(expr=for (i in seq_along(1:2)) plot(i), output=paste0(tempdir(check=TRUE), '/tmp.mp4'))
  }; try(test())
  ffm <- tryCatch({ test() }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
  if (grepl('error|warning', ffm)) return()

  if (!is.null(bar.value.size)) cs.g <- cs.g+theme(axis.text.y=element_text(size=bar.value.size))
  na <- names(gg)
  cat('Video: adjust legend size/rows... \n')
  gg1 <- gg_lgd(gg.all=gg, size.key=lgd.key.size, size.text.key=lgd.text.size, angle.text.key=angle.text.key, position.text.key=position.text.key, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, sub.title.size=sub.title.size, row=lgd.row, sam.dat=sam.uni, tis.trans=tis.trans)
  lay <- rbind(c(NA, NA), c(1, 2), c(NA, NA))
  cat('Saving video... \n')
  res.r=res/144; w.h <- round(as.numeric(strsplit(video.dim, 'x')[[1]])*res.r)
  if (w.h[1] %% 2!=0) w.h[1] <- w.h[1]+1
  if (w.h[2] %% 2!=0) w.h[2] <- w.h[2]+1
  av_capture_graphics(expr=for (i in na) { print(grid.arrange(cs.g, gg1[[i]],widths=unit(c(0.08, width), 'npc'), 
  heights=unit(c(0.05, height, 0.05), 'npc'), layout_matrix=lay)) }, 
  output=paste0(normalizePath(out.dir), "/shm.mp4"), width=w.h[1], height=w.h[2], res=res, vfilter=paste0('framerate=fps=', framerate))

}


library(SummarizedExperiment); library(shiny); library(shinydashboard); library(grImport); library(rsvg); library(ggplot2); library(DT); library(gridExtra); library(ggdendro); library(WGCNA); library(grid); library(xml2); library(plotly); library(data.table); library(genefilter); library(flashClust); library(visNetwork); library(reshape2); library(igraph); library(animation); library(av); library(shinyWidgets)

# Import input matrix.
fread.df <- function(input, isRowGene, header, sep, fill, rep.aggr='mean', check.names=FALSE) {
        
  df0 <- fread(input=input, header=header, sep=sep, fill=fill)
  cna <- make.names(colnames(df0))
  if (cna[1]=='V1') cna <- cna[-1] else cna <- cna[-ncol(df0)] 
  df1 <- as.data.frame(df0); rownames(df1) <- df1[, 1]
  # Subsetting identical column names in a matrix will not trigger appending numbers.
  df1 <- as.matrix(df1[, -1]); colnames(df1) <- cna
  if(isRowGene==FALSE) df1 <- t(df1)
  cna <- colnames(df1); rna <- rownames(df1)
  
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


# enableWGCNAThreads()
shinyServer(function(input, output, session) {

  output$dld.sgl <- downloadHandler(
    filename=function(){"single_aSVG_data.zip" },
    fil.na <- paste0(tempdir(), '/single_aSVG_data.zip'),
    content=function(fil.na){
      zip(fil.na, c('example/arabidopsis_thaliana.root.cross_shm.svg', 'example/expr_arab.txt'))
    }
  )
  output$dld.mul <- downloadHandler(   
    filename=function(){"multiple_aSVG_data.zip" },
    fil.na <- paste0(tempdir(), '/multiple_aSVG_data.zip'),
    content=function(fil.na){
      zip(fil.na, c('example/random_data_multiple_aSVG.txt', 'example/arabidopsis_thaliana.organ_shm1.svg', 'example/arabidopsis_thaliana.organ_shm2.svg'))
    }
  )

  # Instruction.
  output$sum <-renderUI({ includeHTML("file/summary.html") })
  output$input <-renderUI({ includeHTML("file/input.html") })
  output$input1 <-renderUI({ includeHTML("file/input1.html") })
  output$matrix <-renderUI({ includeHTML("file/matrix.html") })
  output$shm.ins <-renderUI({ includeHTML("file/shm.html") })
  output$mhm.ins <-renderUI({ includeHTML("file/mhm.html") })
  output$net.ins <-renderUI({ includeHTML("file/net.html") })
  # Acknowledgement.
  output$ack <-renderUI({ includeHTML("file/acknowledgement.html") })

  # Filter parameters.
  fil <- reactiveValues(P=0, A=0, CV1=-Inf, CV2=Inf)

  observe({

    input$fileIn; input$geneInpath
    updateRadioButtons(session, inputId="dimName", label="Step 4: is column or row gene?", 
    inline=TRUE, choices=c("None", "Row", "Column"), selected="None")
    updateSelectInput(session, 'sep', 'Step 5: separator', c("None", "Tab", "Space", "Comma", "Semicolon"), "None")
    updateRadioButtons(session, inputId='log', label='Log/exp scaling:', choices=c("No", "log2", "exp2"), selected="No", inline=TRUE)
    output$linear.scl <- renderUI({
    numericInput(inputId="lin.scl", label="Fold scaling:", value=1, min=-Inf, max=Inf, step=1, width=95)
    })
    updateRadioButtons(session, inputId='cs.v', label='Color scale based on:', choices=c("Selected rows"="sel.gen", "All rows"="w.mat"), selected="sel.gen", inline=TRUE)
    updateNumericInput(session, inputId="height", label="Overall height:", value=400, min=0.1, max=Inf, step=NA)
    updateNumericInput(session, inputId="width", label="Overall width:", value=760, min=0.1, max=Inf, step=NA)
    updateNumericInput(session, inputId="col.n", label="Columns:", value=2, min=1, max=Inf, step=1)

  })

  observe({

    input$fileIn; input$geneInpath; input$log
    updateNumericInput(session, inputId="A", label="Value (A) to exceed:", value=0) 
    updateNumericInput(session, inputId="P", label="Proportion (P) of samples with values >= A:", value=0, min=0, max=1)
    updateNumericInput(session, inputId="CV1", label="Min coefficient of variation (CV1):", value=-10^4)
    updateNumericInput(session, inputId="CV2", label="Max coefficient of variation (CV2):", value=10^4,) 
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
    if (grepl("_Mustroph$|_Prudencio$|_Merkin$|_Cardoso.Moreira$|_Census$", input$fileIn)) {

        incProgress(0.5, detail="Loading matrix. Please wait.")
        if (input$fileIn=="brain_Prudencio") { df.te <- fread.df(input="example/expr_human.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE) }
        if (input$fileIn=="mouse_Merkin") df.te <- fread.df(input="example/expr_mouse.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="chicken_Cardoso.Moreira") df.te <- fread.df(input="example/expr_chicken.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="shoot_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="organ_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="root_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="shoot_root_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="root_roottip_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="map_Census") df.te <- fread.df(input="example/us_population2018.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE); return(df.te)

    }

    if ((input$fileIn=="custom_computed_data"|input$fileIn=="custom_data") & 
    (is.null(input$geneInpath)|input$dimName=="None"|input$sep=="None")) return(NULL)
    if ((input$fileIn=="custom_computed_data"|input$fileIn=="custom_data") & 
    !is.null(input$geneInpath) & input$dimName!="None" & input$sep!="None") {

      incProgress(0.25, detail="Importing matrix. Please wait.")
      geneInpath <- input$geneInpath; if (input$sep=="Tab") sep <- "\t" else if (input$sep=="Space") sep <- " " else if (input$sep=="Comma") sep <- "," else if (input$sep=="Semicolon") sep <- ";"
      df.upl <- fread.df(input=geneInpath$datapath, isRowGene=(input$dimName=='Row'), header=TRUE, sep=sep, fill=TRUE, rep.aggr='mean'); return(df.upl)
   
    } 

    })

  })

  # Transform data.
  geneIn1 <- reactive({

    if (is.null(geneIn0())|is.null(input$lin.scl)) return(NULL)
    validate(need(try(!is.na(input$lin.scl)), 'Fold should be a numeric!'))
    gene2 <- geneIn0()[['gene2']] 
    if (input$lin.scl==1) { if (input$log=='log2') {
   
      g.min <- min(gene2)
      if (g.min<0) gene2 <- gene2-g.min+1; if (g.min==0) gene2 <- gene2+1; gene2 <- log2(gene2)

    }; if (input$log=='exp2') gene2 <- 2^gene2 } else { gene2 <- gene2*input$lin.scl }

    gene3 <- geneIn0()[['gene3']]; gen.rep <- geneIn0()[['gen.rep']]
    return(list(gene2=gene2, gene3=gene3, gen.rep=gen.rep))

  })

  output$col.order <- renderUI({

    if (is.null(geneIn1())) return()
    col.nas <- colnames(geneIn1()[['gene2']])
    dropdownButton(inputId='dropdown', label='Re-order columns', circle=FALSE, icon=NULL, status='primary',
    actionButton("col.cfm", "Confirm", icon=icon("refresh")), 
    selectizeInput(inputId="col.na", label='', choices=col.nas, selected=col.nas, multiple=TRUE, options= list(plugins=list('remove_button', 'drag_drop')))
    )
  
  })

  geneIn <- reactive({
    
    if (is.null(geneIn1())) return(NULL)    
    gene2 <- geneIn1()[['gene2']]; gene3 <- geneIn1()[['gene3']]; input$fil.but
    if (!all(input$col.na %in% colnames(gene2))) return()
    # Input variables in "isolate" will not triger re-excution, but if the whole reactive object is trigered by "input$fil.but" then code inside "isolate" will re-excute.
    isolate({
  
      se <- SummarizedExperiment(assays=list(expr=as.matrix(gene2)), rowData=gene3)
      if (ncol(gene3)>0) ann.col <- colnames(gene3)[1] else ann.col <- NULL
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=c(fil$CV1, fil$CV2), dir=NULL)
      if (nrow(se)==0) { validate(need(try(nrow(se)>0), 'All rows are filtered out!')); return() }

      gene2 <- as.data.frame(assay(se), stringsAsfactors=FALSE); colnames(gene2) <- make.names(colnames(gene2))
      gene3 <- as.data.frame(rowData(se))[, , drop=FALSE]
      
    })
    cat('Preparing data matrix... \n')
    return(list(gene2=gene2[, input$col.na], gene3=gene3))

  })

  output$dt <- renderDataTable({

    if (is.null(geneIn())) return()
    if (((input$fileIn=="custom_computed_data"|input$fileIn=="custom_data") & is.null(geneIn()))|input$fileIn=="none") return(NULL)

    withProgress(message="Data table: ", value = 0, {

      incProgress(0.5, detail="Displaying. Please wait.")
      if (input$fileIn!="none") {

      gene <- geneIn(); gene.dt <- cbind.data.frame(gene[["gene2"]][, , drop=FALSE], gene[["gene3"]][, , drop=FALSE], stringsAsFactors=FALSE) 

   }; cat('Presenting data matrix... \n')
    datatable(gene.dt, selection=list(mode="multiple", target="row", selected=c(1)),
    filter="top", extensions=c('Scroller'), options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
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
    gID$new <- setdiff(gID$geneID, gID$all); gID$all <- c(gID$all, gID$new)
    if (is.null(r.na)) gID$geneID <- "none"
    })
  observeEvent(input$dt_rows_selected, {
 
    if (is.null(input$dt_rows_selected)) return()
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    gID$new <- setdiff(gID$geneID, gID$all); gID$all <- c(gID$all, gID$new)
    
  })


  geneV <- reactive({

    if (is.null(geneIn())) return(NULL)
    if (input$cs.v=="sel.gen" & is.null(input$dt_rows_selected)) return(NULL)
    if (input$fileIn!="none") { if (input$cs.v=="sel.gen" & !is.null(input$dt_rows_selected)) gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
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
    
    if (is.null(input$col.but)) return()
    if(input$col.but==0) color$col <- colorRampPalette(c('purple', 'yellow', 'blue'))(length(geneV()))

  })
  # As long as a button is used, observeEvent should be used. All variables inside 'observeEvent' trigger code evaluation, not only 'eventExpr'.  
  observeEvent(input$col.but, {

    if (is.null(col.sch())) return (NULL)
    if (input$fileIn!="none") { color$col <- colorRampPalette(col.sch())(length(geneV())) }

  })

  shm.bar <- reactive({

    if (is.null(gID$all)) return(NULL)
    if ((grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn) & !is.null(geneIn()))|((input$fileIn=="custom_computed_data"|input$fileIn=="custom_data") & (!is.null(input$svgInpath)|!is.null(input$svgInpath1)) & !is.null(geneIn()))) {

      if (length(color$col=="none")==0|input$color==""|is.null(geneV())) return(NULL)
      # if(input$col.but==0) color$col <- colorRampPalette(c('purple', 'yellow', 'blue'))(length(geneV()))

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

    if (is.null(gID$all)) return(NULL)
    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    r.na <- rownames(geneIn()[["gene2"]]); gen.sel <- r.na[input$dt_rows_selected]
    updateSelectInput(session, "gen.sel", choices=c("None", gen.sel), selected="None")

  })

  svg.path <- reactive({

    if (input$fileIn=="custom_computed_data"|input$fileIn=="custom_data") { 
      
      if (is.null(input$svgInpath1)) svgIn.df <- input$svgInpath else svgIn.df <- input$svgInpath1
        svg.path <- svgIn.df$datapath; svg.na <- svgIn.df$name
        if (!is.null(input$svgInpath1)) validate(need(try(all(grepl('_shm\\d+', svg.na, perl=TRUE))), "Suffixes of aSVGs should be indexed as '_shm1', '_shm2', '_shm3', ..."))
        ord <- order(gsub('.*_(shm.*)$', '\\1', svg.na))
        svg.path <- svg.path[ord]; svg.na <- svg.na[ord]
    
    } else if (input$fileIn=="brain_Prudencio") { svg.path <- "example/homo_sapiens.brain.svg"; svg.na <- "homo_sapiens.brain.svg" } else if (input$fileIn=="mouse_Merkin") { svg.path <- "example/mus_musculus.male.svg"; svg.na <- "mus_musculus.male.svg" } else if (input$fileIn=="chicken_Cardoso.Moreira") { svg.path <- "example/gallus_gallus.svg"; svg.na <- "gallus_gallus.svg" } else if (input$fileIn=="shoot_Mustroph") { svg.path <- "example/arabidopsis_thaliana.shoot_shm.svg"; svg.na <- "arabidopsis_thaliana.shoot_shm.svg" } else if (input$fileIn=="organ_Mustroph") { svg.path <- "example/arabidopsis_thaliana.organ_shm.svg"; svg.na <- "arabidopsis_thaliana.organ_shm.svg" } else if (input$fileIn=="root_Mustroph") { svg.path <- "example/arabidopsis_thaliana.root.cross_shm.svg"; svg.na <- "arabidopsis_thaliana.root.cross_shm.svg" } else if (input$fileIn=="shoot_root_Mustroph") { svg.path <- "example/arabidopsis_thaliana.shoot.root_shm.svg"; svg.na <- "arabidopsis_thaliana.shoot.root_shm.svg" } else if (input$fileIn=="root_roottip_Mustroph") { svg.path <- "example/arabidopsis_thaliana.root.roottip_shm.svg"; svg.na <- "arabidopsis_thaliana.root.roottip_shm.svg" } else if (input$fileIn=="map_Census") { svg.path <- "example/us_map_shm.svg"; svg.na <- "us_map_shm.svg" } else return(NULL)
    cat('Access aSVG path... \n')
    return(list(svg.path=svg.path, svg.na=svg.na))

  })

  sam <- reactive({ 

    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\1", c.na) else return(NULL) 

  })
  
  svg.df <- reactive({ 

    if (((input$fileIn=="custom_computed_data"|input$fileIn=="custom_data") & 
    (!is.null(input$svgInpath)|!is.null(input$svgInpath1)))|(grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn) & is.null(input$svgInpath))) {

      withProgress(message="Tissue heatmap: ", value=0, {
    
        incProgress(0.5, detail="Extracting coordinates. Please wait.") 
          svg.path <- svg.path()[['svg.path']]
          svg.na <- svg.path()[['svg.na']]; svg.df.lis <- NULL
          # Whether a single or multiple SVGs, all are returned in a list.
          for (i in seq_along(svg.na)) {
         
            cat('Coordinate:', svg.na[i], '\n')
            df_tis <- svg_df(svg.path=svg.path[i], feature=sam())
            validate(need(!is.character(df_tis), paste0(svg.na[i], ': ', df_tis)))
            svg.df.lis <- c(svg.df.lis, list(df_tis))
   
          }; names(svg.df.lis) <- svg.na; 
          return(svg.df.lis)

      })

    }

  })


  observe({

    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; svg.df(); input$hide.lgd 
    tis.tran <- NULL; for (i in seq_along(svg.df())) { tis.tran <- c(tis.tran, svg.df()[[i]][['tis.path']]) }
    updateCheckboxGroupInput(session, inputId="tis", label='Select tissues to be transparent:', choices=intersect(unique(sam()), unique(tis.tran)), selected='', inline=TRUE)

  })

  con <- reactive({ 

    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\3", c.na) else return(NULL) 

  })

  # General selected gene/condition pattern.
  pat.con <- reactive({ con.uni <- unique(con()); if (is.null(con.uni)) return(NULL); paste0(con.uni, collapse='|') })
  pat.gen <- reactive({ if (gID$geneID[1]=='none') return(NULL);  paste0(gID$geneID, collapse='|') })
  pat.all <- reactive({ if (is.null(pat.con())|is.null(pat.gen())) return(NULL); paste0('(', pat.gen(), ')_(', pat.con(), ')') })

  grob <- reactiveValues(all=NULL, gg.all=NULL, lgd.all=NULL); observeEvent(input$fileIn, { grob$all <- grob$gg.all <- grob$lgd.all <- NULL })
  gs.new <- reactive({ 

    if (is.null(svg.df())|is.null(geneIn())|is.null(gID$new)|length(gID$new)==0|is.null(gID$all)|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL); 
    # Avoid repetitive computation.  
    pat.new <- paste0('^', gID$new, '_(', pat.con(), ')_\\d+$')
    if (any(grepl(pat.new, names(grob$all)))) return()
    withProgress(message="Tissue heatmap: ", value=0, {
 
      incProgress(0.25, detail="preparing data.")
      if (input$cs.v=="sel.gen") gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
      if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]]
      svg.df.lis <- svg.df(); grob.lis.all <- w.h.all <- NULL
      # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
      for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      svg.na <- names(svg.df.lis)
      for (i in seq_along(svg.df.lis)) {

        svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]; w.h <- svg.df[['w.h']]
        tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']]
        if (input$pre.scale=='Y') mar <- (1-w.h/w.h.max*0.99)/2 else mar <- NULL
        cat('New grob/ggplot:', gID$new, ' \n')
        grob.lis <- grob_list(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=gID$new, legend.col=fil.cols, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=18, mar.lb=mar, legend.nrow=2, legend.key.size=0.04) # Only gID$new is used.
        validate(need(!is.null(grob.lis), paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')))
        grob.lis.all <- c(grob.lis.all, list(grob.lis))

      }; names(grob.lis.all) <- svg.na
      return(grob.lis.all)

    })

  })

  # Extension of 'observeEvent': any of 'input$log; input$tis; input$col.but; input$cs.v' causes evaluation of all code. 
  # input$tis as an argument in "grob_list" will not cause evaluation of all code, thus it is listed here.
  # Use "observeEvent" to replace "observe" and list events (input$log, input$tis, ...), since if the events are in "observe", every time a new gene is clicked, "input$dt_rows_selected" causes the evaluation of all code in "observe", and the evaluation is duplicated with "gs.new".
  col.reorder <- reactiveValues(col.re='Y')
  observeEvent(input$col.na, { if (input$col.cfm>0) col.reorder$col.re <- 'N' })
  observeEvent(list(log=input$log, tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, pre.scale=input$pre.scale, col.cfm=input$col.cfm), {
    
    grob$all <- grob$gg.all <- grob$lgd.all <- NULL; gs.all <- reactive({ 

      if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
      withProgress(message="Spatial heatmap: ", value=0, {
        incProgress(0.25, detail="preparing data.")
        if (input$cs.v=="sel.gen") gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
        if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]]

      svg.df.lis <- svg.df(); grob.lis.all <- w.h.all <- NULL
      # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
      for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }
      # A set of SHMs are made for each SVG, and all sets of SHMs are placed in a list.
      for (i in seq_along(svg.df.lis)) {

        svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]
        tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']]; w.h <- svg.df[['w.h']]
        if (input$pre.scale=='Y') mar <- (1-w.h/w.h.max*0.99)/2 else mar <- NULL
        cat('All grob/ggplot:', gID$all, ' \n')
        svg.na <- names(svg.df.lis)
        grob.lis <- grob_list(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=gID$all, legend.col=fil.cols, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=18, mar.lb=mar, legend.nrow=2, legend.key.size=0.04) # All gene IDs are used.
        validate(need(!is.null(grob.lis), paste0(svg.na[i], ': no spatial features that have matching sample identifiers in data are detected!')))
        grob.lis.all <- c(grob.lis.all, list(grob.lis))

      }; names(grob.lis.all) <- svg.na
      return(grob.lis.all)

      })

    }); grob.gg.lis <- grob_gg(gs=gs.all())
    grob$all <- grob.gg.lis[['grob']]; grob$gg.all <- grob.gg.lis[['gg']]; grob$lgd.all <- grob.gg.lis[['lgd.all']]

  })

  # when 'color <- reactiveValues(col="none")', upon the app is launched, 'gs.new' is evaluated for 3 time. In the 1st time, 'gID$new'/'gID$all' are NULL, so 'gs.new' is NULL. In the 2nd time, 'color$col[1]=='none'' is TRUE, so NULL is returned to 'gs.new', but 'gID$new'/'gID$all' are 'HRE2'. In the third time, 'color$col[1]=='none'' is FALSE, so 'gs.new' is not NULL, but 'gID$new' is still 'HRE2', so it does not triger evaluation of 'observeEvent' and hence SHMs and legend plot are not returned upon being launched. The solution is to assign colors to 'color$col' in 'observe' upon being launched so that in the 2nd time 'gs.new' is not NULL, and no 3rd time.
  observeEvent(gID$new, { 

    if (is.null(svg.df())|is.null(gID$new)|length(gID$new)==0|is.null(gID$all)|is.null(gs.new())) return(NULL)
    gs.new <- gs.new(); grob.gg.lis <- grob_gg(gs=gs.new())
    grob.all <- c(grob$all, grob.gg.lis[['grob']]); grob$all <- grob.all[unique(names(grob.all))] 
    gg.all <- c(grob$gg.all, grob.gg.lis[['gg']]); grob$gg.all <- gg.all[unique(names(gg.all))]
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
  observeEvent(list(size=input$lgd.key.size, lgd.row=input$lgd.row, tis.trans=input$tis, lgd.label=input$lgd.label, lgd.lab.size=input$lgd.lab.size), {

    lgd.key.size <- input$lgd.key.size; lgd.row <- input$lgd.row
    lgd.label <- input$lgd.label; label.size <- input$lgd.lab.size
    if (is.null(grob$lgd.all)|is.null(lgd.key.size)|is.null(lgd.row)|is.null(lgd.label)) return()
    cat('Adjust legend size/rows... \n')
    grob$lgd.all <- gg_lgd(gg.all=grob$lgd.all, size.key=lgd.key.size, size.text.key=NULL, row=lgd.row, sam.dat=sam(), tis.trans=input$tis, position.text.key='right', label=(lgd.label=='Y'), label.size=label.size)
  
  })
  observeEvent(input$col.cfm, { col.reorder$col.re <- 'Y' })
  # In "observe" and "observeEvent", if one code return (NULL), then all the following code stops. If one code changes, all the code renews.
  observe({
  
    if (is.null(geneIn())|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)
    output$shm <- renderPlot(width=as.numeric(input$width)/2*as.numeric(input$col.n), height=as.numeric(input$height)*length(input$dt_rows_selected), {

    if (col.reorder$col.re=='N') return()
    if (is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)
    if (length(color$col=="none")==0|input$color=="") return(NULL)
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    grob.na <- names(grob$all)
    # Select target grobs.
    # Use definite patterns and avoid using '.*' as much as possible. Try to as specific as possible.
    pat.all <- paste0('^', pat.all(), '(_\\d+$)')
    grob.lis.p <- grob$all[grepl(pat.all, grob.na)] # grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', names(grob.lis.p))); if (length(con)==0) return()
    cat('Plotting spatial heatmaps... \n')
    lay <- input$gen.con; ID <- gID$geneID; ncol <- input$col.n
    shm <- lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, width=input$width, height=input$height, shiny=TRUE)
    if (input$ext!='NA') {
      
      validate(need(try(input$res>0), 'Resolution should be a positive numeric!'))
      validate(need(try(input$lgd.w>=0 & input$lgd.w <1), 'Legend width should be between 0 to 1!'))
      validate(need(try(input$lgd.ratio>0), 'Legend aspect ratio should be a positive numeric!'))
      cs.grob <- ggplotGrob(shm.bar())
      cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc"))
      # Legend size in downloaded SHM is reduced.
      lgd.lis <- grob$lgd.all; lgd.lis <- gg_lgd(gg.all=lgd.lis, sam.dat=sam(), tis.trans=input$tis, label=FALSE)
      lgd.lis <- gg_lgd(gg.all=lgd.lis, size.key=input$lgd.key.size*0.5, size.text.key=NULL, label.size=input$lgd.lab.size, row=input$lgd.row, sam.dat=sam(), tis.trans=input$tis, position.text.key='right', label=(input$lgd.label=='Y'))
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

 output$shm.ui <- renderUI({

      box(title="Spatial Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, width=ifelse(input$hide.lgd=='N', 9, 12), height=NULL, 
      tabBox(title="", width=12, id='shm_all', selected='shm1', side='right', 
      tabPanel(title='Video', value='shm3', 
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '15%', '1%', '7%', '1%', '7%', '1%', '13%'), '', 
      radioButtons(inputId="vdo.but", label="Show video:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE), '',
      numericInput(inputId='vdo.itvl', label='Transition time (s):', value=1, min=0.1, max=Inf, step=1, width=270), '',
      numericInput(inputId='vdo.height', label='Height:', value=0.99, min=0.1, max=0.99, step=0.1, width=270), '',
      numericInput(inputId='vdo.width', label='Width:', value=0.92, min=0.1, max=0.92, step=0.1, width=270), '',
      numericInput(inputId='vdo.res', label='Resolution (dpi):', value=400, min=1, max=1000, step=5, width=270)
      )), uiOutput('video.dim'), textOutput('tran.vdo'), htmlOutput('ffm'), 
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar3"), uiOutput('video'), ""))),
      
      tabPanel(title='Animation', value='shm2', 
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '15%', '1%', '10%', '1%', '10%', '1%', '13%'), '', 
      radioButtons(inputId="ggly.but", label="Show animation:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE), '',
      numericInput(inputId='t', label='Transition time (s):', value=2, min=0.1, max=Inf, step=NA, width=270), '',
      uiOutput('anm.h'), '', uiOutput('anm.w'), '', uiOutput('dld.anm.but')
      )), textOutput('tran'), uiOutput('sld.fm'),
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar2"), htmlOutput("ggly"), ""))),
      
      tabPanel(title="Basic", value='shm1',  
      fluidRow(column(10, splitLayout(cellWidths=c('14%', '1%', '14%', '1%', '9%', '1%', '32%', '1%', '15%'),
      numericInput(inputId='height', label='Overall height:', value=400, min=1, max=Inf, step=NA, width=170), '',
      numericInput(inputId='width', label='Overall width:', value=760, min=1, max=Inf, step=NA, width=170), '',
      numericInput(inputId='col.n', label='Columns:', value=2, min=1, max=Inf, step=1, width=150), '',
      radioButtons(inputId="gen.con", label="Display by:", choices=c("Gene"="gene", "Condition"="con", "None"="none"), selected="gene", inline=TRUE), '', 
     radioButtons(inputId="pre.scale", label="Preserve.scale:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)
      )),
      column(1,
      dropdownButton(inputId='dropdown', label='Color key', circle=FALSE, icon=NULL, status='primary', inline=FALSE, width=250, 
      fluidRow(splitLayout(cellWidths=c('1%', '70%', '25%'), '', textInput("color", "Color scheme:", "purple,yellow,blue", placeholder='Eg: "purple,yellow,blue"', width=200),
      actionButton("col.but", "Go", icon=icon("refresh")))), 
      radioButtons(inputId='cs.v', label='Color scale based on:', choices=c("Selected rows"="sel.gen", "All rows"="w.mat"), selected="sel.gen", inline=TRUE)
      ))
      ), textOutput('h.w.c'), textOutput('msg.col'),

      fluidRow(splitLayout(cellWidths=c('100%'),  checkboxGroupInput(inputId="tis", label="Select tissues to be transparent:", choices='', selected='', inline=TRUE))),
      fluidRow(splitLayout(cellWidths=c('14%', '1%', '24%', '1%', '14%', '1%', '11%', '1%', '15%'), 
      downloadButton("dld.shm", "Download"), '',
      radioButtons(inputId='ext', label='File type:', choices=c('NA'='NA', "png"="png", "jpg"="jpg", "pdf"="pdf"), selected="NA", inline=TRUE), '',
      numericInput(inputId='res', label='Resolustion (dpi):', value=300, min=10, max=Inf, step=10, width=150), '',
      numericInput(inputId='lgd.w', label='Legend width:', value=0.7, min=0, max=1, step=0.1, width=150), '',
      numericInput(inputId='lgd.ratio', label='Legend aspect ratio:', value=1, min=0.0001, max=Inf, step=0.1, width=140)
      )),
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar1"), plotOutput("shm", height='auto'), "")))

      ))
  
  })

  output$shms.o <- renderUI({

    if (is.null(input$svgInpath1)) return(NULL)
    selectInput('shms.in', label='aSVG for legend:', choices=as.list(svg.path()[['svg.na']], selected=svg.path()[['svg.na']][1]))

  })

  output$lgd <- renderPlot(width='auto', height = "auto", {

    validate(need(try(as.integer(input$lgd.row)==input$lgd.row & input$lgd.row>0), 'Legend key rows should be a positive integer!'))
    validate(need(try(input$lgd.key.size>0&input$lgd.key.size<1), 'Legend key size should be between 0 and 1!'))
    svg.path <- svg.path()
    if (is.null(svg.path())|is.null(grob$lgd.all)) return(ggplot())

      # Width and height in original SVG.
      if (!is.null(input$svgInpath1)) svg.na <- input$shms.in else svg.na <- 1
      w.h <- svg.df()[[svg.na]][['w.h']]
      w.h <- as.numeric(gsub("^(\\d+\\.\\d+|\\d+).*", "\\1", w.h)); r <- w.h[1]/w.h[2]
      cat('Plotting legend plot... \n')
      g.lgd <- grob$lgd.all[[svg.na]]; g.lgd <- g.lgd+coord_fixed(ratio=r); return(g.lgd)

  })

  output$lgd.ui <- renderUI({ 
    
    if (is.null(input$hide.lgd)) return(NULL) 
    if (input$hide.lgd=='Y') return(NULL)
    box(title="Original Image", status="primary", solidHeader=TRUE, collapsible=TRUE, uiOutput('shms.o'),
    splitLayout(cellWidths=c("43%", "1%", '43%'),
    numericInput(inputId='lgd.row', label='Legend key rows:', value=2, min=1, max=Inf, step=1, width=150), '',
    numericInput(inputId='lgd.key.size', label='Legend key size:', value=0.04, min=0, max=1, step=0.02, width=150)
    ),
    splitLayout(cellWidths=c("43%", "1%", '43%'),
    radioButtons(inputId="lgd.label", label="Label feature:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE), '',
    numericInput(inputId='lgd.lab.size', label='Label size:', value=2.5, min=0, max=Inf, step=0.5, width=150)
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
    if (input$ggly.but=='N') return()
    if (is.null(geneIn())|is.null(gID$new)|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$gg.all)|input$ggly.but=='N') return(NULL)
    if (length(color$col=="none")==0|input$color=="") return(NULL)

    withProgress(message="Animation: ", value=0, {
    incProgress(0.25, detail="preparing frames...") 
    gg.all <- grob$gg.all; na <- names(gg.all)
    for (i in seq_along(gg.all)) {

      pat <- paste0(na[i], '_', "\\d+\\.html")
      if (length(list.files('www/ggly/', pat))>0) next
      gly <- ggplotly(gg.all[[i]], tooltip='text') %>% layout(showlegend=FALSE)
      gly$sizingPolicy$padding <- 0
      na0 <- paste0(na[i], '_', i, ".html")
      cat('Animation: saving', na0, '\n')
      saveWidget(gly, na0, selfcontained=FALSE, libdir="lib")
      system(paste0('mv ', na0, ' www/ggly/'))

    }; if (!dir.exists('www/ggly/lib')) system('mv lib/ www/ggly/') else if (dir.exists('lib/')) system('rm -rf lib')

    })

  })

  output$sld.fm <- renderUI({
 
    if (is.null(grob$gg.all)|is.null(con())|is.null(gID$geneID)) return(NULL) 
    con.pat <- paste0(unique(con()), collapse='|')
    gen.pat <- paste0(gID$geneID, collapse='|')
    gen.con.pat <- paste0('^(', gen.pat, ')_(', con.pat, ')_\\d+$') 
    sliderInput(inputId='fm', 'Frames', min=1, max=sum(grepl(gen.con.pat, names(grob$gg.all))), step=1, value=1, animate=animationOptions(interval=input$t*10^3, loop=FALSE, playButton=NULL, pauseButton='pause'))
  
  })

  # As long as the variable of 'reactive' is used in the 'ui.R', changes of elements in 'reactive' would cause chain change all the way to 'ui.R'. E.g. the change in "input$ggly.but=='N'" leads to changes in 'output$ggly' and 'ui.R', not necessarily changes in 'output$ggly' call changes in 'gly.url'.
  gly.url <- reactive({ 
    
    if (is.null(input$ggly.but)) return() 
    if (is.null(grob$gg.all)|input$ggly.but=='N'|gID$geneID[1]=='none'|is.null(pat.all())) return(NULL)
    pat <- paste0('^', pat.all(), '_\\d+_(\\d+)\\.html$')
    na <- list.files('www/ggly', pattern=pat)
    na <- na[order(gsub(pat, '\\3', na))][as.integer(input$fm)]
    if (length(na)==0|is.na(na)) return(NULL)
    cat('Animation: access', na, 'path \n')
    paste0('ggly/', na)
  
  })

  # Variables in 'observe' are accessible anywhere in the same 'observe'.
  observe({

    if (is.null(input$ggly.but)) return() 
    if (input$ggly.but=='N'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
    
    gg.all <- grob$gg.all; na <- names(gg.all)
    pat <- paste0(na, collapse='|')
    na <- list.files('www/ggly', pattern=paste0('^(', pat, ')_', input$fm, '\\.html$')); if (length(na)==0) return(NULL)
    pat1 <- paste0('^(', pat, ')_',input$fm, '\\.html$')
    gg.na <- gsub(pat1, '\\1', na); gg <- gg.all[[gg.na]]
    dat <- layer_data(gg); x.max <- max(dat$x); y.max <- max(dat$y)
    w <- 770; h <- y.max/x.max*w
    if (h>550) { h <- 550; w <- x.max/y.max*h }
    
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

    if (input$ggly.but=='N'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
    withProgress(message="Animation: ", value=0, {
    incProgress(0.75, detail="plotting...")
    gly.url <- gly.url(); cat('Animation: plotting', gly.url, '\n')
    tags$iframe(src=gly.url, height=input$height.ly, width=input$width.ly, scrolling='auto')
  
    })

  })

  })

  anm.dld <- reactive({
    
    if (input$ggly.but=='N'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL) 
    withProgress(message="Downloading animation: ", value=0, {
    incProgress(0.1, detail="in progress...")
    gg.all <- grob$gg.all; na <- names(gg.all)
    con.uni <- unique(con()); gen <- gID$geneID
    pat.con <- paste0(con.uni, collapse='|')
    pat.gen <- paste0(gen, collapse='|')
    pat <- paste0('^(', pat.gen, ')_(', pat.con, ')_\\d+$')
    gg <- gg.all[grepl(pat, na)]; gg.na <- names(gg)
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

    selectInput("vdo.dim", label="Fixed dimension:", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected='640x480', width=110)

  })

  output$ffm <- renderText({
    
    ffm <- tryCatch({ system('which ffmpeg', intern=TRUE) }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
    if (!grepl('ffmpeg', ffm)) paste("<span style=\"color:red\">Error: \"ffmpeg\" is not detected!\"</span>")

  })

  observeEvent(list(log=input$log, tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, pre.scale=input$pre.scale, vdo.but=input$vdo.but, vdo.dim=input$vdo.dim, vdo.itvl=input$vdo.itvl, vdo.height=input$vdo.height, vdo.width=input$vdo.width, vdo.res=input$vdo.res), {

    if (is.null(input$vdo.but)) return(NULL) 
    if (input$vdo.but=='N'|is.null(pat.all())) return(NULL)
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
    validate(need(try(!is.na(input$vdo.itvl)&input$vdo.itvl>0), 'Transition time should be a positive numeric!'))
    validate(need(try(!is.na(input$vdo.height)&input$vdo.height>0&input$vdo.height<=0.99), 'Height should be between 0.1 and 0.99!'))
    validate(need(try(!is.na(input$vdo.width)&input$vdo.width>0&input$vdo.width<=0.92), 'Width should be between 0.1 and 0.92!'))
    validate(need(try(!is.na(input$vdo.res)&input$vdo.res>=1&input$vdo.res<=700), 'Resolution should be between 1 and 700!'))
    
    withProgress(message="Video: ", value=0, {
    incProgress(0.75, detail="in progress...")
    gg.all <- grob$gg.all; na <- names(gg.all)
    pat <- paste0(pat.all(), '_\\d+$'); na <- na[grepl(pat, na)]
    gg.all1 <- gg.all[na]
    cat('Making video... \n')
    res <- input$vdo.res; dim <- input$vdo.dim
    if (dim %in% c('1280x800', '1280x1024', '1280x720')&res>450) res <- 450
    if (dim=='1920x1080'&res>300) res <- 300
    selectInput("vdo.dim", label="Fixed dimension:", choices=c('1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480'), selected='640x480', width=110)
    vdo <- video(gg=gg.all1, cs.g=shm.bar(), sam.uni=sam(), tis.trans=input$tis, lgd.key.size=input$lgd.key.size, lgd.text.size=NULL, position.text.key='right', label=(input$lgd.label=='Y'), label.size=input$lgd.lab.size, sub.title.size=8, bar.value.size=6, lgd.row=input$lgd.row, width=input$vdo.width, height=input$vdo.height, video.dim=dim, interval=input$vdo.itvl, res=res, out.dir='./www/video'); if (is.null(vdo)) return()
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
    updateRadioButtons(session, "mat.scale", "Scale: ", c("No", "By column", "By row"), "No", inline=TRUE)

  })


  observe({
  
    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$dt_rows_selected
    updateRadioButtons(session, inputId='ext', label='File type:', choices=c('NA'='NA', "png"="png", "jpg"="jpg", "pdf"="pdf"), selected="NA", inline=TRUE)
    updateRadioButtons(session, inputId="mhm.but", label="Show plot:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)
    updateRadioButtons(session, inputId="ggly.but", label="Show animation:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)
   updateRadioButtons(session, inputId="vdo.but", label="Show video:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)

  })

  observe({
   input$lgd.key.size; input$lgd.row; input$tis; input$lgd.label; input$lgd.lab.size
   updateRadioButtons(session, inputId="vdo.but", label="Show video:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)
  })

  observe({
    input$tis
    updateRadioButtons(session, inputId="ggly.but", label="Show animation:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)
  })
  # Calculate whole correlation or distance matrix.
  cor.dis <- reactive({

    if (is.null(geneIn())|input$mhm.but=='N') return()
    if (((input$fileIn=="custom_computed_data"|input$fileIn=="custom_data") & is.null(geneIn()))|input$fileIn=="none") return(NULL)
    
    withProgress(message="Compute similarity/distance matrix: ", value = 0, {

      incProgress(0.5, detail="Please wait...")
      gene <- geneIn()[['gene2']]
      if (input$measure=='cor') {
      
        m <- cor(x=t(gene))
        if (input$cor.abs==TRUE) { m <- abs(m) }; return(m)

      } else if (input$measure=='dis') { return(-as.matrix(dist(x=gene))) }

    })

  })

  # Subset nearest neighbours for target genes based on correlation or distance matrix.
  submat <- reactive({
    
    if (is.null(cor.dis())|input$mhm.but=='N') return()
    gene <- geneIn()[["gene2"]]; rna <- rownames(gene)
    gen.tar<- rna[input$dt_rows_selected]; mat <- cor.dis()
    
    # Validate filtering parameters in matrix heatmap.<Paste> 
    measure <- input$measure; cor.abs <- input$cor.abs; mhm.v <- input$mhm.v; thr <- input$thr
    if (input$thr=='p') {

      validate(need(try(mhm.v>0 & mhm.v<=1), 'Proportion should be between 0 to 1 !'))

    } else if (input$thr=='n') {

      validate(need(try(mhm.v>=1 & as.integer(mhm.v)==mhm.v & !is.na(mhm.v)), 'Number should be a positive integer !'))

    } else if (input$thr=='v' & measure=='cor') {

      validate(need(try(mhm.v>-1 & mhm.v <1), 'Correlation value should be between -1 to 1 !'))

    } else if (input$thr=='v' & measure=='dis') {

      validate(need(try(mhm.v>=0), 'Distance value should be non-negative !'))

    }

    withProgress(message="Selecting nearest neighbours: ", value = 0, {

      incProgress(0.5, detail="Please wait...")
      arg <- list(p=NULL, n=NULL, v=NULL)
      arg[names(arg) %in% input$thr] <- input$mhm.v
      if (input$measure=='dis' & input$thr=='v') arg['v'] <- -arg[['v']]
      gen.na <- do.call(sub_na, c(mat=list(mat), ID=list(gen.tar), arg)) 
      validate(need(try(length(gen.na)>=2), paste0('Only ', gen.na, ' selected !'))); return(gene[gen.na, ])

    })

  })


  # Plot matrix heatmap.
  output$HMly <- renderPlotly({
    
    if (is.null(submat())|input$mhm.but=='N') return()
    gene <- geneIn()[["gene2"]]; rna <- rownames(gene)
    gen.tar<- rna[input$dt_rows_selected]
    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="Plotting...")
      if (input$mat.scale=="By column") scale.hm <- 'column' else if (input$mat.scale=="By row") scale.hm <- 'row' else scale.hm <- 'no'
      matrix_hm(ID=gen.tar, data=submat(), scale=scale.hm, main='Target Genes and Their Nearest Neighbours', title.size=10, static=FALSE)

    })

  })


  adj.mod <- reactive({ 

    if (input$fileIn=="custom_computed_data") {

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
  
  adj.tree <- reactive({ 
    #gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    if (input$fileIn=="custom_data"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn)) {

      gene <- geneIn()[["gene2"]]; if (is.null(gene)) return()
      type <- input$net.type; sft <- if (type=='distance') 1 else 6

      withProgress(message="Computing: ", value = 0, {
        incProgress(0.3, detail="adjacency matrix.")
        incProgress(0.5, detail="topological overlap matrix.")
        incProgress(0.1, detail="dynamic tree cutting.")
        adjMod <- adj_mod(data=submat(), type=type, minSize=input$min.size, dir=NULL)
        adj <- adjMod[['adj']]; mod4 <- adjMod[['mod']]

      }); return(list(adj=adj, mod4=mod4))

    }

  })

  observe({
    
    input$gen.sel; input$measure; input$cor.abs; input$thr; input$mhm.v
    updateSelectInput(session, 'ds', "Module splitting sensitivity level:", 3:2, selected="3")

  })

  mcol <- reactive({

    if (!is.null(adj.tree()) & (input$fileIn=="custom_data"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn))) { 
      
    withProgress(message="Computing dendrogram:", value=0, {
      incProgress(0.7, detail="hierarchical clustering.")
      mod4 <- adj.tree()[['mod4']]
      
    }); return(mod4) 
    
  }

  })

  observe({
  
    geneIn(); gID$geneID; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    input$gen.sel; input$measure; input$cor.abs; input$thr; input$mhm.v
    updateSelectInput(session, "TOM.in", label="Adjacency threshold:", choices=c("None", sort(seq(0, 1, 0.002), decreasing=TRUE)), selected="None")

  })

  observe({
  
    geneIn(); gID$geneID; input$TOM.in; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    updateRadioButtons(session, inputId="cpt.nw", label="Show plot:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)

  })


  col.sch.net <- reactive({ 

    if(input$color.net=="") return(NULL) 
    col <- gsub(' |\\.|-|;|,|/', '_', input$color.net)
    col <- strsplit(col, '_')[[1]]
    col <- col[col!='']; col1 <- col[!col %in% colors()]
    if (length(col1>0)) validate(need(try(col1 %in% colors()), paste0('Colors not valid: ', col1, ' !'))); col

  }); color.net <- reactiveValues(col.net="none")

  len.cs.net <- 350
  observeEvent(input$col.but.net, {

    if (is.null(col.sch.net())) return (NULL)

      color.net$col.net <- colorRampPalette(col.sch.net())(len.cs.net)

  })

  visNet <- reactive({

    if (input$TOM.in=="None") return(NULL)
    if (input$fileIn=="custom_computed_data") { adj <- adj.mod()[[1]]; mods <- adj.mod()[[2]] } else if (input$fileIn=="custom_data"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn)) { adj <- adj.tree()[[1]]; mods <- mcol() }
    gene <- submat(); if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    lab <- mods[, input$ds][rownames(gene)==input$gen.sel]
    if (lab=="0") { showModal(modalDialog(title="Module", "The selected gene is not assigned to any module. Please select a different gene.")); return() }
    idx.m <- mods[, input$ds]==lab; adj.m <- adj[idx.m, idx.m]; gen.na <- colnames(adj.m) 
    idx.sel <- grep(paste0("^", input$gen.sel, "$"), gen.na); gen.na[idx.sel] <- paste0(input$gen.sel, "_target")
    colnames(adj.m) <- rownames(adj.m) <- gen.na
    withProgress(message="Computing network:", value=0, { 
      incProgress(0.8, detail="making network data frame")
      nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=input$gen.sel, adj.min=input$TOM.in)
      node <- nod.lin[['node']]; colnames(node) <- c('id', 'value')
      link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'
      if (nrow(link1)!=0) { 
        
        link1$title <- link1$value # 'length' is not well indicative of adjacency value, so replaced by 'value'.
        link1$color <- 'lightblue'
        
      }; ann <- geneIn()[[2]]
      if (!is.null(ann)) node <- cbind(node, title=ann[node$id, ], borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      if (is.null(ann)) node <- cbind(node, borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      net.lis <- list(node=node, link=link1)

    }); net.lis

  })

  output$bar.net <- renderPlot({  

    if (input$TOM.in=="None"|input$cpt.nw=="N") return(NULL)
    if (length(color.net$col.net=="none")==0) return(NULL)
    gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    if(input$col.but.net==0) color.net$col.net <- colorRampPalette(c('purple', 'yellow', 'blue'))(len.cs.net) # color.net$col.net is changed alse outside renderPlot, since it is a reactive value.
   
      withProgress(message="Color scale: ", value = 0, {
      incProgress(0.25, detail="Preparing data. Please wait.")
      incProgress(0.75, detail="Plotting. Please wait.")
      node <- visNet()[["node"]]; node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
        cs.net <- col_bar(geneV=v.net, cols=color.net$col.net, width=1, mar=c(3, 0.1, 3, 0.1)); return(cs.net) # '((max(v.net)-min(v.net))/len.cs.net)*0.7' avoids bar overlap.

      })

  })

  output$edge <- renderUI({ 

    if (input$TOM.in=="None") return(NULL)
    if (input$fileIn=="none"|(input$fileIn=="Your own" & is.null(geneIn()))|
    input$gen.sel=="None") return(NULL)
    HTML(paste0("&nbsp&nbsp&nbsp&nbsp Total edges to display (If > 300, the <br/> 
    &nbsp&nbsp&nbsp App can possibly get stuck.): ", dim((visNet()[["link"]]))[1]))

  })

  vis.net <- reactive({ 
    
    if (input$TOM.in=="None"|input$cpt.nw=="N") return(NULL)
    gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.

    withProgress(message="Network:", value=0.5, {
    incProgress(0.3, detail="prepare for plotting.")
    # Match colours with gene connectivity by approximation.
    node <- visNet()[["node"]]; node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
    col.nod <- NULL; for (i in node$value) {

      ab <- abs(i-v.net); col.nod <- c(col.nod, color.net$col.net[which(ab==min(ab))[1]])

    }; node$color <- col.nod
    visNetwork(node, visNet()[["link"]], height="300px", width="100%", background="", main=paste0("Network Module Containing ", input$gen.sel), submain="", footer= "") %>% visIgraphLayout(physics=FALSE, smooth=TRUE) %>% visOptions(highlightNearest=list(enabled=TRUE, hover=TRUE), nodesIdSelection=TRUE)

    })
    
  })

  output$vis <- renderVisNetwork({

    if (input$fileIn=="none"|(input$fileIn=="Your own" & is.null(geneIn()))|input$TOM.in=="None"|input$gen.sel=="None") return(NULL)
    if (input$cpt.nw=="N") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="plotting.")
      vis.net()

    })

  })

  onStop(function() { 

    if (dir.exists('www/ggly/')) {  cat("Removing animation files in 'www/ggly/' ... \n"); system('rm -fr www/ggly/lib/*'); system('rm -f www/ggly/*html') }
    if (dir.exists('html_shm/lib')) {  cat("Removing animation files in 'html_shm/lib/' ... \n"); system('rm -fr html_shm/lib/*'); system('rm -f html_shm/*html') }
    if (dir.exists('www/video/')) {  cat("Removing video file in 'www/video/' ... \n"); system('rm -fr www/video/.mp4*') }
    cat("Session stopped\n")

  })

})



