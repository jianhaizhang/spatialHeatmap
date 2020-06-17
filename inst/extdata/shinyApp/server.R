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

    w <- is.na(as.numeric(gsub('\\D', '', vec)))
    let <- vec[w]; num <- vec[!w]
    let.num <- as.numeric(gsub('\\D', '', vec))
    vec[!w] <- num[order(let.num[!w])]
    vec[w] <- let[order(let.num[w])]; return(vec)
        
  }  

  if (by=='gene') {

    # Sort conditions under each gene.
    con.pat1 <- paste0('.*_(', paste0(con.all, collapse='|'), ')$')
    na.sort <- NULL; for (i in sort_mix(ID.sel)) {
        
      na0 <- na.all[grepl(paste0('^', i, '_'), na.all)]
      con1 <- gsub(con.pat1, '\\1', na0)
      na.sort <- c(na.sort, paste0(i, '_', sort_mix(con1)))

    }

  } else if (by=='con') {

    # Sort conditions and genes.
    gen.pat1 <- paste0('^(', paste0(ID.sel, collapse='|'), ')_.*')
    na.sort <- NULL; for (i in sort_mix(con.all)) {
      
      na0 <- na.all[grepl(paste0('_', i, '$'), na.all)]
      gen1 <- gsub(gen.pat1, '\\1', na0)
      na.sort <- c(na.sort, paste0(sort_mix(gen1), '_', i))

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





col_bar <- function(geneV, cols, width, bar.title.size=10, mar=c(3, 0.1, 3, 0.1)) {        

  color_scale <- y <- NULL
  cs.df <- data.frame(color_scale=geneV, y=1)
  cs.g <- ggplot()+geom_bar(data=cs.df, aes(x=color_scale, y=y), fill=cols, orientation='x', stat="identity", width=((max(geneV)-min(geneV))/length(geneV))*width)+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=margin(t=mar[1], r=mar[2], b=mar[3], l=mar[4], "cm"), panel.grid=element_blank(), panel.background=element_blank(), plot.title= element_text(hjust=0.7, size=bar.title.size))+coord_flip()+labs(title="Value", x=NULL, y=NULL)
  if (max(abs(geneV))<=10000 & max(abs(geneV))>=0.0001) cs.g <- cs.g+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0)) else cs.g <- cs.g+scale_x_continuous(expand=c(0,0), labels=function(x) format(x, scientific=TRUE))+scale_y_continuous(expand=c(0,0))+theme(axis.text.y=element_text(angle=45))
  return(cs.g)

}

lay_shm <- function(lay.shm, con, ncol, ID.sel, grob.list, width, height, shiny) {

  width <- as.numeric(width); height <- as.numeric(height); ncol <- as.numeric(ncol); con <- unique(con)
  grob.all.na <- names(grob.list)
  if (lay.shm=="gene") {

    all.cell <- ceiling(length(con)/ncol)*ncol
    cell.idx <- c(seq_len(length(con)), rep(NA, all.cell-length(con)))
    m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
    lay <- NULL; for (i in seq_len(length(ID.sel))) { lay <- rbind(lay, m+(i-1)*length(con)) }
    # Sort conditions under each gene.
    na.sort <- sort_gen_con(ID.sel=ID.sel, na.all=grob.all.na, con.all=con, by='gene')
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
    
    if (length(intersect(title[dup], feature))>0) return(paste0('Duplicated title text detected: ', paste0(title[dup], collapse=' '), '!')) else {

      w <- title %in% title[dup]
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
  lis <- list(df=df, tis.path=sub('_\\d+$', '', tit), fil.cols=fil.cols); return(lis)

}


grob_list <- function(gene, con.na=TRUE, geneV, coord, ID, cols, tis.path, tis.trans=NULL, sub.title.size, sam.legend='identical', legend.col, legend.title=NULL, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.5, legend.label.size=8, legend.title.size=8, line.size=0.2, line.color='grey70', ...) {

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

    if (lgd==FALSE) scl.fil <- scale_fill_manual(values=g.col, guide=FALSE) else { 

      # Show selected or all samples in legend.
      if (length(sam.legend)==1) if (sam.legend=='identical') sam.legend <- intersect(sam.uni, unique(tis.path)) else if (sam.legend=='all') sam.legend <- unique(tis.path)
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

       }; scl.fil <- scale_fill_manual(values=g.col, breaks=tis.df[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=legend.title, ncol=legend.ncol, nrow=legend.nrow)) 

    }

    g <- ggplot(...)+geom_polygon(data=coord, aes(x=x, y=y, fill=tissue), color=line.color, size=line.size, linetype='solid')+scl.fil+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=sub.title.size))+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))
    if (con.na==FALSE) g.tit <- ggtitle(k) else g.tit <- ggtitle(paste0(k, "_", con)); g <- g+g.tit

    if (lgd==TRUE) {

      g <- g+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.01, 0.01, 0.2, 0, "npc"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=15, face="bold"), legend.position=legend.position, legend.direction=legend.direction, legend.background = element_rect(fill=alpha(NA, 0)), legend.key.size=unit(legend.key.size, "cm"), legend.text=element_text(size=legend.label.size), legend.title=element_text(size=legend.title.size), legend.margin=margin(l=0.1, r=0.1, unit='cm'))+ggtitle('Legend')

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
    tis.col <- gsub("(.*)(__)(.*)", "\\1", c.na); g.lis <- NULL
    grob.na0 <- paste0(k, "_", con.uni); g.lis <- lapply(con.uni, g_list)
    # Repress popups by saving it to a png file, then delete it.
    tmp <- tempfile()
    png(tmp); grob <- lapply(g.lis, ggplotGrob); dev.off(); if (file.exists(tmp)) do.call(file.remove, list(tmp))
    names(g.lis) <- names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob); g.lis.all <- c(g.lis.all, g.lis)

  }; g.lgd <- g_list(con=NULL, lgd=TRUE)
  return(list(grob.lis=grob.lis, g.lgd=g.lgd, g.lis.all=g.lis.all))

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



library(SummarizedExperiment); library(shiny); library(shinydashboard); library(grImport); library(rsvg); library(ggplot2); library(DT); library(gridExtra); library(ggdendro); library(WGCNA); library(grid); library(xml2); library(plotly); library(data.table); library(genefilter); library(flashClust); library(visNetwork); library(reshape2); library(igraph)

# Import input matrix.
fread.df <- function(input, isRowGene, header, sep, fill, rep.aggr='mean') {
        
  df0 <- fread(input=input, header=header, sep=sep, fill=fill)
  cna <- make.names(colnames(df0)[-ncol(df0)])
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
inter.svg <- readLines("example/root_cross_final.svg")
inter.data <- read.table("example/expr_arab.txt", header=TRUE, row.names=1, sep="\t")

shinyServer(function(input, output, session) {

  output$dld.svg <- downloadHandler(

    filename=function(){ "root_cross_final.svg"}, 
    content=function(file){ writeLines(inter.svg, file) }

  )

  output$dld.data <- downloadHandler(

    filename=function(){ "expr_arab.txt" },
    content=function(file){ write.table(inter.data, file, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t") }

  )
  # Instruction.
  output$sum <-renderUI({ includeHTML("file/summary.html") })
  output$input <-renderUI({ includeHTML("file/input.html") })
  output$input1 <-renderUI({ includeHTML("file/input1.html") })
  output$shm.ins <-renderUI({ includeHTML("file/shm.html") })
  output$mhm.ins <-renderUI({ includeHTML("file/mhm.html") })
  output$net.ins <-renderUI({ includeHTML("file/net.html") })
  # Acknowledgement.
  output$ack <-renderUI({ includeHTML("file/acknowledgement.html") })

  # Filter parameters.
  fil <- reactiveValues(P=0, A=-Inf, CV1=-Inf, CV2=Inf)

  observe({

    input$fileIn; input$geneInpath
    updateRadioButtons(session, inputId="dimName", label="Step 3: is column or row gene?", 
    inline=TRUE, choices=c("None", "Row", "Column"), selected="None")
    updateSelectInput(session, 'sep', 'Step 4: separator', c("None", "Tab", "Comma", "Semicolon"), "None")
    updateRadioButtons(session, inputId='log', label='Data transform:', choices=c("No", "log2", "exp2"), selected="No", inline=TRUE)
    updateRadioButtons(session, inputId='cs.v', label='Color scale based on:', choices=c("Selected rows"="sel.gen", "All rows"="w.mat"), selected="sel.gen", inline=TRUE)
    updateNumericInput(session, inputId="height", label="Overall height:", value=400, min=0.1, max=Inf, step=NA)
    updateNumericInput(session, inputId="width", label="Overall width:", value=820, min=0.1, max=Inf, step=NA)
    updateNumericInput(session, inputId="col.n", label="No. of columns:", value=2, min=1, max=Inf, step=1)

  })

  observe({

    input$fileIn; input$geneInpath; input$log
    updateNumericInput(session, inputId="A", label="Value (A) to exceed:", value=0) 
    updateNumericInput(session, inputId="P", label="Proportion (P) of samples with values >= A:", value=0, min=0, max=1)
    updateNumericInput(session, inputId="CV1", label="Min coefficient of variation (CV1):", value=-10^4)
    updateNumericInput(session, inputId="CV2", label="Max coefficient of variation (CV2):", value=10^4,)
    fil$P <- 0; fil$A=-Inf; fil$CV1 <- -Inf; fil$CV2 <- Inf

  })

  observeEvent(input$fil.but, {

    if (input$fileIn=="None") return(NULL)  
    fil$P <- input$P; fil$A <- input$A; fil$CV1 <- input$CV1; fil$CV2 <- input$CV2
  
  })

  output$fil.par <- renderText({
    
    if (input$fileIn=="None") return(NULL)  
    P <- input$P
    validate(need(try(P<=1 & P>=0), 'P should be between 0 to 1 !'))

  })

  geneIn0 <- reactive({

    if (input$fileIn=="None") return(NULL)  
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

    if ((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & 
    (is.null(input$geneInpath)|input$dimName=="None"|input$sep=="None")) return(NULL)
    if ((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & 
    !is.null(input$geneInpath) & input$dimName!="None" & input$sep!="None") {

      incProgress(0.25, detail="Importing matrix. Please wait.")
      geneInpath <- input$geneInpath; if (input$sep=="Tab") sep <- "\t" else if (
      input$sep=="Comma") sep <- "," else if (input$sep=="Semicolon") sep <- ";"
      df.upl <- fread.df(input=geneInpath$datapath, isRowGene=(input$dimName=='Row'), header=TRUE, sep=sep, fill=TRUE, rep.aggr='mean'); return(df.upl)
   
    } 

    })

  })

  # Transform data.
  geneIn1 <- reactive({

    if (is.null(geneIn0())) return(NULL)
    gene2 <- geneIn0()[['gene2']]; if (input$log=='log2') {
   
      g.min <- min(gene2)
      if (g.min<0) gene2 <- gene2-g.min+1; if (g.min==0) gene2 <- gene2+1; gene2 <- log2(gene2)

    }; if (input$log=='exp2') gene2 <- 2^gene2
    gene3 <- geneIn0()[['gene3']]; gen.rep <- geneIn0()[['gen.rep']]
    return(list(gene2=gene2, gene3=gene3, gen.rep=gen.rep))

  })

  geneIn <- reactive({
    
    if (is.null(geneIn1())) return(NULL)    
    gene2 <- geneIn1()[['gene2']]; gene3 <- geneIn1()[['gene3']]; input$fil.but
    # Input variables in "isolate" will not triger re-excution, but if the whole reactive object is trigered by "input$fil.but" then code inside "isolate" will re-excute.
    isolate({
  
      se <- SummarizedExperiment(assays=list(expr=as.matrix(gene2)), rowData=gene3)
      if (ncol(gene3)>0) ann.col <- colnames(gene3)[1] else ann.col <- NULL
      se <- filter_data(data=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=c(fil$CV1, fil$CV2), dir=NULL)
      if (nrow(se)==0) { validate(need(try(nrow(se)>0), 'All rows are filtered out !')); return() }

      gene2 <- as.data.frame(assay(se), stringsAsfactors=FALSE); colnames(gene2) <- make.names(colnames(gene2))
      gene3 <- as.data.frame(rowData(se))[, , drop=FALSE]
      
    }); return(list(gene2=gene2, gene3=gene3))

  })

  output$dt <- renderDataTable({

    if (is.null(geneIn())) return()
    if (((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & is.null(geneIn()))|input$fileIn=="None") return(NULL)

    withProgress(message="Data table: ", value = 0, {

      incProgress(0.5, detail="Displaying. Please wait.")
      if (input$fileIn!="None") {

      gene <- geneIn(); gene.dt <- cbind.data.frame(gene[["gene2"]][, , drop=FALSE], gene[["gene3"]][, , drop=FALSE], stringsAsFactors=FALSE) 

   }

    datatable(gene.dt, selection=list(mode="multiple", target="row", selected=c(1)),
    filter="top", extensions='Scroller', options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
    formatRound(colnames(geneIn()[["gene2"]]), 2)

    })

  })

  gID <- reactiveValues(geneID="none", new=NULL, all=NULL)
  observe({ input$geneInpath; input$fileIn; gID$geneID <- "none" })
  observe({ if (is.null(geneIn())) gID$geneID <- "none" })
  # To make the "gID$new" and "gID$all" updated with the new "input$fileIn", since the selected row is fixed (3rd row), the "gID$new" is not updated when "input$fileIn" is changed, and the downstream is not updated either. The shoot/root examples use the same data matrix, so the "gID$all" is the same (pre-selected 3rd row) when change from the default "shoot" to others like "organ". As a result, the "gene$new" is null and downstream is not updated. Also the "gene$new" is the same when change from shoot to organ, and downstream is not updated, thus "gene$new" and "gene$all" are both set NULL above upon new "input$fileIn".  
  observeEvent(input$fileIn, { gID$all <- gID$new <- NULL })

  observeEvent(input$dt_rows_selected, {
    
    if (is.null(input$dt_rows_selected)) return()
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    gID$new <- setdiff(gID$geneID, gID$all); gID$all <- c(gID$all, gID$new)

    })


  geneV <- reactive({

    if (is.null(geneIn())) return(NULL)
    if (input$cs.v=="sel.gen" & is.null(input$dt_rows_selected)) return(NULL)
    if (input$fileIn!="None") { if (input$cs.v=="sel.gen" & !is.null(input$dt_rows_selected)) gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
    if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]] } 
    seq(min(gene), max(gene), len=1000) # len must be same with that from the function "spatial_hm()". Otherwise the mapping of a gene value to the colour bar is not accurate. 

  })

  col.sch <- reactive({ 

    if(input$color=="") return(NULL)
    col <- gsub(' |\\.|-|;|,|/', '_', input$color)
    col <- strsplit(col, '_')[[1]]
    col <- col[col!='']; col1 <- col[!col %in% colors()]
    if (length(col1>0)) validate(need(try(col1 %in% colors()), paste0('Colors not valid: ', col1, ' !'))); col

  }); color <- reactiveValues(col="none")

  # As long as a button is used, observeEvent should be used. All variables inside 'observeEvent' trigger code evaluation, not only 'eventExpr'.  
  observeEvent(input$col.but, {

    if (is.null(col.sch())) return (NULL)
    if (input$fileIn!="None") { color$col <- colorRampPalette(col.sch())(length(geneV())) }

  })

  output$bar <- renderPlot({

    if (is.null(gID$all)) return(NULL)
    if ((grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn) & !is.null(geneIn()))|((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & !is.null(input$svgInpath) & !is.null(geneIn()))) {

      if (length(color$col=="none")==0|input$color==""|is.null(geneV())) return(NULL)
      if(input$col.but==0) color$col <- colorRampPalette(c('purple', 'yellow', 'blue'))(length(geneV()))

      withProgress(message="Color scale: ", value = 0, {

        incProgress(0.75, detail="Plotting. Please wait.")
        cs.g <- col_bar(geneV=geneV(), cols=color$col, width=1, mar=c(3, 0.1, 3, 0.1)); return(cs.g)

      })

    }

  })


  observe({

    if (is.null(gID$all)) return(NULL)
    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    r.na <- rownames(geneIn()[["gene2"]]); gen.sel <- r.na[input$dt_rows_selected]
    updateSelectInput(session, "gen.sel", choices=c("None", gen.sel), selected="None")

  })

  svg.path <- reactive({

    if (input$fileIn=="Compute locally"|input$fileIn=="Compute online") { svg.path <- input$svgInpath$datapath; svg.na <- input$svgInpath$name } else if (input$fileIn=="brain_Prudencio") { svg.path <- "example/homo_sapiens.brain.svg"; svg.na <- "homo_sapiens.brain.svg" } else if (input$fileIn=="mouse_Merkin") { svg.path <- "example/mus_musculus.male.svg"; svg.na <- "mus_musculus.male.svg" } else if (input$fileIn=="chicken_Cardoso.Moreira") { svg.path <- "example/gallus_gallus.svg"; svg.na <- "gallus_gallus.svg" } else if (input$fileIn=="shoot_Mustroph") { svg.path <- "example/shoot_final.svg"; svg.na <- "shoot_final.svg" } else if (input$fileIn=="organ_Mustroph") { svg.path <- "example/organ_final.svg"; svg.na <- "organ_final.svg" } else if (input$fileIn=="root_Mustroph") { svg.path <- "example/root_cross_final.svg"; svg.na <- "root_cross_final.svg" } else if (input$fileIn=="shoot_root_Mustroph") { svg.path <- "example/shoot_root_final.svg"; svg.na <- "shoot_root_final.svg" } else if (input$fileIn=="root_roottip_Mustroph") { svg.path <- "example/root_roottip_final.svg"; svg.na <- "root_roottip_final.svg" } else if (input$fileIn=="map_Census") { svg.path <- "example/us_map_final.svg"; svg.na <- "us_map_final.svg" } else return(NULL)
    return(list(svg.path=svg.path, svg.na=svg.na))

  })

  sam <- reactive({ 

    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\1", c.na) else return(NULL) 

  })
  
  svg.df <- reactive({ 

    if (is.null(gID$all)) return(NULL)
    if (((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & 
    !is.null(input$svgInpath))|(grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn) & is.null(input$svgInpath))) {

      withProgress(message="Tissue heatmap: ", value=0, {
    
        incProgress(0.5, detail="Extracting coordinates. Please wait.") 
        df_tis <- svg_df(svg.path=svg.path()[['svg.path']], feature=sam())
        validate(need(!is.character(df_tis), df_tis))
        return(df_tis)

      })

    }

  })


  observe({

    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; svg.df()
    updateCheckboxGroupInput(session, inputId="tis", label='Select tissues to be transparent:', choices=intersect(unique(sam()), unique(svg.df()[['tis.path']])), selected='', inline=TRUE)

  })

  con <- reactive({ 

    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\3", c.na) else return(NULL) 

  })

  grob <- reactiveValues(all=NULL, gg.all=NULL); observeEvent(input$fileIn, { grob$all <- grob$gg.all <- NULL })
  gs <- reactive({ 

    if (is.null(svg.df())|is.null(gID$new)|is.null(gID$all)) return(NULL); 
    withProgress(message="Tissue heatmap: ", value=0, {

      incProgress(0.25, detail="preparing data.")
      if (input$cs.v=="sel.gen") gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
      if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]]
      g.df <- svg.df()[["df"]]; tis.path <- svg.df()[["tis.path"]]; fil.cols <- svg.df()[['fil.cols']]
      grob.lis <- grob_list(gene=gene, con.na=geneIn0()[['con.na']], geneV=geneV(), coord=g.df, ID=gID$new, legend.col=fil.cols, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=21) # Only gID$new is used.
      return(grob.lis)

    })

  })

  # Extension of 'observeEvent': any of 'input$log; input$tis; input$col.but; input$cs.v' causes evaluation of all code. 
  observe({
    
    input$log; input$tis; input$col.but; input$cs.v # input$tis as an argument in "grob_list" will not cause evaluation of all code, thus it is listed here.
    grob$all <- grob$gg.all <- NULL; gs <- reactive({ 

      if (is.null(svg.df())) return(NULL); if (input$cs.v=="sel.gen" & is.null(input$dt_rows_selected)) return(NULL)
      withProgress(message="Spatial heatmap: ", value=0, {
        incProgress(0.25, detail="preparing data.")

        g.df <- svg.df()[["df"]]; tis.path <- svg.df()[["tis.path"]]
        if (input$cs.v=="sel.gen") gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
        if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]]
        g.df <- svg.df()[["df"]]; tis.path <- svg.df()[["tis.path"]]; fil.cols <- svg.df()[['fil.cols']]
        grob.lis <- grob_list(gene=gene, geneV=geneV(), coord=g.df, ID=gID$all, legend.col=fil.cols, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=20) # All gene IDs are used.
        return(grob.lis)

      })

    }); grob$all <- gs()[['grob.lis']]; grob$gg.all <- gs()[['g.lis.all']]

  })

  observeEvent(gID$new, { 
                 
    grob.all <- c(grob$all, gs()[['grob.lis']]); grob$all <- grob.all[unique(names(grob.all))] 
    gg.all <- c(grob$gg.all, gs()[['g.lis.all']]); grob$gg.all <- gg.all[unique(names(gg.all))] 
  
  })
  
  output$h.w.c <- renderText({
    
    if (is.null(geneIn())|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)

    height <- input$height; width <- input$width
    col.n <- input$col.n;
    validate(need(height>=0.1 & !is.na(height), 'Height should be a positive numeric !'))
    validate(need(width>=0.1 & !is.na(width), 'Width should be a positive numeric !'))
    validate(need(col.n>=1 & as.integer(col.n)==col.n & !is.na(col.n), 'No. of columns should be a positive integer !'))

  })

  # In "observe" and "observeEvent", if one code return (NULL), then all the following code stops. If one code changes, all the code renews.
  observe({
    
    if (is.null(geneIn())|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)

    output$shm <- renderPlot(width=as.numeric(input$width)/2*as.numeric(input$col.n), height=as.numeric(input$height)*length(input$dt_rows_selected), {

    if (is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)
    if (length(color$col=="none")==0|input$color=="") return(NULL)

    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    grob.na <- names(grob$all); con <- unique(con())
    idx <- NULL; for (i in gID$geneID) { idx <- c(idx, grob.na[grob.na %in% paste0(i, '_', con)]) } 
    grob.lis.p <- grob$all[idx] #grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    lay_shm(lay.shm=input$gen.con, con=con, ncol=input$col.n, ID.sel=gID$geneID, grob.list=grob.lis.p, width=input$width, height=input$height, shiny=TRUE) 

    })

  })

  output$lgd <- renderPlot(width=260, {

    if (is.null(svg.path())|is.null(gs())) return(ggplot())

      svg.path <- svg.path()[['svg.path']]
      # Width and height in original SVG.
      doc <- read_xml(svg.path); w.h <- c(xml_attr(doc, 'width'), xml_attr(doc, 'height'))
      w.h <- as.numeric(gsub("^(\\d+\\.\\d+|\\d+).*", "\\1", w.h)); r <- w.h[1]/w.h[2]
      g.lgd <- gs()[['g.lgd']]; g.lgd <- g.lgd+coord_fixed(ratio =r); return(g.lgd)

  })


  observe({

    geneIn(); input$adj.modInpath; input$A; input$p; input$cv1
    input$cv2; input$min.size; input$net.type
    input$measure; input$cor.abs; input$thr; input$mhm.v
    updateRadioButtons(session, "mat.scale", "Scale: ", c("No", "By column", "By row"), "No", inline=TRUE)

  })


  observe({
  
    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2
    updateRadioButtons(session, inputId="mhm.but", label="Show plot: ", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)

  })


  # Calculate whole correlation or distance matrix.
  cor.dis <- reactive({

    if (is.null(geneIn())|input$mhm.but=='N') return()
    if (((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & is.null(geneIn()))|input$fileIn=="None") return(NULL)
    
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

    if (input$fileIn=="Compute locally") {

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
    if (input$fileIn=="Compute online"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn)) {

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

    if (!is.null(adj.tree()) & (input$fileIn=="Compute online"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn))) { 
      
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
    if (input$fileIn=="Compute locally") { adj <- adj.mod()[[1]]; mods <- adj.mod()[[2]] } else if (input$fileIn=="Compute online"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn)) { adj <- adj.tree()[[1]]; mods <- mcol() }
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
    if (input$fileIn=="None"|(input$fileIn=="Your own" & is.null(geneIn()))|
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

    if (input$fileIn=="None"|(input$fileIn=="Your own" & is.null(geneIn()))|input$TOM.in=="None"|input$gen.sel=="None") return(NULL)
    if (input$cpt.nw=="N") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="plotting.")
      vis.net()

    })

  })

})



