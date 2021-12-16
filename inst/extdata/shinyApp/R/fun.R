library(scater); library(scran); library(BiocSingular); library(SingleCellExperiment); library(Seurat); library(readr); library(Matrix)
library(igraph); library(ggplot2)

# Internal function of spatialHeatmap.
df_is_as <- function(data, fun=is.numeric) {                                                                                     
  dim.na <- dimnames(data)
  # apply(na, 2, is.na): requires much more memory than vapply.                                                                                                     
  if (identical(fun, is.na)) vap <- vapply(seq_len(ncol(data)), function(i) { fun(data[, i]) }, logical(nrow(data))) else if (identical(fun, as.numeric)) vap <- vapply(seq_len(ncol(data)), function(i) { fun(data[, i]) }, numeric(nrow(data))) else if (identical(fun, is.numeric)) vap <- vapply(seq_len(ncol(data)), function(i) { fun(data[, i]) }, logical(1))                        
  # is.numeric always returns one value.
  if (nrow(data)==1 | identical(fun, is.numeric)) vap <- matrix(vap, byrow=TRUE, ncol=ncol(data))                                
  if (!identical(fun, is.numeric)) dimnames(vap) <- dim.na else { rownames(vap) <- dim.na[[1]][1]; colnames(vap) <- dim.na[[2]] }
  return(vap)
}    

# Internal function of spatialHeatmap.
library(SummarizedExperiment)
check_data <- function(data, sam.factor=NULL, con.factor=NULL, usage='other') {
  # save(data, sam.factor, con.factor, usage, file='check.all')
  options(stringsAsFactors=FALSE)
  dat <- fct.cna <- col.meta <- row.meta <- con.na <- NULL
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) {
    if (is(data, 'dgCMatrix')) data <- as.matrix(data) 
    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data))
    if (!identical(cna, colnames(data))) cat('Syntactically valid column names are made! \n')
    # Only in replicate aggregation, data normalization, data filter, replicates are allowed. 
    if (!(usage %in% c('aggr', 'norm', 'filter'))) if (any(duplicated(cna))) stop('Please use function \'aggr_rep\' to aggregate replicates!') 
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    if (nrow(data)==1) na <- matrix(na, byrow=TRUE, ncol=ncol(data))
    na <- as.data.frame(na); rownames(na) <- rna
    # app <- apply(na, 2, is.na): requires much more memory than vapply.
    vap <- df_is_as(na, is.na); rownames(vap) <- rna
    # Exclude rows of all NAs.
    vap <- vap[!rowSums(vap)==ncol(vap), , drop=FALSE]
    idx <- colSums(vap)!=0
    row.meta <- data[idx] # aggr_rep filter_data, submatrix 
    dat <- na[!idx]; colnames(dat) <- fct.cna <- cna[!idx]
    # spatial_hm
    if (usage=='shm') { form <- grepl("__", fct.cna); if (sum(form)==0) { colnames(dat) <- paste0(fct.cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE }

  } else if (is(data, 'SummarizedExperiment')) {

    dat <- assay(data); r.na <- rownames(dat); cna <- make.names(colnames(dat))
    dat <- df_is_as(dat, as.numeric)
    rownames(dat) <- r.na; colnames(dat) <- cna
    if (!identical(cna, colnames(dat))) cat('Syntactically valid column names are made! \n')
    col.meta <- as.data.frame(colData(data))
    # filter_data
    row.meta <- as.data.frame(rowData(data))[, , drop=FALSE]
    # Factors teated by paste0/make.names are vectors.
    if (!is.null(sam.factor) & !is.null(con.factor)) { fct.cna <- colnames(dat) <- paste0(col.meta[, sam.factor], '__', col.meta[, con.factor]); con.na <- TRUE } 

    if (usage=='shm') {

      if (!is.null(sam.factor) & is.null(con.factor)) { sam.na <- as.vector(col.meta[, sam.factor]); fct.cna <- paste0(sam.na, "__", "con"); con.na <- FALSE } else if (is.null(sam.factor)) { form <- grepl("__", cna); if (sum(form)==0) { fct.cna <- paste0(cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE }
      if (!is.null(fct.cna)) { colnames(dat) <- make.names(fct.cna); if (!identical(fct.cna, make.names(fct.cna))) cat('Syntactically valid column names are made! \n') }
      if (any(duplicated(colnames(dat)))) stop('Please use function \'aggr_rep\' to aggregate \'sample__condition\' replicates!')
    
    } else if (usage %in% c('aggr', 'filter')) {
 
      if (!is.null(sam.factor) & is.null(con.factor)) {  fct.cna <- as.vector(col.meta[, sam.factor]) } else if (is.null(sam.factor) & !is.null(con.factor)) { fct.cna <- as.vector(col.meta[, con.factor]) } else fct.cna <- colnames(dat)
      if (!identical(fct.cna, make.names(fct.cna))) cat('Syntactically valid column names are made! \n')
      fct.cna <- colnames(dat) <- make.names(fct.cna) 

    }

  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", "dgCMatrix", or "SummarizedExperiment", except that "spatial_hm" also accepts a "vector".')
 
  }; return(list(dat=dat, fct.cna=fct.cna, col.meta=col.meta, row.meta=row.meta, con.na=con.na)) 

}

# Internal function of spatialHeatmap.
library(genefilter); library(SummarizedExperiment)
filter_data <- function(data, pOA=c(0, 0), CV=c(-Inf, Inf), top.CV=1, ann=NULL, sam.factor, con.factor, dir=NULL, verbose=TRUE) {

  options(stringsAsFactors=FALSE)
  if (top.CV>1|top.CV<0) stop('"top.CV" should be between 0 and 1!')
  # Process data.
  dat.lis <- check_data(data=data, sam.factor=sam.factor, con.factor=con.factor, usage='filter')
  expr <- dat.lis$dat; row.meta <- dat.lis$row.meta; col.meta <- dat.lis$col.meta
  if (verbose==TRUE) { cat('All values before filtering:\n'); print(summary(unlist(as.data.frame(expr)))) }
  expr.t <- as.data.frame(t(expr)); cv.all <- sort(vapply(expr.t, sd, numeric(1))/vapply(expr.t, mean, numeric(1)), decreasing=TRUE)
  if (verbose==TRUE) { cat('All coefficient of variances (CVs) before filtering:\n'); print(summary(cv.all)) }
  if (top.CV<1) {
    cv.min <- cv.all[ceiling(length(cv.all)*top.CV)]
    CV <- c(cv.min, Inf)
  }
  ffun <- filterfun(pOverA(p=pOA[1], A=pOA[2]), cv(CV[1], CV[2]))
  filtered <- genefilter(expr, ffun); expr <- expr[filtered, , drop=FALSE] # Subset one row in a matrix, the result is a numeric vector not a matrix, so drop=FALSE.
  if (verbose==TRUE) { cat('All values after filtering:\n'); print(summary(unlist(as.data.frame(expr)))) }
  expr.t <- as.data.frame(t(expr))
  if (verbose==TRUE & nrow(expr)>0) {
    cv.all <- sort(vapply(expr.t, sd, numeric(1))/vapply(expr.t, mean, numeric(1)), decreasing=TRUE)
    cat('All coefficient of variances (CVs) after filtering:\n'); print(summary(cv.all))
  }
  row.meta <- row.meta[filtered, , drop=FALSE]
  rownames(row.meta) <- NULL # Some rownames could have been automatically mangled
                             # earlier to keep them unique so they don't necessarily
                             # match 'rownames(expr)' anymore. This will cause the
                             # call to 'SummarizedExperiment()' below to fail with
                             # SummarizedExperiment 1.23.2.
  if (!is.null(dir)) {
    dir <- normalizePath(dir, winslash="/", mustWork=FALSE)
    if (!dir.exists(dir)) stop(paste0(dir, ' does not exist!'))
    if (is(data, 'data.frame')|is(data, 'matrix')) {
      expr1 <- cbind.data.frame(expr, row.meta, stringsAsFactors=FALSE)
    }  else if (is(data, 'SummarizedExperiment')) {
      
      if (ncol(row.meta)>0 & !is.null(ann)) {
        expr1 <- cbind.data.frame(expr, row.meta[, ann], stringsAsFactors=FALSE)
        colnames(expr1)[ncol(expr1)] <- ann
      } else expr1 <- expr

    }
    write.table(expr1, paste0(dir, "/customData.txt"), sep="\t", row.names=TRUE, col.names=TRUE)  
  }

  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'dgCMatrix')) { return(cbind(expr, row.meta)) } else if (is(data, 'SummarizedExperiment')) {
  
    rownames(col.meta) <- NULL # If row names present in colData(data), if will become column names of assay(data).
    expr <- SummarizedExperiment(assays=list(expr=expr), rowData=row.meta, colData=col.meta)
    if (!is.null(assayNames(data))) SummarizedExperiment::assayNames(expr) <- assayNames(data)[1]; return(expr)

  }

}

# Filter genes/cells by proportion of min count in row/column.
preprocess_sc <- function(data, gen.rm='^ATCG|^ATCG', min.cnt=1, p.in.cell=0.5, p.in.gen=0.3) {
  # Remove mitochondrial, chloroplast and protoplasting-induced genes.
  if (!is.null(gen.rm)) { 
    rm <- grepl(gen.rm, rownames(data))
    data <- data[!rm, ]; dim(data) 
  }
  if (min.cnt==0|(p.in.cell==0 & p.in.gen==0)|(min.cnt==0 & p.in.cell==0 & p.in.gen==0)) return(data)
  # Filter genes according to min count by genes and by cells. 
  cat('Before filtering :\n'); print(dim(data))
  print(head(sort(rowSums(data)), 5))
  print(head(sort(colSums(data)), 5))
  idx <- data >= min.cnt  
  idx.r <- rowSums(idx) >= p.in.gen*ncol(data)
  idx <- idx[idx.r, , drop=FALSE]
  idx.c <- colSums(idx) >= p.in.cell*nrow(idx)
  data <- data[idx.r, idx.c, drop=FALSE]
  cat('After filtering :\n'); print(dim(data))
  print(head(sort(rowSums(data)), 5))
  print(head(sort(colSums(data)), 5))
  return(data)
}

# Filter single cell data in a list.
filter_sc <- function(lis, blk, gen.rm='^ATCG|^ATCG', min.cnt=1, p.in.cell=0.5, p.in.gen=0.3) {
  gen.na <- list(rownames(blk))
  for (i in seq_along(lis)) {
    lis[[i]] <- preprocess_sc(data=lis[[i]], gen.rm=gen.rm, min.cnt=min.cnt, p.in.cell=p.in.cell, p.in.gen=p.in.gen)
    gen.na <- c(gen.na, list(rownames(lis[[i]])))
  }
  gen.ovl <- Reduce(intersect, gen.na)
  lis <- c(list(bulk=blk[gen.ovl, , drop=FALSE]), lis)
  for (i in seq_along(lis)) {
    lis[[i]] <- lis[[i]][gen.ovl, , drop=FALSE]; print(dim(lis[[i]]))
  }; return(lis)
}

# Normalize by CPM. The log2 values are transformed to power of 2 and then CPM. Lastly, CPM-values are transformed back to log2.
# sce.nor: output of "logNormCounts", where library size factors are applied already.
cpm_sc <- function(sce.nor) {
  log.cnt <- logcounts(sce.nor)
  cnt.cpm <- calculateCPM(2^log.cnt-1)
  logcounts(sce.nor) <- log2(cnt.cpm+1)
  return(sce.nor)
}

# Normalize one or more data sets.
# dat.lis: if multiple data are contained, they are combined, and normalized together. The normalized data is split by the original data and the original data is replaced by normalized data in the returned list. 
norm_sc <- function(dat.lis, cpm=FALSE) {
  library(scran); library(SingleCellExperiment)
  nas <- names(dat.lis)
  if (any(nas=='')) stop('The list should be named!')
  # Anchor sce.
  dat0 <- dat.lis[[1]]
  if (is(dat0, 'SingleCellExperiment')) dat0 <- assays(dat0)[[1]]
  sce.sc <- SingleCellExperiment(assays=list(counts=as.matrix(dat0[, 1, drop=FALSE])), colData=DataFrame(set=0))
  for (i in seq_along(dat.lis)) {
    dat0 <- dat.lis[[i]]
    if (is(dat0, 'SingleCellExperiment')) sce0 <- dat0
    if (is(dat0, 'dgCMatrix')|is(dat0, 'matrix')|is(dat0, 'data.frame')) sce0 <- SingleCellExperiment(assays=list(counts=as.matrix(dat0)))
    colData(sce0)$set <- nas[i]; sce.sc <- cbind(sce.sc, sce0)
    }; sce.sc <- sce.sc[, -1]

    min.size <- 100
    if (min.size > ncol(sce.sc)) { print('Fewer cells than min size in quickCluster!'); return()
    }
    set.seed(1000); clusters <- quickCluster(sce.sc, min.size=min.size)
    sce.sc <- computeSumFactors(sce.sc, cluster=clusters, min.mean=1)
    sce.sc.nor <- logNormCounts(sce.sc)
    # CPM.
    if (cpm==TRUE) sce.sc.nor <- cpm_sc(sce.sc.nor)
    for (i in nas) {
      dat.lis[[i]] <- subset(sce.sc.nor, , set==i)
    }
    return(dat.lis)
}


# Cluster single cells or mixture of single cells and bulk.
# data: log2 data (normalized) in one of SingleCellExperiment, dgCMatrix, matrix, data.frame.
cluster_sc <- function(data, min.dim=5, max.dim=55, pca=FALSE, method='snn', dimred='PCA') {
  library(scran); library(SingleCellExperiment); library(scater)
  if (is(data, 'SingleCellExperiment')) sce <- data
  if (is(data, 'dgCMatrix')|is(data, 'matrix')|is(data, 'data.frame')) {
    if (all(round(data)==data)) stop('The data to cluster should be in log2 scale!')
    sce <- SingleCellExperiment(list(logcounts=as.matrix(data)))
  }
  df.var.sc <- modelGeneVar(sce) # Use logcounts by default.
  top.hvgs.sc <- getTopHVGs(df.var.sc, prop=0.1)

  library(BiocSingular)
  set.seed(101011001)
  sce.dimred <- denoisePCA(sce, assay.type = "logcounts", technical=df.var.sc, subset.row=top.hvgs.sc, min.rank=min.dim, max.rank=max.dim)
  if (pca==TRUE) return(sce.dimred)
  library(scater)
  # Other argument: n_dimred, ntop. By default only 2 dimensions are returned by runTSNE/runUMAP.
  # runUMAP returns different results before/after runTSNE.
  sce.dimred <- runUMAP(sce.dimred, dimred="PCA", ncomponents=min.dim)
  sce.dimred <- runTSNE(sce.dimred, dimred="PCA", ncomponents=2)
  
  # dims <- list(pca='PCA', tsne='TSNE', umap='UMAP') 
  # Only one is accepted: assay.type = "logcounts" or use.dimred="PCA".
  if (method=='snn') gr.sc <- buildSNNGraph(sce.dimred, use.dimred=dimred)
  if (method=='knn') gr.sc <- buildKNNGraph(sce.dimred, use.dimred=dimred)
  library(igraph)
  colLabels(sce.dimred) <- as.character(cluster_walktrap(gr.sc)$membership)
  cdat.sc <- colData(sce.dimred); cdat.sc$cell <- rownames(cdat.sc)
  colData(sce.dimred) <- cdat.sc
  return(sce.dimred)
}

# Assign true bulk to cells in colData slot. 
true_bulk <- function(sce, df.match) {
  # match.lis <- lapply(match.lis, function(x) { paste0(x, collapse=',') })
  true.bulk <- df.match$trueBulk; names(true.bulk) <- df.match$cell
  svg.bulk <- df.match$SVGBulk; names(svg.bulk) <- df.match$cell
  colData(sce)$trueBulk <- true.bulk[colData(sce)$cell]
  colData(sce)$SVGBulk <- svg.bulk[colData(sce)$cell]
  return(sce)
}

# Plot target cells/bulk in top layers with target colors.
# sce.dim: output of plotTSNE
# tar.col: a vector of colors, named with target cells.
dimred <- function(sce.dim, tar.col) {
  cdat <- colData(sce.dim); selected <- cdat$id <- rownames(cdat)
  selected[!selected %in% names(tar.col)] <- 'other'
  cdat$selected <- selected
  colData(sce.dim) <- cdat
  set.seed(10)
  g.dim <- plotTSNE(sce.dim, colour_by="selected")
  df1 <- g.dim$data
  df1$colour_by <- factor(df1$colour_by, levels=c('other', names(tar.col))) 
  df1 <- df1[order(df1$colour_by), ]
  # size <- seq_along(levels(df1$colour_by))
  size <- rep(2, length(levels(df1$colour_by)))
  names(size) <- levels(df1$colour_by)
  col.o <- 'grey70'; names(col.o) <- 'other'
  cols <- c(tar.col, col.o)
  g <- ggplot(df1, aes(x=X, y=Y, shape=colour_by, colour=colour_by, size=colour_by))+geom_point()+scale_size_manual(values=size)+ scale_colour_manual(values=cols)+scale_shape_manual(values=seq_along(levels(df1$colour_by)))+labs(colour='', size='', shape='', x='TSNE1', y='TSNE2')+theme_classic()
  return(g)

}


# Filter cells according to sim values in each cell cluster. If dimred is NULL, log-scale values are used and the output is converted to counts. If dimred='PCA', the sims are calculated on PCs.
filter_cell <- function(sce.tsne, dimred=NULL, sim=0.5, sim.p=0.5, sim.meth='spearman') {
  colData(sce.tsne)$idx.all <- seq_len(ncol(sce.tsne))
  lab <- as.character(colData(sce.tsne)$label); lab.uni <- sort(unique(lab))
  dat <- dat.dim <- cdat <- NULL
  if (is.null(dimred)) { # Use expression values. 
  dat.aggr <- data.frame(matrix(NA, nrow=nrow(sce.tsne), ncol=1))
  idx.kp <- NULL; for (i in lab.uni) {
    sce <- sce.tsne[, lab==i]
    # Calculate PCCs on log-values. 
    clus <- logcounts(sce)
    sims <- cor(clus, method=sim.meth); idx <- sims > sim
    w <- which(colSums(idx) > ncol(clus)*sim.p)
    idx.kp <- c(idx.kp, colData(sce)$idx.all[w])
    # In practice, cell ids are not known.
    sce <- sce[, w]; cat(paste0('Cluster', i, ': ', length(w)), '\n')
    print(table(colnames(sce)))

  }; sce.kp <- sce.tsne[, idx.kp]; return(sce.kp)
  } else { # Use top dims.
    if (!dimred %in% reducedDimNames(sce.tsne)) stop('The name of reduced dimensionality is not found!')
    cdat.all <- colData(sce.tsne)
    if (log2==FALSE) { dat.all <- logcounts(sce.tsne)^2-1; dat.all[dat.all<0] <- 0; dat.all <- round(dat.all) } else dat.all <- logcounts(sce.tsne) 
    reddim <- t(reducedDim(sce.tsne, dimred))
    if (nrow(reddim)<5) stop('Dimensions should be >= 5!')
    dat.aggr.dim <- data.frame(matrix(NA, nrow=nrow(reddim), ncol=1))
    dat.aggr <- data.frame(matrix(NA, nrow=nrow(sce.tsne), ncol=1))
    for (i in lab.uni) {	    
      # Use PC-based similarity to filter.
      dim0 <- reddim[, lab==i]; dat.all0 <- dat.all[, lab==i]      
      sims.dim <- cor(dim0, method=sim.meth); idx.dim <- sims.dim > sim
      w.dim <- which(colSums(idx.dim) > ncol(dim0)*sim.p)
      # Remaining cell data frame of PCs.
      dim0 <- dim0[, w.dim, drop=FALSE]; cat(paste0('Cluster', i, ': ', length(w.dim)), '\n')
      print(table(colnames(dim0)))
      dat.dim <- cbind(dat.dim, dim0); 
      # Remaining cell data frame of log counts.
      dat.all0 <- dat.all0[, w.dim, drop=FALSE] 
      dat <- cbind(dat, dat.all0)
      # Meta data.
      cdat <- rbind(cdat, cdat.all[lab==i, ][w.dim, , drop=FALSE])
    if (aggr==TRUE & length(w.dim)>0) {
      # Aggregate remaining cells by PCs.
      df.ag.dim <- data.frame(i=rowMeans(dim0))
      rownames(df.ag.dim) <- rownames(dim0)
      colnames(df.ag.dim) <- i
      dat.aggr.dim <- cbind(dat.aggr.dim, df.ag.dim)
      # Aggregate remaining cells by log counts.
      df.ag <- data.frame(i=rowMeans(dat.all0))
      rownames(df.ag) <- rownames(dat.all0)
      colnames(df.ag) <- i
      dat.aggr <- cbind(dat.aggr, df.ag)
    }
  }; dat.aggr.dim <- dat.aggr.dim[, -1, drop=FALSE]
  dat.aggr <- dat.aggr[, -1, drop=FALSE]
  }; print(dim(dat.dim))
  return(list(data.dim=dat.dim, data=dat, cdat=cdat, dat.aggr.dim=dat.aggr.dim, dat.aggr=dat.aggr))
}


# Create pseudo bulk by randomly sampling cells in each cell type, and sum these sampled cells together by cell type. 
sd_bulk <- function(sc, p) {
  idx.sc.sel <- NULL
  cna.sc <- colnames(sc); cna.sc.uni <- unique(cna.sc)
  # Randomly sample cells in each cell type
  for (i in cna.sc.uni) {
    idx.sc <- seq_len(ncol(sc)); idx.sc0 <- idx.sc[cna.sc==i]
    # If length(x)==1, sample(x, n) is sample(1:x, n).
    if (length(idx.sc0)==1) idx.sel <- idx.sc0 else {
      # Size to sample in each cell type.
      size <- ceiling(length(idx.sc0)*p)
      set.seed(size); idx.sel <- sample(idx.sc0, size) 
    }
    idx.sc.sel <- c(idx.sc.sel, idx.sel)
  }
  sc.sam <- sc[, sort(idx.sc.sel)]
  sc.sam.cna <- colnames(sc.sam)
  # Sum each sampled cells in each cell type.
  blk.sd <- vapply(sort(unique(sc.sam.cna)), function(x) rowSums(sc.sam[, sc.sam.cna==x, drop=FALSE]), numeric(nrow(sc.sam)))
  rownames(blk.sd) <- rownames(sc.sam)
  colnames(blk.sd) <- toupper(colnames(blk.sd))
  blk <- blk.sd
  # Match list.
  blk.sd.na.uni <- unique(colnames(blk.sd))
  match.lis <- as.list(toupper(blk.sd.na.uni))
  na <- tolower(blk.sd.na.uni)
  # substr(na, 1, 1) <- toupper(substr(na, 1, 1))
  names(match.lis) <- na; return(list(bulk.sd=blk.sd, match.lis=match.lis))
}


# Calculate ROC for the combination strategy.
com_roc <- function(dat.com, cdat.com, dat.blk, df.match, sim.meth='spearman') { 
  if (nrow(dat.com)<5) stop('Dimensions should be >= 5!')
  # Co-cluster labels.
  lab.com <- cdat.com$label
  labs.com <- as.character(sort(lab.com))
  lab.com.uni <- unique(labs.com)
  # Co-clusters in a list.
  coclus <- lapply(lab.com.uni, function(x) sort(table(rownames(cdat.com)[lab.com==x])))
  names(coclus) <- lab.com.uni; bulk.na <- unique(colnames(dat.blk))
  # Index all bulk and cells for convenience to selected cells with a bulk tissue assignment.
  cna.all <- colnames(dat.com)
  names(cna.all) <- seq_along(cna.all)

  # In each co-cluster, compute PCC between cells and bulk.
  df.roc <- NULL; for (i in lab.com.uni) {
    # Cell and bulk names in a co-cluster.
    co.na <- names(coclus[[i]])
    blk.na0 <- unique(bulk.na[bulk.na %in% co.na])
    if (length(blk.na0)==0) next
    cat('Bulk in cocluster ', i, ':\n', sep='')
    cat(blk.na0, '\n', sep=' ')

    # Subset a cocluster in the combined data frame.
    clus0 <- dat.com[, lab.com==i, drop=FALSE]
    cna <- colnames(clus0); cna.all0 <- cna.all[lab.com==i]
    # Aggregate bulk in the subsetted co-cluster.
    mat <- vapply(blk.na0, function(x) rowMeans(clus0[, cna==x, drop=FALSE]), numeric(nrow(clus0)))

    # Combine aggregated bulk back with cells. In practice, cell ids are not known.
    sc0 <- clus0[, !colnames(clus0) %in% blk.na0, drop=FALSE]
    cna.sc0 <- cna.all0[!colnames(clus0) %in% blk.na0] 
    sc0.na.uni <- unique(colnames(sc0))
    blk.ag.sc <- cbind(mat, sc0)
    # Correlation between aggregated bulk and cells.
    sim <- cor(blk.ag.sc, method=sim.meth)
    sim.sub <- sim[rownames(sim) %in% blk.na0, !colnames(sim) %in% blk.na0, drop=FALSE]
    sc0.na <- colnames(sim.sub)
    # In each cocluster, compare each cell with each bulk according to sims, TRUE/FALSE is decided by the bulk with the largest sim with the cell.
   for (j in seq_len(ncol(sim.sub))) {
     # All sims of one cell between each bulk. Select the largest.
     blk.sel <- sort(sim.sub[, j], decreasing=TRUE)[1]
     if (length(blk.na0)==1) names(blk.sel) <- blk.na0

     # Matching list between the cells and bulk in the co-cluster.
     # match.lis0 <- match.lis[sc0.na[j]]
     true.bulk.pas <- df.match$trueBulk[df.match$cell==sc0.na[j]]
     true.bulk <- strsplit(true.bulk.pas, ',|;| ')[[1]]
     true.bulk <- true.bulk[true.bulk!='']
     svg.bulk <- df.match$SVGBulk[df.match$cell==sc0.na[j]]
     # Check if the cell matches with the true bulk.
     # idx is based on bulk+cell.
     df0 <- data.frame(assignedBulk=names(blk.sel), cell=sc0.na[j], response=names(blk.sel) %in% true.bulk, predictor=blk.sel, idx=names(cna.sc0[j]), trueBulk=true.bulk.pas, SVGBulk=svg.bulk)
     df.roc <- rbind(df.roc, df0)
   }

  }
  if (is.null(df.roc)) return()
  df.tr <- subset(df.roc, response==TRUE)
  cat('TRUE and FALSE:\n'); print(dim(df.roc))
  cat('TRUE:\n'); print(dim(df.tr))
  if (length(unique(df.roc$response))!=2) return()
  # idx is based on cells only.
  # idx.sc.asg <- seq_along(dat.com)[as.numeric(df.roc$idx)]
  # idx.all <- c(rep(0, ncol(dat.blk)), seq_along(colnames(dat.com)[!colnames(dat.com) %in% bulk.na]))
  # df.roc$idx <- idx.all[idx.sc.asg]
  df.roc$idx <- as.numeric(df.roc$idx) - ncol(dat.blk)
  library(pROC)
  roc.obj <- roc(df.roc$response, df.roc$predictor, smoothed = TRUE, ci=TRUE, ci.alpha=0.9, stratified=FALSE, plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=TRUE, show.thres=TRUE, direction='<', levels=c('FALSE', 'TRUE'))
  return(list(roc.obj=roc.obj, df.roc=df.roc))
}


# Co-cluster bulk and filtered cells, and compute AUCs.
# cell.kp.lis: returned by filter_cell.
coclus_roc <- function(blk=NULL, cell.kp, df.match, min.dim=10, max.dim=55, method='snn', dimred='PCA', sim.meth='spearman') {
  # Combine bulk and single cell data.
  dat.kp.sc <- logcounts(cell.kp)
  cat('Single cell data: \n'); print(dim(dat.kp.sc))
  # Two or more cells are needed.
  if (length(unique(colnames(dat.kp.sc)))<=1) {
    cat('At least two cells are needed!\n'); return()
  }
  if (is(blk, 'SingleCellExperiment')) blk <- logcounts(blk)
  if (all(round(blk)==blk)) stop('The bulk data should be in log2 scale!')
  inter <- intersect(rownames(blk), rownames(dat.kp.sc))
  blk.kp <- as.matrix(blk[inter, ]); sc.kp <- dat.kp.sc[inter, ]
  
  # cell.aggr <- cell.aggr[inter, , drop=FALSE]
  # Check correlations between cell clusters and bulk. 
  # cna.blk <- colnames(blk.kp)
  # Aggregate every bulk.
  # blk.aggr <- vapply(unique(cna.blk), function(x) rowMeans(blk.kp[, cna.blk==x, drop=FALSE]), numeric(nrow(blk.kp)))
  # cor.blk.sc <- cor(cbind(blk.aggr, cell.aggr), method=sim.meth)
  # cna.co <- colnames(cor.blk.sc)
  # mean(cor.blk.sc[cna.co %in% cna.blk, !cna.co %in% cna.blk])

  com.kp <- as.matrix(cbind(blk.kp, sc.kp)); dim(com.kp)
  table(colnames(blk.kp)); table(colnames(sc.kp))

  # Cluster cells+bulk.
  sce.tsne.com <- cluster_sc(data=com.kp, min.dim=min.dim, max.dim=max.dim, method=method, dimred=dimred)
  # sce.tsne.com <- clus.com$sce.tsne; sce.com.nor <- clus.com$sce.nor
  cdat.com <- colData(sce.tsne.com); lab.com <- cdat.com$label
  # data.com: logcounts(sce.com.nor), reducedDim(sce.tsne.com, 'PCA')
  # roc <- com_roc(dat.com=logcounts(sce.com.nor), cdat.com, dat.blk=blk)
  roc.lis <- com_roc(dat.com=t(reducedDim(sce.tsne.com, dimred)), cdat.com=cdat.com, dat.blk=blk, df.match=df.match, sim.meth=sim.meth)
  return(roc.lis)
}


# Assign bulk to cell clusters, and compare cells with bulk by PCC in each cluster.
asg_blk_roc <- function(blk=NULL, cell.kp.lis, match.lis, cpm=FALSE, sim.p.blk=0.9, min.dim=10, max.dim=55, sim.meth='spearman', dimred='PCA') {
  oat.kp.sc <- cell.kp.lis$data
  cat('Single cell data: \n'); print(dim(dat.kp.sc))
  # Two or more cells are needed.
  if (length(unique(colnames(dat.kp.sc)))<=1) {
    cat('At least two cells are needed!\n'); return()
  }
  inter <- intersect(rownames(blk), rownames(dat.kp.sc))
  blk.kp <- as.matrix(blk[inter, ]); sc.kp <- dat.kp.sc[inter, ]
  dat.aggr <- cell.kp.lis$dat.aggr[inter, , drop=FALSE]

  # Check correlations between cell clusters and bulk. 
  cna.blk <- colnames(blk.kp); cna.blk.uni <- unique(cna.blk)
  # Aggregate every bulk.
  blk.aggr <- vapply(unique(cna.blk), function(x) rowMeans(blk.kp[, cna.blk==x, drop=FALSE]), numeric(nrow(blk.kp)))
  cor.blk.sc <- cor(cbind(blk.aggr, dat.aggr), method=sim.meth)
  cna.co <- colnames(cor.blk.sc)
  mean(cor.blk.sc[cna.co %in% cna.blk, !cna.co %in% cna.blk])

  # Combine bulk and single cell data.
  com.kp <- as.matrix(cbind(blk.aggr, sc.kp))
  table(colnames(blk.kp)); table(colnames(sc.kp))

  # Top PCs of cells+bulk.
  com <- cluster_sc(data=com.kp, cpm=cpm, min.dim=min.dim, max.dim=max.dim, pca=TRUE, dimred=dimred)
  com.pca <- t(reducedDim(com, dimred))
  if (nrow(com.pca)<5) stop('Dimensions should be >= 5!')
  cna.com <- colnames(com.pca)
  blk.com <- com.pca[, cna.blk.uni, drop=FALSE]
  sc.com <- com.pca[, !cna.com %in% cna.blk.uni, drop=FALSE]
  cdat.sc <- cell.kp.lis$cdat; lab.sc <- cdat.sc$label
  lab.sc.uni <- unique(lab.sc); table(lab.sc)
  # In each cell cluster, assign bulk and compare cell with bulk by PCC.
  df.roc <- NULL; for (i in lab.sc.uni) {
    # Assign bulk to each cell cluster.
    cell0 <- sc.com[, lab.sc==i, drop=FALSE]; dim(cell0)
    table(colnames(cell0))
    co.blk.cell <- cor(cbind(blk.com, cell0), method=sim.meth)
    cna0 <- colnames(co.blk.cell)
    sim.sub <- co.blk.cell[cna0 %in% cna.blk.uni, !cna0 %in% cna.blk.uni, drop=FALSE]
    rsum.sim <- rowSums(sim.sub)
    if (all(rsum.sim<0)) {
      blk.sel0 <- names(rsum.sim[rsum.sim>=max(rsum.sim)/sim.p.blk])
    } else {
      blk.sel0 <- names(rsum.sim[rsum.sim>=max(rsum.sim)*sim.p.blk])
    }
    sc0.na <- colnames(sim.sub)
    sim.sub <- sim.sub[blk.sel0, , drop=FALSE]

    # In each cocluster, compare each cell with each bulk according to sims, TRUE/FALSE is decided by the bulk with the largest sim with the cell.
   for (j in seq_len(ncol(sim.sub))) {
     # All sims of one cell between each bulk. Select the largest.
     blk.max <- sort(sim.sub[, j], decreasing=TRUE)[1]
     if (length(blk.sel0)==1) names(blk.max) <- blk.sel0
     # Matching list between the cells and bulk in the co-cluster.
     match.lis0 <- match.lis[sc0.na[j]]
     # Check if the cell matches with the true bulk.
     df0 <- data.frame(bulk=names(blk.max), cell=sc0.na[j], response=names(blk.max) %in% match.lis0[[1]], predictor=blk.max)
     df.roc <- rbind(df.roc, df0)
   }
  }

  if (is.null(df.roc)) return()
  df.tr <- subset(df.roc, response==TRUE)
  cat('TRUE and FALSE:\n'); print(dim(df.roc))
  cat('TRUE:\n'); print(dim(df.tr))
  if (length(unique(df.roc$response))!=2) return()
  library(pROC)
  roc.obj <- roc(df.roc$response, df.roc$predictor, smoothed = TRUE, ci=TRUE, ci.alpha=0.9, stratified=FALSE, plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=TRUE, show.thres=TRUE, direction='<', levels=c('FALSE', 'TRUE'))
  return(list(roc.obj=roc.obj, df.roc=df.roc))

}

# Plot AUCs in grouped bars.
# df.par: output of coclus_roc
bar_group <- function(df.par, sim='sim', sim.p='sim.p', auc='auc', dim='dim', bar.lab='dim', bar.width=0.8, dodge=0.8, lab.size=3,   title=NULL, lgd.key.size=0.05, ylim=c(0, 1)) {    
  cna <- colnames(df.par)                                         
  cna[cna==sim] <- 'sim'; cna[cna==sim.p] <- 'sim.p'                                                                               
  cna[cna==auc] <- 'auc'; cna[cna==dim] <- 'dim'                                                                                     
  colnames(df.par) <- cna                                                                                                          
  df.par$sim <- as.character(df.par$sim)                                                                                           
  df.par$sim.p <- as.character(df.par$sim.p)                                                                                       
  df.par$auc <- round(df.par$auc, 2)                                                                                               
  gg <- ggplot(df.par, aes(x=sim, y=auc, fill=sim.p)) +                                                                            
  geom_bar(color='black', width=bar.width, position=position_dodge2(width=dodge, preserve="single"), stat="identity") +            
  labs(title=title, x="sim threshold", y='AUC', fill='sim.p')+theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.key.size = unit(lgd.key.size, 'npc')) + ylim(ylim[1], ylim[2])
                                                                                                                                   
  if (bar.lab=='dim') gg <- gg + geom_text(aes(label=dim), vjust=-0.2, size=lab.size,  position=position_dodge2(width=dodge,         preserve="single"))                                                                                                                
  if (bar.lab=='auc') gg <- gg + geom_text(aes(label=auc), vjust=-0.2, size=lab.size,  position=position_dodge2(width=dodge,       preserve="single"))                                                                                                                
  return(gg)                                                                                                                       
}      

# Test the co-clustering workflow with different parameter combinations.
# bulk, scell: normalized externally, in log2 scale.
tests <- function(bulk, scell, df.match, sc.dim.min=10, max.dim=50, sim=round(seq(0.2, 0.8, by=0.1), 1), sim.p=round(seq(0.2, 0.8, by=0.1), 1), dim=seq(5, 40, by=1), method='snn', dimred='PCA', sim.meth='spearman', return.all=FALSE, file='df.par') {
  # save(bulk, scell, df.match, sc.dim.min, max.dim, sim, sim.p, dim, method, dimred, sim.meth, return.all, file='tests.par')
  df.par <- data.frame(sim=numeric(), sim.p=numeric(), dim=numeric(), auc=numeric(), total=numeric(), true=numeric(), thr=numeric(), sens=numeric(), spec=numeric())
  for (i in sim) {
    for (j in sim.p) {
      df.par <- rbind(df.par, data.frame(sim=i, sim.p=j, dim=dim, auc=0, total=0, true=0, thr=0, sens=0, spec=0))
    }
  }; print(head(df.par, 3))
  if (is(bulk, 'SingleCellExperiment')|is(bulk, 'SummarizedExperiment')) bulk <- logcounts(bulk)
  if (all(round(bulk)==bulk)) stop('The bulk data should be in log2 scale!')
  inter <- intersect(rownames(bulk), rownames((scell)))
  bulk <- bulk[inter, , drop=FALSE]; scell <- scell[inter, , drop=FALSE]
  # dat.com <- cbind(bulk, scell)

  if (!is.null(sc.dim.min)) {
    # dat.com.nor <- norm_sc(data=dat.com, cpm=cpm)
    # scell <- dat.com.nor[, !colnames(bulk) %in% colnames(dat.com.nor)]
    clus.sc <- cluster_sc(data=scell, min.dim=sc.dim.min, max.dim=max.dim, pca=FALSE, method=method, dimred=dimred)
  }
  for (i in seq_len(nrow(df.par))) {
    print(df.par[i, , drop=FALSE])
    if (is.null(sc.dim.min)) clus.sc <- cluster_sc(data=scell, min.dim=df.par$dim[i], max.dim=max.dim, method=method, dimred=dimred)
    clus.sc <- true_bulk(clus.sc, df.match)
    # Filter cells by PCC and proportion p.
    cell.kp <- filter_cell(clus.sc, dimred=NULL, sim=df.par$sim[i], sim.p=df.par$sim.p[i], sim.meth=sim.meth)
    # Cocluster.
    roc.lis <- coclus_roc(blk=bulk, cell.kp=cell.kp, df.match=df.match, min.dim=df.par$dim[i], max.dim=max.dim, method=method, dimred=dimred, sim.meth=sim.meth); roc.lis[[1]]
    if (return.all==TRUE) return(c(list(sce=cell.kp), roc.lis))
    roc.obj <- roc.lis$roc.obj; if (is.null(roc.obj)) next
    df.roc <- roc.lis$df.roc
    df.par$auc[i] <- round(auc(roc.obj), 3)
    df.par$total[i] <- nrow(df.roc)
    df.par$true[i] <- sum(df.roc$response)
    best <- round(coords(roc.obj, x='best'), 3)
    if (nrow(best)==1) { 
      df.par$thr[i] <- round(best$threshold, 3)
      df.par$sens[i] <- round(best$sensitivity, 3)
      df.par$spec[i] <- round(best$specificity, 3)
    } else if (unique(abs(best$threshold))==Inf) { 
      df.par$thr[i] <- df.par$sens[i] <- df.par$spec[i] <- Inf
    }
    print(df.par[seq_len(i), , drop=FALSE])
    if (!is.null(file)) saveRDS(df.par, file=file)
  }; return(df.par)

}

# Find optimal combinations of coclustering parameters.
# lis: list of parameter tests results. A slot is one parameter data frame returned by tests.
opt_par <- function(lis, auc=0.7, total=500) {
  rows <- NULL; nas <- names(lis)
  for (i in seq_along(lis)) {
    df.par <- get(load(lis[[i]]))
    lis[[i]] <- df.par <- df.par[df.par$auc >= auc & df.par$total >= total, , drop=F]
    rows <- c(rows, rownames(df.par))
  }; names(lis) <- nas
  tab <- sort(table(rows)); tab[tab==max(tab)]
  row.sel <- names(tab[tab==max(tab)])
  # Extract optimal parameters from parameter data frame.
  opt <- lapply(seq_along(lis), function(i) { 
    df0 <- lis[[i]][row.sel, ]
    df0[!rowSums(is.na(df0)) == ncol(df0), , drop=F]
    }
  )
  names(opt) <- nas
  return(list(lis=lis, rows=sort(table(rows)), opt=opt))
}        

# Plot optimal parameters.
# lis: list of extracted data frames of optimal parameter combinations.
bar_opt <- function(lis, row) {
  gg.lis <- NULL
  for (i in seq_along(lis)) {
    df.opt <- lis[[i]][as.character(row), , drop = FALSE]
    if (all(is.na(df.opt))) next
    gg <- bar_group(df.par=df.opt, sim='sim', sim.p='sim.p', auc='auc', dim='dim', bar.width=0.8, dodge=0.8, lab.size=3,
title=names(lis)[i], lgd.key.size=0.02, ylim=c(0, 1)) 
    gg.lis <- c(gg.lis, list(gg)) 
  }; return(gg.lis) 
}

# Statistics of test results.

stat <- function(path, auc.min=0.000001, total.min=300, true.min=100) {
  df.par <- readRDS(path)
  df.par <- subset(df.par, auc >= auc.min & total >= total.min & true >= true.min)
  print(dim(df.par)); print(mean(df.par$auc))
  return(df.par)
}




