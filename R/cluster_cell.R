#' Cluster single cells or combination of single cells and bulk
#'
#' Cluster only single cell data or combination of single cell and bulk data. Clusters are created by first building a graph, where nodes are cells and edges represent connections between nearest neighbors, then partitioning the graph. The cluster labels are stored in the \code{cluster} column of \code{colData} slot of \code{SingleCellExperiment}.

#' @param sce The single cell data or combination of single cell and bulk data at log2 scale after dimensionality reduction in form of \code{SingleCellExperiment}.
#' @param graph.meth Method to build a nearest-neighbor graph, \code{snn} (see \code{\link[bluster]{makeSNNGraph}}) or \code{knn} (default, see \code{\link[bluster]{makeKNNGraph}}). The clusters are detected by first creating a nearest neighbor graph using \code{snn} or \code{knn} then partitioning the graph. 
#' @param dimred A string of \code{PCA} (default) or \code{UMAP} specifying which reduced dimensions to use for creating a nearest neighbor graph. 
#' @param dims The number of dimensions in `dimred` to use.
#' @param assay.type If `dimred=NULL`, the assay in `assayNames(sce)` to use.
#' @param knn.gr Additional arguments in a named list passed to \code{\link[bluster]{makeKNNGraph}}.
#' @param snn.gr Additional arguments in a named list passed to \code{\link[bluster]{makeSNNGraph}}.
#' @param cluster The clustering method. One of \code{wt} (\code{\link[igraph]{cluster_walktrap}}, default), \code{fg} (\code{\link[igraph]{cluster_fast_greedy}}), \code{le} (\code{\link[igraph]{cluster_leading_eigen}}), \code{sl} (\code{\link[igraph]{cluster_fast_greedy}}), \code{eb} (\code{\link[igraph]{cluster_edge_betweenness}}).  
#' @param wt.arg,fg.arg,sl.arg,le.arg,eb.arg A named list of arguments passed to \code{wt}, \code{fg}, \code{le}, \code{sl},       \code{eb} respectively.

#' @return A \code{SingleCellExperiment} object. 

#' @examples

#' library(SingleCellExperiment)
#' set.seed(123)  # for reproducibility

#' # random counts: 1000 genes × 100 cells
#' counts_mat <- matrix(
#'   rpois(1000 * 100, lambda = 5), nrow = 1000, ncol = 100,
#'   dimnames = list(
#'     paste0("Gene", 1:1000),
#'     paste0("Cell", 1:100)
#'   )
#' )

#' sce <- SingleCellExperiment(
#'   assays = list(counts = counts_mat)
#' )

#' sce.qc = qc_cell(sce)
#' sce.norm <- norm_cell(sce.qc)
#' sce.dimred <- reduce_dim(sce.norm)

#' \donttest{
#' sce.clus <- cluster_cell(sce=sce.dimred, graph.meth='snn', dimred='PCA')
#' # Clusters.
#' table(colData(sce.clus)$cluster)
#'
#' }
#' # See details in function "coclus_meta" by running "?coclus_meta".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. https://igraph.org
#' Lun A (2026). _bluster: Clustering Algorithms for Bioconductor_. doi:10.18129/B9.bioc.bluster <https://doi.org/10.18129/B9.bioc.bluster>. R package version 1.22.0, <https://bioconductor.org/packages/bluster>.

#' @export 
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment SingleCellExperiment 

cluster_cell <- function(sce, graph.meth='knn', dimred='PCA', dims=50, assay.type = NULL, knn.gr=list(k = 10, directed = FALSE), 
                         snn.gr=list(k = 10), cluster='wt', wt.arg=list(steps = 4), fg.arg=list(), sl.arg=list(spins = 25), 
                         le.arg=list(), eb.arg=list()) {
  # save(sce, graph.meth, dimred, dims, assay.type, knn.gr, snn.gr, cluster, wt.arg, fg.arg, sl.arg, le.arg, eb.arg, file='cluster_cell.arg')
  pkg <- check_pkg('bluster'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  if (!is(sce, 'SingleCellExperiment')) stop('The "sce" should be an "SingleCellExperiment" object!')
  if (!dimred %in% reducedDimNames(sce)) stop('The "dimred" is not detected in "reducedDimNames"!')
  if ('cluster' %in% colnames(colData(sce))) stop('The "cluster" is a reserved colname in "colData" to store cluster assignments in this function!')
  # Only one is accepted: assay.type = "logcounts" or use.dimred="PCA".
  if (!is.null(dimred) && dimred %in% c('PCA', 'UMAP', 'TSNE')) { 
    dat.gr = reducedDim(sce, dimred)
    if (dims < ncol(dat.gr)) dat.gr = dat.gr[, seq_len(dims), drop=FALSE]
    if (dims == 1) stop('The "dims" should be at least 2!')
  } else {
    if (!assay.type %in% assayNames(sce)) stop('The "assay.type" should be one of "assayNames(sce)"!')
    dat.gr=t(assays(sce)[[assay.type]])
  }  
  
  if (graph.meth=='knn') { 
    gr.sc <- do.call(bluster::makeKNNGraph, c(list(x=dat.gr), knn.gr))
  }
  if (graph.meth=='snn') {
    gr.sc <- do.call(bluster::makeSNNGraph, c(list(x=dat.gr), snn.gr))
  }
  # cluster: detected cell clusters. label: customer clusters.

  clus.all <- detect_cluster(graph=gr.sc, clustering=cluster, wt.arg=wt.arg, fg.arg=fg.arg, sl.arg=sl.arg, le.arg=le.arg, eb.arg=eb.arg)
  if (is.null(clus.all)) return()
  clus <- as.character(clus.all$membership)

  # clus <- as.character(do.call(cluster_walktrap, c(list(graph=gr.sc), cluster.wk))$membership)
  clus <- paste0('clus', clus)
  cdat.sc <- colData(sce); rna <- rownames(cdat.sc)
  lab.lgc <- 'label' %in% make.names(colnames(cdat.sc))
  if (lab.lgc) {
    cdat.sc <- cbind(cluster=clus, cdat.sc)
    idx <- colnames(cdat.sc) %in% c('cluster', 'label')
    cdat.sc <- cdat.sc[, c(which(idx), which(!idx))]
  } else cdat.sc <- cbind(cluster=clus, cdat.sc)
  # "cbind" removes row names in "cdat.sc". If "cdat.sc" has no row names, the existing column names in "sce" are erased.
  rownames(cdat.sc) <- rna
  colnames(cdat.sc) <- make.names(colnames(cdat.sc))
  colData(sce) <- cdat.sc; return(sce)
}
