#' Build a nearest-neighbor graph 
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param The method to build the graph, \code{SNN} (\code{\link[bluster]{makeSNNGraph}}) or \code{KNN} (\code{\link[bluster]{makeKNNGraph}}). 
#' @param use.dimred The reduced dimentionality to use, such as \code{PCA}, \code{TSNE}, \code{UMAP}.
#' @param snn.arg A list of basic arguments passed to \code{\link[bluster]{makeSNNGraph}}.
#' @param snn.arg.more A list of additional arguments passed to \code{\link[bluster]{makeSNNGraph}}.
#' @param knn.arg A list of basic arguments passed to \code{\link[bluster]{makeKNNGraph}}.
#' @param knn.arg.more A list of additional arguments passed to \code{\link[bluster]{makeKNNGraph}}.
#' @inheritParams cluster_cell

#' @return A graph where nodes are cells and edges represent connections between nearest neighbors.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Lun A (2026). _bluster: Clustering Algorithms for Bioconductor_. doi:10.18129/B9.bioc.bluster <https://doi.org/10.18129/B9.bioc.bluster>. R package version 1.22.0, <https://bioconductor.org/packages/bluster>.

nn_graph <- function(sce, method='SNN', use.dimred=NULL, dims=50, assay.type = NULL, snn.arg=list(k=10, type=c("rank", "number", "jaccard")), 
                     snn.arg.more=list(), knn.arg=list(k=10), knn.arg.more=list()) {
  # assays(sce)[['logcounts']] <- as.matrix(assays(sce)[['logcounts']])
  cat('Scell: nearest neighbor graph ... \n')
  pkg <- check_pkg('bluster'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
 
  if (!is.null(use.dimred) && use.dimred %in% c('PCA', 'UMAP', 'TSNE')) { 
    dat.gr = reducedDim(sce, use.dimred)
    if (dims < ncol(dat.gr)) dat.gr = dat.gr[, seq_len(dims), drop=FALSE]
    if (dims == 1) stop('The "dims" should be at least 2!')
  } else {
    if (!assay.type %in% assayNames(sce)) stop('The "assay.type" should be one of "assayNames(sce)"!')
    dat.gr=t(assays(sce)[[assay.type]])
  } 
  
  if (method=='SNN') {
    # Only one is accepted: assay.type = "logcounts" or use.dimred="PCA".
    g <- do.call(bluster::makeSNNGraph, c(list(x=dat.gr), snn.arg, snn.arg.more))
  } else if (method=='KNN') {
    g <- do.call(bluster::makeKNNGraph, c(list(x=dat.gr), knn.arg, knn.arg.more))
  }; return(g)
}
