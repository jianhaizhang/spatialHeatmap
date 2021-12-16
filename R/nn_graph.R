#' Build a nearest-neighbor graph 
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param The method to build the graph, \code{SNN} or \code{KNN}. 
#' @param use.dimred The reduced dimentionality to use, such as \code{PCA}, \code{TSNE}, \code{UMAP}.
#' @param snn.arg A list of basic arguments passed to \code{buildSNNGraph}.
#' @param snn.arg.more A list of additional arguments passed to \code{buildSNNGraph}.
#' @param knn.arg A list of basic arguments passed to \code{buildKNNGraph}.
#' @param knn.arg.more A list of additional arguments passed to \code{buildKNNGraph}.

#' @return A graph where nodes are cells and edges represent connections between nearest neighbors.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2

#' @importFrom scran buildSNNGraph buildKNNGraph 

nn_graph <- function(sce, method='SNN', use.dimred=NULL, snn.arg=list(k=10, type=c("rank", "number", "jaccard")), snn.arg.more=list(), knn.arg=list(k=10), knn.arg.more=list()) {
  # assays(sce)[['logcounts']] <- as.matrix(assays(sce)[['logcounts']])
  cat('Scell: nearest neighbor graph ... \n')
  if (method=='SNN') {
    # Only one is accepted: assay.type = "logcounts" or use.dimred="PCA".
    g <- do.call('buildSNNGraph', c(list(x=sce, use.dimred=use.dimred), snn.arg, snn.arg.more))
  } else if (method=='KNN') {
    g <- do.call('buildKNNGraph', c(list(x=sce, use.dimred=use.dimred), knn.arg, knn.arg.more)) 
  }; return(g)
}
