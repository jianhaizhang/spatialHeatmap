#' Cluster single cells or combination of single cells and bulk
#'
#' Cluster only single cell data or combination of single cell and bulk data. The cluster assignments are stored in the \code{cluster} column of \code{colData} slot of \code{SingleCellExperiment}.
#' @param data The normalized single cell data or nomalized combination of single cell and bulk data at log2 scale in form of \code{SingleCellExperiment}, \code{dgCMatrix}, \code{matrix}, or \code{data.frame}.
#' @param prop Numeric scalar specifying the proportion of genes to report as highly variable genes (HVGs) in \code{\link[scran]{getTopHVGs}}. The default is \code{0.1}.
#' @param min.dim,max.dim Integer scalars specifying the minimum (\code{min.dim}) and maximum (\code{max.dim}) number of (principle components) PCs to retain respectively in \code{\link[scran]{denoisePCA}}. The default is \code{min.dim=13}, \code{max.dim=50}.
#' @param pca Logical, if \code{TRUE} only the data with reduced dimentionality by PCA is returned and no clustering is performed. The default is \code{FALSE} and clustering is performed after dimensionality reduction.
#' @param graph.meth Method to build a nearest-neighbor graph, \code{snn} (see \code{\link[scran]{buildSNNGraph}}) or \code{knn} (default, see \code{\link[scran]{buildKNNGraph}}). The clusters are detected by first creating a nearest neighbor graph using \code{snn} or \code{kn} then partitioning the graph. 
#' @param dimred A string of \code{PCA} (default) or \code{UMAP} specifying which reduced dimensionality to use in creating a nearest neighbor graph. Internally, before building a nearest neighbor graph the data dimensionalities are reduced by PCA and UMAP respectively.
#' @param model.var Additional arguments in a named list passed to \code{\link[scran]{modelGeneVar}}.
#' @param top.hvg Additional arguments in a named list passed to \code{\link[scran]{getTopHVGs}}, such as \code{top.hvg=list(n = 3000)}. 
#' @param de.pca Additional arguments in a named list passed to \code{\link[scran]{denoisePCA}}, such as \code{de.pca=list(assay.type = "logcounts")}.
#' @param tsne Additional arguments in a named list passed to \code{\link[scater]{runTSNE}}, such as \code{tsne=list(dimred="PCA", ncomponents=2)}.
#' @param umap Additional arguments in a named list passed to \code{\link[scater]{runUMAP}}, such as \code{umap=list(dimred="PCA")}.
#' @param knn.gr Additional arguments in a named list passed to \code{\link[scran]{buildKNNGraph}}.
#' @param snn.gr Additional arguments in a named list passed to \code{\link[scran]{buildSNNGraph}}.
#' @param cluster.wk Additional arguments in a named list passed to \code{\link[igraph]{cluster_walktrap}}, such as \code{cluster.wk=list(steps = 4)}.

#' @return A list of normalized single cell and bulk data. 

#' @examples

#' library(scran); library(scuttle) 
#' sce <- mockSCE(); sce <- logNormCounts(sce)
#' # Modelling the variance.
#' var.stats <- modelGeneVar(sce)
#' sce <- denoisePCA(sce, technical=var.stats, subset.row=rownames(var.stats)) 
#' \donttest{
#' sce.clus <- cluster_cell(data=sce, prop=0.1, min.dim=13, max.dim=50, graph.meth='snn', dimred='PCA')
#' # Clusters.
#' table(colData(sce.clus)$label)
#'
#' }
#' # See details in function "cocluster" by running "?cocluster".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.
#' Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. https://igraph.org


#' @export cluster_cell 
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment SingleCellExperiment 
#' @importFrom scran buildSNNGraph buildSNNGraph 
#' @importFrom igraph cluster_walktrap

cluster_cell <- function(data, prop=0.1, min.dim=13, max.dim=50, pca=FALSE, graph.meth='knn', dimred='PCA', model.var=list(), top.hvg=list(n = 3000), de.pca=list(assay.type = "logcounts"), tsne=list(dimred="PCA", ncomponents=2), umap=list(dimred="PCA"), knn.gr=list(), snn.gr=list(), cluster.wk=list(steps = 4)) {
  if (is(data, 'SingleCellExperiment')) sce <- data
  if (is(data, 'dgCMatrix')|is(data, 'matrix')|is(data, 'data.frame')) {
    if (all(round(data)==data)) stop('The data to cluster should be in log2 scale!')
    sce <- SingleCellExperiment(list(logcounts=as.matrix(data)))
  }
  if ('cluster' %in% colnames(colData(sce))) stop('The "cluster" is a reserved colname in "colData" to store cluster assignments in this function!')
  
  sce.dimred <- reduce_dim(sce=sce, prop=prop, min.dim=min.dim, max.dim=max.dim, model.var=model.var, top.hvg=top.hvg, de.pca=de.pca, pca=pca, tsne=tsne, umap=umap)

  # dims <- list(pca='PCA', tsne='TSNE', umap='UMAP') 
  # Only one is accepted: assay.type = "logcounts" or use.dimred="PCA".
  if (graph.meth=='knn') { 
    gr.sc <- do.call(buildKNNGraph, c(list(x=sce.dimred, use.dimred=dimred), knn.gr))
  }
  if (graph.meth=='snn') {
    gr.sc <- do.call(buildSNNGraph, c(list(x=sce.dimred, use.dimred=dimred), snn.gr))
  }
  # cluster: detected cell clusters. label: customer clusters.
  clus <- as.character(do.call(cluster_walktrap, c(list(graph=gr.sc), cluster.wk))$membership)
  clus <- paste0('clus', clus)
  cdat.sc <- colData(sce.dimred); cdat.sc$cell <- rownames(cdat.sc)
  lab.lgc <- 'label' %in% make.names(colnames(cdat.sc))
  if (lab.lgc) {
    cdat.sc <- cbind(cluster=clus, cdat.sc)
    idx <- colnames(cdat.sc) %in% c('cluster', 'label')
    cdat.sc <- cdat.sc[, c(which(idx), which(!idx))]
  } else cdat.sc <- cbind(cluster=clus, cdat.sc)
  # "cbind" removes row names in "cdat.sc". If "cdat.sc" has no row names, the existing column names in "sce.dimred" are erased.
  rownames(cdat.sc) <- cdat.sc$cell
  colnames(cdat.sc) <- make.names(colnames(cdat.sc))
  colData(sce.dimred) <- cdat.sc; return(sce.dimred)
}
