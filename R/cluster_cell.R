#' Cluster single cells or combination of single cells and bulk
#'
#' In co-clustering, cluster only single cell data or combination of single cell and bulk data. The cluster assignments are stored in the \code{label} column of \code{colData} slot of \code{SingleCellExperiment}.
#' @param data The normalized single cell data or nomalized combination of single cell and bulk data at log2 scale in form of \code{SingleCellExperiment}, \code{dgCMatrix}, \code{matrix}, or \code{data.frame}.
#' @param prop Numeric scalar specifying the proportion of genes to report as highly variable genes (HVGs). The default is \code{0.1}.
#' @param min.dim,max.dim Integer scalars specifying the minimum (\code{min.dim}) and maximum (\code{max.dim}) number of (principle components) PCs to retain respectively in \code{\link[scran]{denoisePCA}}. The default is \code{min.dim=5}, \code{max.dim=50}.
#' @param pca Logical, if \code{TRUE} only the data with reduced dimentionality by PCA is returned and no clustering is performed. The default is \code{FALSE} and clustering is performed after dimensionality reduction.
#' @param graph.meth Method to build a nearest-neighbor graph, \code{snn} (see \code{\link[scran]{buildSNNGraph}}) or \code{knn} (default, see \code{\link[scran]{buildKNNGraph}}). The clusters are detected by first creating a nearest neighbor graph using \code{snn} or \code{kn} then partitioning the graph. 
#' @param dimred A string of \code{PCA} (default) or \code{UMAP} specifying which reduced dimensionality to use in creating a nearest neighbor graph. Internally, before building a nearest neighbor graph the data dimensionalities are reduced by PCA and UMAP respectively.

#' @return A list of normalized single cell and bulk data. 

#' @examples

#' library(scran); library(scuttle) 
#' sce <- mockSCE(); sce <- logNormCounts(sce)
#' # Modelling the variance.
#' var.stats <- modelGeneVar(sce)
#' sce <- denoisePCA(sce, technical=var.stats, subset.row=rownames(var.stats)) 
#' \donttest{
#' sce.clus <- cluster_cell(data=sce, prop=0.1, min.dim=5, max.dim=50, graph.meth='snn', dimred='PCA')
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
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.
#' Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. https://igraph.org


#' @export cluster_cell 
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment SingleCellExperiment colLabels colLabels<-
#' @importFrom scran modelGeneVar getTopHVGs denoisePCA buildSNNGraph buildSNNGraph 
#' @importFrom scater runUMAP runTSNE
#' @importFrom igraph cluster_walktrap

cluster_cell <- function(data, prop=0.1, min.dim=5, max.dim=50, pca=FALSE, graph.meth='knn', dimred='PCA') {
  if (is(data, 'SingleCellExperiment')) sce <- data
  if (is(data, 'dgCMatrix')|is(data, 'matrix')|is(data, 'data.frame')) {
    if (all(round(data)==data)) stop('The data to cluster should be in log2 scale!')
    sce <- SingleCellExperiment(list(logcounts=as.matrix(data)))
  }
  df.var.sc <- modelGeneVar(sce) # Use logcounts by default.
  top.hvgs.sc <- getTopHVGs(df.var.sc, prop=prop)
  max.pc <- (length(top.hvgs.sc) - 1) / 2
  # Avoid warning (the max returned PCS should be <= 50% of length(top.hvgs.sc) - 1).
  if (max.pc < max.dim | max.pc < min.dim) {
    top.hvgs.sc <- getTopHVGs(df.var.sc, prop=1)
    cat('"prop" is set 1 in "getTopHVGs" due to too less genes. \n')
    max.pc <- (length(top.hvgs.sc) - 1) / 2
    if (max.pc < max.dim | max.pc < min.dim) {
      warning('The input single cell data has too less genes!'); return()
    }
  }
  # prop=1 cannot avoid all warnings if the real max.dim does not increase even though prop=1.
  sce.dimred <- tryCatch(
    expr = {
      denoisePCA(sce, assay.type = "logcounts", technical=df.var.sc, subset.row=top.hvgs.sc, min.rank=min.dim, max.rank=max.dim)
    },
    warning = function(w){ 'w' }, error = function(e){ 'e' } 
  ) 
  if (!is(sce.dimred, 'SingleCellExperiment')) {
    message('min.dim = ', min.dim, ' and ', 'max.dim = ', max.dim, ' cannot be satisfied!'); return()
  }
  if (pca==TRUE) return(sce.dimred)
  # Other argument: n_dimred, ntop. By default only 2 dimensions are returned by runTSNE/runUMAP.
  # runUMAP returns different results before/after runTSNE.
  # Avoid warnings due to duplicated column names.
  cna <- colnames(sce.dimred); colnames(sce.dimred) <- seq_len(ncol(sce.dimred)) 
  sce.dimred <- runUMAP(sce.dimred, dimred="PCA", ncomponents=min.dim)
  # Row names in colData, reducedDim change accordingly.
  colnames(sce.dimred) <- cna
  sce.dimred <- runTSNE(sce.dimred, dimred="PCA", ncomponents=2)
  
  # dims <- list(pca='PCA', tsne='TSNE', umap='UMAP') 
  # Only one is accepted: assay.type = "logcounts" or use.dimred="PCA".
  if (graph.meth=='snn') gr.sc <- buildSNNGraph(sce.dimred, use.dimred=dimred)
  if (graph.meth=='knn') gr.sc <- buildKNNGraph(sce.dimred, use.dimred=dimred)
  colLabels(sce.dimred) <- as.character(cluster_walktrap(gr.sc)$membership)
  cdat.sc <- colData(sce.dimred); cdat.sc$cell <- rownames(cdat.sc)
  colData(sce.dimred) <- cdat.sc
  return(sce.dimred)
}


