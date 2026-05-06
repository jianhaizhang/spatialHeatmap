#' Reducing dimensionality in count data
#'
#' A meta function for reducing dimensionality in count data.

#' @param sce Normalized single cell data in \code{SingleCellExperiment} returned by \code{norm_cell}. Alternative forms include \code{dgCMatrix}, \code{matrix}, \code{data.frame}.

#' @param choose.hvg Arguments in a named list passed to \code{\link[scrapper]{chooseRnaHvgs.se}}. 
#' @param min.dim,max.dim Integer scalars specifying the minimum (\code{min.dim}) and maximum (\code{max.dim}) number of (principle components) PCs to retain respectively in \code{\link[scrapper]{runPca.se}}. 
#' @param pca.arg Arguments in a named list passed to \code{\link[scrapper]{runPca.se}}.
#' @param umap.arg Arguments in a named list passed to \code{\link[scrapper]{runUmap.se}}. 
#' @param tsne.arg Arguments in a named list passed to \code{\link[scrapper]{runTsne.se}}. 
#' @param pca Logical, if \code{TRUE} only the data with reduced dimentionality of PCA is returned. The default is \code{FALSE}, and UMAP and TSNE are also returned. 

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

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods,    17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Lun ATL, Kancherla J (2023). "Powering single-cell analyses in the browser with WebAssembly." _Journal of Open Source Software_, *8*(89), 5603. doi:10.21105/joss.05603 <https://doi.org/10.21105/joss.05603>.

#' @export reduce_dim 
#' @importFrom SingleCellExperiment SingleCellExperiment 

reduce_dim <- function(sce, choose.hvg=list(assay.type = "logcounts", block = NULL,
                                            more.var.args=list(transform = TRUE),
                                            top = 4000, more.choose.args=list(keep.ties = TRUE)),
                       min.dim=11, max.dim=50,
                       pca.arg = list(assay.type = "logcounts", block = NULL),
                       umap.arg = list(num.dim = 2, reddim.type = "PCA"),
                       tsne.arg = list(perplexity = 30, reddim.type = "PCA", output.name = "TSNE"), 
                       pca=FALSE) {
  # getTopHVGs: prop=0.1 is super more important than n=3000 in co-clustering.
  # save(sce, choose.hvg, min.dim, max.dim, pca.arg, umap.arg, tsne.arg, pca, file='reduce_dim.arg')
  pkg <- check_pkg('scrapper'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  if (is(sce, 'dgCMatrix')|is(sce, 'matrix')|is(sce, 'data.frame')) {
    if (all(round(sce)==sce)) stop('The "sce" should be in log2 scale!')
    sce <- SingleCellExperiment(list(logcounts=as.matrix(sce)))
  }
  # Use logcounts by default.
  message('Log-expression values are expected.')
  # stats <- do.call(scrapper::modelGeneVariances, c(list(x=logcounts(sce)), model.var))
  # stats <- modelGeneVariances(logcounts(sce))
  # out <- fitVarianceTrend(stats$statistics$means, stats$statistics$variances)
  # hvg <- do.call(scrapper::chooseHighlyVariableGenes, c(list(stats=stats$statistics$residuals), choose.hvg))
  # hvg <- chooseHighlyVariableGenes(stats$statistics$residuals, top = 4000)
  
  sce <- do.call(scrapper::chooseRnaHvgs.se, c(list(x=sce), choose.hvg))
  if (max.dim < min.dim) max.dim = min.dim
  if (min.dim < 2) stop('"min.dim" should be at least 2!')
  pca.arg$number=min.dim
  sce.dimred <- do.call(scrapper::runPca.se, c(list(x=sce, features=rowData(sce)$hvg), pca.arg))
  
  if (pca==TRUE) return(sce.dimred)
  # Other argument: n_dimred, ntop. By default only 2 dimensions are returned by runTSNE/runUMAP.
  # runUMAP returns different results before/after runTSNE.
  # Avoid warnings due to duplicated column names.
  cna <- colnames(sce.dimred); colnames(sce.dimred) <- seq_len(ncol(sce.dimred)) 
  sce.dimred <- do.call(scrapper::runUmap.se, c(list(x=sce.dimred), umap.arg))
  # Row names in colData, reducedDim change accordingly.
  colnames(sce.dimred) <- cna
  sce.dimred <- do.call(scrapper::runTsne.se, c(list(x=sce.dimred), tsne.arg))
  return(sce.dimred)
}
