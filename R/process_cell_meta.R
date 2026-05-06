#' Processing single cell RNA-seq count data 
#'
#' A meta function for processing single cell RNA-seq count data, including quality control, normalization, dimensionality reduction.

#' @param sce Single cell RNA-seq count data in \code{SingleCellExperiment}.
#' @inheritParams qc_cell
#' @inheritParams norm_cell
#' @inheritParams reduce_dim

#' @return A \code{SingleCellExperiment} object. 
#' @details
#' In the QC, compute per-cell QC metrics from an initialized matrix of RNA counts, and use the metrics to suggest filter thresholds to retain high-quality cells. Refer to \code{filterRnaQcMetrics} in the `scrapper` package for more details. 
#' In the normalization, compute log-normalized expression values after performing scaling normalization of an RNA count matrix. See more details in \code{normalizeRnaCounts.se} from the scrapper package.
#' In dimensionality reduction, the high-dimensional gene expression data are embedded into a low dimensional space using PCA, tSNE and UMAP. All three embedding result sets are stored in a \code{SingleCellExperiment} object. Details are seen in \code{runPca.se}, \code{runUmap.se}, and \code{runTsne.se} from `scrapper`. 

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

#' sce.dimred <- process_cell_meta(sce)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods,    17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation  of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.

#' @export process_cell_meta 

process_cell_meta <- function(sce, qc.metric=list(assay.type = "counts", subsets=list(), altexp.proportions=NULL, block = NULL), 
                              center.sf=list(block = NULL, mode = c("lowest", "per-block")),
                              log.norm=list(assay.type = "counts", log = TRUE, pseudo.count = 1),
                              choose.hvg=list(assay.type = "logcounts", block = NULL,
                                              more.var.args=list(transform = TRUE),
                                              top = 4000, more.choose.args=list(keep.ties = TRUE)),
                              min.dim=5, max.dim=25,
                              pca.arg = list(assay.type = "logcounts", number = 25, block = NULL),
                              umap.arg = list(num.dim = 2, reddim.type = "PCA"),
                              tsne.arg = list(perplexity = 30, reddim.type = "PCA", output.name = "TSNE"),
                              pca=FALSE
                              ) {
  # save(sce, qc.metric, center.sf, log.norm, choose.hvg, min.dim, max.dim, pca.arg, umap.arg, tsne.arg, pca, file='process_cell_meta.arg')
  # Quality control.
  sce.qc <- qc_cell(sce=sce, qc.metric=qc.metric)
  # Normalization.
  sce.norm <- norm_cell(sce=sce.qc, center.sf=center.sf, log.norm=log.norm, com=TRUE)
  # Dimensionality reduction.
  sce.dimred <- reduce_dim(sce=sce.norm, choose.hvg=choose.hvg, min.dim=min.dim, max.dim=max.dim, pca.arg=pca.arg, umap.arg=umap.arg, tsne.arg=tsne.arg, pca=pca)
  return(sce.dimred)
}
