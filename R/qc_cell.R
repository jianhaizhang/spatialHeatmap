#' Quality control in single cell data
#'
#' A meta function for quality control in single-cell RNA-seq data.
#' @param sce Raw single cell count data in form of \code{SingleCellExperiment}.
#' @param qc.metric Quality control arguments in a named list passed to \code{\link[scrapper]{quickRnaQc.se}}. Eg: `list(assay.type = "counts", subsets=list(mito=grepl("^mt", rownames(sce))), altexp.proportions="ERCC")`.
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


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.
#' Lun ATL, Kancherla J (2023). "Powering single-cell analyses in the browser with WebAssembly." _Journal of Open Source Software_, *8*(89), 5603. doi:10.21105/joss.05603 <https://doi.org/10.21105/joss.05603>.

#' @export qc_cell
#' @importFrom SingleCellExperiment altExpNames 

qc_cell <- function(sce, qc.metric=list(assay.type = "counts", subsets=list(), altexp.proportions=NULL, block = NULL)) {
  # save(sce, qc.metric, file='qc_cell.arg')
  pkg <- check_pkg('scrapper'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  
  do.call(scrapper::quickRnaQc.se, c(list(x=sce), qc.metric))

  # Combine main dataset in assay and all datasets in altExpNames, and return index.
  # lis <- combine_sce_with_altexps(sce=sce, assay_name = NULL)
  # Combined datasets
  # dat = lis$combined_matrix
  # Index for each dataset.
  # idx = lis$dataset_index; idx = idx[!names(idx) %in% 'main']
  # QC
  # qc=do.call(scrapper::computeRnaQcMetrics, c(list(x=dat, subsets=idx), qc.metric))
  # filt=do.call(scrapper::suggestRnaQcThresholds, c(list(metrics=qc), qc.sug))
  # keep=do.call(scrapper::filterRnaQcMetrics, c(list(thresholds=filt, metrics=qc), qc.filter))
  # Filtered sce
  # sce[, keep]
}
