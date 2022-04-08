#' Subset the SingleCellExperiment
#'
#' @param sce The \code{SingleCellExperiment} to subset.
#' @param mat The processed matrix (filtered, aggregated) in \code{assays} slot. Its column names will replace column names in \code{sce}, so it is better to be in "sample__condition" format.
#' @param cna The column names of \code{assay} slot. It is required if \code{mat} or \code{cdat} is assigned a value, since when assign \code{mat} to \code{assay} slot the two items should have the same column names, and \code{cdat} will erase existing column names of the \code{assay} slot.
#' @param col.idx A vector used to subset \code{sce} by columns. If \code{mat} is aggregated by replicates in columns, this argument is the vector indicating the first non-duplicate replicate.
#' @param row.idx A vector to subset \code{sce} by rows.
#' @param cdat The data frame in \code{colData} after subsetting, which will be added to the final \code{sce}.
#' @param rdat The data frame in \code{rowData} after subsetting, which will be added to the final \code{sce}.
#' @inheritParams filter_data


#' @return A \code{SingleCellExperiment} object.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment

#' @importFrom SummarizedExperiment assays colData rowData rowData<-

sce_sub <- function(sce, assay.na=NULL, mat=NULL, cna=NULL, col.idx=NULL, row.idx=NULL, cdat=NULL, rdat=NULL) {
  if (!is.null(mat) | !is.null(cdat)) if (is.null(cna)) stop("The 'cna' is required if 'cdat' is assigned a value!")
  if (is.null(assay.na)) assay.na <- 1
  if (!is.null(col.idx)) sce <- sce[, col.idx]
  if (!is.null(row.idx)) sce <- sce[row.idx, ]

  if (!is.null(mat)) {
    colnames(sce) <- colnames(mat) <- cna
    # assays(sce) and mat should have same column names.
    assays(sce)[[assay.na]] <- mat
  }
  if (!is.null(cdat)) {
    # Erase colnames in assay.
    colData(sce) <- as(cdat, 'DataFrame'); colnames(sce) <- cna
  }
  if (!is.null(rdat)) {
    # Rownames are not erased.
    rowData(sce) <- as(rdat, 'DataFrame')
  }; return(sce)
}




