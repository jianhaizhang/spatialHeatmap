#' Quality control in single cell data
#'
#' A meta function for quality control in single-cell RNA-seq data.
#' @param sce Raw single cell count data in form of \code{SingleCellExperiment}.
#' @param qc.metric Quality control arguments in a named list passed to \code{\link[scuttle]{perCellQCMetrics}}, such as \code{qc.metric=list(threshold=1)}.
#' @param qc.filter Quality control filtering arguments in a named list passed to \code{\link[scuttle]{perCellQCFilters}}, such as \code{qc.filter=list(nmads=3)}. 

#' @return A \code{SingleCellExperiment} object. 

#' @examples

#' library(scran); library(scuttle) 
#' sce <- mockSCE()
#' qc_cell(sce, qc.metric=list(subsets=list(Mt=rowData(sce)$featureType=='mito'), threshold=1))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.


#' @export qc_cell
#' @importFrom SingleCellExperiment altExpNames 
#' @importFrom scuttle perCellQCMetrics perCellQCFilters

qc_cell <- function(sce, qc.metric=list(threshold=1), qc.filter=list(nmads=3)) {
  # Quality control.
  stats <- do.call(perCellQCMetrics, c(list(x=sce), qc.metric))
  # Discard unreliable cells.
  sub.fields <- 'subsets_Mt_percent'
  ercc <- 'ERCC' %in% altExpNames(sce)
  if (ercc) sub.fields <- c('altexps_ERCC_percent', sub.fields)
  qc <- do.call(perCellQCFilters, c(list(x=stats, sub.fields=sub.fields), qc.filter))
  sce <- sce[, !qc$discard]; return(sce)
}
