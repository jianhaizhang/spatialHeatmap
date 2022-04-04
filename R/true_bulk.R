#' Assign true bulk to cells in \code{colData} slot. 
#'
#' In co-clustering, assign true bulk to cells in \code{colData} slot. 

#' @param sce A \code{SingleCellExperiment} of clustered single cell data.
#' @param df.match The matching table between cells and true bulk.

#' @return A \code{SingleCellExperiment} object.

#' @examples

#' See the example in the "cocluster" function by running "?cocluster".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.  0, https://bioconductor.org/packages/SummarizedExperiment.

#' @export true_bulk
#' @importFrom SummarizedExperiment colData

true_bulk <- function(sce, df.match) {
  true.bulk <- df.match$trueBulk; names(true.bulk) <- df.match$cell
  cdat <- colData(sce) 
  svg.bulk <- df.match$SVGBulk; names(svg.bulk) <- df.match$cell
  colData(sce)$trueBulk <- true.bulk[cdat$cell]
  colData(sce)$SVGBulk <- svg.bulk[cdat$cell]
  cdat <- colData(sce) 
  idx <- is.na(cdat$trueBulk) & is.na(cdat$SVGBulk)
  colData(sce)[idx, c('trueBulk', 'SVGBulk')] <- 'none'
  return(sce)
}




