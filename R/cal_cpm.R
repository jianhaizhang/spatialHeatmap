#' Normalize by CPM. 
#'
#' The log2 values are transformed to power of 2 (counts) and then to CPM. To maintain log2-scale values, the CPM values are then transformed back to log2. The returned \code{SingleCellExperiment} contains values of log2-transformed CPM. 

#' @param sce.nor The output of \code{\link[scrapper]{normalizeRnaCounts.se}} in form of \code{SingleCellExperiment}, where raw counts can be accessed by `counts(sce.nor)`.
#' @param sf A numeric vector containing size factors to adjust the library sizes. If NULL, the library sizes are used directly. See \code{\link[scrapper]{centerSizeFactors}}
#' @return A \code{SingleCellExperiment} object.
#' @keywords Internal
#' @noRd

#' @details
#' If size.factors are provided or available in x, they are used to define the effective library sizes. This is done by scaling all size factors such that the mean factor is equal to the mean sum of counts across all features. The effective library sizes are then used as the denominator of the CPM calculation.
#' 
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Douglas Bates and Martin Maechler (2021). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.4-0. https://CRAN.R-project.org/package=Matrix

#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom Matrix Matrix 
#' @importFrom SingleCellExperiment counts

cal_cpm <- function(sce.nor, sf = NULL) {
  
  ## x: matrix-like object, rows = genes/features, cols = samples/cells
  if (is(sce.nor, 'SummarizedExperiment')|is(sce.nor, 'SingleCellExperiment')) x <- as.matrix(counts(sce.nor))
  ## library sizes
  lib_sizes <- colSums(x)
  
  ## use raw library sizes if sf is NULL
  if (is.null(sf)) {
    effective_lib <- lib_sizes
  } else {
    if (length(sf) != ncol(x)) {
      stop("'sf' must have length equal to ncol(x)")
    }
    if (any(sf <= 0, na.rm = TRUE)) {
      stop("'sf' must contain positive values")
    }
    
    ## mimic scuttle behavior:
    ## scale size factors so that:
    ## mean(sf_scaled) == mean(library sizes)
    
    sf_scaled <- sf / mean(sf, na.rm = TRUE) *
      mean(lib_sizes, na.rm = TRUE)
    
    effective_lib <- sf_scaled
  }
  
  ## CPM
  cpm <- t(
    t(x) / effective_lib
  ) * 1e6
  
  logcounts(sce.nor) <- Matrix(as.matrix(log2(cpm+1)), sparse=TRUE)
  return(sce.nor)
}
