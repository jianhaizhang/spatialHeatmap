#' Assign true bulk to cells in \code{colData} slot. 
#'
#' In co-clustering, assign true bulk to cells in \code{colData} slot. 

#' @param sce A \code{SingleCellExperiment} of clustered single cell data.
#' @param df.match The matching table between cells and true bulk.

#' @return A \code{SingleCellExperiment} object.

#' @examples

#' # Matching table.
#' match.mus.brain.pa <- system.file("extdata/shinyApp/example", "match_mouse_brain_cocluster.txt", package="spatialHeatmap")
#' df.match.mus.brain <- read.table(match.mus.brain.pa, header=TRUE, row.names=1, sep='\t')
#' df.match.mus.brain
#' 
#' # Create random data matrix.
#' df.random <- matrix(rexp(30), nrow=5)
#' dimnames(df.random) <- list(paste0('gene', seq_len(nrow(df.random))), c('cere', 'cere', 'hipp', 'hipp', 'corti.sub', 'corti.sub'))
#' 
#' library(SingleCellExperiment); library(S4Vectors)
#' cell.refined <- SingleCellExperiment(assays=list(logcounts=df.random), colData=DataFrame(cell=colnames(df.random)))
#' 
#' #cell.refined <- true_bulk(cell.refined, df.match.mus.brain)
#' #colData(cell.refined)
#' 
#' # See detailed example in the "coclus_meta" function by running "?coclus_meta".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.  0, https://bioconductor.org/packages/SummarizedExperiment.

#' @export true_bulk
#' @importFrom SummarizedExperiment colData

true_bulk <- function(sce, df.match) {
  if (!'cell' %in% colnames(df.match)) stop('The "cell" column is missing in the matching table!')
  if (!'SVGBulk' %in% colnames(df.match)) df.match$SVGBulk <- 'none'
  true.bulk <- df.match$dataBulk; names(true.bulk) <- df.match$cell
  cdat <- colData(sce) 
  svg.bulk <- df.match$SVGBulk; names(svg.bulk) <- df.match$cell
  colData(sce)$dataBulk <- as.character(true.bulk[cdat$cell])
  colData(sce)$SVGBulk <- as.character(svg.bulk[cdat$cell])
  cdat <- colData(sce) 
  idx <- is.na(cdat$dataBulk) & is.na(cdat$SVGBulk)
  colData(sce)[idx, c('dataBulk', 'SVGBulk')] <- 'none'
  return(sce)
}




