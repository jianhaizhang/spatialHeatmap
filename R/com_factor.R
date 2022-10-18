#' Combine Factors in Targets File 
#'
#' This is a helper function for data/aSVGs involving three or more factors such as sample, time, condition. It combine factors in targets file to make composite factors.

#' @param se A \code{SummarizedExperiment} object.
#' @param target A \code{data.frame} object of targets file. 
#' @param factors2com A character vector of column names or a numeric vector of column indeces in the targets file. Entries in these columns are combined. 
#' @param sep The separator in the combined factors. One of \code{_}, and \code{.} (default).
#' @param factor.new The column name of the new combined factors.
#' @return If \code{se} is provided, a \code{SummarizedExperiment} object is returned, where the \code{colData} slot contains the new column of combined factors. Otherwise, a\code{data.frame} object is returned, where the new column of combined factors is appended.

#' @examples

#' library(SummarizedExperiment)
#' mus.se.pa <- system.file('extdata/shinyApp/example/mus_brain_vars_se.rds', package='spatialHeatmap')
#' mus.se <- readRDS(mus.se.pa); targets.info <- colData(mus.se)
#' targets.new <- com_factor(target=targets.info, factors2com=c('time', 'treatment', 'injury'), factor.new='comDim')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Attilio, Peter J, Dustin M Snapper, Milan Rusnak, Akira Isaac, Anthony R Soltis, Matthew D Wilkerson, Clifton L Dalgard, and Aviva J Symes. 2021. “Transcriptomic Analysis of Mouse Brain After Traumatic Brain Injury Reveals That the Angiotensin Receptor Blocker Candesartan Acts Through Novel Pathways.” Front. Neurosci. 15 (March): 636259

#' @export com_factor
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame

com_factor <- function(se, target, factors2com, sep='.', factor.new) {
  target <- as.data.frame(target)
  target$new <- as.character(lapply(seq_len(nrow(target)), function(x) paste0(target[x, factors2com], collapse=sep)))
  colnames(target)[colnames(target)=='new'] <- factor.new
  cat('New combined factors:', unique(target[, factor.new]), '\n')
  if (missing(se)) return(target)
  # Keep columns present in colData but not in target.
  cdat <- colData(se); cdat <- cdat[, !colnames(cdat) %in% colnames(target), drop=FALSE]
  colData(se) <- cbind(DataFrame(target), cdat)
  return(se)
}
