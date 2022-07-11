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

#' clp.tar <- system.file('extdata/shinyApp/example/target_coleoptile.txt', package='spatialHeatmap')
#' target.clp <- read_fr(clp.tar)
#' target.clp <- com_factor(target=target.clp, factors2com=c('organism_part', 'age'), factor.new='samTime')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Narsai, Reena, David Secco, Matthew D Schultz, Joseph R Ecker, Ryan Lister, and James Whelan. 2017. "Dynamic and Rapid Changes in the Transcriptome and Epigenome During Germination and in Developing Rice (Oryza Sativa) Coleoptiles Under Anoxia and Re-Oxygenation." Plant J. 89 (4): 805–24

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
