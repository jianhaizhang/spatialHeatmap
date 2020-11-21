#' Edit Targets Files
#'
#' Replace existing entries in a chosen column of a targets file with desired ones.

#' @param df.tar The data frame of a targets file.
#' @param column The column to edit, either the column name or an integer of the column index. 
#' @param old A vector of existing entries to replace, where the length must be the same with \code{new}.
#' @param new A vector of desired entries to replace that in \code{old}, where each entry corresponds to a counterpart in \code{old} respectively.
#' @param sub.row A vector of integers corresponding to target rows for editing, or a vector of TRUE and FALSE corresponding to each row. Default is all rows in the targets file.
#' @return A data frame.

#' @examples

#' sh.tar <- system.file('extdata/shinyApp/example/target_arab.txt', package='spatialHeatmap')
#' target.sh <- read_fr(sh.tar)
#' target.sh.new <- edit_tar(df.tar=target.sh, column='condition', old=c('control', 'hypoxia'), new=c('C', 'H'), sub.row=c(1:12))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @export edit_tar

edit_tar <- function(df.tar, column, old, new, sub.row) {

  if (length(old)!=length(new)) stop('The number of entries in "old" must be the same with "new"!')
  
  if (missing(sub.row)) sub.row <- seq_len(nrow(df.tar))
  df.tar[, column] <- as.vector(df.tar[, column])
  clm <- df.tar[sub.row, column]
  for (i in seq_along(old)) { clm[clm %in% old[i]] <- new[i] }
  df.tar[sub.row, column] <- clm; return(df.tar)

}

