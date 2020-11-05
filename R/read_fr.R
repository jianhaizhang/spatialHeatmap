#' Import Data from Tabular Files
#'
#' This function reads data from a tabular file, which is a wrapper of \link[data.table]{fread}. If the tabular file contains both character and numeric columns, it is able to maintain the character or numeric attribute for each column in the returned data frame. In addition, it is able to detect separators automatically.

#' @param input The file path.
#' @param header One of TRUE, FALSE, or "auto". Default is TRUE. Does the first data line contain column names, according to whether every non-empty field on the first data line is type character? If "auto" or TRUE is supplied, any empty column names are given a default name. 
#' @inheritParams data.table::fread
#' @param fill Logical (default is TRUE). If TRUE then in case the rows have unequal length, blank fields are implicitly filled.

#' @return A data frame.

#' @examples

#' sh.tar <- system.file('extdata/shinyApp/example/target_arab.txt', package='spatialHeatmap')
#' target.sh <- read_fr(sh.tar); target.sh[60:63, ]

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Matt Dowle and Arun Srinivasan (2019). data.table: Extension of `data.frame`. R package version 1.12.8. https://CRAN.R-project.org/package=data.table

#' @export read_fr
#' @importFrom data.table fread

read_fr <- function(input, header=TRUE, sep='auto', fill=TRUE, check.names=FALSE) {
  
  options(stringsAsFactors=FALSE)
  df0 <- fread(input=input, header=header, sep=sep, fill=fill, check.names=check.names)
  cna <- make.names(colnames(df0))
  if (cna[1]=='V1') cna <- cna[-1] else cna <- cna[-ncol(df0)] 
  df1 <- as.data.frame(df0); rownames(df1) <- df1[, 1]
  # Subsetting identical column names in a matrix will not trigger numbers to be appened.
  df1 <- df1[, -1]; colnames(df1) <- cna; rna <- rownames(df1)
  # Applies to data frame of numeric-character mixture.
  na <- vapply(seq_len(ncol(df1)), function(i) { tryCatch({ as.numeric(df1[, i]) }, warning=function(w) { return(rep(NA, nrow(df1)))
    }) }, FUN.VALUE=numeric(nrow(df1)) )
  na <- as.data.frame(na); rownames(na) <- rna
  idx <- colSums(apply(na, 2, is.na))!=0
  row.meta <- df1[idx]; expr <- na[!idx]
  colnames(expr) <- cna[!idx]; return(cbind(expr, row.meta))

}




