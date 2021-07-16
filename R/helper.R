#' Apply basic functions on data frame/matrix columns.
#'
#' @param data A data frame or matrix.
#' @param fun The function to apply on the column: \code{is.numeric}, \code{is.na}, \code{as.numeric}.
#' @return A matrix.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

df_is_as <- function(data, fun=is.numeric) {                                                                                     
  dim.na <- dimnames(data)
  # apply(na, 2, is.na): requires much more memory than vapply.                                                                                                     
  if (identical(fun, is.na)) vap <- vapply(seq_len(ncol(data)), function(i) { fun(data[, i]) }, logical(nrow(data))) else if (identical(fun, as.numeric)) vap <- vapply(seq_len(ncol(data)), function(i) { fun(data[, i]) }, numeric(nrow(data))) else if (identical(fun, is.numeric)) vap <- vapply(seq_len(ncol(data)), function(i) { fun(data[, i]) }, logical(1))                        
  # is.numeric always returns one value.
  if (nrow(data)==1 | identical(fun, is.numeric)) vap <- matrix(vap, byrow=TRUE, ncol=ncol(data))                                
  if (!identical(fun, is.numeric)) dimnames(vap) <- dim.na else { rownames(vap) <- dim.na[[1]][1]; colnames(vap) <- dim.na[[2]] }
  return(vap)
}    

#' Check validity of a function on data
#'
#' @param data A vector of single value.
#' @param fun The function to apply on the data: \code{is.numeric}, \code{is.na}, \code{as.numeric}.
#' @return A vector, 'warning', or 'error'.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

check <- function(data, fun) {                                                                                                  
  er.wa <- tryCatch({ fun(data) }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )            
  return(er.wa)                                                                                                                  
}
