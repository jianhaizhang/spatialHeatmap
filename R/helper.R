
#' Check availability of packages.
#' @keywords Internal
#' @noRd

check_pkg <- function(x) { 
  tryCatch(
    utils::packageVersion(x), error = function(e){'e'},
    warning = function(w){'w'} 
  )  
} 

#' Check validity of an object
#'
#' @param x An object.

#' @return TRUE or FALSE
#' @keywords Internal
#' @noRd

check_obj <- function(x) {
  if (length(x)==0) return(FALSE); if (is.na(x)) return(FALSE)
  return(TRUE)
}

#' Shown popup window
#'
#' @param msg The main content to show.

#' @return A pop-up window in Shiny app.
#' @keywords Internal
#' @noRd

#' @importFrom shiny modalDialog span tagList modalButton

modal <- function(title = NULL, msg, easyClose=FALSE) {
  modalDialog(title = title, span(msg),
    footer = tagList(modalButton("Dismiss")), size = c("m"), easyClose=easyClose
  )
}

#' Check if SVGs are paired with templates of raster images
#'
#' @param svg.na SVG names or paths.
#' @param raster.ext A vector of raster image file extensions such as \code{c('.jpg', '.JPG', '.png', '.PNG')}.
#' @param shiny If used in Shiny app, the value should be \code{TRUE}.

#' @return The raster path.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @importFrom shiny showModal

svg_raster <- function(svg.na, raster.ext, shiny=TRUE) {                                                                                           
  ext <- paste0('\\', c('.svg', raster.ext), '$', collapse='|')
  # Ensure the raster and paired SVG names are same except for extensions.
  paired <- all(table(sub(ext, '', svg.na))==2)
  if (!paired) {
    msg <- 'If raster images are provided as templates, ensure they have the same names with corresponding aSVG images. E.g. test.png, test.svg.'
    if (shiny==TRUE) {
      showModal(modal(msg=msg)); validate(need(try(paired), ''))                                                                                            
    } else return(msg)
  } else return()
}

#' Extract tmplate path for a selected SVG
#'
#' @param svg.paths,svg.nas All svg/raster paths and names, which are in the same order.
#' @param svg.na Selected SVG name.
#' @return The raster path.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

raster_path <- function(svg.paths, svg.nas, svg.na, raster.ext){
  raster.pat <- sub('\\.svg', paste0('(', paste0('\\', raster.ext, '$', collapse='|'), ')'), svg.na)
  raster.pa <- svg.paths[grepl(raster.pat, svg.nas)]; return(raster.pa)
}

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

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

dat_fun <- function(data, fun) {                                                                                                  
  er.wa <- tryCatch({ fun(data) }, error=function(e){ return('e') }, warning=function(w) { return('w') } )            
  return(er.wa)                                                                                                                  
}
