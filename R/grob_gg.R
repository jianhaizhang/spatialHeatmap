#' Separate SHMs of grobs from ggplot
#'
#' Separate SHMs of grobs from ggplot. Different SHMs of same 'gene_condition' are indexed with suffixed of '_1', '_2', ...
#' @param gs A list of SHMs of grobs, legend plot, and SHMs of ggplot, returned by \code{\link{grob_list}}.

#' @return A list of grob SHMs and ggplot SHMs.  
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

grob_gg <- function(gs) {
 
  # One SVG and multiple SVGs are treated same way. The list of gs slots are named with SVG file names.
  grob.all <- gg.all <- lgd.all <- NULL; idx <- seq_along(gs)
  for (i in idx) {
    
    grob0 <- gs[[i]][['grob.lis']]; names(grob0) <- paste0(names(grob0), '_', i)
    grob.all <- c(grob.all, grob0); gg0 <- gs[[i]][['g.lis.all']]
    names(gg0) <- names(grob0); gg.all <- c(gg.all, gg0)
    lgd.all <- c(lgd.all, list(gs[[i]][['g.lgd']]))

  }; names(lgd.all) <- names(gs)
  return(list(grob=grob.all, gg=gg.all, lgd.all=lgd.all))

}
