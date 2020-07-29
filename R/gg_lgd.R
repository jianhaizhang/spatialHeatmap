#' Adjust Legend Key Size and Rows in Ggplot.
#'
#' @param gg.all A list of spatial heatmaps of ggplot.
#' @param size.key A numeric of legend key size. If \code{size.text} is NULL, it also applies to legend text size.
#' @param size.text A numeric of legend text size.
#' @param sub.title.size The title size of ggplot.
#' @param row An integer of rows in legend key.
#' @param sam.dat A vector of samples in the data matrix.
#' @param tis.trans A vector of tissues to be transparent.
#' @return A list of ggplots.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @importFrom ggplot2 theme layer_data scale_fill_manual unit element_text guide_legend

gg_lgd <- function(gg.all, size.key=0.02, size.text=8, sub.title.size=NULL, row=2, sam.dat, tis.trans=NULL) {

  for (i in seq_along(gg.all)) {
  
    g <- gg.all[[i]] 
    if (!is.null(size.key)) g <- g+theme(legend.key.size=unit(size.key, "npc"), legend.text=element_text(size=ifelse(is.null(size.text), 8*size.key*33, size.text)))
    if (!is.null(sub.title.size)) g <- g+theme(plot.title=element_text(hjust=0.5, size=sub.title.size))
    if (!is.null(row)) {

      lay.dat <- layer_data(g)
      dat <- g$data; g.col <- lay.dat$fill; names(g.col) <- dat$tissue
      g.col <- g.col[!duplicated(names(g.col))]; tis.path <- dat$feature
      sam.legend <- intersect(unique(sam.dat), unique(tis.path))
      sam.legend <- setdiff(sam.legend, tis.trans) 
      leg.idx <- !duplicated(tis.path) & (tis.path %in% sam.legend)
      gg.all[[i]] <- g+scale_fill_manual(values=g.col, breaks=as.vector(dat$tissue)[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=NULL, nrow=row))

    }

  }; return(gg.all)

}
