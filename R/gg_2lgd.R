#' Add Secondary Legend
#'
#' @param gg.all A list of spatial heatmaps of ggplot.
#' @param sam.dat A vector of samples in the data matrix.
#' @param ft.trans A vector of tissues to be transparent.
#' @param position.2nd The position of the secondary legend. One of "top", "right", "bottom", "left", or a two-component numeric vector. The default is "bottom". Applies to the static image and video.
#' @param legend.nrow.2nd An integer of rows of the secondary legend keys. Applies to the static image and video.
#' @param legend.ncol.2nd An integer of columns of the secondary legend keys. Applies to the static image and video.
#' @param legend.key.size.2nd A numeric of legend key size. The default is 0.03. Applies to the static image and video.
#' @param add.feature.2nd Logical TRUE or FALSE. Add feature identifiers to the secondary legend or not. The default is FALSE. Applies to the static image.
#' @param legend.text.size.2nd A numeric of the secondary legend text size. The default is 10. Applies to the static image and video.
#' @param angle.text.key.2nd A value of angle of key text in the secondary legend. Default is 0. Applies to the static image and video.
#' @param position.text.key.2nd The position of key text in the secondary legend, one of "top", "right", "bottom", "left". Default is "right". Applies to the static image and video.

#' @return A list of ggplots.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @importFrom ggplot2 theme layer_data scale_fill_manual unit element_text guide_legend margin

gg_2lgd <- function(gg.all, sam.dat, ft.trans, position.2nd='bottom', legend.nrow.2nd=NULL, legend.ncol.2nd=NULL, legend.key.size.2nd=0.03, add.feature.2nd=FALSE, legend.text.size.2nd=10, angle.text.key.2nd=0, position.text.key.2nd='right') {

  for (i in seq_along(gg.all)) {

    g <- gg.all[[i]]; lay.dat <- layer_data(g)
    dat <- g$data; g.col <- lay.dat$fill
    names(g.col) <- dat$tissue; df.val <- round(dat$value, 2)
    df.tis <- as.vector(dat$tissue); tis.path <- dat$feature
    g.col <- g.col[!duplicated(names(g.col))]; tis.path <- dat$feature
    ft.legend <- intersect(unique(sam.dat), unique(tis.path))
    ft.legend <- setdiff(ft.legend, ft.trans) 
    leg.idx <- !duplicated(tis.path) & (tis.path %in% ft.legend)
    df.tar <- df.tis[leg.idx]; lab <- val.tar <- df.val[leg.idx]
    path.tar <- tis.path[leg.idx]
    if (add.feature.2nd==TRUE) lab <- paste0(path.tar, ' (', val.tar, ')') 
    gde <- guide_legend(title=NULL, nrow=legend.nrow.2nd, ncol=legend.ncol.2nd, label.theme=element_text(angle=angle.text.key.2nd, size=legend.text.size.2nd), label.position=position.text.key.2nd)
    g <- g+scale_fill_manual(values=g.col, breaks=df.tar, labels=lab, guide=gde)+theme(legend.box.margin=margin(-20, 0, 2, 0, unit='pt'))
    if (!is.null(legend.key.size.2nd)) g <- g+theme(legend.key.size=unit(legend.key.size.2nd, "npc"))
    if (position.2nd!='bottom') g <- g+theme(legend.position=position.2nd)
    gg.all[[i]] <- g
 
  }; return(gg.all)

}
