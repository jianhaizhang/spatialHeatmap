#' Adjust Legend Key Size and Rows in Ggplot.
#'
#' @param gg.all A list of spatial heatmaps of ggplot.
#' @param size.key A numeric of legend key size. If \code{size.text} is NULL, it also applies to legend text size.
#' @param size.text.key A numeric of legend text size.
#' @param angle.text.key A value of key text angle in legend plot. The default is NULL, equivalent to 0.
#' @param position.text.key The position of key text in legend plot, one of "top", "right", "bottom", "left". Default is NULL, equivalent to "right".
#' @param legend.value.vdo Logical TRUE or FALSE. If TRUE, the numeric values of matching spatial features are added to video legend. The default is NULL.
#' @param sub.title.size The title size of ggplot.
#' @param row An integer of rows in legend key.
#' @param col An integer of columns in legend key.
#' @param label Logical. If TRUE, spatial features having matching samples are labeled by feature identifiers. The default is FALSE. It is useful when spatial features are labeled by similar colors. 
#' @param label.size The size of spatial feature labels in legend plot. The default is 4.
#' @param label.angle The angle of spatial feature labels in legend plot. Default is 0.
#' @param hjust The value to horizontally adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param vjust The value to vertically adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param opacity The transparency of colored spatial features in legend plot. Default is 1. If 0, features are totally transparent.
#' @param key Logical. The default is TRUE and keys are added in legend plot. If \code{label} is TRUE, the keys could be removed. 
#' @param sam.dat A vector of samples in the data matrix.
#' @param ft.trans A vector of tissues to be transparent.
#' @return A list of ggplots.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @importFrom ggplot2 theme layer_data scale_fill_manual unit element_text guide_legend

gg_lgd <- function(gg.all, size.key=NULL, size.text.key=8, angle.text.key=NULL, position.text.key=NULL, legend.value.vdo=NULL, sub.title.size=NULL, row=NULL, col=NULL, label=FALSE, label.size=3, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, sam.dat, ft.trans=NULL) {

  tissue <- x0 <- y0 <- NULL
  # Function to remove feature labels. 
  rm_label <- function(g) {
        
    g.layer <- g$layer; if (length(g.layer)==1) return(g) 
    for (k in rev(seq_along(g.layer))) {

      na.lay <- unique(names(as.list(g.layer[[k]])$geom_params))
      if (all(c('check_overlap', 'angle', 'size') %in% na.lay)) g$layers[[k]] <- NULL

    }; return(g)

  }
  for (i in seq_along(gg.all)) {
  
    g <- gg.all[[i]] 
    if (!is.null(size.key)) g <- g+theme(legend.key.size=unit(size.key, "npc"), legend.text=element_text(size=ifelse(is.null(size.text.key), 8*size.key*33, size.text.key)))
    if (!is.null(sub.title.size)) g <- g+theme(plot.title=element_text(hjust=0.5, size=sub.title.size))
    if (!is.null(row)|!is.null(col)|label==TRUE|opacity!=1|!is.null(angle.text.key)|!is.null(position.text.key)|!is.null(legend.value.vdo)) {

      lay.dat <- layer_data(g)
      dat <- g$data; g.col <- lay.dat$fill
      names(g.col) <- dat$tissue; df.tis <- as.vector(dat$tissue); df.val <- round(dat$value, 2)
      g.col <- g.col[!duplicated(names(g.col))]; tis.path <- dat$feature
      ft.legend <- intersect(unique(sam.dat), unique(tis.path))
      ft.legend <- setdiff(ft.legend, ft.trans) 
      leg.idx <- !duplicated(tis.path) & (tis.path %in% ft.legend)
      df.tar <- df.tis[leg.idx]; lab <- path.tar <- tis.path[leg.idx]; val.tar <- df.val[leg.idx]
      if (sum(legend.value.vdo)==1) lab <- paste0(path.tar, ' (', val.tar, ')') 
      if (opacity!=1) g.col <- alpha(g.col, opacity)
      if (key==TRUE) gde <- guide_legend(title=NULL, nrow=row, ncol=col, label.theme=element_text(angle=angle.text.key, size=g$theme$legend.text$size), label.position=position.text.key)
      if (key==FALSE) gde <- FALSE
      if (!is.null(row)|!is.null(col)|opacity!=1|key==FALSE|!is.null(angle.text.key)|!is.null(position.text.key)|!is.null(legend.value.vdo)) g <- g+scale_fill_manual(values=g.col, breaks=df.tar, labels=lab, guide=gde)
      if (label==TRUE) {

        dat$x0 <- dat$y0 <- dat$label <- NA
        lab.idx <- dat$feature %in% path.tar
        dat1 <- dat[lab.idx, ]; dat1$label <- dat1$feature
        df.lab <- data.frame() 
        for (j in unique(dat1$tissue)) {

          df0 <- subset(dat1, tissue==j)
          x <- mean(df0$x); y <- mean(df0$y)
          df0$x0 <- x; df0$y0 <- y
          df.lab <- rbind(df.lab, df0)
       
         } 
         g <- rm_label(g)+geom_text(data=df.lab, aes(label=label, x=x0, y=y0), check_overlap=TRUE, size=label.size, angle=label.angle, hjust=hjust, vjust=vjust)

      }; gg.all[[i]] <- rm_label(g)

    }; if (label==FALSE) { gg.all[[i]] <- rm_label(g) }

  }; return(gg.all)

}





