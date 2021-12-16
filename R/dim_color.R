#' Map colours in SHMs to embedding plots 
#'
#' @param sce.aggr A \code{SingleCellExperiment} containing the (un-)aggregated cells that have source tissue assignments. The \code{lis.match} will be built from \code{colLabels} internally.
#' @param gg.dim The ggplot of embedding plot.
#' @param gg.shm.all The list of SHM ggplot.
#' @param grob.shm.all The list of SHM grob.
#' @param col.shm.all The list of SHM colours.
#' @param color.by A column name in the \code{colData} slot such as \code{label}.  
#' @param lgd.all.dim Logical. The default is \code{FALSE}, and only cells with source tissue assignment have legends.
#' @param con.na Logical, TRUE or FALSE. Default is TRUE, meaning conditions are available.
#' @param lis.match The remaching list of spatial features in data and aSVGs.
#' @param sub.title.size The title size of embedding plots. The default is 11.

#' @return A nested list of embedding and SHM plots.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016

#' @importFrom ggplot2 layer_data ggplot geom_point theme_classic theme element_text element_blank labs 

dim_color <- function(sce.aggr=NULL, gg.dim, gg.shm.all, grob.shm.all, col.shm.all, color.by, lgd.all.dim=FALSE, con.na=TRUE, lis.match=NULL, sub.title.size=11) {
  # save(sce.aggr, gg.dim, gg.shm.all, grob.shm.all, col.shm.all, color.by, lgd.all.dim, con.na, lis.match, sub.title.size, file='dim.color.all')
  if (!is.null(sce.aggr) & !is.null(lis.match)) stop("Only one of 'sce.aggr' and 'lis.match' is required!")
  if (!is.null(sce.aggr)) {
    # The matching list between aggregated cells and aSVG spatial features. The former are cells with a source tissue assignment in co-clustering.   
    blk.uni <- unique(colLabels(sce.aggr))
    lis.match <- as(blk.uni, 'list'); names(lis.match) <- blk.uni
  }
  lis.match <- lis.match[!unlist(lapply(lis.match, is.null))]
  # Ggplots of all reduced dim.
  n <- length(grob.shm.all); gg.dim.all <- rep(list(gg.dim), n)
  names(gg.dim.all) <- paste0('dim_', names(grob.shm.all))
  # Match colors in SHMs to dim plots.
  for (i in seq_along(gg.dim.all)) {
    gg.dim <- gg.dim.all[i]; na <- sub('^dim_', '', names(gg.dim))
    g.col <- col.shm.all[[paste0('col_', na)]]
    dat.ft.na <- names(lis.match)
    # 'gray80' is a reserved color.
    dim.col <- rep('gray80', length(lis.match))
    names(dim.col) <- dat.ft.na
    for (j in dat.ft.na) {
    # Matched svg fts have the same color, so only the 1st is taken. 
      ft.svg <- lis.match[[j]][1]; matched <- g.col[ft.svg]
      if (length(matched)==0) next else if (is.na(matched)) {
        matched <- g.col[sub('__\\d+', '', names(g.col))==ft.svg][1]
        if (!is.na(matched)) if (matched!='NA') dim.col[j] <- matched
      } else if (length(matched)>0) { 
        if (matched!='NA') dim.col[j] <- matched }
    }
    # Max shapes: 128.
    # sp <- seq_along(dim.col); names(sp) <- names(dim.col)
    # Merge colour and shape legend: dim.col and sp have the same names.
    # save(gg.dim, dim.col, sp, color.by, file='gdsc')

    gg.dim0 <- gg.dim[[1]]
    dat <- gg.dim0$data; lay.dat <- layer_data(gg.dim0) 
    dat.all <- cbind(lay.dat, colour_by=dat$colour_by)
    dat.all$fill <- dat.all$colour <- dim.col[dat.all$colour_by]
    dat.all <- dat.all[, !colnames(dat.all) %in% 'colour']
    colnames(dat.all)[colnames(dat.all)=='colour_by'] <- 'feature'
    dat.all$tissue <- dat.all$feature
    ft.all <- unique(dat.all$feature)
    ft.o <- ft.all[!ft.all %in% names(dim.col)]
    col.o <- rep('gray80', length(ft.o))
    names(col.o) <- ft.o; dim.col.all <- c(dim.col, col.o)
    dat.all$fill <- dim.col.all[dat.all$feature]
    # sp <- seq_along(dim.col.all); names(sp) <- names(dim.col.all)
    dat.all <- rbind(subset(dat.all, fill == 'gray80'), subset(dat.all, fill != 'gray80'))

    # gg.dim.all[[i]] <- ggplot(dat.all, aes(x=x, y=y, shape=feature, text=dat.all$feature)) + geom_point(size=2, alpha=1, aes(colour=feature)) + scale_color_manual(values=dim.col) + scale_shape_manual(values=sp) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank()) + labs(title=gsub('^dim_(.*)_\\d+$', '\\1', names(gg.dim)), x=gg.dim0$labels$x, y=gg.dim0$labels$y, colour=color.by, shape=color.by)
  if (con.na==TRUE) tit <- gsub('^dim_(.*)_\\d+$', '\\1', names(gg.dim)) else tit <- gsub('^dim_(.*)_con_\\d+$', '\\1', names(gg.dim))
  # Re-plot dimensionlaity plot.
  gg <- ggplot(dat.all, aes(x=x, y=y, text=dat.all$feature)) + geom_point(size=2, alpha=1, aes(colour=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank()) + labs(title=tit, x=gg.dim0$labels$x, y=gg.dim0$labels$y, colour=color.by, shape=color.by)
  if (lgd.all.dim==FALSE) gg <- gg + scale_color_manual(values=dim.col.all, breaks=names(dim.col)) else gg <- gg + scale_color_manual(values=dim.col.all, breaks=names(dim.col.all))
  gg.dim.all[[i]] <- gg

  }
  # Convert all reduced dim of ggplots to grobs.
  grob.dim.all <- grob_shm(gg.dim.all, lgd.pos='right')
  # Empty list of all reduced dim and SHMs. 
  dim.shm.gg.lis <- dim.shm.grob.lis <- rep(list(NULL), 2*n)

  # Assign all reduced dims to the empty list.
  dim.shm.gg.lis[seq(1, 2*n, 2)] <- gg.dim.all
  dim.shm.grob.lis[seq(1, 2*n, 2)] <- grob.dim.all
  names(dim.shm.gg.lis)[seq(1, 2*n, 2)] <- names(dim.shm.grob.lis)[seq(1, 2*n, 2)] <- names(grob.dim.all)
  # Assign all SHMs to the empty list.
  dim.shm.gg.lis[seq(2, 2*n, 2)] <- gg.shm.all
  dim.shm.grob.lis[seq(2, 2*n, 2)] <- grob.shm.all
  names(dim.shm.gg.lis)[seq(2, 2*n, 2)] <- names(dim.shm.grob.lis)[seq(2, 2*n, 2)] <- names(grob.shm.all)
  return(list(dim.shm.gg.lis=dim.shm.gg.lis, dim.shm.grob.lis=dim.shm.grob.lis))
}
