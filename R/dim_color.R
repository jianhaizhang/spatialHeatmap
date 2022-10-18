#' Map colours in SHMs to embedding plots through manually-created matching list 
#'
#' @param gg.dim The ggplot of embedding plot.
#' @param gg.shm.all The list of SHM ggplots.
#' @param grob.shm.all The list of SHM grobs.
#' @param col.shm.all The list of SHM colours.
#' @param cell.group A column name in \code{colData} such as \code{cluster} (auto-generated), \code{label}. Cells are divided into clusters by this column name and these clusters are matched to bulk tissues. It is also the legend title in the embedding plots.
#' @param tar.cell The names of target cell clusters to show in embedding plot. The default is \code{matched} and only matching cell clusters have legends in the embedding plot.
#' @param con.na Logical, \code{TRUE} or \code{FALSE}. Default is \code{TRUE}, meaning conditions are available.
#' @param lis.match The maching list of spatial features between data and aSVGs.
#' @param sub.title.size The title size of embedding plots. The default is 11.
#' @param dim.lgd.pos The legend position. The default is \code{bottom}.
#' @param dim.lgd.nrow The number of legend rows. The default is \code{2}.
#' @param dim.lgd.text.size The size of legend text. The default is \code{8}.

#' @return A nested list of embedding and SHM plots.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016

#' @importFrom ggplot2 layer_data ggplot geom_point theme_classic theme element_text element_blank labs scale_color_manual scale_shape_manual margin guide_legend 

dim_color <- function(gg.dim, gg.shm.all, grob.shm.all, col.shm.all, gg.lgd.all, col.lgd.all, grob.lgd.all, profile=TRUE, cell.group, tar.cell, con.na=TRUE, lis.match=NULL, sub.title.size=11, dim.lgd.pos='bottom', dim.lgd.nrow=2, dim.lgd.key.size=4, dim.lgd.text.size=13) {
  # save(gg.dim, gg.shm.all, grob.shm.all, col.shm.all, cell.group, gg.lgd.all, col.lgd.all, grob.lgd.all, profile, tar.cell, con.na, lis.match, sub.title.size, dim.lgd.pos, dim.lgd.nrow, dim.lgd.key.size, dim.lgd.text.size, file='dim.color.arg')
  x <- y <- fill <- feature <- NULL
  lis.match <- lis.match[!unlist(lapply(lis.match, is.null))]
  # if (tar.cell[1]=='matched') tar.cell <- unique(names(lis.match))
  dat.ft.all <- unique(gg.dim$data$colour_by)

  if (any(!tar.cell %in% dat.ft.all)) stop("Make sure all entries in 'tar.cell' are in 'names(lis.match))'!")
  # Ggplots of all reduced dim. 
   if (profile==TRUE) {
    n <- length(grob.shm.all); gg.dim.all <- rep(list(gg.dim), n)
    names(gg.dim.all) <- paste0('dim_', names(grob.shm.all))
  } else if (profile==FALSE) {
    n <- length(gg.lgd.all); gg.dim.all <- rep(list(gg.dim), n)
    names(gg.dim.all) <- paste0('dim_', names(gg.lgd.all))
  }
 
  # Match colors in SHMs to dim plots. Colour order: data -> svg feature -> embedding plot.
  for (i in seq_along(gg.dim.all)) {
    gg.dim <- gg.dim.all[i]
    if (profile==TRUE) dim.col <- col_dim_toblk(gg.dim, gcol.all=col.shm.all, lis.match)
    if (profile==FALSE) dim.col <- col_dim_toblk(gg.dim, gcol.all=col.lgd.all, lis.match)
    # Max shapes: 128.
    # sp <- seq_along(dim.col); names(sp) <- names(dim.col)
    # Merge colour and shape legend: dim.col and sp have the same names.
    # save(gg.dim, dim.col, sp, cell.group, file='gdsc')
    gg.dim0 <- gg.dim[[1]]
    dat <- gg.dim0$data; lay.dat <- layer_data(gg.dim0) 
    dat.all <- cbind(lay.dat, colour_by=dat$colour_by)
    dat.all$fill <- dat.all$colour <- dim.col[dat.all$colour_by]
    dat.all <- dat.all[, !colnames(dat.all) %in% 'colour']
    colnames(dat.all)[colnames(dat.all)=='colour_by'] <- 'feature'
    # dat.all$tissue <- dat.all$feature
    ft.all <- unique(dat.all$feature)
    ft.o <- ft.all[!ft.all %in% names(dim.col)]
    col.o <- rep('gray80', length(ft.o))
    names(col.o) <- ft.o; dim.col.all <- c(dim.col, col.o)
    dat.all$fill <- dim.col.all[dat.all$feature]
    # sp <- seq_along(dim.col.all); names(sp) <- names(dim.col.all)
    dat.all <- rbind(subset(dat.all, fill == 'gray80'), subset(dat.all, fill != 'gray80'))

    # gg.dim.all[[i]] <- ggplot(dat.all, aes(x=x, y=y, shape=feature, text=dat.all$feature)) + geom_point(size=2, alpha=1, aes(colour=feature)) + scale_color_manual(values=dim.col) + scale_shape_manual(values=sp) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank()) + labs(title=gsub('^dim_(.*)_\\d+$', '\\1', names(gg.dim)), x=gg.dim0$labels$x, y=gg.dim0$labels$y, colour=cell.group, shape=cell.group)
  if (con.na==TRUE) tit <- gsub('^dim_(.*)_\\d+$', '\\1', names(gg.dim)) else tit <- gsub('^dim_(.*)_con_\\d+$', '\\1', names(gg.dim))
  if (profile==FALSE) tit <- NULL
  # Non-target cell clusters have colour of 'gray80'.
  # tar.cell is not required, since aSVG features corresponding to non-target cells are already transparent in SHM.
  non.tar <- setdiff(dat.ft.all, tar.cell)
  dim.col.all[non.tar] <- 'gray80'
  # Legal shapes: c(0:25, 32:127)
  sp.sel <- c(15:18, 7:14)
  sp.all <- c(0, 2:25, 32:127)
  sp.all <- c(sp.sel, setdiff(sp.all, sp.sel))
  # Cell cluster shapes.
  sp <- sp.all[seq_along(dat.ft.all)]
  names(sp) <- dat.ft.all
  if (length(non.tar) > 0) br <- tar.cell else br <- dat.ft.all
  # Re-plot dimensionlaity plot.
  gg <- ggplot(dat.all, aes(x=x, y=y, text=dat.all$feature)) + geom_point(size=2, alpha=1, aes(colour=feature, shape=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank(), legend.position=dim.lgd.pos, legend.text=element_text(size=dim.lgd.text.size), legend.margin = margin(t=-0.02, l=0.05, r=0.1, unit='npc')) + scale_color_manual(values=dim.col.all, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow)) + scale_shape_manual(values=sp, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size))) + labs(title=tit, x=gg.dim0$labels$x, y=gg.dim0$labels$y, colour=cell.group, shape=cell.group) 
  gg.dim.all[[i]] <- gg

  }
  # Convert all reduced dim of ggplots to grobs.
  grob.dim.all <- grob_shm(gg.dim.all, lgd.pos=NULL)
  # Empty list of all reduced dim and SHMs. 
  dim.shm.gg.lis <- dim.shm.grob.lis <- rep(list(NULL), 2*n)

  # Assign all reduced dims to the empty list.
  dim.shm.gg.lis[seq(1, 2*n, 2)] <- gg.dim.all
  dim.shm.grob.lis[seq(1, 2*n, 2)] <- grob.dim.all
  names(dim.shm.gg.lis)[seq(1, 2*n, 2)] <- names(dim.shm.grob.lis)[seq(1, 2*n, 2)] <- names(grob.dim.all)
  # Assign all SHMs to the empty list.
  if (profile==TRUE) {
    dim.shm.gg.lis[seq(2, 2*n, 2)] <- gg.shm.all
    dim.shm.grob.lis[seq(2, 2*n, 2)] <- grob.shm.all
    names(dim.shm.gg.lis)[seq(2, 2*n, 2)] <- names(dim.shm.grob.lis)[seq(2, 2*n, 2)] <- names(grob.shm.all)
  } else {
    dim.shm.gg.lis[seq(2, 2*n, 2)] <- gg.lgd.all
    dim.shm.grob.lis[seq(2, 2*n, 2)] <- grob.lgd.all
    names(dim.shm.gg.lis)[seq(2, 2*n, 2)] <- names(dim.shm.grob.lis)[seq(2, 2*n, 2)] <- names(grob.lgd.all)
  }
  return(list(dim.shm.gg.lis=dim.shm.gg.lis, dim.shm.grob.lis=dim.shm.grob.lis))
}


#' Assign colors from SVG features to cell labels in embedding plot through a matching list when mapping cells to bulk.
#' 
#' @return A vector.
#' @keywords Internal
#' @noRd 

col_dim_toblk <- function(gg.dim, gcol.all, lis.match) { 
  na <- sub('^dim_', '', names(gg.dim))
  g.col <- gcol.all[[paste0('col_', na)]]
  dat.ft.na <- names(lis.match)
  # 'gray80' is a reserved color.
  dim.col <- rep('gray80', length(lis.match))
  names(dim.col) <- dat.ft.na
  for (j in dat.ft.na) {
  # Matched svg fts have the same color, so only the 1st is taken. 
    ft.svg <- lis.match[[j]][1]; matched <- g.col[ft.svg]
    # If "matched" is valid, there is no "sub('__\\d+', '', names(g.col))", so speed is faster.
    if (length(matched)==0) next else if (is.na(matched)) {
      matched <- g.col[sub('__\\d+', '', names(g.col))==ft.svg][1]
      if (!is.na(matched)) if (matched!='NA') dim.col[j] <- matched
    } else if (length(matched)>0) { 
      if (matched!='NA') dim.col[j] <- matched }
    }
    return(dim.col)
}


