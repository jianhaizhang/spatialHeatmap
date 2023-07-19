#' Make Spatial Heatmap Video
#'
#' @param gg A list of spatial heatmaps of ggplot.
#' @param cs.g The color key of ggplot.
#' @param lgd.key.size The size of legend key (including text). Default is 0.02.
#' @param lgd.text.size The size of legend text. Default is 8.
#' @param lgd.row An integer of legend rows.
#' @param lgd.col An integer of legend columns.
#' @inheritParams spatial_hm

#' @param angle.text.key A value of key text angle in legend plot. The default is NULL, equivalent to 0.
#' @param position.text.key The position of key text in legend plot, one of "top", "right", "bottom", "left". Default is NULL, equivalent to "right".
#' @param legend.value.vdo Logical TRUE or FALSE. If TRUE, the numeric values of matching spatial features are added to video legend. The default is NULL.
#' @param sub.title.size The title size of ggplot.
#' @param bar.width The color bar width, between 0 and 1.
#' @param label Logical. If TRUE, spatial features having matching samples are labeled by feature identifiers. The default is FALSE. It is useful when spatial features are labeled by similar colors. 
#' @param label.size The size of spatial feature labels. The default is 4.
#' @param label.angle The angle of spatial feature labels in legend plot. Default is 0.
#' @param hjust The value to horizontally adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param vjust The value to vertically adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param opacity The transparency of colored spatial features in legend plot. Default is 1. If 0, features are totally transparent.
#' @param key Logical. The default is TRUE and keys are added in legend plot. If \code{label} is TRUE, the keys could be removed. 


#' @return A video is saved in \code{out.dir}.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Jeroen Ooms (2020). av: Working with Audio and Video in R. R package version 0.5.0. https://CRAN.R-project.org/package=av
#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra

#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 layer_data margin
#' @importFrom grid grobTree unit

video <- function(gg, cs.g, lgd, gcol.lgd=NULL, lgd.r=0.4, h=0.7, type='shm', lgd.title=NULL, sub.title.size=10, bar.width=0.12, bar.axis.title.size=5, bar.value.size=4, lgd.key.size=0.02, lgd.text.size=8, lgd.nrow.2nd=2, lgd.key.size.2nd=0.01, lgd.text.size.2nd=5, angle.text.key=NULL, position.text.key=NULL, lgd.row=2, lgd.row.2nd=2, lgd.col=NULL, legend.value.vdo=NULL, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, dim.lgd.text.size=13, dim.lgd.key.size=4, dim.lgd.nrow=2, dim.lgd.plot.margin=NULL, video.dim='640x480', res=500, interval=1, framerate=1, out.dir) {  
  # save(gg, cs.g, lgd, gcol.lgd, lgd.r, h, type, lgd.title, sub.title.size, bar.width, bar.value.size, lgd.key.size, lgd.text.size, lgd.key.size.2nd, lgd.text.size.2nd, angle.text.key, position.text.key, lgd.row, lgd.row.2nd, lgd.col, legend.value.vdo, label, label.size, label.angle, hjust, vjust, opacity, key, dim.lgd.text.size, dim.lgd.key.size, dim.lgd.nrow, video.dim, res, interval, framerate, out.dir, file='video.arg')

  x <- y <- feature <- NULL
  pkg <- check_pkg('av'); if (is(pkg, 'character')) stop(pkg)
  ffm <- check_exp(test_ffm())   
  if ('w' %in% ffm | 'e' %in% ffm) return()

  if (!is.null(bar.value.size)) cs.g <- cs.g+theme(axis.text.y=element_text(size=bar.value.size))
  if (!is.null(bar.axis.title.size)) cs.g <- cs.g+theme(axis.title.y=element_text(size=bar.axis.title.size))
  if (!'col.idp' %in% type) {
    v.mar <- (1-h)/2
    cs.g <- cs.g+theme(plot.margin=margin(t=v.mar, r=0.01, b=v.mar, l=0.01, unit="npc"))
  }
  na <- names(gg)
  na.shm <- grep('^dim_', na, invert=TRUE, value=TRUE) 
  na.dim <- grep('^dim_', na, value=TRUE) 

  cat('Video: adjust legend size/rows in SHMs ... \n')
  gg[na.shm] <- gg_lgd(gg.all=gg[na.shm], gcol.lgd=gcol.lgd, size.key=lgd.key.size.2nd, size.text.key=lgd.text.size.2nd, angle.text.key=angle.text.key, position.text.key=position.text.key, legend.value.vdo=legend.value.vdo, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, sub.title.size=sub.title.size, row=lgd.row.2nd, col=lgd.col)
  
  cat('Saving video... \n')
  res.r=res/144; w.h <- round(as.numeric(strsplit(video.dim, 'x')[[1]])*res.r)
  if (w.h[1] %% 2!=0) w.h[1] <- w.h[1]+1
  if (w.h[2] %% 2!=0) w.h[2] <- w.h[2]+1
  dir <- normalizePath(out.dir, winslash="/", mustWork=FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  tmp <- normalizePath(tempfile(), winslash='/', mustWork=FALSE)

  if ('shm' %in% type) { # SHM.
    lay <- rbind(c(NA, NA), c(1, 2), c(NA, NA))
    av::av_capture_graphics(expr=for (i in na.shm) { print(grid.arrange(cs.g, gg[[i]], widths=unit(c(bar.width, 1-bar.width), 'npc'), 
  heights=unit(c(0.05, 0.99, 0.05), 'npc'), layout_matrix=lay)) }, 
  output=file.path(dir, 'shm.mp4'), width=w.h[1], height=w.h[2], res=res, vfilter=paste0('framerate=fps=', framerate))
  } else { 
    cat('Video: adjust legend size/rows in dimension reduction plots ... \n')
    w <- (1-bar.width)/2
    for (i in na.dim) { # Refine dims.
      g0 <- gg[[i]]; lgd.com <- lgd_com(g0); dat <- lgd.com$dat
      gg[[i]] <- ggplot(dat, aes(x=x, y=y, text=dat$feature)) + geom_point(size=1, alpha=0.7, colour=dat$fill, aes(colour=feature, shape=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), legend.position='none', legend.text=element_text(size=dim.lgd.text.size), legend.background = element_rect(fill='transparent'), axis.text = element_blank(), axis.ticks = element_blank(), axis.title=element_text(size=4), axis.title.y = element_text(margin = margin(r=-0.2, unit='npc')), axis.title.x = element_text(margin = margin(t=0, unit='npc'))) + scale_shape_manual(values=lgd.com$shape, breaks=lgd.com$br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size))) + labs(title=g0$labels$title, x=g0$labels$x, y=g0$labels$y, colour='', shape='')
    } 
    if ('col.grp' %in% type) { # Covis: coloring by group.
    lay <- rbind(c(NA, NA, NA), c(1, 2, 3), c(NA, NA, NA))
    av::av_capture_graphics(expr=for (i in na.shm) { 
    print(grid.arrange(cs.g, gg[[paste0('dim_', i)]], gg[[i]], widths=unit(c(bar.width, w, w), 'npc'), 
    heights=unit(c(0.05, 0.7, 0.05), 'npc'), layout_matrix=lay)) }, 
    output=file.path(dir, 'shm.mp4'), width=w.h[1], height=w.h[2], res=res, vfilter=paste0('framerate=fps=', framerate))
  } else if ('col.idp' %in% type) { # Covis: coloring by cells.
    png(tmp); cs.grob <- ggplotGrob(cs.g); dev.off()
    if (!is(h, 'numeric')) { msg <- '"Height" should be between 0 and 1!'; warning(msg); return() }
    if (h > 1) h <- 1 else if (h < 0) h <- 0
    # Legend
    cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc"), heights=unit(0.95 * h, "npc"))
    w <- (1-bar.width)/3

    # Dim and SHM.
    png(tmp); grob.lis <- lapply(gg, function(x) { ggplotGrob(x+theme(legend.position="none"))}); dev.off()
    tr <- lapply(grob.lis, grobTree)
    # SHM in legend.
    na.lgd.shm <- grep('\\.svg$', names(lgd), value=TRUE)
    lgd[na.lgd.shm] <- gg_lgd(gg.all=lgd[na.lgd.shm], gcol.lgd=gcol.lgd, size.key=lgd.key.size, size.text.key=lgd.text.size, angle.text.key=angle.text.key, position.text.key=position.text.key, legend.value.vdo=legend.value.vdo, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, sub.title.size=sub.title.size, row=lgd.row, col=lgd.col, title=lgd.title, lgd.space.x=0.01)
    na.lgd.dim <- grep('dim\\.lgd$', names(lgd), value=TRUE)
    # Dim in legend.
    g0 <- lgd[[na.lgd.dim]]; lgd.com <- lgd_com(g0)
    dat <- lgd.com$dat
    if (is.null(dim.lgd.plot.margin)) dim.lgd.plot.margin <- margin(t=0.01, r=0.15, b=0.01, l=0.15, unit="npc")
    lgd[[na.lgd.dim]] <- ggplot(dat, aes(x=x, y=y, text=dat$feature)) + geom_point(size=1, alpha=0.7, aes(colour=feature, shape=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), legend.position='bottom', legend.text=element_text(size=dim.lgd.text.size), legend.background = element_rect(fill='transparent'), axis.text = element_blank(), axis.ticks = element_blank(), axis.title=element_text(size=4), axis.title.y = element_text(margin = margin(r=-0.2, unit='npc')), axis.title.x = element_text(margin = margin(t=0, unit='npc')), legend.spacing.x = unit(-0.005, 'npc'), legend.spacing.y = unit(-0.025, 'npc'), plot.margin=dim.lgd.plot.margin) + scale_color_manual(values=lgd.com$color, breaks=lgd.com$br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, byrow = TRUE)) + scale_shape_manual(values=lgd.com$shape, breaks=lgd.com$br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, byrow = TRUE, override.aes = list(size=dim.lgd.key.size))) + labs(title=g0$labels$title, x=g0$labels$x, y=g0$labels$y, colour='', shape='')

    png(tmp); grob.lgd <- lapply(lgd, ggplotGrob); dev.off()
    lgd.tr <- lapply(grob.lgd, grobTree)

    # If legend.r = 0, legend plot size is a square.
    lgd.arr <- arrangeGrob(grobs=lgd.tr, layout_matrix=matrix(seq_along(lgd), ncol=1), widths=unit(0.99, "npc"), heights=unit(rep(w + (0.99 - w) * lgd.r, length(lgd)), "npc"))
    lay <- rbind(c(NA, NA, NA, NA), c(1, 2, 3, 4), c(NA, NA, NA, NA))

    av::av_capture_graphics(expr=for (i in na.shm) {    
    print(grid.arrange(cs.arr, tr[[paste0('dim_', i)]], tr[[i]], lgd.arr, widths=unit(c(bar.width, w, w, w), 'npc'), 
    heights=unit(c(0.05, h, 0.05), 'npc'), layout_matrix=lay)) 
    }, 
    output=file.path(dir, 'shm.mp4'), width=w.h[1], height=w.h[2], res=res, vfilter=paste0('framerate=fps=', framerate))
  }
  }
}

#' # Retrieve colors, shapes, and breaks in dimension reduction plots.
#'
#' @keywords Internal
#' @noRd

lgd_com <- function(gg) {
  feature <- NULL
  dat <- gg$data; dat.lay <- layer_data(gg)
  dat.lay <- cbind(dat.lay, dat[, c('key', 'feature')])
  df.sub <- subset(dat.lay, !duplicated(feature))
  dim.col.all <- setNames(df.sub$colour, df.sub$feature)
  br <- df.sub$feature[!grepl('^gray80$', df.sub$colour)]
  sp <- setNames(df.sub$shape, df.sub$feature)
  return(list(color=dim.col.all, shape=sp, br=br, dat=dat))
}


#' # Test if "av" works
#'
#' @keywords Internal
#' @noRd

test_ffm <- function() {
  if (any(c('e', 'w') %in% check_pkg('av'))) stop('The package "av" is not detected!')
  av::av_capture_graphics(expr=for (i in seq_len(2)) plot(i), output=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/tmp.mp4'))
}



