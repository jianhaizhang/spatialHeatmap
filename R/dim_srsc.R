#' Embedding plots and overlaid plots 
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
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.  0, https://bioconductor.org/packages/SummarizedExperiment.
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 layer_data ggplot geom_point theme_classic theme element_text element_blank labs scale_color_manual scale_shape_manual margin guide_legend element_rect 
#' @importFrom stats setNames

dim_srsc <- function(cell, ID, geneV, cols, tar.cell=NULL, tar.bulk=NULL, con.na.cell, dimred, size.pt=1, alpha.pt=0.8, shape=NULL, svg.all, profile, gg.shm.all, grob.shm.all, gg.lgd.all, col.lgd.all, grob.lgd.all, cell.group, lis.match=NULL, sub.title.size=11, sub.title.vjust=2, dim.lgd.pos='bottom', dim.lgd.nrow=2, dim.lgd.key.size=4, dim.lgd.text.size=13, dim.axis.font.size=10, linewidth=0.1, line.color='gray50', dim.lgd.direc=NULL, size.r=1, verbose=TRUE) {
  # align <- spatialHeatmap:::align
  # shp <- spatialHeatmap:::shp
  # col_dim_tocell <- spatialHeatmap:::col_dim_tocell
  # col_dim_toblk <- spatialHeatmap:::col_dim_toblk
  # svg_separ <- spatialHeatmap:::svg_separ
  # value2color <- spatialHeatmap:::value2color
  # grob_shm <- spatialHeatmap:::grob_shm
  # save(cell, ID, geneV, cols, tar.cell, tar.bulk, con.na.cell, dimred, size.pt, alpha.pt, shape, svg.all, profile, gg.shm.all, grob.shm.all, gg.lgd.all, col.lgd.all, grob.lgd.all, cell.group, lis.match, sub.title.size, sub.title.vjust, dim.lgd.pos, dim.lgd.nrow, dim.lgd.key.size, dim.lgd.text.size, dim.axis.font.size, linewidth, line.color, dim.lgd.direc, size.r, verbose, file='dim.srsc.arg')
  if (verbose==TRUE) message('Spatial single cells: embedding/overlay plots ...')
  feature <- x <- y <- group <- X <- Y <- NULL
  cna <- colnames(cell)
  if (con.na.cell==FALSE & length(grep('__', cna))==0) {
    colnames(cell) <- paste0(cna, '__con')
    cell$variable <- 'con'; vars <- 'con'
  } else vars <- unique(sub('.*__', '', cna))

  cord.pt <- data.frame(colData(cell)[, c('X', 'Y', cell.group)])
  cordn <- coordinate(svg.all)[[1]]; agl <- angle(svg.all)[[1]]
  df.ovl <- subset(cordn, feature=='overlay')
  cord.pt <- align(df.pt=cord.pt, df.ovl=df.ovl, angle=agl, size.r=size.r)
  # Target cell groups.
  grp <- unique(cord.pt[, cell.group])
  if (is.null(tar.cell) & !is.null(tar.bulk)) { 
    tar.cells <- as.character(unname(unlist(lis.match[tar.bulk])))
  } else if (!is.null(tar.cell) & is.null(tar.bulk)) {
    tar.cells <- tar.cell
    if (tar.cell[1] %in% 'all') tar.cells <- grp
    if (any(!tar.cells %in% grp)) { 
      msg <- 'Ensure the target cells are in the "cell.group"!'
      warning(msg); return(msg)
    }
  }
  br <- grp[grp %in% tar.cells]
  dat <- data.frame(reducedDim(cell, dimred)[, 1:2])
  xlab <- paste0(dimred, '1'); ylab <- paste0(dimred, '2')
  colnames(dat) <- c('x', 'y')
  dat$group <- colData(cell)[, cell.group] 
  sp <- shp(shape, grp) 
  dim.lgd.lis <- ovl.lgd.lis <- dim.lis <- ovl.lis <- NULL 
  # Dim and scSHM lgd plots under each variable.
  for (i in vars) {
    cell.col <- setNames(rep('gray80', length(grp)), grp)
    dim.na <- names(gg.lgd.all)[1]
    # Copy colors from SHMs to cells.
    if (!is.null(tar.bulk)) {
      tar.cell.col <- col_dim_tocell(dim.na=dim.na, gcol.all=col.lgd.all, lis.match=lis.match[tar.bulk])   
    } else if (!is.null(tar.cell)) {
      tar.cell.col <- col_dim_toblk(dim.na=dim.na, gcol.all=col.lgd.all, lis.match=lis.match[tar.cell])
    }; cell.col[names(tar.cell.col)] <- tar.cell.col
    # colnames(cell): can be duplicated, not coerced by appending with numbers. 
    idx <- grepl(paste0('__', i, '$'), colnames(cell))
    if (con.na.cell==TRUE) tit <- i else tit <- NULL
    # Dim lgd plots under each variable.
    dat0 <- dat[idx, , drop=FALSE]
    idx.tar <- dat0$group %in% tar.cells
    dat0 <- rbind(dat0[!idx.tar, , drop=FALSE], dat0[idx.tar, , drop=FALSE])
    g.dim.fix <- ggplot(dat0) + geom_point(mapping=aes(x=x, y=y, colour=group, shape=group), size=size.pt, alpha=alpha.pt) + scale_color_manual(values=cell.col, labels=br, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow)) + labs(title=tit)
    # Dim lgd plots are organized under each variable in a nested list.
    dim.lgd.lis <- c(dim.lgd.lis, setNames(list(setNames(list(g.dim.fix), 'dim.lgd')), i))
    # scSHM lgd plot under each variable.
    g.ovl <- ggplot(cord.pt[idx, , drop=FALSE], aes(text=rownames(cord.pt))) + geom_point(mapping=aes(x=X, y=Y, colour=cord.pt[, cell.group], shape=cord.pt[, cell.group]), size=size.pt, alpha=alpha.pt) + scale_color_manual(values=cell.col, labels=br, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow)) + labs(title=tit)
    # scSHM lgd plots are organized under each variable in a nested list.
    ovl.lgd.lis <- c(ovl.lgd.lis, setNames(list(setNames(list(g.ovl), 'ovl.lgd')), i))
  }
  # Dim and scSHM plots under each gene and each variable.
  dat.sc <- assays(cell)$logcounts
  for (i in ID) {
    for (j in vars) {
      # One gene under one variable.
      idx <- grepl(paste0('__', j, '$'), colnames(cell))
      dat0 <- dat[idx, , drop=FALSE]; dat0$colour <- 'gray80'
      idx.tar <- dat0$group %in% tar.cells & idx
      color.dat <- value2color(dat.sc[i, idx.tar, drop=FALSE], geneV, cols)
      dat0$colour[idx.tar] <- color.dat
      # Target cells on top.
      dat0 <- rbind(dat0[!idx.tar, , drop=FALSE], dat0[idx.tar, , drop=FALSE])
      tit <- ifelse(con.na.cell==TRUE, paste0(i, ' ', j), i)
      # Dim plots with cell-by-value coloring.  
      g.dim <- ggplot(dat0) + geom_point(mapping=aes(x=x, y=y, shape=group), size=size.pt, alpha=alpha.pt, colour=dat0$colour) + labs(title=tit)
      na0 <- paste0('dim_', i, '_', j, '_1')
      dim.lis <- c(dim.lis, setNames(list(g.dim), na0))

      cord.pt0 <- cord.pt[idx, , drop=FALSE]
      cord.pt0$colour <- 'gray80'
      cord.pt0$colour[idx.tar] <- color.dat
      # Section plots with cell-by-value coloring. 
      g.ovl <- ggplot(cord.pt0, aes(text=rownames(cord.pt0))) + geom_point(mapping=aes(x=X, y=Y, colour=cord.pt0[, cell.group], shape=cord.pt0[, cell.group]), size=size.pt, alpha=alpha.pt, colour=cord.pt0$colour) + labs(title=tit)
      na0 <- paste0('ovl_', i, '_', j, '_1')
      ovl.lis <- c(ovl.lis, setNames(list(g.ovl), na0))
    }
  }
  dim.all <- NULL
  for (i in vars) dim.all <- c(dim.all, dim.lgd.lis[[i]])
  dim.all <- c(dim.all, dim.lis)
  # Show target cells in dim plots.
  for (i in seq_along(dim.all)) {
    dim.all[[i]] <- dim.all[[i]] + scale_shape_manual(values=sp, labels=br, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size))) + labs(x=xlab, y=ylab, colour=grp, shape=grp) + theme_classic() + theme(plot.title=element_text(hjust=0.5, vjust=sub.title.vjust, size=sub.title.size), plot.margin=margin(0.005, 0.005, 0.005, 0.005, "npc"), legend.box.margin=margin(-3, 0, 2, 0, unit='pt'), legend.background=element_rect(color=NA, fill='transparent'), legend.position=dim.lgd.pos, legend.direction=dim.lgd.direc, legend.text=element_text(size=dim.lgd.text.size), legend.margin=margin(l=0.01, r=0.01, unit='npc'), axis.text = element_blank(), axis.ticks = element_blank(), axis.title=element_text(size=dim.axis.font.size), aspect.ratio=1)
  }
  for (i in seq_along(vars)) dim.lgd.lis[[i]] <- dim.all[i]
  dim.lis <- dim.all[names(dim.lis)]
  for (i in seq_along(dim.lis)) dim.lis[[i]] <- dim.lis[[i]] + guides(colour="none", shape='none')
  
  ovl.all <- NULL
  for (i in vars) ovl.all <- c(ovl.all, ovl.lgd.lis[[i]])
  ovl.all <- c(ovl.all, ovl.lis) 
  aspect.ratio <- svg_separ(svg.all)$aspect.r
  # Show target cells in scSHM.
  for (i in seq_along(ovl.all)) { 
    ovl.all[[i]] <- ovl.all[[i]] + scale_shape_manual(values=sp, labels=br, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size))) + labs(x='', y='', colour=grp, shape=grp) + geom_polygon(data=cordn, mapping=aes(x=x, y=y, group=feature), fill=NA, color=line.color, linewidth=linewidth, linetype='solid', alpha=1, inherit.aes = FALSE)+ theme_void() + theme(plot.title=element_text(hjust=0.5, vjust=sub.title.vjust+3, size=sub.title.size), plot.margin=margin(0.005, 0.005, 0.005, 0.005, "npc"), legend.box.margin=margin(-3, 0, 2, 0, unit='pt'), legend.background=element_rect(color=NA, fill='transparent'), aspect.ratio = 1/aspect.ratio, legend.position=dim.lgd.pos, legend.direction=dim.lgd.direc, legend.text=element_text(size=dim.lgd.text.size), legend.margin=margin(l=0.01, r=0.01, unit='npc'))+scale_y_continuous(expand=expansion(mult=c(0, 0)))+scale_x_continuous(expand=expansion(mult=c(0, 0)))
  }
  for (i in seq_along(vars)) ovl.lgd.lis[[i]] <- ovl.all[i]
  ovl.lis <- ovl.all[names(ovl.lis)]
  for (i in seq_along(ovl.lis)) ovl.lis[[i]] <- ovl.lis[[i]] + guides(colour="none", shape='none')
  # Grobs of dim and scSHM lgd plots.
  for (i in vars) dim.lgd.lis[[i]]$dim.lgd.grob <-  grob_shm(list(dim.lgd.lis[[i]]$dim.lgd), lgd.pos=NULL, verbose=verbose)[[1]]
  for (i in vars) ovl.lgd.lis[[i]]$ovl.lgd.grob <-  grob_shm(list(ovl.lgd.lis[[i]]$ovl.lgd), lgd.pos=NULL, verbose=verbose)[[1]]

  # Grobs of dim and scSHM plots with cell-by-value coloring.
  grob.dim.lis <- grob_shm(dim.lis, lgd.pos=NULL, verbose=verbose)
  grob.ovl.lis <- grob_shm(ovl.lis, lgd.pos=NULL, verbose=verbose)
  # Empty list of all dim, scSHM, and SHM plots.
  if (profile==FALSE) n <- length(vars)
  if (profile==TRUE) n <- length(gg.shm.all)
  n.all <- 3 * n
  dim.ovl.shm.gg <- dim.ovl.shm.grob <- rep(list(NULL), n.all)

  # Assign all reduced dims to the empty list.
  dim.ovl.shm.gg[seq(1, n.all, 3)] <- dim.lis
  dim.ovl.shm.grob[seq(1, n.all, 3)] <- grob.dim.lis
  names(dim.ovl.shm.gg)[seq(1, n.all, 3)] <- names(dim.ovl.shm.grob)[seq(1, n.all, 3)] <- names(grob.dim.lis)
  # Assign all scSHMs to the empty list.
  dim.ovl.shm.gg[seq(2, n.all, 3)] <- ovl.lis
  dim.ovl.shm.grob[seq(2, n.all, 3)] <- grob.ovl.lis
  names(dim.ovl.shm.gg)[seq(2, n.all, 3)] <- names(dim.ovl.shm.grob)[seq(2, n.all, 3)] <- names(grob.ovl.lis)
  # Assign all SHMs to the empty list.
  if (profile==TRUE) {
    dim.ovl.shm.gg[seq(3, n.all, 3)] <- gg.shm.all
    dim.ovl.shm.grob[seq(3, n.all, 3)] <- grob.shm.all
    names(dim.ovl.shm.gg)[seq(3, n.all, 3)] <- names(dim.ovl.shm.grob)[seq(3, n.all, 3)] <- names(grob.shm.all)
  } else { # Assign dim and scSHM lgd plots under each variable to the empty list.
    dim.ovl.shm.gg[seq(1, n.all, 3)] <- lapply(vars, function(i) dim.lgd.lis[[i]]$dim.lgd)[1]
    dim.ovl.shm.grob[seq(1, n.all, 3)] <- lapply(vars, function(i) dim.lgd.lis[[i]]$dim.lgd.grob)[1]
    dim.ovl.shm.gg[seq(2, n.all, 3)] <- lapply(vars, function(i) ovl.lgd.lis[[i]]$ovl.lgd)[1]
    dim.ovl.shm.grob[seq(2, n.all, 3)] <- lapply(vars, function(i) ovl.lgd.lis[[i]]$ovl.lgd.grob)[1]
    # The first legend SHM is considered.
    dim.ovl.shm.gg[seq(3, n.all, 3)] <- gg.lgd.all[1]
    dim.ovl.shm.grob[seq(3, n.all, 3)] <- grob.lgd.all[1]

    names(dim.ovl.shm.gg)[seq(1, n.all, 3)] <- names(dim.ovl.shm.grob)[seq(1, n.all, 3)] <- paste0('dim_', names(dim.lgd.lis))  
    names(dim.ovl.shm.gg)[seq(2, n.all, 3)] <- names(dim.ovl.shm.grob)[seq(2, n.all, 3)] <- paste0('ovl_', names(ovl.lgd.lis))
    names(dim.ovl.shm.gg)[seq(3, n.all, 3)] <- names(dim.ovl.shm.grob)[seq(3, n.all, 3)] <- names(grob.lgd.all)
  }
  return(list(dim.ovl.shm.gg=dim.ovl.shm.gg, dim.ovl.shm.grob=dim.ovl.shm.grob, dim.lgd.lis=dim.lgd.lis, ovl.lgd.lis=ovl.lgd.lis))
}
 

#' Rotating or moving single cell coordinates.
#'
#' # https://github.com/tibo31/GeoXp/blob/master/R/rotation.R  
#' # rotation <- function(coords, angle) { 
#' #  radian <- (angle * pi) / 180 
#' #  x <- c(cos(radian), -sin(radian)) 
#' #  y <- c(sin(radian), cos(radian))
#' #  nlecoord <- coords %*% cbind(x, y); return(nlecoord)
#' # } 

#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' https://github.com/ErasmusOIC/SMoLR/blob/master/R/rotate.R

rotate_cord <- function(x, y, angle, type=c("degrees","radial"), method=c("transform","polar","polar_extended"), center=c(0,0), translate=NULL, stretch=NULL, flip=FALSE){
  type <- match.arg(type); method <- match.arg(method)
  if(!(length(translate)==2 || is.null(translate))){stop("translation coordinates should be a vector of length 2")}
  if(!(is.logical(flip))){stop("Flip should be TRUE or FALSE")}
  if(flip){ x <- -x }
  if(!is.null(stretch)){
    x <- x*stretch; y <- y*stretch; center <- center*stretch
    if(!is.null(translate)){translate<- translate*stretch}
  }
  x <- x-center[1]; y <- y-center[2]
  if(type=="degrees"){angle <- angle*pi/180}
  if(type=="radial" && angle>(2*pi)){warning("Angle is bigger than 2pi are you sure it's in rads", call. = FALSE)}
  if(method=="polar" || method=="polar_extended"){
    r <-sqrt(x^2+y^2); phi <- atan2(x,y)
    new_x <- r*sin(phi+angle); new_y <- r*cos(phi+angle)
    xy <- cbind(new_x, new_y)
  }
  if(method=="polar_extended"){
    switch(type, degrees={ phi <- (phi+angle)*180/pi}, radial={phi <- phi+angle}
    )
  ext_list <- list(Coordinates=xy, Angles=phi, Distance_from_center=r); return(invisible(ext_list))
  }
  if(method=="transform"){
    conversionmatrix <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), ncol=2, nrow=2)
    xy <- cbind(x,y)%*%conversionmatrix
  }
  xy[, 1] <- xy[, 1]+center[1]; xy[, 2] <- xy[, 2]+center[2]
  if(!is.null(translate)){
    xy[, 1] <- xy[, 1]+translate[1]
    xy[, 2] <- xy[, 2]+translate[2]
  }; return(xy)
}


#' Convert seurat objects to sce objects, including image coordinates. 
#'
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.  0, https://bioconductor.org/packages/SummarizedExperiment.
#' Hao and Hao et al. Integrated analysis of multimodal single-cell data. Cell (2021) [Seurat V4]
#' Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019) [Seurat V3]
#' Butler et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat Biotechnol (2018) [Seurat V2]
#' Satija and Farrell et al. Spatial reconstruction of single-cell gene expression data. Nat Biotechnol (2015) [Seurat V1]
#' Pagès H, Lawrence M, Aboyoun P (2022). _S4Vectors: Foundation of vector-like and list-like containers in Bioconductor_. R package version 0.36.1, <https://bioconductor.org/packages/S4Vectors>.

#' @importFrom SummarizedExperiment colData<- assays<- 
#' @importFrom S4Vectors DataFrame 

srt2sce <- function(cell, assay, image, x, y, cell.group, var.cell) {
  pkg <- check_pkg('Seurat'); if (is(pkg, 'character')) stop(pkg)
  cord.pt <- cbind(cell@images[[image]]@coordinates, cell@meta.data)
  cna <- colnames(cord.pt)
  cna[cna %in% x] <- 'X'; cna[cna %in% y] <- 'Y'
  colnames(cord.pt) <- cna
  cord.pt[, cell.group] <- as.vector(cord.pt[, cell.group])
  sce.sp <- Seurat::as.SingleCellExperiment(cell, assay=assay)
  assays(sce.sp)$counts <- assays(sce.sp)$scaledata <- NULL
  colData(sce.sp) <- cdat <- DataFrame(cord.pt); 
  colnames(sce.sp) <- paste0(cdat[, cell.group], '__', cdat[, var.cell])
  sce.sp
}

#' Aligning single cells on a section with target region in SHM.
#'
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

align <- function(df.pt, df.ovl, angle, size.r=1, move.x=0, move.y=0) {
  # Rotate section.
  df.pt[, c('X', 'Y')] <- rotate_cord(x=df.pt$X, y=df.pt$Y, angle=angle, type="degrees", method="transform", stretch=1)
  # Scale section.
  x.scl <- (max(df.pt$X) - min(df.pt$X)) / (max(df.ovl$x) - min(df.ovl$x))
  y.scl <- (max(df.pt$Y) - min(df.pt$Y)) / (max(df.ovl$y) - min(df.ovl$y))
  df.pt$X <- df.pt$X/x.scl * size.r
  df.pt$Y <- df.pt$Y/y.scl * size.r
  # Move section.
  df.pt$X <- (max(df.ovl$x) + min(df.ovl$x))/2 - (max(df.pt$X) + min(df.pt$X))/2 + df.pt$X + move.x
  df.pt$Y <- (max(df.ovl$y) + min(df.ovl$y))/2 - (max(df.pt$Y) + min(df.pt$Y))/2 + df.pt$Y + move.y
  return(df.pt)
}

#' Assigning shapes to cell groups. 
#'
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

shp <- function(shape, group) { 
  # Legal shapes: c(0:25, 32:127)
  sp.all <- c(0, 2:25, 32:127)
  if (is.null(shape)) if (length(group) <= 122) {
    sp.sel <- c(15:18, 7:14, 3:4)
    sp.all <- c(sp.sel, setdiff(sp.all, sp.sel))
    # Cell cluster shapes.
    sp <- sp.all[seq_along(group)]; names(sp) <- group
  } else sp <- setNames(rep(19, length(group)), group)
  if (!is.null(shape)) {
    if (!all(names(shape) %in% group)) return("All shapes should be named with cell group labels!")
    if (!all(shape %in% sp.all)) return("Legal shapes are c(0, 2:25, 32:127)!")
    sp <- shape
  }; return(sp)
}




