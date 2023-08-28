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
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.  0, https://bioconductor.org/packages/SummarizedExperiment.
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 layer_data ggplot geom_point theme_classic theme element_text element_blank labs scale_color_manual scale_shape_manual margin guide_legend element_rect
#' @importFrom stats setNames 
     
dim_decon <- function(sce=NULL, ID, data, geneV, cols, tar.cell=NULL, con.na.cell, dimred, profile, size.pt=1.5, gg.shm.all, grob.shm.all, col.shm.all, gg.lgd.all, col.lgd.all, grob.lgd.all, cell.group, lis.match=NULL, sub.title.size=11, dim.lgd.pos='bottom', dim.lgd.nrow=2, dim.lgd.key.size=4, dim.lgd.text.size=13, dim.axis.font.size=10, alpha.pt=0.8, size.lab.pt=5, hjust.lab.pt=0.5, vjust.lab.pt=1.5, shape=NULL) {
  # save(sce, ID, data, geneV, cols, tar.cell, con.na.cell, dimred, profile, size.pt, gg.shm.all, grob.shm.all, col.shm.all, gg.lgd.all, col.lgd.all, grob.lgd.all, cell.group, lis.match, sub.title.size, dim.lgd.pos, dim.lgd.nrow, dim.lgd.key.size, dim.lgd.text.size, dim.axis.font.size, alpha.pt, size.lab.pt, hjust.lab.pt, vjust.lab.pt, shape, file='dim.decon.arg')
  x <- y <- fill <- feature <- variable <- label <- NULL
  cdat <- colData(sce); blk <- unique(cdat$bulk)
  lis.match <- NULL; for (i in blk) {
    cdat0 <- subset(cdat, bulk==i)
    cell.uni <- unique(cdat0[, cell.group])
    if (!is.null(tar.cell)) {
      tar0 <- cell.uni[cell.uni %in% tar.cell]
    }
    lis.match <- c(lis.match, setNames(list(tar0), i))
  }
  # Legal shapes: c(0:25, 32:127)
  # sp.sel <- c(15:18, 7:14); sp.all <- c(0, 2:25, 32:127)
  # sp.all <- c(sp.sel, setdiff(sp.all, sp.sel))
  # Cell cluster shapes.
  cell.all <- unique(cdat[, cell.group]) 
  # sp <- sp.all[seq_along(cell.all)]; names(sp) <- cell.all
  sp <- shp(shape, cell.all)
  vars <- unique(cdat$variable) 
   dim.lgd.lis <- lapply(seq_along(vars), function(i) { 
     # Transfer colors from SHM to dim plots. The 1st legend SHM is considered.
     dim.col <- col2cell(col.lgd.all[[1]], lis.match)
     # Non-target cells have color gray80.
     cell.all <- unique(unlist(lis.match))    
     dim.col.all <- setNames(rep('gray80', length(cell.all)), cell.all)
     dim.col.all[names(dim.col)] <- dim.col 

     cdat0 <- subset(cdat, variable==vars[i])
     grps <- cdat0[, cell.group]; grp <- unique(grps)
     br <- grp[grp %in% tar.cell]
     dat <- data.frame(reducedDim(sce, dimred)[, 1:2])
     xlab <- paste0(dimred, '1'); ylab <- paste0(dimred, '2')
     colnames(dat) <- c('x', 'y')
     dat$feature <- grps
     dat$fill <- dim.col.all[dat$feature]
     dat$fill[is.na(dat$fill)] <- 'gray80'
     dat$label <- sce$prop
     if (is.null(dat$label)) dat$label <- ''
     dat <- rbind(subset(dat, fill == 'gray80'), subset(dat, fill != 'gray80'))
     tit <- NULL; if (TRUE %in% con.na.cell) tit <- vars[i] 
     # Re-plot dimensionlaity plot.

     dim.lgd <- ggplot(dat, aes(x=x, y=y, text=dat$feature, label=label)) + geom_point(size=size.pt, alpha=alpha.pt, aes(shape=feature), colour=dat$fill) + geom_text(size=size.lab.pt, hjust=hjust.lab.pt, vjust=vjust.lab.pt, check_overlap = TRUE) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank(), legend.position=dim.lgd.pos, legend.text=element_text(size=dim.lgd.text.size), legend.margin = margin(t=-0.02, l=0.05, r=0.1, unit='npc'), legend.background = element_rect(fill='transparent'), axis.title=element_text(size=dim.axis.font.size), aspect.ratio=1) + scale_shape_manual(values=sp, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size))) + labs(title=tit, x=xlab, y=ylab, shape=cell.group); list(dim.lgd=dim.lgd, dat=dat, br=br)
   }); names(dim.lgd.lis) <- vars
  # Dim plots per ID per variable.
  gg.dim.all <- NULL; for (id in ID) {
    for (vari in names(dim.lgd.lis)) {
      # Re-use data extracted from original dim plots
      dat <- dim.lgd.lis[[vari]]$dat
      # Re-use shapes in dim legend plots.
      dim.lgd <- dim.lgd.lis[[vari]]$dim.lgd
      dat$fill <- 'gray80'
      # Non-target cells have color gray80.
      idx.mat <- dat$feature %in% tar.cell
      # Isolate assay data by ID and variable.
      gene0 <- data[, grepl(paste0('__', vari, '$'), colnames(data)), drop=FALSE]
      # Order in gene0 is the same with dat.
      color.dat <- value2color(gene0[id, idx.mat, drop=FALSE], geneV, cols)
      # Colors of target cells.
      dat$fill[idx.mat] <- color.dat
      dat <- rbind(subset(dat, fill == 'gray80'), subset(dat, fill != 'gray80'))
      tit <- ifelse(con.na.cell==TRUE, paste0(id, '_', vari), id)
      br <- dim.lgd.lis[[vari]]$br
      # Re-plot dimensionlaity plot.
      gg <- ggplot(dat, aes(x=x, y=y, text=dat$feature, label=label)) + geom_point(size=size.pt, alpha=alpha.pt, colour=dat$fill, aes(colour=feature, shape=feature)) + geom_text(size=size.lab.pt, hjust=hjust.lab.pt, vjust=vjust.lab.pt, check_overlap = TRUE) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank(), legend.position=dim.lgd.pos, legend.text=element_text(size=dim.lgd.text.size), legend.margin = margin(t=-0.02, l=0.05, r=0.1, unit='npc'), legend.background = element_rect(fill='transparent'), axis.title=element_text(size=dim.axis.font.size), aspect.ratio=1) + labs(title=tit, x=dim.lgd$labels$x, y=dim.lgd$labels$y) + scale_shape_manual(values=sp, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size))) + guides(colour="none", shape='none') 
    gg.dim.all <- c(gg.dim.all, setNames(list(gg), paste0('dim_', id, '_', vari, '_1')) )
    } 
  } 
  # Convert all reduced dim of ggplots to grobs.
  grob.dim.all <- grob_shm(gg.dim.all, lgd.pos=NULL)
  # Empty list of all reduced dim and SHMs.
  if (profile==FALSE) n.all <- length(vars) + length(gg.lgd.all[1]) 
  if (profile==TRUE) n.all <- length(gg.shm.all) * 2

  dim.shm.gg.lis <- dim.shm.grob.lis <- rep(list(NULL), n.all)
  # Add grobs of dim legend plots.
  for (i in vars) dim.lgd.lis[[i]]$dim.lgd.grob <-  grob_shm(list(dim.lgd.lis[[i]]$dim.lgd), lgd.pos=NULL)[[1]]

  # Assign all SHMs to the empty list.
  if (profile==TRUE) {
    # Assign all reduced dims to the empty list.
    dim.shm.gg.lis[seq(1, n.all, 2)] <- gg.dim.all
    dim.shm.grob.lis[seq(1, n.all, 2)] <- grob.dim.all
    names(dim.shm.gg.lis)[seq(1, n.all, 2)] <- names(dim.shm.grob.lis)[seq(1, n.all, 2)] <- names(grob.dim.all)
    
    dim.shm.gg.lis[seq(2, n.all, 2)] <- gg.shm.all
    dim.shm.grob.lis[seq(2, n.all, 2)] <- grob.shm.all
    names(dim.shm.gg.lis)[seq(2, n.all, 2)] <- names(dim.shm.grob.lis)[seq(2, n.all, 2)] <- names(grob.shm.all)
  } else {
    # Assign all reduced dims to the empty list.
    dim.shm.gg.lis[seq(1, n.all, 2)] <- lapply(vars, function(i) dim.lgd.lis[[i]]$dim.lgd)[1]
    dim.shm.grob.lis[seq(1, n.all, 2)] <- lapply(vars, function(i) dim.lgd.lis[[i]]$dim.lgd.grob)[1]
    names(dim.shm.gg.lis)[seq(1, n.all, 2)] <- names(dim.shm.grob.lis)[seq(1, n.all, 2)] <- names(dim.lgd.lis)
    
    dim.shm.gg.lis[seq(2, n.all, 2)] <- gg.lgd.all[1]
    dim.shm.grob.lis[seq(2, n.all, 2)] <- grob.lgd.all[1]
    names(dim.shm.gg.lis)[seq(2, n.all, 2)] <- names(dim.shm.grob.lis)[seq(2, n.all, 2)] <- names(gg.lgd.all[1])
  }
  return(list(dim.shm.gg.lis=dim.shm.gg.lis, dim.shm.grob.lis=dim.shm.grob.lis, dim.lgd.lis=dim.lgd.lis))
}

#' Assign colors from SVG features to cell labels in embedding plot through a matching list when mapping bulk to cells.
#' 
#' @return A vector.
#' @keywords Internal
#' @noRd

col2cell <- function(shm.col, lis.match) {
  svg.ft.na <- names(lis.match)
  # 'gray80' is a reserved color. 
  cell.labs <- unlist(lis.match)
  dim.col <- rep('gray80', length(cell.labs))
  names(dim.col) <- cell.labs
  for (j in svg.ft.na) { # Colors: svg to cell.
    matched.col <- shm.col[sub('__\\d+', '', names(shm.col))==j][1]
    if (matched.col!='NA') dim.col[lis.match[[j]]] <- matched.col
  }; return(dim.col)
}

