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
#' @importFrom ggplot2 layer_data ggplot geom_point theme_classic theme element_text element_blank labs scale_color_manual scale_shape_manual margin guide_legend element_rect
#' @importFrom stats setNames

dim_color_idp <- function(sce=NULL, row.sel=NULL, covis.type=NULL, targ=NULL, ID, gene, geneV, cols, profile, tar.cell=NULL, tar.bulk=NULL, con.na.cell, gg.dim, gg.shm.all, grob.shm.all, col.shm.all, gg.lgd.all, col.lgd.all, grob.lgd.all, cell.group, lis.match=NULL, sub.title.size=11, dim.lgd.pos='bottom', dim.lgd.nrow=2, dim.lgd.key.size=4, dim.lgd.text.size=13, dim.axis.font.size=10, alpha.pt=0.8, shape=NULL, lgd.plot.margin=margin(t=0.01, r=0.01, b=0.01, l=0.01, unit="npc")) {
  # save(sce, row.sel, covis.type, targ, ID, gene, geneV, cols, profile, tar.cell, tar.bulk, con.na.cell, gg.dim, gg.shm.all, grob.shm.all, col.shm.all, gg.lgd.all, col.lgd.all, grob.lgd.all, cell.group, lis.match, sub.title.size, dim.lgd.pos, dim.lgd.nrow, dim.lgd.key.size, dim.lgd.text.size, alpha.pt, shape, dim.axis.font.size, lgd.plot.margin, file='dim.color.idp.arg')
  x <- y <- fill <- feature <- NULL
  if (is.null(lgd.plot.margin)) lgd.plot.margin <- margin(t=0.01, r=0.01, b=0.01, l=0.01, unit="npc")
  auto <- covis.type %in% c('toBulkAuto', 'toCellAuto')
  if (auto) {
    cdat <- colData(sce)
    # The matching list between aggregated cells and aSVG spatial features. The former are cells with a source tissue assignment in co-clustering.   
    blk.uni <- unique(cdat$assignedBulk)
    blk.uni <- setdiff(blk.uni, 'none')
    lis.match <- as(blk.uni, 'list'); names(lis.match) <- blk.uni
    targ <- setdiff(targ, 'none')
    if (any(!targ %in% blk.uni)) stop("Make sure all entries in 'targ' are in 'assignedBulk'!")
  } 
  lis.match <- lis.match[!unlist(lapply(lis.match, is.null))]
  # Non-target cells.
  if (auto) tar.cells <- targ else { 
    if (!is.null(tar.cell)) tar.cells <- unique(names(lis.match)) else if (!is.null(tar.bulk)) tar.cells <- unique(unlist(lis.match))
  }
   vars.cell <- names(gg.dim)
   # Transfer colors from SHM to dim plots according to match list.
   dim.lgd.lis <- lapply(seq_along(gg.dim), function(i) { 
     gg.dim0 <- gg.dim[i]
     names(gg.dim0) <- paste0('dim_', names(gg.lgd.all))
     # Extract data from ggplots.
     gg.dat <- gg_dat(gg.dim0[[1]]); lis <- list(gg.dat=gg.dat)
     # Transfer colors from SHM to dim plots.
     if (!is.null(tar.bulk)|auto) { 
       dim.col <- col_dim_tocell(gg.dim0, gcol.all=col.lgd.all, lis.match) 
     } else if (!is.null(tar.cell)) { 
       dim.col <- col_dim_toblk(gg.dim0, gcol.all=col.lgd.all, lis.match)
     }
     # Non-target cells have color gray80.
     cell.all <- unique(gg.dat$feature)
     dim.col.all <- setNames(rep('gray80', length(cell.all)), cell.all)
     dim.col.all[names(dim.col)] <- dim.col 
     gg.dat$fill <- dim.col.all[gg.dat$feature] 
     gg.dat <- rbind(subset(gg.dat, fill == 'gray80'), subset(gg.dat, fill != 'gray80'))

     # Legal shapes: c(0:25, 32:127)
     #sp.sel <- c(15:18, 7:14); sp.all <- c(0, 2:25, 32:127)
     #sp.all <- c(sp.sel, setdiff(sp.all, sp.sel))
     # Cell cluster shapes.
     #sp <- sp.all[seq_along(cell.all)]; names(sp) <- cell.all
     sp <- shp(shape, cell.all)
     non.tar <- setdiff(cell.all, tar.cells) 
     if (length(non.tar) > 0) br <- tar.cells else br <- cell.all
     tit <- NULL; if (TRUE %in% con.na.cell) tit <- vars.cell[i] 
     # Re-plot dimensionlaity plot.
     dim.lgd <- ggplot(gg.dat, aes(x=x, y=y, text=gg.dat$feature)) + geom_point(size=2, alpha=alpha.pt, aes(colour=feature, shape=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), legend.position=dim.lgd.pos, legend.text=element_text(size=dim.lgd.text.size), legend.margin = margin(t=-0.02, l=0.05, r=0.1, unit='npc'), plot.margin = lgd.plot.margin, legend.background = element_rect(fill='transparent'), axis.text = element_blank(), axis.ticks = element_blank(), axis.title=element_text(size=dim.axis.font.size)) + scale_color_manual(values=dim.col.all, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow)) + scale_shape_manual(values=sp, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size))) + labs(title=tit, x=gg.dim0[[1]]$labels$x, y=gg.dim0[[1]]$labels$y, colour=cell.group, shape=cell.group) 
     lis <- c(lis, list(sp=sp), list(dim.lgd=dim.lgd)); lis
   }); names(dim.lgd.lis) <- names(gg.dim)
  # Dim plots per ID per variable.
  gg.dim.all <- NULL; for (id in ID) {
    for (vari in names(dim.lgd.lis)) {
      # Re-use data extracted from original dim plots
      gg.dat <- dim.lgd.lis[[vari]]$gg.dat 
      # Re-use shapes in dim legend plots.
      sp <- dim.lgd.lis[[vari]]$sp
      dim.lgd <- dim.lgd.lis[[vari]]$dim.lgd
      gg.dat$fill <- 'gray80'
      # Non-target cells have color gray80.
      idx.mat <- gg.dat$feature %in% tar.cells
      # Isolate assay data by ID and variable.
      gene0 <- gene[, grepl(paste0('__', vari, '$'), colnames(gene)), drop=FALSE]
      # Order in gene0 is the same with gg.dat.
      color.dat <- value2color(gene0[id, idx.mat, drop=FALSE], geneV, cols)
      # Colors of target cells.
      gg.dat$fill[idx.mat] <- color.dat
      gg.dat <- rbind(subset(gg.dat, fill == 'gray80'), subset(gg.dat, fill != 'gray80'))
      tit <- ifelse(con.na.cell==TRUE, paste0(id, '_', vari), id)
      cell.all <- unique(gg.dat$feature)
      non.tar <- setdiff(cell.all, tar.cells) 
      if (length(non.tar) > 0) br <- tar.cells else br <- cell.all
      br  <- unlist(ifelse(length(non.tar) > 0, list(tar.cells), list(cell.all)))
      # Re-plot dimensionlaity plot.
      gg <- ggplot(gg.dat, aes(x=x, y=y, text=gg.dat$feature)) + geom_point(size=2, alpha=alpha.pt, colour=gg.dat$fill, aes(colour=feature, shape=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank(), legend.position=dim.lgd.pos, legend.text=element_text(size=dim.lgd.text.size), legend.margin = margin(t=-0.02, l=0.05, r=0.1, unit='npc'), legend.background = element_rect(fill='transparent'), axis.title=element_text(size=dim.axis.font.size)) + labs(title=tit, x=dim.lgd$labels$x, y=dim.lgd$labels$y) + scale_shape_manual(values=sp, breaks=br, guide=guide_legend(title=NULL, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size))) + guides(colour="none", shape='none') 
    gg.dim.all <- c(gg.dim.all, setNames(list(gg), paste0('dim_', id, '_', vari, '_1')) )
    } 
  }
 # Convert all reduced dim of ggplots to grobs.
  grob.dim.all <- grob_shm(gg.dim.all, lgd.pos=NULL)
  # Empty list of all reduced dim and SHMs. 
  if (profile==FALSE) n.all <- length(vars.cell) + length(gg.lgd.all[1]) 
  if (profile==TRUE) n.all <- length(gg.shm.all) * 2 
  dim.shm.gg.lis <- dim.shm.grob.lis <- rep(list(NULL), n.all) 
  # Add grobs of dim legend plots.
  for (i in vars.cell) dim.lgd.lis[[i]]$dim.lgd.grob <-  grob_shm(list(dim.lgd.lis[[i]]$dim.lgd), lgd.pos=NULL)[[1]] 
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
    dim.shm.gg.lis[seq(1, n.all, 2)] <- lapply(vars.cell, function(i) dim.lgd.lis[[i]]$dim.lgd)[1]
    dim.shm.grob.lis[seq(1, n.all, 2)] <- lapply(vars.cell, function(i) dim.lgd.lis[[i]]$dim.lgd.grob)[1]
    names(dim.shm.gg.lis)[seq(1, n.all, 2)] <- names(dim.shm.grob.lis)[seq(1, n.all, 2)] <- names(dim.lgd.lis)
    dim.shm.gg.lis[seq(2, n.all, 2)] <- gg.lgd.all[1]
    dim.shm.grob.lis[seq(2, n.all, 2)] <- grob.lgd.all[1]
    names(dim.shm.gg.lis)[seq(2, n.all, 2)] <- names(dim.shm.grob.lis)[seq(2, n.all, 2)] <- names(gg.lgd.all[1]) 
  }
  return(list(dim.shm.gg.lis=dim.shm.gg.lis, dim.shm.grob.lis=dim.shm.grob.lis, dim.lgd.lis=dim.lgd.lis))
}

#' Extract data from ggplot.
#' 
#' @return A vector.
#' @keywords Internal
#' @noRd 

#' @importFrom ggplot2 layer_data

gg_dat <- function(gg, shm=TRUE) {
  dat <- gg$data; lay.dat <- layer_data(gg) 
  if (shm==TRUE) {
    dat.all <- cbind(lay.dat, colour_by=dat$colour_by)
    dat.all <- dat.all[, !colnames(dat.all) %in% 'colour']
    colnames(dat.all)[colnames(dat.all)=='colour_by'] <- 'feature'
  } else dat.all <- cbind(dat, lay.dat)
  return(dat.all)
}



