#' Embedding plots of single cells/bulk tissues after co-clustering
#'
#' @param sce A \code{SingleCellExperiment} object with reduced dimensions seen by \code{reducedDimNames(sce)}. 
#' @param dim One of \code{PCA}, \code{UMAP}, \code{TSNE}, the method for reducing dimensionality.
#' @param color.by One of the column names in the \code{colData} slot of \code{sce}.
#' @param group.sel An entry in the \code{color.by} column. All cells under this entry are selected as a group to show.
#' @param row.sel A numeric vector of row numbers in the \code{colData} slot of \code{sce}. The cells corresponding to these rows are highlighted and plotted on top of other cells.
#' @param cocluster.only Logical, if \code{TRUE} (default), only coclusters (including bulk and cells) are colored and the rest are in gray.
#' @param x.break,y.break Two numeric vectors for x, y axis breaks respectively. E.g. \code{seq(-10, 10, 2)}. The default is \code{NULL}.
#' @param panel.grid Logical. If \code{TRUE}, the panel grid will be shown.
#' @param lgd.title.size,lgd.key.size,lgd.text.size The size of legend plot title, legend key, legend text respectively.
#' @param point.size,bulk.size The size of cells and bulk tissues respectively.
#' @param alpha The transparency of cells and bulk tissues. The default is 0.6.
#' @param stroke The line width of cells and bulk tissues.
#' @param axis.text.size,axis.title.size The size of axis text and title respectively.

#' @return An object of ggplot.

#' @examples

#' library(scran); library(scuttle) 
#' sce <- mockSCE(); sce <- logNormCounts(sce) 
#' # Modelling the variance.
#' var.stats <- modelGeneVar(sce) 
#' sce <- denoisePCA(sce, technical=var.stats, subset.row=rownames(var.stats)) 
#' plot_dim(sce, dim='PCA', color.by='Cell_Cycle')

#' # See function "coclus_meta" by running "?coclus_meta".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.


#' @export 
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom ggplot2 ggplot scale_color_manual geom_point aes theme_classic theme element_text labs scale_x_continuous scale_y_continuous

plot_dim <- function(sce, dim=NULL, color.by, group.sel=NULL, row.sel=NULL, cocluster.only=TRUE, x.break=NULL, y.break=NULL, panel.grid=FALSE, lgd.title.size=10, lgd.key.size=0.03, lgd.text.size=10, point.size=1.5, bulk.size=3, alpha=0.6, stroke=0.5, axis.text.size=10, axis.title.size=11) {
  # All scenarios: 1. Only cell, select by row, group, or all. 2. Bulk and cell, select by row, group, or all.
  x <- y <- key <- colour_by <- NULL
  cdat <- colData(sce); blk.cell <- FALSE
  if (all(c('bulk', 'cell') %in% cdat$bulkCell)) blk.cell <- TRUE

  # Construct data for plotting.
  if (is.null(dim)) dim <- reducedDimNames(sce)[1]
  dim.all <- reducedDim(sce, dim)
  # First two dims.
  df.all <- dim.all[, c(1, 2)]; # rna <- rownames(df.all)
  df.all <- as.data.frame(df.all)
  colnames(df.all) <- c('x', 'y') 
  df.all$colour_by <- cdat[, color.by]
  df.all$sample <- cdat$sample
  if (blk.cell) df.all$bulkCell <- cdat$bulkCell
  # Used in event_data/ggplotly.
  df.all$key <- seq_len(nrow(df.all))
  labs <- paste0(dim, c(1, 2))
  # Percentage of variance explained.
  per.var <- attr(dim.all, 'percentVar')
  if (!is.null(per.var)) { # Axis labels.
    per <- paste(round(per.var[seq_len(2)], 0), "%", sep="") 
    labs <- paste0(labs, ' (', per, ')')
  }

  if (!is.null(group.sel) & !is.null(row.sel)) stop('At least one of "group.sel" and "row.sel" should be "NULL"!')
  if (is.null(group.sel) & is.null(row.sel)) row.sel <- seq_len(ncol(sce))
  if (!is.null(group.sel)) row.sel <- which(cdat[[color.by]]==group.sel)
  ae <- aes(colour=colour_by)
  lab <- labs(x=labs[1], y=labs[2], colour=color.by)

  # Select by rows or groups.
  if (!is.null(row.sel) & length(row.sel) < ncol(sce)) {
    df.sel <- df.all[row.sel, ]
    if (blk.cell & !is.null(group.sel)) { # Select by group in the scenario of bulk + cell.
      df.sel$colour_by <- paste0(df.sel$colour_by, '.selected.', df.sel$bulkCell)
      sel <- unique(df.sel$colour_by)
      sam <- unique(df.sel$bulkCell)
      # Bulk: red, cell: blue
      if (all(c('bulk', 'cell') %in% sam)) {
        col.sel <- c('red', 'blue')
        names(col.sel)[1] <- grep('bulk', sel, value=TRUE)
        names(col.sel)[2] <- grep('cell', sel, value=TRUE)
      } else if ('bulk' %in% sam) { 
        col.sel <- 'red'; names(col.sel) <- sel
      } else if ('cell' %in% sam) {
        col.sel <- 'blue'; names(col.sel) <- sel 
      }
    } else if (blk.cell & is.null(group.sel)) { # Select by row in the scenario of bulk + cell
      df.sel$colour_by <- paste0(df.sel$colour_by, '.selected.', df.sel$bulkCell)
      sel <- unique(df.sel$colour_by)
      sam <- unique(df.sel$bulkCell)
      # Bulk: red, cell: blue
      if (all(c('bulk', 'cell') %in% sam)) {
        blk <- grep('bulk', sel, value=TRUE)
        blk.col <- rep('red', length(blk))
        names(blk.col) <- blk
        cell <- grep('cell', sel, value=TRUE)
        cell.col <- rep('blue', length(cell))
        names(cell.col) <- cell
        col.sel <- c(blk.col, cell.col)
      } else if ('bulk' %in% sam) {
        col.sel <- rep('red', length(sel))
        names(col.sel) <- sel
      } else if ('cell' %in% sam) {
        col.sel <- rep('blue', length(sel))
        names(col.sel) <- sel
      } 
    } else { # Select by group in the scenario of only cell
      df.sel$colour_by <- paste0(df.sel$colour_by, '.selected')
      sel <- unique(df.sel$colour_by)
      col.sel <- rep('blue', length(sel)); names(col.sel) <- sel 
    }
    # Non-target cell/bulk.
    df.other <- df.all[setdiff(df.all$key, row.sel), ] 
    other <- unique(df.other$colour_by)
    col.other <- rep('gray80', length(other))
    names(col.other) <- other
    col.all <- c(col.other, col.sel)
    df.all <- rbind(df.other, df.sel)
    if (blk.cell & !is.null(group.sel)) { # Select by group in the scenario of bulk + cell
      scl.col <- scale_color_manual(values=col.all, breaks=sel, labels=sub('^.*selected\\.', '', sel)) 
    } else if (blk.cell & is.null(group.sel)) {# Select by row in the scenario of bulk + cell
      sel0 <- c(grep('bulk', sel, value=TRUE)[1], grep('cell', sel, value=TRUE)[1])
      sel0 <- sel0[!is.na(sel0)]
      scl.col <- scale_color_manual(values=col.all, breaks=sel0, labels=sub('^.*selected\\.', '', sel0)) 
    } else { 
      if (!is.null(group.sel)) { # Select by group in the scenario of only cell  
        scl.col <- scale_color_manual(name='Selected', values=col.all, breaks=sel, labels=sub('\\.selected$', '', sel))
      } else { # Select by row in the scenario of bulk + cell
        scl.col <- scale_color_manual(name='Selected', values=col.all, breaks=sel[1], labels=sub('.*', '', sel[1]))
      }
    }
  } else scl.col <- NULL

  if (blk.cell) { # In the scenario of bulk + cell
    df.all <- rbind(subset(df.all, bulkCell=='cell'), subset(df.all, bulkCell=='bulk')) 
    if (is.null(group.sel) & length(row.sel) == ncol(sce)) { # All bulk and cells.     
      # Shape and size refer to the same variable.
      df.all$sample.shape <- paste0(df.all$sample, '.', df.all$bulkCell)
      sam.sp <- unique(df.all$sample.shape)
      blk <- grep('bulk', sam.sp, value=TRUE)  
      cell <- grep('cell', sam.sp, value=TRUE)  
      # Shape
      sp.blk <- rep(17, length(blk)); names(sp.blk) <- blk
      sp.cell <- rep(19, length(cell)); names(sp.cell) <- cell
      sp.all <- c(sp.blk, sp.cell)
      scl.sp <- scale_shape_manual(name='sample', values=sp.all, breaks=c(blk[1], cell[1]), labels=sub('^.*\\.', '', c(blk[1], cell[1])))
      # Size
      sz.all <- rep(point.size, length(sp.all))
      sp.all.na <- names(sp.all) 
      names(sz.all) <- sp.all.na
      sz.all[grep('bulk$', sp.all.na)] <- bulk.size
      scl.sz <- scale_size_manual(values=sz.all, breaks=c(blk[1], cell[1]), labels=sub('^.*\\.', '', c(blk[1], cell[1])))
      ae <- aes(colour=colour_by, shape=sample.shape, size=sample.shape)
      lab <- labs(x=labs[1], y=labs[2], colour=color.by, shape='sample', size='sample')
   } else if (!is.null(group.sel)) { # Select a co-cluster.
      # Shape
      col.all.na <- names(col.all)
      sp.all <- rep(19, length(col.all))
      names(sp.all) <- col.all.na
      sp.all[grep('bulk$', col.all.na)] <- 17 
      scl.sp <- scale_shape_manual(values=sp.all, breaks=sel, labels=sub('^.*selected\\.', '', sel))
      # Size
      sz.all <- rep(point.size, length(col.all)) 
      names(sz.all) <- col.all.na
      sz.all[grep('bulk$', col.all.na)] <- bulk.size
      scl.sz <- scale_size_manual(values=sz.all, breaks=sel, labels=sub('^.*selected\\.', '', sel))
      ae <- aes(colour=colour_by, shape=colour_by, size=colour_by)
      lab <- labs(x=labs[1], y=labs[2], colour=group.sel, shape=group.sel, size=group.sel)
    } else if (is.null(group.sel) & length(row.sel) < ncol(sce)) { # Select by row.
      # Shape
      col.all.na <- names(col.all)
      sp.all <- rep(19, length(col.all))
      names(sp.all) <- col.all.na
      sp.all[grep('bulk$', col.all.na)] <- 17
      sel0 <- c(grep('bulk', sel, value=TRUE)[1], grep('cell', sel, value=TRUE)[1])
      sel0 <- sel0[!is.na(sel0)]
      scl.sp <- scale_shape_manual(values=sp.all, breaks=sel0, labels=sub('^.*selected\\.', '', sel0))
      # Size
      sz.all <- rep(point.size, length(col.all))
      names(sz.all) <- col.all.na
      sz.all[grep('bulk$', col.all.na)] <- bulk.size
      scl.sz <- scale_size_manual(values=sz.all, breaks=sel0, labels=sub('^.*selected\\.', '', sel0))
      ae <- aes(colour=colour_by, shape=colour_by, size=colour_by)
      lab <- labs(x=labs[1], y=labs[2], colour='Selected', shape='Selected', size='Selected')
    }
  } else { scl.sp <- NULL; scl.sz <- NULL }

  thm <- theme(plot.title=element_text(hjust=0.5, size=14), legend.title=element_text(size=lgd.title.size), legend.key.size=unit(lgd.key.size, "npc"), legend.text=element_text(size=lgd.text.size), axis.text=element_text(size=axis.text.size), axis.title=element_text(size=axis.title.size, face="plain"))
  gg <- ggplot(df.all, aes(x=x, y=y, key=key)) + geom_point(alpha=alpha, stroke=stroke, ae) + thm + lab + scl.col + scl.sp + scl.sz

  if (cocluster.only==TRUE & is.null(group.sel) & length(row.sel) == ncol(sce) & 'bulk' %in% sce$bulkCell) {
    df1 <- gg$data; df2 <- layer_data(gg)
    int.na <- intersect(colnames(df1), colnames(df2))
    df.all <- cbind(df1[, !colnames(df1) %in% int.na], layer_data(gg))
    # All default colors.
    df.col <- subset(df.all, !duplicated(colour_by))
    clus.all <- df.col$colour_by
    col.all <- df.col$colour; names(col.all) <- clus.all
    clus.cell <- clus.all[!clus.all %in% subset(df.all, bulkCell=='bulk')$colour_by]
    col.all[clus.cell] <- 'gray80'
    coclus <- setdiff(clus.all, clus.cell)
    scl.col <- scale_color_manual(values=col.all, breaks=coclus, labels=coclus) 

    sam.sp <- unique(df.all$sample.shape)
    blk <- grep('bulk', sam.sp, value=TRUE)  
    cell <- grep('cell', sam.sp, value=TRUE)  
    # Shape
    sp.blk <- rep(17, length(blk)); names(sp.blk) <- blk
    sp.cell <- rep(19, length(cell)); names(sp.cell) <- cell
    sp.all <- c(sp.blk, sp.cell)
    scl.sp <- scale_shape_manual(name='sample', values=sp.all, breaks=c(blk[1], cell[1]), labels=sub('^.*\\.', '', c(blk[1], cell[1])))
    # Size
    sz.all <- rep(point.size, length(sp.all))
    sp.all.na <- names(sp.all) 
    names(sz.all) <- sp.all.na
    sz.all[grep('bulk$', sp.all.na)] <- bulk.size
    scl.sz <- scale_size_manual(values=sz.all, breaks=c(blk[1], cell[1]), labels=sub('^.*\\.', '', c(blk[1], cell[1])))
    ae <- aes(colour=colour_by, shape=sample.shape, size=sample.shape)
    lab <- labs(x=labs[1], y=labs[2], colour=color.by, shape='sample', size='sample')
    df.all <- rbind(subset(df.all, colour_by %in% clus.cell), subset(df.all, colour_by %in% coclus))
    df.all <- rbind(subset(df.all, bulkCell=='cell'), subset(df.all, bulkCell=='bulk'))

    gg <- ggplot(df.all, aes(x=x, y=y, key=key)) + geom_point(alpha=alpha, stroke=stroke, ae) + thm + lab + scl.col + scl.sp + scl.sz

  }
  # Overwrite theme.
  if (panel.grid==FALSE) gg <- gg + theme_classic() + thm + lab + scl.col + scl.sp + scl.sz 
  if (!is.null(x.break) & is.null(y.break)) gg <- gg + scale_x_continuous(breaks=x.break)
  if (!is.null(y.break) & is.null(x.break)) gg <- gg + scale_y_continuous(breaks=y.break)
  if (!is.null(x.break) & !is.null(y.break)) gg <- gg + scale_x_continuous(breaks=x.break) + scale_y_continuous(breaks=y.break)
  gg
}
