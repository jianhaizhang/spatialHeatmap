#' Embedding plots of single cells/bulk tissues after co-clustering
#'
#' @param sce A \code{SingleCellExperiment} object with reduced dimensions seen by \code{reducedDimNames(sce)}. 
#' @param dim One of \code{PCA}, \code{UMAP}, \code{TSNE}, the method for reducing dimensionality.
#' @param color.by One of the column names in the \code{colData} slot of \code{sce}.
#' @param group.sel An entry in the \code{color.by} column. All cells under this entry are selected as a group to show.
#' @param row.sel A numeric vector of row numbers in the \code{colData} slot of \code{sce}. The cells corresponding to these rows are highlighted and plotted on top of other cells.
#' @param cocluster.only Logical, only applicable when \code{color.by='cluster'}. If \code{TRUE} (default), only coclusters (including bulk and cells) are colored and the rest are in gray.
#' @param x.break,y.break Two numeric vectors for x, y axis breaks respectively. E.g. \code{seq(-10, 10, 2)}. The default is \code{NULL}.
#' @param panel.grid Logical. If \code{TRUE}, the panel grid will be shown.
#' @param lgd.title.size,lgd.key.size,lgd.text.size The size of legend plot title, legend key, legend text respectively.
#' @param point.size,bulk.size The size of cells and bulk tissues respectively.
#' @param alpha The transparency of cells and bulk tissues. The default is 0.6.
#' @param stroke,bulk.stroke The line width of cells and bulk tissues respectively.
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
#' @importFrom ggplot2 ggplot scale_fill_manual scale_size_manual geom_point aes theme_classic theme element_text labs scale_x_continuous scale_y_continuous guides guide_legend

plot_dim <- function(sce, dim=NULL, color.by, group.sel=NULL, row.sel=NULL, cocluster.only=TRUE, x.break=NULL, y.break=NULL, panel.grid=FALSE, lgd.title.size=13, lgd.key.size=0.03, lgd.text.size=12, point.size=3, bulk.size=5, alpha=0.7, stroke=0.2, bulk.stroke=1, axis.text.size=10, axis.title.size=11) {
  # All scenarios: 1. Only cell, select by row, group, or all. 2. Bulk and cell, select by row, group, or all.
  x <- y <- key <- colour_by <- bulkCell <- NULL
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
  df.all$assignedBulk <- cdat$assignedBulk
  
  if (blk.cell) {
    df.all$bulkCell <- cdat$bulkCell
    clus.all <- unique(df.all$colour_by)
    clus.cell <- clus.all[!clus.all %in% subset(df.all, bulkCell=='bulk')$colour_by]
    coclus <- setdiff(clus.all, clus.cell)  
  }
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
  ae <- aes(fill=colour_by)
  lab <- labs(x=labs[1], y=labs[2], fill=color.by)
  row.only <- is.null(group.sel) & length(row.sel) < ncol(sce)
  row.grp <- (!is.null(row.sel) | !is.null(group.sel)) & length(row.sel) < ncol(sce)

  # Select by rows or groups.
  if (row.grp) {
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
    } else { # Select by row/group in the scenario of only cell
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
      scl.col <- scale_fill_manual(values=col.all, breaks=sel, labels=sub('^.*selected\\.', '', sel))
    } else if (blk.cell & is.null(group.sel)) {# Select by row in the scenario of bulk + cell
      sel0 <- c(grep('bulk', sel, value=TRUE)[1], grep('cell', sel, value=TRUE)[1])
      sel0 <- sel0[!is.na(sel0)]
      scl.col <- scale_fill_manual(values=col.all, breaks=sel0, labels=sub('^.*selected\\.', '', sel0)) 
    } else { 
      if (!is.null(group.sel)) { # Select by row/group in the scenario of only cell  
        scl.col <- scale_fill_manual(name='Selected', values=col.all, breaks=sel, labels=sub('\\.selected$', '', sel))
      } else { # Select by row in the scenario of bulk + cell
        scl.col <- scale_fill_manual(name='Selected', values=col.all, breaks=sel[1], labels=sub('.*', '', sel[1]))
      }
    }
  } else scl.col <- NULL

  if (blk.cell) { # In the scenario of bulk + cell
    df.all <- rbind(subset(df.all, bulkCell=='cell'), subset(df.all, bulkCell=='bulk')) 
    if (is.null(group.sel) & length(row.sel) == ncol(sce)) { # All bulk and cells.     
      # Size
      sz.all <- c(bulk=bulk.size, cell=point.size)
      scl.sz <- scale_size_manual(values=sz.all, breaks=names(sz.all), labels=names(sz.all))
      ae <- aes(fill=colour_by, size=bulkCell)
      lab <- labs(x=labs[1], y=labs[2], fill=ifelse(color.by=='bulkCell', 'sample', color.by), shape='sample', size='sample')
   } else if (!is.null(group.sel)) { # Select a co-cluster.
      # Size
      sz.all <- c(bulk=bulk.size, cell=point.size)
      scl.sz <- scale_size_manual(values=sz.all, breaks=names(sz.all), labels=names(sz.all))
      ae <- aes(fill=colour_by, size=bulkCell)
      if (all(c('bulk', 'cell') %in% df.sel$bulkCell)) { 
        # Fill -> colour_by, size -> bulkCell: fill and size do not conflict on the same point, so they can be merged in legend though they are referring to different varialbles.
        lab.fil.sz <- group.sel
        df.coclus.blk.sel <- subset(df.all, colour_by %in% sel & bulkCell=='bulk')
        rep0 <- table(paste0(df.coclus.blk.sel$colour_by, ': ', df.coclus.blk.sel$sample))
        if (all(rep0==1)) { # No bulk replicates.
         lab.fil.sz <- paste0(group.sel, ': ', df.coclus.blk.sel$sample)
        }
        lab <- labs(x=labs[1], y=labs[2], fill=lab.fil.sz, size=lab.fil.sz)
      } else if ('cell' %in% df.sel$bulkCell) {
        lab <- labs(x=labs[1], y=labs[2], fill=group.sel, size='sample') 
      } else if ('bulk' %in% df.sel$bulkCell) {
        lab <- labs(x=labs[1], y=labs[2], fill=group.sel, size='sample')
      }
    } else if (row.only) { # Select by row.
      # Size
      sz.all <- c(bulk=bulk.size, cell=point.size)
      ae <- aes(fill=colour_by, shape=colour_by, size=bulkCell)
      scl.sz <- scale_size_manual(values=sz.all, breaks=names(sz.all), labels=names(sz.all))
      if (all(c('bulk', 'cell') %in% df.sel$bulkCell)) { 
        lab <- labs(x=labs[1], y=labs[2], fill='Selected', size='Selected')
      } else if ('cell' %in% df.sel$bulkCell) {
        lab <- labs(x=labs[1], y=labs[2], fill='Selected', size='sample') 
      } else if ('bulk' %in% df.sel$bulkCell) {
        lab <- labs(x=labs[1], y=labs[2], fill='Selected', size='sample')
      }
    }
  } else { scl.sz <- NULL }

  thm <- theme(plot.title=element_text(hjust=0.5, size=14), legend.title=element_text(size=lgd.title.size), legend.key.size=unit(lgd.key.size, "npc"), legend.text=element_text(size=lgd.text.size), axis.text=element_text(size=axis.text.size), axis.title=element_text(size=axis.title.size, face="plain"))
  
  # Stroke
  sk <- rep(stroke, nrow(df.all))
  sk[df.all$bulkCell == 'bulk'] <- bulk.stroke
  df.all$stroke <- sk
  gg <- ggplot(df.all, aes(x=x, y=y, key=key)) + geom_point(colour='black', shape=21, alpha=alpha, stroke=df.all$stroke, ae) + thm + lab + scl.col + scl.sz

  if (cocluster.only==TRUE & is.null(group.sel) & length(row.sel) == ncol(sce) & blk.cell & color.by=='cluster') {
    df1 <- gg$data; df2 <- layer_data(gg)
    int.na <- intersect(colnames(df1), colnames(df2))
    df.all <- cbind(df1[, !colnames(df1) %in% int.na], layer_data(gg))
    # All default colors.
    df.col <- subset(df.all, !duplicated(colour_by))
    clus.all <- df.col$colour_by
    col.all <- df.col$fill; names(col.all) <- clus.all
    col.all[clus.cell] <- 'gray80'
    scl.col <- scale_fill_manual(values=col.all, breaks=coclus, labels=coclus)

    df.coclus.blk <- subset(df.all, colour_by %in% coclus & bulkCell=='bulk')
    lab.fil <- table(paste0(df.coclus.blk$colour_by, ': ', df.coclus.blk$sample))
    if (all(lab.fil==1)) { # No bulk replicates.
      # If multiple bulk are co-clustered in the same cluster, only one bulk is indicated in the legend.
      scl.col <- scale_fill_manual(values=col.all, breaks=sub(': .*$', '', names(lab.fil)), labels=names(lab.fil))
    }
    # Size
    sz.all <- c(bulk=bulk.size, cell=point.size)
    scl.sz <- scale_size_manual(values=sz.all, breaks=names(sz.all), labels=names(sz.all))
    ae <- aes(fill=colour_by, size=bulkCell)
    lab <- labs(x=labs[1], y=labs[2], fill=color.by, size='sample')
    df.all <- rbind(subset(df.all, colour_by %in% clus.cell), subset(df.all, colour_by %in% coclus))
    df.all <- rbind(subset(df.all, bulkCell=='cell'), subset(df.all, bulkCell=='bulk'))

    gg <- ggplot(df.all, aes(x=x, y=y, key=key)) + geom_point(colour='black', alpha=alpha, stroke=df.all$stroke, shape=21, ae) + thm + lab + scl.col + scl.sz
  }
  gg <- gg + guides(fill = guide_legend(override.aes = list(size=point.size)))
  # Overwrite theme.
  if (panel.grid==FALSE) gg <- gg + theme_classic() + thm + lab + scl.col + scl.sz 
  if (!is.null(x.break) & is.null(y.break)) gg <- gg + scale_x_continuous(breaks=x.break)
  if (!is.null(y.break) & is.null(x.break)) gg <- gg + scale_y_continuous(breaks=y.break)
  if (!is.null(x.break) & !is.null(y.break)) gg <- gg + scale_x_continuous(breaks=x.break) + scale_y_continuous(breaks=y.break)
  gg
}
