#' Plotting single cells in reduced dimensionalities 
#'
#' @param sce A \code{SingleCellExperiment} object with reduced dimensions seen by \code{reducedDimNames(sce)}. 
#' @param dim One of \code{PCA}, \code{UMAP}, \code{TSNE}, the method for reducing dimensionality.
#' @param color.by One of the column names in the \code{colData} slot of \code{sce}.
#' @param group.sel An entry in the \code{color.by} column. All cells under this entry are selected as a group to show.
#' @param row.sel A numeric vector of row numbers in the \code{colData} slot of \code{sce}. The cells corresponding to these rows are highlighted and plotted on top of other cells.
#' @param x.break,y.break Two numeric vectors for x, y axis breaks respectively. E.g. \code{seq(-10, 10, 2)}. The default is \code{NULL}.
#' @param panel.grid Logical. If \code{TRUE}, the panel grid will be shown.

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

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom ggplot2 ggplot scale_color_manual geom_point aes theme_classic theme element_text labs scale_x_continuous scale_y_continuous

#' @export plot_dim

plot_dim <- function(sce, dim=NULL, color.by, group.sel=NULL, row.sel=NULL, x.break=NULL, y.break=NULL, panel.grid=FALSE) {
  x <- y <- key <- colour_by <- NULL
  if (is.null(dim)) dim <- reducedDimNames(sce)[1]
  dim.all <- reducedDim(sce, dim)
  # First two dims.
  df.all <- dim.all[, c(1, 2)]; # rna <- rownames(df.all)
  df.all <- as.data.frame(df.all)
  colnames(df.all) <- c('x', 'y') 
  df.all$colour_by <- colData(sce)[, color.by]
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
  if (!is.null(group.sel)) row.sel <- colData(sce)[[color.by]]==group.sel

  if (!is.null(row.sel)) {
    df.sel <- df.all[row.sel, ]
    df.sel$colour_by <- paste0(df.sel$colour_by, '.selected')
    sel <- unique(df.sel$colour_by)
    col.sel <- rep('blue', length(sel))
    names(col.sel) <- sel
    df.other <- df.all[setdiff(df.all$key, row.sel), ] 
    other <- unique(df.other$colour_by)
    col.other <- rep('gray80', length(other))
    names(col.other) <- other
    col.all <- c(col.other, col.sel)
    df.all <- rbind(df.other, df.sel)
    scl.col <- scale_color_manual(name='Selected', values=col.all, breaks=sel, labels=sub('\\.selected$', '', sel))
  } else scl.col <- NULL
  gg <- ggplot(df.all, aes(x=x, y=y, key=key)) + geom_point(size=1.5, alpha=0.6, shape=19, stroke=0.5, aes(colour=colour_by)) + theme(plot.title=element_text(hjust=0.5, size=14)) + labs(x=labs[1], y=labs[2], colour=color.by) + scl.col
  if (panel.grid==FALSE) gg <- gg + theme_classic() 
  if (!is.null(x.break) & is.null(y.break)) gg <- gg + scale_x_continuous(breaks=x.break)
  if (!is.null(y.break) & is.null(x.break)) gg <- gg + scale_y_continuous(breaks=y.break)
  if (!is.null(x.break) & !is.null(y.break)) gg <- gg + scale_x_continuous(breaks=x.break) + scale_y_continuous(breaks=y.break)
  gg
}
