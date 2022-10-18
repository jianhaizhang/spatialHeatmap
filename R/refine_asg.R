#' Refine the bulk-cell assignments
#' 
#' Refine the bulk-cell assignments by including custom bulk-cell assignments. 

#' @param sce.all The coclustering results returned by \code{cocluster}. 
#' @param df.desired.bulk A "data.frame" of desired bulk for some cells. The cells could be specified by providing x-y axis ranges in an embedding plot ("UMAP", "PCA", "TSNE") returned by \code{plot_dim}. E.g. \code{df.desired.bulk <- data.frame(x.min=c(4, -6), x.max=c(5, -5), y.min=c(-2.5, 2), y.max=c(-2, 2.5), desiredSVGBulk=c('CORT', 'STELE'), dimred='UMAP')}, where columns \code{x.min}, \code{x.max}, \code{y.min}, \code{y.max}, \code{desiredSVGBulk}, \code{dimred} are required. In this example, cells located in 4 <= x <= 5 and -2.5 <= y <= -2 in the "UMAP" plot are assigned "STELE", and cells located in -6 <= x <= -5 and 2 <= y <= 2.5 in the "UMAP" plot are assigned "CORT". \cr Alternatively, the "data.frame" could be downloaded from the Shiny app launched by \code{desired_bulk_shiny}. \cr \code{df.desired.bulk} is used to tailor the co-clustering results. That is to say additional true bulk-cell assignments are created and included in the final assignments. If these assignments conflict with the co-colustering results the latter would be overwritten.

#' @return A \code{SingleCellExperiment} of remaining bulk-cell assignments. 

#' @rdname cocluster

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.

#' @export
#' @importFrom SingleCellExperiment reducedDimNames 
#' @importFrom SummarizedExperiment colData 

refine_asg <- function(sce.all, df.desired.bulk=NULL) {
  # save(sce.all, df.desired.bulk, file='sub.asg.arg')
  key <- bulkCell <- NULL
  if (is.null(df.desired.bulk)) return(sce.all)
  if ('key' %in% colnames(df.desired.bulk)) df.desired.bulk <- subset(df.desired.bulk, !key %in% subset(sce.all, , bulkCell=='bulk')$index)
  x <- y <- desiredBulk <- index <- predictor <- total <- true <- NULL
  # Validate labels.
  #sce <- res.lis$cell.refined; 
  labs <- colData(sce.all)$cluster
  er.wa <- dat_fun(as.vector(labs), as.numeric)
  if (is.numeric(er.wa) & length(er.wa) > 0) colData(sce.all)$cluster <- as.character(labs)
  # sce <- subset(sce.all, , bulkCell=='cell')
  # If desired bulk is provided with x-y coordinate range, convert "df.desired.bulk" to the form downloaded from Shiny app.
  xy.ran.cna <- c('x.min', 'x.max', 'y.min', 'y.max') 
  if (any(xy.ran.cna %in% colnames(df.desired.bulk))) {
    if (!all(c(xy.ran.cna, 'desiredBulk', 'dimred') %in% colnames(df.desired.bulk))) stop(paste0('Make sure the "df.desired.bulk" has these columns: ', paste0(xy.ran.cna, collapse=', '), '!'))
    x.min <- df.desired.bulk$x.min
    x.max <- df.desired.bulk$x.max
    y.min <- df.desired.bulk$y.min
    y.max <- df.desired.bulk$y.max
    dimred <- unique(df.desired.bulk$dimred)
    if (!(all(x.min <= x.max) & all(y.min <= y.max))) stop('Make sure x.min <= x.max and y.min <= y.max !')
    if (! dimred %in% reducedDimNames(sce.all)) stop(paste0(dimred, ' is not found!'))
    # Cells are selected in single cell embedding plot, so if bulk tissues are included in df.desired.bulk, they are not counted.  
    gg <- plot_dim(sce.all, dim=dimred, color.by=colnames(colData(sce.all))[1])
    dat <- gg$data; df.desire <- NULL
    dat <- subset(dat, !key %in% subset(sce.all, , bulkCell=='bulk')$index)
    blk <- df.desired.bulk$desiredBulk
    # Assign each bulk to cells in corresponding x-y range.
    for (i in seq_len(nrow(df.desired.bulk))) { 
      df.desired.bulk0 <- df.desired.bulk[i, ]
      x.min0 <- df.desired.bulk0$x.min
      x.max0 <- df.desired.bulk0$x.max
      y.min0 <- df.desired.bulk0$y.min
      y.max0 <- df.desired.bulk0$y.max
      df0 <- subset(dat, x >= x.min0 & x <= x.max0 & y >= y.min0 & y <= y.max0)
      if (nrow(df0)==0) {
        cat('No cells selected for row: ', i, '!\n', sep=''); next
      }
      # Same form as downloaded from Shiny app.
      df.desire0 <- data.frame(x=df0$x, y=df0$y, key=df0$key, desiredBulk=blk[i], dimred=dimred) 
      df.desire <- rbind(df.desire, df.desire0)   
    }
    if (is.null(df.desire)) { 
      warning('No cells selected!'); return(sce.all)
    }
    df.desired.bulk <- df.desire
  }

  # sce$index.cell <- seq_len(ncol(sce))
  # Check "key" and "index.cell" follow the same order in all cells.
  # all(colData(sce)[df.desired.bulk$key, ]$index.cell==df.desired.bulk$key)
  colData(sce.all)[df.desired.bulk$key, ]$assignedBulk <- df.desired.bulk$desiredBulk  
  colData(sce.all)[df.desired.bulk$key, ]$similarity <- 1
  # sce$index.cell <- NULL; sce.all[, sce$index] <- sce
  sce.all$dimred <- unique(df.desired.bulk$dimred)
  return(sce.all)
}
