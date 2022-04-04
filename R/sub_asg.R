#' Subset the bulk-cell assignments
#' 
#' Subset the bulk-cell assignments according to a threshold, which is a similarity value between bulk and cells. 

#' @param res.lis The result list of coclustering, which is the output of \code{tests} and comprises three slots \code{sce}, \code{roc.obj}, \code{df.roc}.
#' @param thr The threshold for subsetting bulk-cell assignments, which is a similarity value (Pearson's or Spearman's correlation coefficient) between bulk and cells. Only bulk-cell assignments with similarity values above the thresold would remain. The default is 0.
#' @param df.desired.bulk A "data.frame" of desired bulk for some cells. The cells could be specified by providing x-y axis ranges in an embedding plot ("UMAP", "PCA", "TSNE") returned by \code{plot_dim}. E.g. \code{df.desired.bulk <- data.frame(x.min=c(4, -6), x.max=c(5, -5), y.min=c(-2.5, 2), y.max=c(-2, 2.5), desiredSVGBulk=c('CORT', 'STELE'), dimred='UMAP')}, where columns \code{x.min}, \code{x.max}, \code{y.min}, \code{y.max}, \code{desiredSVGBulk}, \code{dimred} are required. In this example, cells located in 4 <= x <= 5 and -2.5 <= y <= -2 in the "UMAP" plot are assigned "STELE", and cells located in -6 <= x <= -5 and 2 <= y <= 2.5 in the "UMAP" plot are assigned "CORT". \cr Alternatively, the "data.frame" could be downloaded from the Shiny app launched by \code{desired_bulk_shiny}. \cr \code{df.desired.bulk} and \code{df.match} together are used to tailor the co-clustering results. That is to say additional true bulk-cell assignments are created and included in the final assignments. If these assignments conflict with the co-colustering results the latter would be overwritten.
#' @param df.match The ground-truth matching between cells and bulk. See the example of \code{data(df.match)}.
#' @param true.only Logical. If \code{TRUE}, only the true assignments are returned in the subset \code{SingleCellExperiment}. This argument affects the values for plottting SHMs.

#' @return A \code{SingleCellExperiment} of remaining bulk-cell assignments. 

#' @examples

#' See function "cocluster" by running "?cocluster".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.

#' @export sub_asg
#' @importFrom SingleCellExperiment colLabels  
#' @importFrom SummarizedExperiment colData 

sub_asg <- function(res.lis, thr=0, df.desired.bulk=NULL, df.match=NULL, true.only=TRUE) {
  # save(res.lis, thr, df.desired.bulk, df.match, true.only, file='sub.asg.arg')
  # Validate labels.
  sce <- res.lis$cell.refined; labs <- colLabels(sce)
  er.wa <- check(as.vector(labs), as.numeric)
  if (is.numeric(er.wa) & length(er.wa) > 0) colLabels(sce) <- as.character(labs)
  # If desired bulk is provided with x-y coordinate range, convert "df.desired.bulk" to the form downloaded from Shiny app.
  xy.ran.cna <- c('x.min', 'x.max', 'y.min', 'y.max') 
  if (any(xy.ran.cna %in% colnames(df.desired.bulk))) {
    if (!all(c(xy.ran.cna, 'desiredSVGBulk', 'dimred') %in% colnames(df.desired.bulk))) stop(paste0('Make sure the "df.desired.bulk" has these columns: ', paste0(xy.ran.cna, collapse=', '), '!'))
    x.min <- df.desired.bulk$x.min
    x.max <- df.desired.bulk$x.max
    y.min <- df.desired.bulk$y.min
    y.max <- df.desired.bulk$y.max
    dimred <- unique(df.desired.bulk$dimred)
    if (!(all(x.min <= x.max) & all(y.min <= y.max))) stop('Make sure x.min <= x.max and y.min <= y.max !')
    if (! dimred %in% reducedDimNames(sce)) stop(paste0(dimred, ' is not found!'))
    gg <- plot_dim(sce, dim=dimred, color.by=colnames(colData(sce))[1])
    dat <- gg$data; df.desire <- NULL
    blk <- df.desired.bulk$desiredSVGBulk
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
      df.desire0 <- data.frame(x=df0$x, y=df0$y, key=df0$key, desiredSVGBulk=blk[i], dimred=dimred) 
      df.desire <- rbind(df.desire, df.desire0)   
    }; df.desired.bulk <- df.desire
  }

  df.roc <- res.lis$df.roc
  # Include desired bulk to "df.roc".
  if (!is.null(df.desired.bulk)) {
    xy.cna <- c('x', 'y', 'key', 'desiredSVGBulk', 'dimred') 
    if (!all(xy.cna %in% colnames(df.desired.bulk))) stop(paste0('Make sure "df.desired.bulk" has these columns: ', paste0(xy.cna, collapse=', '), '!'))
    key <- as.numeric(df.desired.bulk$key)
    # "none" indicates no change.
    desire.blk <- setdiff(unique(df.desired.bulk$desiredSVGBulk), 'none')
    if (length(desire.blk)>0) {
      for (i in desire.blk) {
      df.desired.bulk0 <- subset(df.desired.bulk, desiredSVGBulk==i)
      df.match0 <- subset(df.match, SVGBulk==i)[1, ]
      # Use "df.match" and "desire.blk" as gold standard to construct true coclustering results "df.roc0" for each "desire.blk", then append "df.roc0" to "df.roc".
      df.roc0 <- data.frame(matrix(ncol=ncol(df.roc), nrow=nrow(df.desired.bulk0), dimnames=list(NULL, colnames(df.roc))))
      df.roc0[, colnames(df.match0)] <- df.match0 
      df.roc0$assignedBulk <- df.roc0$trueBulk
      df.roc0$response <- TRUE; df.roc0$predictor <- 1
      df.roc0$index <- df.desired.bulk0$key
      df.roc <- rbind(subset(df.roc, !index %in% df.roc0$index), df.roc0)
      }; res.lis$df.roc <- df.roc
    }
  }
 
  # Subset bulk-cell assignment according to thr.
  # False assignments are aslo included and will be aggregated for plotting SHMs.
  df.roc <- subset(df.roc, predictor >= thr)
  # Cells without bulk assigments have original labels.
  cdat <- colData(sce)
  cdat$assignedBulk <- cdat$response <- 'none'
  cdat$index <- seq_len(ncol(sce))
  cdat$assignedBulk[df.roc$index] <- df.roc$assignedBulk
  cdat$response[df.roc$index] <- df.roc$response
  cdat$index[df.roc$index] <- df.roc$index
  # In consideration of rematching in coclustering.
  cdat$cell[df.roc$index] <- df.roc$cell
  cdat$SVGBulk[df.roc$index] <- df.roc$SVGBulk
  cdat$trueBulk[df.roc$index] <- df.roc$trueBulk

  cdat <- cdat[, c('label', 'cell', 'assignedBulk', 'response', 'SVGBulk', 'trueBulk', 'index')]
  colData(sce) <- cdat
  sce.sub <- sce[, df.roc$index]
  if (true.only==TRUE) sce.sub <- sce.sub[, colData(sce.sub)$response==TRUE]
  # Indicate the dim plot through which desired bulk is assigned, so that in SHMs the same dim plot is generated.
  if (!is.null(df.desired.bulk)) colData(sce.sub)$dimred <- unique(df.desired.bulk$dimred)
  return(list(cell.refined=sce, cell.sub=sce.sub, df.roc=df.roc))
}
