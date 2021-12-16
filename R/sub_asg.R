#' Subset the bulk-cell assignments
#' 
#' Subset the bulk-cell assignments according to a threshold, which is a similarity value between bulk and cells. 

#' @param res.lis The result list of coclustering, which is the output of \code{tests} and comprises three slots \code{sce}, \code{roc.obj}, \code{df.roc}.
#' @param thr The threshold for subsetting bulk-cell assignments, which is a similarity value (Pearson's or Spearman's correlation coefficient) between bulk and cells. Only bulk-cell assignments with similarity values above the thresold would remain.
#' @param true.only Logical. If \code{TRUE}, only the true assignments are returned in the subset \code{SingleCellExperiment}. This argument affects the values for plottting SHMs.

#' @return A \code{SingleCellExperiment} of remaining bulk-cell assignments. 

#' @examples



#' @export sub_asg

sub_asg <- function(res.lis, thr, true.only=TRUE) {
  # Validate labels.
  sce <- res.lis$sce; labs <- colLabels(sce)
  er.wa <- check(as.vector(labs), as.numeric)
  if (is.numeric(er.wa) & length(er.wa) > 0) colLabels(sce) <- as.character(labs)
  # Subset bulk-cell assignment according to thr.
  # False assignments are aslo included and will be aggregated for plotting SHMs.
  df.roc <- res.lis$df.roc
  df.roc <- subset(df.roc, predictor >= thr)
  # Cells without bulk assigments have original labels.
  cdat <- colData(sce)
  cdat$assignedBulk <- cdat$response <- 'none'
  cdat$idx <- seq_len(ncol(sce))
  cdat$assignedBulk[df.roc$idx] <- df.roc$assignedBulk
  cdat$response[df.roc$idx] <- df.roc$response
  cdat$idx[df.roc$idx] <- df.roc$idx
  cdat <- cdat[, c('label', 'cell', 'assignedBulk', 'response', 'SVGBulk', 'trueBulk', 'idx')]
  colData(sce) <- cdat
  sce.sub <- sce[, df.roc$idx]
  if (true.only==TRUE) sce.sub <- sce.sub[, colData(sce.sub)$response==TRUE]
  return(list(sce=sce, sce.sub=sce.sub))
}
