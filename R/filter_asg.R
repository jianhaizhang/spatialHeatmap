#' Refine the bulk-cell assignments
#' 
#' Refine the bulk-cell assignments by subsetting the assignments according to a threshold, which is a similarity value between bulk and cells. 

#' @param res The coclustering results returned by \code{cocluster}.
#' @param min.sim The similarity cutoff for filterig bulk-cell assignments, which is a Pearson's or Spearman's correlation coefficient between bulk and cells. Only bulk-cell assignments with similarity values above the thresold would remain. The default is 0.

#' @return A \code{SingleCellExperiment} of remaining bulk-cell assignments. 
#' @rdname cocluster


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.

#' @export
#' @importFrom SingleCellExperiment reducedDimNames 
#' @importFrom SummarizedExperiment colData colData<- 

filter_asg <- function(res, min.sim=0) {
  bulkCell <- similarity <- NULL
  if (is(res, 'SingleCellExperiment')) sce.all <- res else if (is(res, 'list')) sce.all <- res$sce.all
  # Set 'none' to assignments filtered out.
  cdat <- colData(sce.all)
  cdat0 <-  subset(cdat, bulkCell=='cell' & similarity!='none')
  cdat0$similarity <- as.numeric(cdat0$similarity)
  idx <- cdat0$similarity < min.sim
  cdat0$assignedBulk[idx] <- cdat0$similarity[idx] <- 'none'
  # Optimization.
  if (is(res, 'list')) cdat0$response[idx] <- cdat0$trueBulk[idx] <- 'none'
  cdat[cdat0$index, ] <- cdat0; colData(sce.all) <- cdat

  if (is(res, 'SingleCellExperiment')) res <- sce.all else if (is(res, 'list')) {
    if (any(c('e', 'w') %in% check_pkg('pROC'))) stop('The package "pROC" is not detected!')
    cdat <- colData(sce.all)
    cdat0 <- subset(cdat, similarity!='none')
    roc.obj <- pROC::roc(cdat0$response, as.numeric(cdat0$similarity), smoothed = TRUE, ci=TRUE, ci.alpha=0.9, stratified=FALSE, plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=TRUE, show.thres=TRUE, direction='<', levels=c('FALSE', 'TRUE'))
    res$sce.all <- sce.all; res$roc.obj <- roc.obj
  }
  return(res)
}
