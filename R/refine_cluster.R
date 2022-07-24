#' Refine single cell clusters
#'
#' In each cell cluster, the pairwise Spearman or Pearson correlation coefficients (similarities) are calculated between cells. Cells having similarities over \code{sim} with other cells in the same cluster at proportion over \code{sim.p} remain, and other cells are filtered out. The resulting clusters are more homogeneous. 

#' @param sce.clus The single cell data in form of \code{SummarizedExperiment}, where cluster assignments are stored in the \code{label} column in \code{colData} slot.  
#' @param sim,sim.p Both are numeric scalars, ranging from 0 to 1. \code{sim} is a similarity (Spearman or Pearson correlation coefficient) cutoff between cells and \code{sim.p} is a proportion cutoff. In a certain cell cluster, cells having similarity >= \code{sim} with other cells in the same cluster at proportion >= \code{sim.p} would remain. Otherwise, they are discarded.
#' @param sim.meth Method to compute similarities between cells, \code{spearman} or \code{pearson}. The \code{logcount} values in \code{sce.clus} are used. 
#' @param verbose Logical. If \code{TRUE} (default), intermediate messages are printed.

#' @return A \code{SummarizedExperiment} object with some cells discarded. 

#' @examples

#' library(scran); library(scuttle) 
#' sce <- mockSCE(); sce <- logNormCounts(sce) 
#' # Modelling the variance.
#' var.stats <- modelGeneVar(sce)
#' sce.dimred <- denoisePCA(sce, technical=var.stats, subset.row=rownames(var.stats))
#' \donttest{
#' sce.clus <- cluster_cell(data=sce.dimred, graph.meth='snn', dimred='PCA')
#' # Clusters. 
#' table(colData(sce.clus)$label)
#' 
#' cell.refined <- refine_cluster(sce.clus, sim=0.5, sim.p=0.8, sim.meth='spearman', verbose=TRUE)
#' }
#'
#' # See details in function "coclus_meta" by running "?coclus_meta".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.
#' mezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with       Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation  of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.

#' @export refine_cluster 
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment logcounts

refine_cluster <- function(sce.clus, sim=0.2, sim.p=0.8, sim.meth='spearman', verbose=TRUE) {
  colData(sce.clus)$index.all <- seq_len(ncol(sce.clus))
  lab <- as.character(colData(sce.clus)$cluster); lab.uni <- sort(unique(lab))
  dat <- dat.dim <- cdat <- NULL
  # Use expression values to calculate similarity. 
  dat.aggr <- data.frame(matrix(NA, nrow=nrow(sce.clus), ncol=1))
  idx.kp <- NULL; for (i in lab.uni) {
    sce <- sce.clus[, lab==i]
    # Calculate PCCs on log-values. 
    clus <- logcounts(sce)
    sims <- cor(clus, method=sim.meth); idx <- sims > sim
    w <- which(colSums(idx) > ncol(clus)*sim.p)
    idx.kp <- c(idx.kp, colData(sce)$index.all[w])
    # In practice, cell ids are not known.
    sce <- sce[, w]; if (verbose==TRUE) cat(paste0('Cluster', i, ': ', length(w)), '\n')
    if (verbose==TRUE) print(table(colnames(sce)))
  }; sce.kp <- sce.clus[, idx.kp]; return(sce.kp)
}
