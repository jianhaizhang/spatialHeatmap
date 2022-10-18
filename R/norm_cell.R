#' Normalizing single cell data
#'
#' A meta function for normalizing single-cell RNA-seq data.
#' @param sce Single cell count data in form of \code{SingleCellExperiment} after quality control, which is returned by \code{qc_cell}.
#' @param bulk Bulk tissue count data in form of \code{SingleCellExperiment}, \code{SummarizedExperiment}, or \code{data.frame}.
#' @param cpm Logical. If \code{FALSE} (default), the count data are only normalized by \code{\link[scran]{computeSumFactors}}. If \code{TRUE}, the data are first normalized by \code{\link[scran]{computeSumFactors}} then transformed to counts per million by \code{\link[scuttle]{calculateCPM}}. 
#' @param count.kp Logical. If \code{FALSE} (default), the count data is discarded and only log2-scale data are kept.
#' @param quick.clus Arguments in a named list passed to \code{\link[scran]{quickCluster}}, such as \code{quick.clus=list(min.size = 100)}. 
#' @param com.sum.fct Arguments in a named list passed to \code{\link[scran]{computeSumFactors}}, such as \code{com.sum.fct=list(max.cluster.size = 3000))}. 
#' @param log.norm Arguments in a named list passed to \code{\link[scuttle]{logNormCounts}}. 
#' @param com Logical, if \code{TRUE} the returned cell and bulk data are column-wise combined, otherwise they are separated in a \code{list}. 
#' @param wk.dir The directory path to save normalized data. 
#' @return A \code{SingleCellExperiment} object. 

#' @examples

#' library(scran); library(scuttle); library(SummarizedExperiment) 
#' sce <- mockSCE()
#' sce.qc <- qc_cell(sce, qc.metric=list(subsets=list(Mt=rowData(sce)$featureType=='mito'), threshold=1))
#' sce.norm <- norm_cell(sce.qc)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.
#' Morgan M, Obenchain V, Hester J, Pagès H (2022). SummarizedExperiment: SummarizedExperiment container. R package version 1.26.1, https://bioconductor.org/packages/SummarizedExperiment


#' @export
#' @importFrom SingleCellExperiment altExpNames 
#' @importFrom SummarizedExperiment colData<-  
#' @importFrom scuttle logNormCounts
#' @importFrom scran quickCluster computeSumFactors 

norm_cell <- function(sce, bulk=NULL, cpm=FALSE, count.kp=FALSE, quick.clus=list(min.size = 100), com.sum.fct=list(max.cluster.size = 3000, min.mean=1), log.norm=list(), com=FALSE, wk.dir=NULL) {
  bulkCell <- NULL
  if (!is.null(wk.dir)) norm.dir <- file.path(wk.dir, 'norm_res') else norm.dir <- NULL
  if (!is.null(norm.dir)) if (!dir.exists(norm.dir)) dir.create(norm.dir, recursive = TRUE)
  # if (!is(sce, 'list')) sce <- list(sce=sce); nas <- names(sce)
  # if (any(nas=='')) stop('The input data list should be named!')
  #for (i in nas) { 
    # sce0 <- sce[[i]]
  if (!is(sce, 'SingleCellExperiment')) sce <- SingleCellExperiment(assays=list(counts=as.matrix(sce)))
  if (!is.null(bulk)) {
    if (!is(bulk, 'SummarizedExperiment') & !is(bulk, 'SingleCellExperiment')) bulk <- SingleCellExperiment(assays=list(counts=as.matrix(bulk)))
    if (is(bulk, 'SummarizedExperiment')) bulk <- as(bulk, 'SingleCellExperiment') 
    colData(sce) <- colData(bulk) <- NULL
    bulk$bulkCell <- 'bulk'; sce$bulkCell <- 'cell'
    bulk$sample <- colnames(bulk); sce$sample <- colnames(sce)
    int <- intersect(rownames(bulk), rownames(sce)) 
    sce <- cbind(bulk[int, ], sce[int, ])
  }
    if (quick.clus$min.size > ncol(sce)) { message('fewer cells than min size in quickCluster!'); return() }
    # Normalization.
    clusters <- do.call(quickCluster, c(list(x=sce), quick.clus))
    sce <- do.call(computeSumFactors, c(list(x=sce, cluster=clusters), com.sum.fct))
    sce <- do.call(logNormCounts, c(list(x=sce), log.norm))
    # CPM.
    if (cpm==TRUE) sce <- cal_cpm(sce)
    if (count.kp==FALSE) assays(sce)$counts <- NULL
    # sce[[i]] <- sce0 
  # } 
  res <- sce 
  if (!is.null(bulk) & com==FALSE) {
    bulk <- subset(res, , bulkCell=='bulk')
    cell <- subset(res, , bulkCell=='cell')
    res <- list(bulk=bulk, cell=cell)
  }
  if (!is.null(norm.dir)) saveRDS(res, file=paste0(norm.dir, '/', ifelse(cpm==TRUE, 'cpm', 'fct'), '.rds'))
  return(res)
}       

#' Normalize by CPM. 
#'
#' The log2 values are transformed to power of 2 (counts) and then to CPM. To maintain log2-scale values, the CPM values are transformed back to log2. The returned \code{SingleCellExperiment} contains values at both CPM and log2.
#' @param sce.nor The output of \code{\link[scuttle]{logNormCounts}} in form of \code{SingleCellExperiment}, where library size factors are applied already.

#' @return A \code{SingleCellExperiment} object.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.
#' Douglas Bates and Martin Maechler (2021). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.4-0. https://CRAN.R-project.org/package=Matrix

#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom scuttle calculateCPM
#' @importFrom Matrix Matrix 

cal_cpm <- function(sce.nor) {
  log.cnt <- logcounts(sce.nor)
  cnt.cpm <- calculateCPM(2^log.cnt-1)
  logcounts(sce.nor) <- Matrix(log2(cnt.cpm+1), sparse=TRUE)
  return(sce.nor)
}
