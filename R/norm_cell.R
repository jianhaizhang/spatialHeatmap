#' Normalizing single cell data
#'
#' A meta function for normalizing single-cell RNA-seq data.
#' @param sce Single cell count data in form of \code{SingleCellExperiment} after quality control, which is returned by \code{qc_cell}.
#' @param bulk Bulk tissue count data in form of \code{SingleCellExperiment}, \code{SummarizedExperiment}, or \code{data.frame}.
#' @param cpm Logical. If \code{FALSE} (default), the count data are only normalized by \code{\link[scrapper]{normalizeRnaCounts.se}}. If \code{TRUE}, the data are further transformed to counts per million. 
#' @param count.kp Logical. If \code{FALSE} (default), the count data is discarded and only log2-scale data are kept.
#' @param center.sf Arguments in a named list passed to \code{\link[scrapper]{centerSizeFactors}}. 
#' @param log.norm Arguments in a named list passed to \code{\link[scrapper]{normalizeRnaCounts.se}}. 
#' @param com Logical, if \code{TRUE} the returned cell and bulk data are column-wise combined, otherwise they are separated in a \code{list}. 
#' @param wk.dir The directory path to save normalized data. 
#' @return A \code{SingleCellExperiment} object. 

#' @examples

#' library(SingleCellExperiment)
#' set.seed(123)  # for reproducibility

#' # random counts: 1000 genes × 100 cells
#' counts_mat <- matrix(
#'   rpois(1000 * 100, lambda = 5), nrow = 1000, ncol = 100,
#'   dimnames = list(
#'     paste0("Gene", 1:1000),
#'     paste0("Cell", 1:100)
#'   )
#' )

#' sce <- SingleCellExperiment(
#'   assays = list(counts = counts_mat)
#' )

#' sce.qc = qc_cell(sce)
#' sce.norm <- norm_cell(sce.qc)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Morgan M, Obenchain V, Hester J, Pagès H (2022). SummarizedExperiment: SummarizedExperiment container. R package version 1.26.1, https://bioconductor.org/packages/SummarizedExperiment
#' Lun ATL, Kancherla J (2023). "Powering single-cell analyses in the browser with WebAssembly." _Journal of Open Source Software_, *8*(89), 5603. doi:10.21105/joss.05603 <https://doi.org/10.21105/joss.05603>.

#' @export
#' @importFrom SingleCellExperiment altExpNames 
#' @importFrom SummarizedExperiment colData<- colData

norm_cell <- function(sce, bulk=NULL, cpm=FALSE, count.kp=FALSE,
                      center.sf=list(block = NULL, mode = c("lowest", "per-block")),
                      log.norm=list(assay.type = "counts", log = TRUE, pseudo.count = 1), com=FALSE, wk.dir=NULL) {
  # save(sce, bulk, cpm, count.kp, center.sf, log.norm, com, wk.dir, file='norm_cell.arg')
  pkg <- check_pkg('scrapper'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  bulkCell <- NULL
  if (!is.null(wk.dir)) norm.dir <- file.path(wk.dir, 'norm_res') else norm.dir <- NULL
  if (!is.null(norm.dir)) if (!dir.exists(norm.dir)) dir.create(norm.dir, recursive = TRUE)
  # if (!is(sce, 'list')) sce <- list(sce=sce); nas <- names(sce)
  # if (any(nas=='')) stop('The input data list should be named!')
  #for (i in nas) { 
    # sce0 <- sce[[i]]
  if (!is(sce, 'SingleCellExperiment')) sce <- SingleCellExperiment(assays=list(counts=as.matrix(sce)))
  cdat.sc <- colData(sce)
  if (!is.null(bulk)) {
    if (!is(bulk, 'SummarizedExperiment') & !is(bulk, 'SingleCellExperiment')) bulk <- SingleCellExperiment(assays=list(counts=as.matrix(bulk)))
    if (is(bulk, 'SummarizedExperiment')) bulk <- as(bulk, 'SingleCellExperiment')
    cdat.blk <- colData(bulk) 
  
    # Erase input metadata in colData if input bulk and cell data do not have the same colnames in colData.
    #if (ncol(cdat.sc)==0 | ncol(cdat.blk)== 0 | ncol(cdat.sc)!=ncol(cdat.blk)) { 
    #  colData(sce) <- colData(bulk) <- NULL
    #}
    #if (ncol(cdat.sc) > 0 & ncol(cdat.blk) > 0 & ncol(cdat.sc)==ncol(cdat.blk)) {
    #  if (!all(colnames(cdat.sc)==colnames(cdat.blk))) colData(sce) <- colData(bulk) <- NULL
    #}
    bulk$bulkCell <- 'bulk'; sce$bulkCell <- 'cell'
    bulk$sample <- colnames(bulk); sce$sample <- colnames(sce)
    # int <- intersect(rownames(bulk), rownames(sce))
    # pkg <- check_pkg('BiocGenerics'); if (is(pkg, 'character')) stop(pkg)
    sce <- cbind_se(bulk, sce)
  }
    # Normalization
    sf <- do.call(scrapper::centerSizeFactors, c(list(size.factors=colSums(assay(sce))), center.sf))
    sce <- do.call(scrapper::normalizeRnaCounts.se, c(list(x=sce, size.factors=sf), log.norm))
    
    # CPM.
    if (cpm==TRUE) sce <- cal_cpm(sce.nor=sce, sf)
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
