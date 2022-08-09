#' Normalize one or multiple count data sets.
#'
#' Normalize count data of single cell and bulk provided in a \code{list} for co-clustering. The single cell and bulk data are combined, normalized and subsequently separated. The input single cell and bulk data are replaced by normalized data respectively.
#' @param dat.lis A named \code{list} containing count data of single cell and bulk, returned by \code{filter_iter}.
#' @inheritParams norm_cell

#' @return A list of normalized single cell and bulk data. 

#' @examples

#' # Example bulk data of mouse brain for coclustering (Vacher et al 2021).   
#' blk.mus.pa <- system.file("extdata/shinyApp/example", "bulk_mouse_cocluster.txt", package="spatialHeatmap")
#' blk.mus <- as.matrix(read.table(blk.mus.pa, header=TRUE, row.names=1, sep='\t', check.names=FALSE))
#' blk.mus[1:3, 1:5]
#' 
#' # Example single cell data for coclustering (Ortiz et al 2020).
#' sc.mus.pa <- system.file("extdata/shinyApp/example", "cell_mouse_cocluster.txt", package="spatialHeatmap")
#' sc.mus <- as.matrix(read.table(sc.mus.pa, header=TRUE, row.names=1, sep='\t', check.names=FALSE))
#' sc.mus[1:3, 1:5]  
#'
#' # Initial filtering.
#' blk.mus <- filter_data(data=blk.mus, sam.factor=NULL, con.factor=NULL, pOA=c(0.1, 5), CV=c(0.2, 100), dir=NULL)
#' dim(blk.mus)  
#' mus.lis <- filter_cell(lis=list(sc.mus=sc.mus), bulk=blk.mus, gen.rm=NULL, min.cnt=1, p.in.cell=0.5, p.in.gen=0.1)
   
#' \donttest{
#' # Normalization: bulk and single cell are combined and normalized, then separated.
#' mus.lis.nor <- norm_multi(dat.lis=mus.lis, cpm=FALSE)
#' }

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.
#' Douglas Bates and Martin Maechler (2021). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.4-0. https://CRAN.R-project.org/package=Matrix
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.  0, https://bioconductor.org/packages/SummarizedExperiment.
#' Vacher CM, Lacaille H, O'Reilly JJ, Salzbank J et al. Placental endocrine function shapes cerebellar development and social behavior. Nat Neurosci 2021 Oct;24(10):1392-1401. PMID: 34400844.
#' Ortiz C, Navarro JF, Jurek A, Märtin A et al. Molecular atlas of the adult mouse brain. Sci Adv 2020 Jun;6(26):eabb3446. PMID:  32637622

#' @export norm_multi
#' @importFrom SingleCellExperiment SingleCellExperiment cbind 
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom scran quickCluster computeSumFactors 
#' @importFrom scuttle logNormCounts


norm_multi <- function(dat.lis, cpm=FALSE, count.kp=FALSE, quick.clus=list(min.size = 100), com.sum.fct=list(max.cluster.size = 3000), log.norm=list(), wk.dir=NULL) {
  if (!is.null(wk.dir)) norm.dir <- file.path(wk.dir, 'norm_res') else norm.dir <- NULL 
  if (!is.null(norm.dir)) if (!dir.exists(norm.dir)) dir.create(norm.dir)
  set <- NULL; nas <- names(dat.lis)
  if (any(nas=='')) stop('The list should be named!')
  # Anchor sce.
  dat0 <- dat.lis[[1]]
  if (is(dat0, 'SingleCellExperiment') | is(dat0, 'SummarizedExperiment')) dat0 <- assays(dat0)[[1]]
  sce.sc <- SingleCellExperiment(assays=list(counts=as.matrix(dat0[, 1, drop=FALSE])), colData=DataFrame(set=0))
  for (i in seq_along(dat.lis)) {
    dat0 <- dat.lis[[i]]
    if (is(dat0, 'SummarizedExperiment')) dat0 <- as(dat0, 'SingleCellExperiment')
    if (is(dat0, 'SingleCellExperiment') | is(dat0, 'SummarizedExperiment')) sce0 <- dat0
    if (is(dat0, 'dgCMatrix')|is(dat0, 'matrix')|is(dat0, 'data.frame')) sce0 <- SingleCellExperiment(assays=list(counts=as.matrix(dat0)))
    colData(sce0)$set <- nas[i]; sce.sc <- cbind(sce.sc, sce0)
  }; sce.sc <- sce.sc[, -1]
  if (quick.clus$min.size > ncol(sce.sc)) { print('Fewer cells than min size in quickCluster!'); return() }
  sce.sc.nor <- norm_cell(sce.sc, quick.clus=quick.clus, com.sum.fct=com.sum.fct, log.norm=log.norm, cpm=cpm, count.kp=count.kp)
  for (i in nas) {
    sce0 <- subset(sce.sc.nor, , set==i)
    colData(sce0)$set <- NULL; dat.lis[[i]] <- sce0
  }
  if (!is.null(norm.dir)) saveRDS(dat.lis, file=paste0(norm.dir, '/', ifelse(cpm==TRUE, 'cpm', 'fct'), '.rds'))
  return(dat.lis)
}
