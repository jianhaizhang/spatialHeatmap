#' Normalize one or multiple count data sets.
#'
#' Normalize count data of single cell and bulk provided in a \code{list} in co-clustering. The single cell and bulk data are combined, normalized and subsequently separated. The input single cell and bulk data are replaced by normalized data respectively.
#' @param dat.lis A named \code{list} containing count data of single cell and bulk, which are in form of \code{matrix}, \code{data.frame}, \code{dgCMatrix}, or \code{SingleCellExperiment}.
#' @param cpm Logical. The count data are first normalized by \code{\link[scran]{computeSumFactors}}. If \code{TRUE}, the data is further normalized by counts per million (cpm). The default is \code{FALSE}.
#' @param count.kp Logical. If \code{FALSE} (default), the count data is discarded and only log2-scale data are kept.

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

norm_multi <- function(dat.lis, cpm=FALSE, count.kp=FALSE) {
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
  min.size <- 100
  if (min.size > ncol(sce.sc)) { print('Fewer cells than min size in quickCluster!'); return() }
  clusters <- quickCluster(sce.sc, min.size=min.size)
  sce.sc <- computeSumFactors(sce.sc, cluster=clusters, min.mean=1)
  sce.sc.nor <- logNormCounts(sce.sc)
  # CPM.
  if (cpm==TRUE) sce.sc.nor <- cal_cpm(sce.sc.nor)
  if (count.kp==FALSE) assays(sce.sc.nor)$counts <- NULL
  for (i in nas) {
    sce0 <- subset(sce.sc.nor, , set==i)
    colData(sce0)$set <- NULL; dat.lis[[i]] <- sce0
  }; return(dat.lis)
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



