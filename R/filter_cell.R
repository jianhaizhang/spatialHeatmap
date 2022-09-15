#' Filter single cell data separately in a list
#'
#' Filter single cell data separately in a list and take overlap genes between all single cell and bulk data. The bulk data are not filtered and are only used to obtain overlap genes.

#' @param sce A \code{SingleCellExperiment} of single cell data.  
#' @param bulk The bulk data in form of \code{data.frame}, \code{SummarizedExperiment}, or \code{SingleCellExperiment}. They are only used to obtain overlapping genes with single cell data and not filtered. The default is \code{NULL}.
#' @param gen.rm A regular expression of gene identifiers in single cell data to remove before filtering. E.g. mitochondrial, chloroplast and protoplasting-induced genes (\code{^ATCG|^ATCG}). The default is \code{NULL}.
#' @param cutoff The minmun count of gene expression. The default is \code{1}.
#' @param p.in.cell The proportion cutoff of counts above \code{cutoff} in a cell. The default is \code{0.4}.
#' @param p.in.gen The proportion cutoff of counts above \code{cutoff} in a gene. The default is \code{0.2}.
#' @param verbose Logical. If \code{TRUE} (default), intermediate messages are printed.
#' @inheritParams norm_cell

#' @return A list of filtered single cell data and bulk data, which have common genes.

#' @examples

#' # Example bulk data of mouse brain for coclustering (Vacher et al 2021). 
#' blk.mus.pa <- system.file("extdata/shinyApp/example", "bulk_mouse_cocluster.txt", package="spatialHeatmap")
#' blk.mus <- as.matrix(read.table(blk.mus.pa, header=TRUE, row.names=1, sep='\t', check.names=FALSE))  
#' blk.mus[1:3, 1:5] 
#'# Example single cell data for coclustering (Ortiz et al 2020).
#' sc.mus.pa <- system.file("extdata/shinyApp/example", "cell_mouse_cocluster.txt", package="spatialHeatmap")  
#' sc.mus <- as.matrix(read.table(sc.mus.pa, header=TRUE, row.names=1, sep='\t', check.names=FALSE)) 
#'sc.mus[1:3, 1:5]  
 
#' # Initial filtering. 
#' blk.mus <- filter_data(data=blk.mus, sam.factor=NULL, con.factor=NULL, pOA=c(0.1, 5), CV=c(0.2, 100), dir=NULL)
#' dim(blk.mus) 
#' mus.lis <- filter_cell(lis=list(sc.mus=sc.mus), bulk=blk.mus, gen.rm=NULL, cutoff=1, p.in.cell=0.5, p.in.gen=0.1)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0. https://bioconductor.org/packages/SummarizedExperiment
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). "Orchestrating single-cell analysis with Bioconductor." _Nature Methods_, *17*, 137-145. <URL: https://www.nature.com/articles/s41592-019-0654-x>
#' Douglas Bates and Martin Maechler (2021). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.4-0. https://CRAN.R-project.org/package=Matrix
#' Vacher CM, Lacaille H, O'Reilly JJ, Salzbank J et al. Placental endocrine function shapes cerebellar development and social     behavior. Nat Neurosci 2021 Oct;24(10):1392-1401. PMID: 34400844.                                                                  
#' Ortiz C, Navarro JF, Jurek A, Märtin A et al. Molecular atlas of the adult mouse brain. Sci Adv 2020 Jun;6(26):eabb3446. PMID:  32637622
# Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x

#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment

filter_cell <- function(sce, bulk=NULL, gen.rm=NULL, cutoff=1, p.in.cell=0.4, p.in.gen=0.2, com=FALSE, verbose=TRUE) {
  # if (any(names(lis)=='')) stop('The list should be named!')
  # gen.na <- NULL                                
  # if (!is.null(bulk)) gen.na <- list(rownames(bulk))
  # for (i in seq_along(lis)) {
  #  if (verbose==TRUE) cat(names(lis[i]), '\n')
  res <- preprocess_sc(data=sce, gen.rm=gen.rm, cutoff=cutoff, p.in.cell=p.in.cell, p.in.gen=p.in.gen, verbose=verbose)
    # gen.na <- c(gen.na, list(rownames(lis[[i]]))) 
  # } 
  # gen.ovl <- Reduce(intersect, gen.na)
  if (!is.null(bulk)) {
    if (!is(bulk, 'SummarizedExperiment') & !is(bulk, 'SingleCellExperiment')) bulk <- SingleCellExperiment(assays=list(counts=as.matrix(bulk)))
    int <- intersect(rownames(bulk), rownames(res))
    bulk <- bulk[int, ]; res <- res[int, ]
    if (com==TRUE) {  
      bulk$bulkCell <- 'bulk'; res$bulkCell <- 'cell'
      res <- cbind(bulk, res)
    } else if (com==FALSE) { res <- list(bulk=bulk, cell=res) }
    # for (i in seq_along(lis)) {
    #  lis0 <- lis[[i]]; int <- intersect(rownames(lis0), rownames(bulk))
    #  lis[[i]] <- cbind(bulk[int, ], lis0[int, ])
    # }
  } 
  #names(lis) <- paste0('bulk.', names(lis)); return(lis) 
  return(res)
} 


#' Filter genes/cells by proportion of min count in row/column.
#'
#' @param data The single cell data in form of \code{data.frame}, \code{SingleCellExperiment}, or \code{SummarizedExperiment}.
#' @param gen.rm The pattern of gene identifiers to remove such as mitochondrial, chloroplast and protoplasting-induced genes (\code{^ATCG|^ATCG}).
#' @param cutoff The min count of gene expression. The default is \code{1}.
#' @param p.in.cell The proportion of counts above \code{cutoff} in a cell.
#' @param p.in.gene The proportion of counts above \code{cutoff} in a gene.
#' @param verbose Logical. If \code{TRUE} (default), intermediate messages are printed.

#' @return The filtered single cell data.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0. https://bioconductor.org/packages/SummarizedExperiment
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). "Orchestrating single-cell analysis with Bioconductor." _Nature Methods_, *17*, 137-145. <URL: https://www.nature.com/articles/s41592-019-0654-x>
#' Douglas Bates and Martin Maechler (2021). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.4-0. https://CRAN.R-project.org/package=Matrix
#' R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
# Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x

#' @importFrom SummarizedExperiment assay 
#' @importFrom Matrix rowSums colSums
#' @importFrom utils head 
#' @importFrom SingleCellExperiment SingleCellExperiment

preprocess_sc <- function(data, gen.rm=NULL, cutoff=1, p.in.cell=0.4, p.in.gen=0.2, verbose=TRUE) { 
  if (is(data, 'SummarizedExperiment') | is(data, 'SingeCellExperiment')) {
    dat <- assay(data)
    if (!is(dat, 'dgCMatrix')) dat <- as(dat, 'dgCMatrix')
  } else dat <- data 
  # Remove mitochondrial, chloroplast and protoplasting-induced genes. E.g. gen.rm='^ATCG|^ATCG'
  if (!is.null(gen.rm)) {
    rm <- grepl(gen.rm, rownames(dat))
    dat <- dat[!rm, ]; dim(dat); data <- data[!rm, ] 
  }
  if (cutoff==0|(p.in.cell==0 & p.in.gen==0)|(cutoff==0 & p.in.cell==0 & p.in.gen==0)) return(data)
  # Filter genes according to min count by genes and by cells.                                       
  if (verbose==TRUE) {
    cat('Before filtering :\n'); print(dim(dat))
    print(head(sort(rowSums(dat)), 5)) 
    print(head(sort(colSums(dat)), 5)) 
  }
  idx <- dat >= cutoff 
  idx.r <- rowSums(idx) >= p.in.gen * ncol(dat) 
  idx <- idx[idx.r, , drop=FALSE] 
  idx.c <- colSums(idx) >= p.in.cell*nrow(idx) 
  dat <- dat[idx.r, idx.c, drop=FALSE] 
  data <- data[idx.r, idx.c] 
  if (verbose==TRUE) {
  cat('After filtering :\n'); print(dim(dat)) 
    print(head(sort(rowSums(dat)), 5))  
    print(head(sort(colSums(dat)), 5))
  }
  if (!is(data, 'SummarizedExperiment') & !is(data, 'SingeCellExperiment')) data <- SingleCellExperiment(assays=list(counts=as.matrix(data)))
  return(data)
} 



