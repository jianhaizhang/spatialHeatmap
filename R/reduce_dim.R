#' Reducing dimensionality in count data
#'
#' A meta function for reducing dimensionality in count data.

#' @param sce Normalized single cell data in \code{SingleCellExperiment} returned by \code{norm_cell}. Alternative forms include \code{dgCMatrix}, \code{matrix}, \code{data.frame}.
#' @param prop Numeric scalar specifying the proportion of genes to report as highly variable genes (HVGs) in \code{\link[scran]{getTopHVGs}}. The default is \code{0.1}. 
#' @param min.dim,max.dim Integer scalars specifying the minimum (\code{min.dim}) and maximum (\code{max.dim}) number of (principle components) PCs to retain respectively in \code{\link[scran]{denoisePCA}}. The default is \code{min.dim=11}, \code{max. dim=50}.  
#' @param model.var Additional arguments in a named list passed to \code{\link[scran]{modelGeneVar}}. 
#' @param top.hvg Additional arguments in a named list passed to \code{\link[scran]{getTopHVGs}}, such as \code{top.hvg=list(n =   3000)}. 
#' @param de.pca Additional arguments in a named list passed to \code{\link[scran]{denoisePCA}}, such as \code{de.pca=list(assay.  type = "logcounts")}.
#' @param pca Logical, if \code{TRUE} only the data with reduced dimentionality by PCA is returned and no clustering is performed. The default is \code{FALSE} and clustering is performed after dimensionality reduction. 
#' @param tsne Additional arguments in a named list passed to \code{\link[scater]{runTSNE}}, such as \code{tsne=list(dimred="PCA", ncomponents=2)}. 
#' @param umap Additional arguments in a named list passed to \code{\link[scater]{runUMAP}}, such as \code{umap=list(dimred="PCA")}. 

#' @return A \code{SingleCellExperiment} object. 

#' @examples

#' library(scran); library(scuttle) 
#' sce <- mockSCE()
#' sce.qc <- qc_cell(sce, qc.metric=list(subsets=list(Mt=rowData(sce)$featureType=='mito'), threshold=1))
#' sce.norm <- norm_cell(sce.qc)
#' sce.dimred <- reduce_dim(sce.norm)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods,    17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation  of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.


#' @export reduce_dim 
#' @importFrom SingleCellExperiment SingleCellExperiment 
#' @importFrom scran modelGeneVar getTopHVGs denoisePCA 
#' @importFrom scater runUMAP runTSNE

reduce_dim <- function(sce, prop=0.1, min.dim=13, max.dim=50, model.var=list(assay.type="logcounts"), top.hvg=list(), de.pca=list(assay.type = "logcounts"), pca=FALSE, tsne=list(dimred="PCA", ncomponents=2), umap=list(dimred="PCA")) {
  # getTopHVGs: prop=0.1 is super more important than n=3000 in co-clustering.
  #  save(sce, prop, min.dim, max.dim, model.var, top.hvg, de.pca, pca, tsne, umap, file='reduce.dim.arg')
  min.dim <- as.numeric(min.dim); max.dim <- as.numeric(min.dim)
  if (is(sce, 'dgCMatrix')|is(sce, 'matrix')|is(sce, 'data.frame')) {
    if (all(round(sce)==sce)) stop('The "sce" should be in log2 scale!')
    sce <- SingleCellExperiment(list(logcounts=as.matrix(sce)))
  }
  # Use logcounts by default.
  df.var.sc <- do.call(modelGeneVar, c(list(x=sce), model.var))
  top.hvgs.sc <- do.call(getTopHVGs, c(list(stats=df.var.sc, prop=prop), top.hvg))
  max.pc <- (length(top.hvgs.sc) - 1) / 2
  # Avoid warning (the max returned PCS should be <= 50% of length(top.hvgs.sc) - 1).
  if (max.pc < max.dim | max.pc < min.dim) {
    top.hvgs.sc <- getTopHVGs(df.var.sc, prop=1)
    cat('"prop" is set 1 in "getTopHVGs" due to too less genes. \n')
    max.pc <- (length(top.hvgs.sc) - 1) / 2
    if (max.pc < max.dim | max.pc < min.dim) {
      warning('The input single cell data has too less genes!'); return()
    }
  }
  # prop=1 cannot avoid all warnings if the real max.dim does not increase even though prop=1.
  sce.dimred <- tryCatch(
    expr = {
      do.call(denoisePCA, c(list(x=sce, technical=df.var.sc, subset.row=top.hvgs.sc, min.rank=min.dim, max.rank=max.dim), de.pca))
    },
    warning = function(w){ 'w' }, error = function(e){ 'e' } 
  ) 
  if (!is(sce.dimred, 'SingleCellExperiment')) {
    message('min.dim = ', min.dim, ' and ', 'max.dim = ', max.dim, ' cannot be satisfied!'); return()
  }
  if (pca==TRUE) return(sce.dimred)
  # Other argument: n_dimred, ntop. By default only 2 dimensions are returned by runTSNE/runUMAP.
  # runUMAP returns different results before/after runTSNE.
  # Avoid warnings due to duplicated column names.
  cna <- colnames(sce.dimred); colnames(sce.dimred) <- seq_len(ncol(sce.dimred)) 
  sce.dimred <- do.call(runUMAP, c(list(x=sce.dimred, ncomponents=min.dim), umap))
  # Row names in colData, reducedDim change accordingly.
  colnames(sce.dimred) <- cna
  sce.dimred <- do.call(runTSNE, c(list(x=sce.dimred), tsne))
  return(sce.dimred)
}
