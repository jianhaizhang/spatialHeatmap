#' Processing single cell RNA-seq count data 
#'
#' A meta function for processing single cell RNA-seq count data, including quality control, normalization, dimensionality reduction.

#' @param sce Single cell RNA-seq count data in \code{SingleCellExperiment}.
#' @inheritParams qc_cell
#' @inheritParams norm_cell
#' @inheritParams reduce_dim

#' @return A \code{SingleCellExperiment} object. 
#' @details
#' In the QC, frequently used per-cell metrics are calculated for identifying problematic cells, such as library size, number of detected features above a threshold, mitochodrial gene percentage, etc. Then these metrics are used to determine outlier cells     based on median-absolute-deviation (MAD). Refer to \code{perCellQCMetrics} and \code{perCellQCFilters} in the scuttle package for more details. 
#' In the normalization, a quick-clustering method is applied to divide cells into clusters. Then a scaling normalization method is performed to obtain per-cluster size factors. Next, the size factor in each cluster is decomposed into per-cell size factors by a deconvolution strategy. Finally, all cells are normalized by per-cell size factors. See more details in \code{quickCluster}, \code{computeSumFactors} from the scran package, and \code{logNormCounts} from the scuttle package.   
#' In dimensionality reduction, the high-dimensional gene expression data are embedded into a 2-3 dimensional space using PCA, tSNE and UMAP. All three embedding result sets are stored in a \code{SingleCellExperiment} object. Details are seen in \code{denoisePCA} from scran, and \code{runUMAP}, \code{runTSNE} from scater. 

#' @examples

#' library(scran); library(scuttle) 
#' sce <- mockSCE()
#' sce.dimred <- process_cell_meta(sce, qc.metric=list(subsets=list(Mt=rowData(sce)$featureType=='mito'), threshold=1))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods,    17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation  of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186. doi: 10.1093/bioinformatics/btw777.
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.

#' @export process_cell_meta 

process_cell_meta <- function(sce, qc.metric=list(threshold=1), qc.filter=list(nmads=3), quick.clus=list(min.size = 100), com.sum.fct=list(max.cluster.size = 3000), log.norm=list(), prop=0.1, min.dim=13, max.dim=50, model.var=list(), top.hvg=list(n = 3000), de.pca=list(assay.type ="logcounts"), pca=FALSE, tsne=list(dimred="PCA", ncomponents=2), umap=list(dimred="PCA")) {
  # Quality control.
  sce.qc <- qc_cell(sce=sce, qc.metric=qc.metric, qc.filter=qc.filter)
  # Normalization.
  sce.norm <- norm_cell(sce=sce.qc, quick.clus=quick.clus, com.sum.fct=com.sum.fct, log.norm=log.norm, com=TRUE)
  # Dimensionality reduction.
  sce.dimred <- reduce_dim(sce=sce.norm, prop=prop, min.dim=min.dim, max.dim=max.dim, model.var=model.var, top.hvg=top.hvg, de.pca=de.pca, pca=pca, tsne=tsne, umap=umap)
  return(sce.dimred)
}
