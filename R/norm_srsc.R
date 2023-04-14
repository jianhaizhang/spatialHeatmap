
#' Jointl normalization of spatially resolved single cell data and bulk data
#'
#' @param cell A \code{Seurat} object.
#' @param assay The assay to use for normalization in the spatial single-cell data.
#' @param bulk The bulk assay data. 

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @examples

#' library(Seurat); library(SeuratData); library(SummarizedExperiment)
#' if (!'stxBrain' %in% InstalledData()$Dataset) InstallData("stxBrain")
#' # Bulk data. 
#' blk.mus <- readRDS(system.file("extdata/shinyApp/data", "bulk_mouse_cocluster.rds", package="spatialHeatmap"))
#' assay(blk.mus)[1:3, 1:5] 
#' # Spatial single-cell data.
#' brain <- LoadData("stxBrain", type = "anterior1") 
#' # Joint normalization. 
#' # nor.lis <- norm_srsc(cell=brain, assay='Spatial', bulk=blk.mus) 
#' # Separate bulk and cell data. 
#' # srt.sc <- nor.lis$cell; bulk <- nor.lis$bulk 

#' @references
#' Hao and Hao et al. Integrated analysis of multimodal single-cell data. Cell (2021) [Seurat V4]
#' Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019) [Seurat V3]
#' Butler et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat Biotechnol (2018) [Seurat V2]
#' Satija and Farrell et al. Spatial reconstruction of single-cell gene expression data. Nat Biotechnol (2015) [Seurat V1]
#' Pagès H, Lawrence M, Aboyoun P (2022). _S4Vectors: Foundation of vector-like and list-like containers in Bioconductor_. R package version 0.36.1, <https://bioconductor.org/packages/S4Vectors>.
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.  0, https://bioconductor.org/packages/SummarizedExperiment.
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” _Nature Methods_, *17*, 137-145. <https://www.nature.com/articles/s41592-019-0654-x>.

#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays<- 
#' @importFrom S4Vectors DataFrame 

norm_srsc <- function(cell, assay, bulk) {
  pkg <- check_pkg('Seurat'); if (is(pkg, 'character')) stop(pkg)
  if (!is(cell, 'Seurat')) stop('"cell" should be a "Seurat" object!')
  if (!is(bulk, 'SingleCellExperiment') & !is(bulk, 'SummarizedExperiment')) stop('"bulk" should be a "SingleCellExperiment" or "SummarizedExperiment" object!')

  sce.sp <- Seurat::as.SingleCellExperiment(cell, assay=assay)
  sce.sp$sample <- make.names(colnames(sce.sp))
  bulk$sample <- make.names(colnames(bulk))
  # Take overlap genes.
  inter <- intersect(rownames(sce.sp), rownames(bulk))
  if (length(inter)==0) stop('"cell" and "bulk" do not have common biomolecules!')
  sce.sp <- sce.sp[inter, ]; bulk <- bulk[inter, ]
  idx.blk <- seq_len(ncol(bulk))
  idx.sc <- setdiff(seq_len(ncol(sce.sp)+ncol(bulk)), idx.blk)
  colnames(bulk) <- idx.blk; colnames(sce.sp) <- idx.sc
  # Joint normalization.
  sce <- SingleCellExperiment(
    list(counts=cbind(assay(bulk), assay(sce.sp))), 
    colData=DataFrame(sample=c(bulk$sample, sce.sp$sample))
  )
  srt <- Seurat::as.Seurat(sce, counts = "counts", data=NULL)
  srt <- Seurat::SCTransform(srt, assay = "originalexp", verbose = FALSE)
  # Separate bulk.
  bulk <- subset(srt, cells=idx.blk)
  bulk <- Seurat::as.SingleCellExperiment(bulk, assay='SCT')
  assays(bulk)$scaledata <- NULL; colnames(bulk) <- bulk$sample
  # Separate cells.
  srt.sc <- subset(srt, cells=idx.sc)
  srt.sc <- Seurat::RenameCells(srt.sc, new.names=unname(srt.sc$sample))
  cell <- Seurat::RenameCells(cell, new.names=unname(srt.sc$sample))
  cell <- subset(cell, features=rownames(srt.sc))
  cell@assays[[assay]] <- srt.sc@assays$originalexp
  cell@assays$SCT <- srt.sc@assays$SCT
  Seurat::DefaultAssay(cell) <- 'SCT'; cell@meta.data <- srt.sc@meta.data
  return(list(cell=cell, bulk=bulk))
}
