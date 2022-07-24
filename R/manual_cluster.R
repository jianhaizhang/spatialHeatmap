#' Add manual cell cluster labels to SingleCellExperiment 
#'
#' Add manually created cell cluster labels in a \code{data.frame} to \code{SingleCellExperiment}.

#' @param sce A \code{SingleCellExperiment}. 
#' @param df.clus A \code{data.frame} of manually created single cell clusters. The \code{cell} column contains cell identifiers that are present in \code{rownames(colData(sce))} and the \code{cluster} column contains the manually created cluster labels. 

#' @return An object of \code{SingleCellExperiment}.

#' @examples

#' set.seed(10)
#' # Read single cell data.
#' sce.pa <- system.file("extdata/shinyApp/example", "sce_manual_mouse.rds", package="spatialHeatmap")
#' sce <- readRDS(sce.pa)
#' # Quality control, normalization, dimensionality reduction on the single cell data.
#' sce.dimred <- process_cell_meta(sce.manual, qc.metric=list(subsets=list(Mt=rowData(sce.manual)$featureType=='mito'), threshold=1))
#' # Read manual cell cluster labels.
#' manual.clus.mus.sc.pa <- system.file("extdata/shinyApp/example", "manual_cluster_mouse_brain.txt", package="spatialHeatmap") 
#' manual.clus.mus.sc <- read.table(manual.clus.mus.sc.pa, header=TRUE, sep='\t')
#' # Include manual cell cluster labels in "SingleCellExperiment".
#' sce.clus <- manual_cluster(sce=sce.dimred, df.clus=df.clus)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Morgan M, Obenchain V, Hester J, Pagès H (2022). SummarizedExperiment: SummarizedExperiment container. R package version 1.26.1, https://bioconductor.org/packages/SummarizedExperiment

#' @importFrom SummarizedExperiment colData colData<-  

#' @export manual_cluster

manual_cluster <- function(sce, df.clus) { 
  if (!all(c('cell', 'cluster') %in% colnames(df.clus))) stop('Both "cell" and "cluster" columns are required in "df.clus"!') 
  if (any(duplicated(df.clus$cell))) stop('Cell identifiers in "df.clus" should be unique!')
  cdat <- colData(sce); rna <- rownames(cdat) 
  if ('cluster' %in% colnames(cdat)) message('The existing column "cluster" in "colDdata(sce)" will be overwritten!')  
  int <- intersect(df.clus$cell, rna) 
  if (length(int)==0) stop('No cells in "df.clus" are not detected in "sce"!') 
  miss <- setdiff(df.clus$cell, rna)
  if (length(miss) > 0) { 
    msg <- paste0(length(miss), ' cells in "df.clus" are not detected in "sce"!'); message(msg) 
  } 
  sce <- sce[, rownames(cdat) %in% int]; cdat <- colData(sce)  
  df.clus <- subset(df.clus, cell %in% int)
  clus <- df.clus$cluster; names(clus) <- df.clus$cell 
  cdat$cluster <- clus[rownames(cdat)]
  cdat <- cdat[, c('cluster', setdiff(colnames(cdat), 'cluster'))] 
  colData(sce) <- cdat; return(sce) 
}     

