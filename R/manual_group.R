#' Add manual cell group labels to SingleCellExperiment 
#'
#' Add manually created cell group labels in a \code{data.frame} to \code{SingleCellExperiment}.

#' @param sce A \code{SingleCellExperiment}. 
#' @param df.group A \code{data.frame} of manually created single cell groups. At least two columns are required, corresponding to cell identifiers that are present in \code{rownames(colData(sce))} and the manually created group labels respectively.
#' @param cell The column name in \code{df.group} indicating cells.
#' @param cell.group The column name in \code{df.group} indicating cell group labels.

#' @return An object of \code{SingleCellExperiment}.

#' @examples

#' set.seed(10)
#' # Read single cell data.
#' sce.pa <- system.file("extdata/shinyApp/example", "sce_manual_mouse.rds", package="spatialHeatmap")
#' sce <- readRDS(sce.pa)
#' # Quality control, normalization, dimensionality reduction on the single cell data.
#' sce.dimred <- process_cell_meta(sce.manual, qc.metric=list(subsets=list(Mt=rowData(sce.manual)$featureType=='mito'), threshold=1))
#' # Read manual cell group labels.
#' manual.clus.mus.sc.pa <- system.file("extdata/shinyApp/example", "manual_cluster_mouse_brain.txt", package="spatialHeatmap") 
#' manual.clus.mus.sc <- read.table(manual.clus.mus.sc.pa, header=TRUE, sep='\t')
#' # Include manual cell group labels in "SingleCellExperiment".
#' sce.clus <- manual_group(sce=sce.dimred, df.group=manual.clus.mus.sc, cell='cell', cell.group='cell.group')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Morgan M, Obenchain V, Hester J, Pagès H (2022). SummarizedExperiment: SummarizedExperiment container. R package version 1.26.1, https://bioconductor.org/packages/SummarizedExperiment

#' @importFrom SummarizedExperiment colData colData<-  

#' @export manual_group

manual_group <- function(sce, df.group, cell, cell.group) {                                                                        
  # if (!all(c('cell', 'cluster') %in% colnames(df.group))) stop('Both "cell" and "cluster" columns are required in "df.group"!')  
  if (any(duplicated(df.group[, cell]))) stop('Cell identifiers in "df.group" should be unique!')                                  
  cdat <- colData(sce); rna <- rownames(cdat)                                                                                      
  if (cell.group %in% colnames(cdat)) {                                                                                            
    msg <- paste0("The existing column '", cell.group, "' in  'colDdata(sce)' will be overwritten!")                               
    message(msg); cdat[, cell.group] <- NULL; colData(sce) <- cdat                                                                 
  }                                                                                                                                
  int <- intersect(df.group[, cell], rna)                                                                                          
  if (length(int)==0) stop('No cells in "df.group" are detected in "sce"!')                                                        
  miss <- setdiff(df.group[, cell], rna)                                                                                           
  if (length(miss) > 0) {                                                                                                          
    msg <- paste0(length(miss), ' cells in "df.group" are not detected in "sce"!'); message(msg)                                   
  }                                                                                                                                
  sce <- sce[, rownames(cdat) %in% int]; cdat <- colData(sce)                                                                      
  df.group <- subset(df.group, cell %in% int)                                                                                      
  clus <- df.group[, cell.group]; names(clus) <- df.group[, cell]                                                                  
  df.clus <- DataFrame(clus=clus[rownames(cdat)])                                                                                  
  colnames(df.clus) <- cell.group; cdat <- cbind(cdat, df.clus)                                                                    
  cdat <- cdat[, c(cell.group, setdiff(colnames(cdat), cell.group))]                                                               
  colData(sce) <- cdat; return(sce)                                                                                                
}    


