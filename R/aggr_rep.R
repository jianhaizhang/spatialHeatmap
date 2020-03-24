#' Aggregate Sample Replicates in Data Matrix
#'
#' This function averages sample replicates. It inputs the gene expression matrix and sample metadata in form of "SummarizedExperiment". In "colData" slot, at least "sample" and "condition" information should be included.

#' @param se A "SummarizedExperiment" storing a gene expression matrix and metadata. The "assays" slot stores a expression matrix with row and column names being gene IDs and sample/conditions, respectively. The "rowData" can store a data frame of row (gene) anntation, but is optional. \cr The "colData" slot is required and contains a data frame with at least 2 columns corresponding sample and condition information respectively. It is important that only letters, digits, single underscore, dots, single space are allowed in the 2 columns. The 2 columns are ultimately concatenated by double underscore "__" to replace the original column names in the expression matrix. Thus the syntax of the column names in the output expression matrix is "sample__condition". E.g. "shoot_pGL2__hypoxia" (Mustroph et al. 2009), where "shoot_pGL2" is the sample and "hypoxia" is the condition. If the original column names in the expression matrix are already formatted as "sample__condition", then the "colData" slot is optional.  

#' @param samples A character. The column name corresponding to samples in the "colData" of "se" argument. Can be NULL if column names of expression matrix in "se" argument are already formatted as "sample__condition".

#' @param conditions A character. The column name corresponding to conditions in the "colData" of "se" argument. Can be NULL if column names of expression matrix in "se" argument are already formatted as "sample__condition".

#' @param aggr Aggregate the replicates by "mean" or "median". Default is "mean".

#' @return A "SummarizedExperiment" object containing averaged expression matrix and metadata. The column names of the expression matrix are formatted as "sample__condition". 

#' @examples
#' ## The toy data below is truncated from GEO dataset GSE14502.
#' # Targets file.
#' path1 <- system.file('extdata/example_data/target_geo.txt', package='spatialHeatmap')
#' target.geo <- read.table(path1, header=TRUE, row.names=1, sep='\t',  stringsAsFactors=F); target.geo[1:3, ]

#' # Gene expression matrix.
#' path2 <- system.file('extdata/example_data/arab_geo.txt', package='spatialHeatmap')
#' expr <- read.table(path2, header=TRUE, row.names=1, sep='\t'); expr[1:3, 1:5]

#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays=list(expr=expr), colData=target.geo)
#' se <- aggr_rep(se=se, samples='samples', conditions='conditions')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @export aggr_rep
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment


aggr_rep <- function(se, samples, conditions, aggr='mean') {

  mat <- assay(se); col.meta <- as.data.frame(colData(se))
  fct <- colnames(mat) <- paste0(make.names(col.meta[, samples]), '__', make.names(col.meta[, conditions]))
  # To keep colnames, "X" should be a character, not a factor.
  if (aggr=='mean') mat <- sapply(X=unique(fct), function(x) rowMeans(mat[, fct==x, drop=FALSE]))
  if (aggr=='median') { 
  
    mat <- sapply(X=unique(fct), function(x) Biobase::rowMedians(mat[, fct==x, drop=FALSE]))
    rownames(mat) <- rownames(se)

  }
  col.meta <- col.meta[!duplicated(fct), ]; rownames(col.meta) <- NULL
  se <- SummarizedExperiment(assays=list(expr=mat), rowData=rowData(se), colData=col.meta); return(se)

}




