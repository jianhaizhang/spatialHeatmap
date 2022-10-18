#' Subset Target Data for Spatial Enrichment
#'
#' This function subsets the target spatial features (\emph{e.g.} cells, tissues, organs) and variables (\emph{e.g.} experimental treatments, time points) for the subsequent spatial enrichment.   

#' @param data A \code{SummarizedExperiment} object. The \code{colData} slot is required to contain at least two columns of "feature" and "variable" respectively. The \code{rowData} slot can optionally contain a column of discriptions of each gene and the column name should be \code{metadata}.  
#' @param feature The column name of "features" in the \code{colData} slot.  
#' @param ft.sel A vector of at least two selected features for spatial enrichment, which come from the \code{feature} column. The default is \code{NULL} and the first two features will be selected. If \code{all}, then all features will be selected.  
#' @param variable The column name of "variables" in the \code{colData} slot. 
#' @param var.sel A vector of at least two selected variables for spatial enrichment, which come from the \code{variable} column. The default is \code{NULL} and the first two variables will be selected. If \code{all}, then all variables will be selected. 
#' @param com.by One of \code{feature}, \code{variable}. If \code{feature}, pairwise comparisons will be perfomed between the selected features (\code{ft.sel}) and the \code{var.sel} will be treated as replicates. If \code{variable}, pairwise comparisons will be perfomed between the selected variables (\code{var.sub}) and the features will be treated as replicates. The default is \code{feature}. The corresponding column will be moved to the first in the \code{colData} slot and be recognized in the spatial enrichment process.  
#' @param target A single-component vector of the target for spatial enrichment. If \code{com.by='feature'}, the target will be one of the entries in \code{ft.sel}. If \code{com.by='variable'}, the target will be one of the entries in \code{var.sub}. The default is \code{NULL}, and the first entity in \code{ft.sel} is selected, since the default \code{com.by} is \code{feature}. A \code{target} column will be included in the \code{colData} slot and will be recognized in spatial enrichment.   

#' @return A subsetted \code{SummarizedExperiment} object, where the \code{com.by} is placed in the first column in \code{colData} slot.

#' @examples 

#' ## In the following examples, the toy data come from an RNA-seq analysis on development of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, it is
#' ## included in this package. The complete raw count data are downloaded using the R package
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769".   
#'
#' ## Set up toy data.
#' 
#' # Access toy data. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing samples and conditions is required for toy data. It should be made
#' # based on the experiment design, which is accessible through the accession number 
#' # "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in this
#' # package and accessed below. 

#' # Access the count table. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]

#' # Access the example targets file. 
#' tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
#' target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
#' # Every column in toy data corresponds with a row in targets file. 
#' target.chk[1:5, ]
#' # Store toy data in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene metadata, but not required. Only the 
#' # column named "metadata" will be recognized. 
#' # Pseudo row metadata.
#' metadata <- paste0('meta', seq_len(nrow(count.chk))); metadata[1:3]
#' rowData(se.chk) <- DataFrame(metadata=metadata)
#'
#' ## As conventions, raw sequencing count data should be normalized and filtered to
#' ## reduce noise. Since normalization will be performed in spatial enrichment, only filtering
#' ## is required before subsetting the data.  
#'
#' # Filter out genes with low counts and low variance. Genes with counts over 5 in
#' # at least 10% samples (pOA), and coefficient of variance (CV) between 3.5 and 100 are 
#' # retained.
#' se.fil.chk <- filter_data(data=se.chk, sam.factor='organism_part', con.factor='age',
#' pOA=c(0.1, 5), CV=c(3.5, 100), dir=NULL)
#' # Subset the data.
#' data.sub <- tar_ref(data=se.fil.chk, feature='organism_part', ft.sel=c('brain', 'heart',
#' 'kidney'), variable='age', var.sel=c('day10', 'day12'), com.by='feature', target='brain')


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 

#' @export tar_ref
#' @importFrom SummarizedExperiment colData assayNames colData<- assayNames<- 

tar_ref <- function(data, feature, ft.sel=NULL, variable, var.sel=NULL, com.by='feature', target=NULL) {
  cdat <- colData(data)
  cdat[, feature] <- make.names(cdat[, feature])
  cdat[, variable] <- make.names(cdat[, variable])
  fts <- cdat[, feature]; fcts <- cdat[, variable]
  if (is.null(ft.sel)) ft.sel <- unique(fts)[seq_len(2)] else if (ft.sel[1]=='all') ft.sel <- unique(fts)
  if (is.null(var.sel)) var.sel <- unique(fcts)[seq_len(2)] else if (var.sel[1]=='all') var.sel <- unique(fcts)

  idx <- cdat[, feature] %in% ft.sel & cdat[, variable] %in% var.sel
  data <- data[, idx]; cdat <- colData(data)
  rna <- rownames(cdat); tar <- rep('yes', nrow(cdat))
  
  if (com.by=='feature') {
    cdat <- cdat[, c(feature, variable, setdiff(colnames(cdat), c(feature, variable)))]
  } else if (com.by=='variable') { 
    cdat <- cdat[, c(variable, feature, setdiff(colnames(cdat), c(feature, variable)))]
  } else if (com.by=='feature.variable') { # Combine features and variables.
    ft.fct <- paste0(cdat[, feature], '__', cdat[, variable])
    # Irrelevant metadata.
    cdat1 <- cdat[, setdiff(colnames(cdat), c(feature, variable))]
    cdat <- cbind(feature.variable=ft.fct, cdat[, c(feature, variable)], cdat1)
  }
  cdat <- cbind(target=tar, cdat); cna <- colnames(cdat)
  cdat <- cdat[, c(cna[c(2, 1)], setdiff(cna, cna[seq_len(2)]))]
  if (is.null(target)) target <- ft.sel[1]
  cdat$target[!cdat[, 1] %in% target] <- 'no'
  rownames(cdat) <- rna; colData(data) <- cdat
  colnames(data) <- paste0(cdat[, feature], '__', cdat[, variable])
  # Name the assay: required in distinct.
  if (is.null(assayNames(data))) assayNames(data) <- 'count'
  return(data)
}


