#' Aggregate Sample Replicates in Data Matrix
#'
#' This function aggregate sample/condition replicates by mean or median. The input data is a gene expression matrix in form of "SummarizedExperiment". The "colData" slot contains at least 2 columns corresponding to replicates of samples and conditions respectively. 

#' @param aggr Aggregate concatenated "sample__condition" replicates by "mean" or "median". Default is "mean". The 2 columns describing sample and condition replicates in "colData" slot of "se" parameter are concatenated by double underscore "__" to form "sample__condition" replicates. E.g. In "cerebellum__normal" (Prudencio et al. 2015), "cerebellum" is the sample and "normal" is the condition. The concatenated replicates are used for aggregating and then replacing the original column names of the data matrix in "assay" slot.  If the original column names in the data matrix are already formatted in the syntax "sample__condition", then the "colData" slot is not required in the "se" parameter. 

#' @inheritParams filter_data

#' @param aggr One of "mean", "median", and "none", methods of aggregating the replicates. Default is "mean". If "none", replicates are not aggregated.

#' @return A "SummarizedExperiment" object containing aggregated data matrix and metadata. The column names of the data matrix are formatted as "sample__condition". 

#' @examples

#' # The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
#' library(ExpressionAtlas); library(SummarizedExperiment)
#' rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]; assay(rse.hum)[1:3, 1:3]
#'
#' # A targets file describing replicates of samples and conditions is required, which should be made based on the "colData" slot in "RangedSummarizedExperiment". See the "se" parameter for details. This targets file is available in spatialHeatmap.
#' brain.pa <- system.file('extdata/shinyApp/example/target_brain.txt', package='spatialHeatmap')
#' target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
#' # The "organism_part" and "disease" column describes tissue and condition replicates respectively. Note that the replicates of the same tissue or condition should have the identical name.
#' target.hum[c(1:3, 41:42), 4:5]
#' # Place the targets file into "colData" slot. 
#' colData(rse.hum) <- DataFrame(target.hum)
#' 
#' # For users with little R expertise, if the gene expression matrix comes as a data frame, it should be placed into "SummarizedExperiment" before proceeding to next step. An example is shown below by borrowing a data frame from the brain data.
#' # Borrow a data matrix.
#' df <- assay(rse.hum); df[1:2, 1:3]
#' # Place the data matrix and targets file (target.hum) into "SummarizedExperiment".
#' rse.hum <- SummarizedExperiment(assay=df, colData=target.hum, rowData=NULL)
#' 
#' # The count matrix is normalised with estimateSizeFactors (type=‘ratio’).
#' se.nor.hum <- norm_data(data=rse.hum, method.norm='CNF', data.trans='log2')
#'
#' # Average replicates of concatenated sample__condition.
#' se.aggr.hum <- aggr_rep(data=se.nor.hum, sam.factor='organism_part', con.factor='disease', aggr='mean')
#' assay(se.aggr.hum)[49939:49942, ] # The concatenated tissue__conditions are the column names of the output data matrix.

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Prudencio, Mercedes, Veronique V Belzil, Ranjan Batra, Christian A Ross, Tania F Gendron, Luc J Pregent, Melissa E Murray, et al. 2015. "Distinct Brain Transcriptome Profiles in C9orf72-Associated and Sporadic ALS." Nat. Neurosci. 18 (8): 1175–82
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97

#' @export aggr_rep
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment

aggr_rep <- function(data, sam.factor, con.factor, aggr='mean') {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- colnames(data) 
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    ann <- data[idx]; mat <- na[!idx]; fct <- colnames(mat) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) {

    mat <- assay(data); col.meta <- as.data.frame(colData(data))
    # Factors teated by paste0/make.names are vecters.
    if (!is.null(sam.factor) & !is.null(con.factor)) { fct <- colnames(mat) <- paste0(make.names(col.meta[, sam.factor]), '__', make.names(col.meta[, con.factor])) } else if (!is.null(sam.factor) & is.null(con.factor)) {  fct <- colnames(mat) <- make.names(col.meta[, sam.factor]) } else if (is.null(sam.factor) & !is.null(con.factor)) { fct <- colnames(mat) <- make.names(col.meta[, con.factor]) }
 
  }
  # To keep colnames, "X" should be a character, not a factor.
  if (aggr=='mean') mat <- sapply(X=unique(fct), function(x) rowMeans(mat[, fct==x, drop=FALSE]))
  if (aggr=='median') {
  
    mat <- sapply(X=unique(fct), function(x) Biobase::rowMedians(mat[, fct==x, drop=FALSE]))
    rownames(mat) <- rownames(data)

  }
  
  if (is(data, 'data.frame')|is(data, 'matrix')) { return(cbind(mat, ann)) } else if (is(data, 'SummarizedExperiment')) { 
  
    col.meta <- col.meta[!duplicated(fct), ]; rownames(col.meta) <- NULL
    data <- SummarizedExperiment(assays=list(expr=mat), rowData=rowData(data), colData=col.meta); return(data)

  } 

}


