#' Aggregate "Sample__Condition" Replicates in Data Matrix
#'
#' This function aggregates "sample__condition" (see \code{data} argument) replicates by mean or median. The input data is either a \code{data.frame} or \code{SummarizedExperiment}. 

#' @param aggr Aggregate "sample__condition" replicates by "mean" or "median". The default is "mean". If the \code{data} argument is a \code{SummarizedExperiment}, the "sample__condition" replicates are internally formed by connecting samples and conditions with "__" in \code{colData} slot, and are subsequently replace the original column names in \code{assay} slot. If no condition specified to \code{con.factor}, the data are aggregated by sample replicates. If "none", no aggregation is applied.  

#' @inheritParams filter_data

#' @return The returned value is the same class with the input data, a \code{data.frame} or \code{SummarizedExperiment}. In either case, the column names of the data matrix follows the "sample__condition" scheme.

#' @examples

#' ## In the following examples, the 2 toy data come from an RNA-seq analysis on developments of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are
#' ## included in this package. The complete raw count data are downloaded using the R package 
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as a
#' ## "data frame" input to exemplify data with simple samples/conditions, while toy data2 as 
#' ## "SummarizedExperiment" to illustrate data involving complex samples/conditions.   
#'
#' ## Set up toy data.
#' 
#' # Access toy data1.
#' cnt.chk.simple <- system.file('extdata/shinyApp/example/count_chicken_simple.txt', 
#' package='spatialHeatmap')
#' df.chk <- read.table(cnt.chk.simple, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
#' # Columns follow the namig scheme "sample__condition", where "sample" and "condition" stands
#' # for organs and time points respectively.
#' df.chk[1:3, ]
#'
#' # A column of gene annotation can be appended to the data frame, but is not required.  
#' ann <- paste0('ann', seq_len(nrow(df.chk))); ann[1:3]
#' df.chk <- cbind(df.chk, ann=ann)
#' df.chk[1:3, ]
#'
#' # Access toy data2. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing samples and conditions is required for toy data2. It should be made
#' # based on the experiment design, which is accessible through the accession number "E-MTAB-6769"
#' # in the R package ExpressionAtlas. An example targets file is included in this package and 
#' # accessed below. 

#' # Access the example targets file. 
#' tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
#' target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
#' # Every column in toy data2 corresponds with a row in targets file. 
#' target.chk[1:5, ]
#' # Store toy data2 in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene annotation, but not required.
#' rowData(se.chk) <- DataFrame(ann=ann)
#'
#' # Aggregate "sample_condition" replicates in toy data1.
#' df.aggr.chk <- aggr_rep(data=df.chk, aggr='mean')
#' df.aggr.chk[1:3, ]
#'
#' # Aggregate "sample_condition" replicates in toy data2, where "sample" is "organism_part" and
#' # "condition" is "age". 
#' se.aggr.chk <- aggr_rep(data=se.chk, sam.factor='organism_part', con.factor='age', aggr='mean')
#' assay(se.aggr.chk)[1:3, 1:3]

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' \cr McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export aggr_rep
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment

aggr_rep <- function(data, sam.factor, con.factor, aggr='mean') {

  options(stringsAsFactors=FALSE)
  # Process data.
  dat.lis <- check_data(data=data, sam.factor=sam.factor, con.factor=con.factor, usage='aggr')
  mat <- dat.lis$dat; fct.cna <- dat.lis$fct.cna; row.meta <- dat.lis$row.meta; col.meta <- dat.lis$col.meta

  # To keep colnames, "X" should be a character, not a factor.
  if (aggr=='mean') mat <- vapply(unique(fct.cna), function(x) rowMeans(mat[, fct.cna==x, drop=FALSE]), numeric(nrow(mat)))
  if (aggr=='median') {
  
    mat <- vapply(unique(fct.cna), function(x) Biobase::rowMedians(mat[, fct.cna==x, drop=FALSE]), numeric(nrow(mat)))
    rownames(mat) <- rownames(data)

  }
  
  if (is(data, 'data.frame')|is(data, 'matrix')) { return(cbind(mat, row.meta)) } else if (is(data, 'SummarizedExperiment')) { 
  
    col.meta <- col.meta[!duplicated(fct.cna), ]; rownames(col.meta) <- NULL
    data <- SummarizedExperiment(assays=list(expr=mat), rowData=rowData(data), colData=col.meta); return(data)

  }

}


