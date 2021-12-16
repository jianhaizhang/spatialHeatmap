#' Normalize Sequencing Count Matrix
#'
#' This function normalizes sequencing count data. It accepts the count matrix and sample metadata (optional) in form of \code{SummarizedExperiment} or \code{data.frame}. In either class, the columns and rows of the count matix should be sample/conditions and genes respectively.  

#' @inheritParams filter_data   

#' @param norm.fun One of the normalizing functions: "CNF", "ESF", "VST", "rlog", "none". Specifically, "CNF" stands for \code{\link[edgeR]{calcNormFactors}} from edgeR (McCarthy et al. 2012), and "EST", "VST", and "rlog" is equivalent to \code{\link[DESeq2]{estimateSizeFactors}}, \cr \code{\link[DESeq2]{varianceStabilizingTransformation}}, and \code{\link[DESeq2]{rlog}} from DESeq2 respectively (Love, Huber, and Anders 2014). If "none", no normalization is applied. The default is "CNF" and the output data is processed by \code{\link[edgeR]{cpm}} (Counts Per Million). The parameters of each normalization function are provided through \code{parameter.list}.

#' @param parameter.list A list of parameters for each normalizing function assigned in \code{norm.fun}. The default is NULL and \code{list(method='TMM')}, \code{list(type='ratio')}, \cr \code{list(fitType='parametric', blind=TRUE)}, \cr \code{list(fitType='parametric', blind=TRUE)} is internally set for "CNF", "ESF", "VST", "rlog" respectively. Note the slot name of each element in the list is required. \emph{E.g.} \code{list(method='TMM')} is expected while \code{list('TMM')} would cause errors. \cr Complete parameters of "CNF": https://www.rdocumentation.org/packages/edgeR/ \cr versions/3.14.0/topics/calcNormFactors \cr Complete parameters of "ESF": https://www.rdocumentation.org/packages/ \cr DESeq2/versions/1.12.3/topics/estimateSizeFactors \cr Complete parameters of "VST": https://www.rdocumentation.org/packages/ \cr DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation \cr Complete parameters of "rlog": https://www.rdocumentation.org/packages/ \cr DESeq2/versions/1.12.3/topics/rlog

#' @param log2.trans Logical, TRUE or FALSE. If TRUE (default) and the selected normalization method does not use log2 scale by default ("ESF"), the output data is log2-transformed after normalization. If FALSE and the selected normalization method uses log2 scale by default ("VST", "rlog"), the output data is 2-exponent transformed after normalization. 
#' @param data.trans This argument is deprecated and replaced by \code{log2.trans}. One of "log2", "exp2", and "none", corresponding to transform the count matrix by "log2", "2-based exponent", and "no transformation" respecitvely. The default is "none".

#' @return If the input data is \code{SummarizedExperiment}, the retured value is also a \code{SummarizedExperiment} containing normalized data matrix and metadata (optional). If the input data is a \code{data.frame}, the returned value is a \code{data.frame} of normalized data and metadata (optional). 

#' @seealso \code{\link[edgeR]{calcNormFactors}} in edgeR, and \code{\link[DESeq2]{estimateSizeFactors}}, \code{\link[DESeq2]{varianceStabilizingTransformation}}, \code{\link[DESeq2]{rlog}} in DESeq2.

#' @examples

#' ## In the following examples, the 2 toy data come from an RNA-seq analysis on development of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are
#' ## included in this package. The complete raw count data are downloaded using the R package
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as
#' ## a "data frame" input to exemplify data of simple samples/conditions, while toy data2 as
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
#' # Store toy data2 in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk)
#'
#' # Normalize raw count data. The normalizing function "calcNormFactors" (McCarthy et al. 2012)
#' # with default settings is used.  
#' df.nor.chk <- norm_data(data=df.chk, norm.fun='CNF', log2.trans=TRUE)
#' se.nor.chk <- norm_data(data=se.chk, norm.fun='CNF', log2.trans=TRUE)


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' \cr McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' \cr McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export norm_data
#' @importFrom SummarizedExperiment assays
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts varianceStabilizingTransformation rlog


norm_data <- function(data, assay.na=NULL, norm.fun='CNF', parameter.list=NULL, log2.trans=TRUE, data.trans) {
  calls <- names(vapply(match.call(), deparse, character(1))[-1]) 
  if("data.trans" %in% calls) { 
    if (data.trans=='log2') log2.trans <- TRUE else log2.trans <- FALSE
    warning('"data.trans" is deprecated and replaced by "log2.trans"! \n') 
  }
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) {
    dat.lis <- check_data(data=data, assay.na=assay.na, usage='norm'); expr <- dat.lis$dat; ann <- dat.lis$row.meta
  } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
    if (is.null(assay.na)) {
      if (length(assays(data)) > 1) stop("Please specify which assay to use by assigning the assay name to 'assay.na'!") else if (length(assays(data)) == 1) assay.na <- 1
    }
    expr <- SummarizedExperiment::assays(data)[[assay.na]] 
  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", "dgCMatrix", "SummarizedExperiment", or "SingleCellExperiment", except that "spatial_hm" also accepts a "vector".') }
  
  if (is.null(norm.fun)) norm.fun <- 'none'
  if (!norm.fun %in% c("CNF", "ESF", "VST", "rlog", 'none')) stop('"norm.fun" should be one of "CNF", "ESF", "VST", "rlog", or "none"!')
  if (norm.fun=='none') {
    if (log2.trans==TRUE) { 
      v.min <- min(expr); if (v.min<0) expr <- expr-v.min
      expr <- log2(expr+1) 
    } 
  }
  
  if (min(expr)>=0 & all(round(expr)==expr)) {

    if (norm.fun=='CNF') {

      na <- names(parameter.list); if (!('method' %in% na)|is.null(parameter.list)) {
        parameter.list <- c(list(method='TMM'), parameter.list)
      }
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      y <- DGEList(counts=expr); 
      y <- do.call(calcNormFactors, c(list(object=y), parameter.list))
      cat('Computing CPM ... \n')
      expr <- cpm(y, normalized.lib.sizes=TRUE, log=(log2.trans==TRUE))

    } else dds <- DESeqDataSetFromMatrix(countData=expr, colData=data.frame(col.dat=colnames(expr)), design=~1) # "design" does not affect "rlog" and "varianceStabilizingTransformation".

    if (norm.fun=='ESF') {

      na <- names(parameter.list); if (!('type' %in% na)|is.null(parameter.list)) { parameter.list <- c(list(type='ratio'), parameter.list) }
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      dds <- do.call(estimateSizeFactors, c(list(object=dds), parameter.list))
      expr <- counts(dds, normalized=TRUE)
      if (log2.trans==TRUE) expr <- log2(expr+1)

    } else if (norm.fun=='VST') { 
    
      # Apply A Variance Stabilizing Transformation (VST) To The Count Data.
      if (is.null(parameter.list)) parameter.list <- list(fitType='parametric', blind=TRUE)
      na <- names(parameter.list); if (!('fitType' %in% na)) { parameter.list <- c(list(fitType='parametric'), parameter.list) }
      if (!('blind' %in% na)) { parameter.list <- c(list(blind=TRUE), parameter.list) }    
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      # Returns log2-scale data. 
      vsd <- do.call(varianceStabilizingTransformation, c(list(object=dds), parameter.list))
      expr <- SummarizedExperiment::assay(vsd)
      if (log2.trans==FALSE) expr <- 2^expr

    } else if (norm.fun=='rlog') { 
  
      if (is.null(parameter.list)) parameter.list <- list(fitType='parametric', blind=TRUE)
      na <- names(parameter.list); if (!('fitType' %in% na)) { parameter.list <- c(list(fitType='parametric'), parameter.list) }
      if (!('blind' %in% na)) { parameter.list <- c(list(blind=TRUE), parameter.list) }    
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      # Apply A 'Regularized Log' Transformation. 
      rld <- do.call(rlog, c(list(object=dds), parameter.list))
      expr <- SummarizedExperiment::assay(rld) 
      if (log2.trans==FALSE) expr <- 2^expr  

    }

  } else if (min(expr)<0 | !all(round(expr)==expr)) stop('Nornalization only applies to data matrix of all non-negative integers! \n')
  if (log2.trans == FALSE) expr <- round(expr)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) { return(cbind(expr, ann)) } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) { SummarizedExperiment::assays(data)[[assay.na]] <- expr; return(data) }

} 



