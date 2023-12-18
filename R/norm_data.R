#' Normalize Sequencing Count Matrix
#'
#' This function normalizes sequencing count data in form of \code{SummarizedExperiment} or \code{data.frame}. In either class, the columns and rows of the count matix should be samples/conditions and genes respectively.  

#' @inheritParams filter_data   

#' @param norm.fun Normalizing functions, one of "CNF", "ESF", "VST", "rlog", "none". Specifically, "CNF" stands for \code{\link[edgeR]{calcNormFactors}} from edgeR (McCarthy et al. 2012), and "EST", "VST", and "rlog" is equivalent to \code{\link[DESeq2]{estimateSizeFactors}}, \cr \code{\link[DESeq2]{varianceStabilizingTransformation}}, and \code{\link[DESeq2]{rlog}} from DESeq2 respectively (Love, Huber, and Anders 2014). If "none", no normalization is applied. The default is "CNF" and the output data is processed by \code{\link[edgeR]{cpm}} (Counts Per Million). The parameters of each normalization function are provided through \code{par.list}.

#' @param par.list A list of parameters for each normalizing function assigned in \code{norm.fun}. The default is NULL and \code{list(method='TMM')}, \code{list(type='ratio')}, \cr \code{list(fitType='parametric', blind=TRUE)}, \cr \code{list(fitType='parametric', blind=TRUE)} is internally set for "CNF", "ESF", "VST", "rlog" respectively. Note the slot name of each element in the \code{list} is required, \emph{e.g.} \code{list(method='TMM')} rather than \code{list('TMM')}. \cr Complete parameters of "CNF": https://www.rdocumentation.org/packages/edgeR/ \cr versions/3.14.0/topics/calcNormFactors \cr Complete parameters of "ESF": https://www.rdocumentation.org/packages/ \cr DESeq2/versions/1.12.3/topics/estimateSizeFactors \cr Complete parameters of "VST": https://www.rdocumentation.org/packages/ \cr DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation \cr Complete parameters of "rlog": https://www.rdocumentation.org/packages/ \cr DESeq2/versions/1.12.3/topics/rlog

#' @param log2.trans Logical. If \code{TRUE} (default) and the selected normalization method does not use log2 scale by default ("ESF"), the output data is log2-transformed after normalization. If \code{FALSE} and the selected normalization method uses log2 scale by default ("VST", "rlog"), the output data is 2-exponent transformed after normalization. 
#' @param data.trans This argument is deprecated and replaced by \code{log2.trans}. One of "log2", "exp2", and "none", corresponding to transform the count matrix by "log2", "2-based exponent", and "no transformation" respecitvely. The default is "none".

#' @return An object of \code{SummarizedExperiment} or \code{data.frame}, depending on the input data.   

#' @seealso \code{\link[edgeR]{calcNormFactors}} in edgeR, and \code{\link[DESeq2]{estimateSizeFactors}}, \code{\link[DESeq2]{varianceStabilizingTransformation}}, \code{\link[DESeq2]{rlog}} in DESeq2.

#' @inherit filter_data examples

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export 
#' @importFrom SummarizedExperiment assays assay assays<-
#' @importFrom edgeR DGEList calcNormFactors cpm

norm_data <- function(data, assay.na=NULL, norm.fun='CNF', par.list=NULL, log2.trans=TRUE, data.trans) {
  if (!missing(data.trans)) {
    if ('log2' %in% data.trans) log2.trans <- TRUE else log2.trans <- FALSE
    warning('"data.trans" is deprecated and replaced by "log2.trans"! \n')
  }
  #calls <- names(vapply(match.call(), deparse, character(1))[-1]) 
  #if("data.trans" %in% calls) { 
  #  if (data.trans=='log2') log2.trans <- TRUE else log2.trans <- FALSE
  #  warning('"data.trans" is deprecated and replaced by "log2.trans"! \n') 
  #}
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) {
    dat.lis <- check_data(data=data, assay.na=assay.na, usage='norm'); expr <- dat.lis$dat; ann <- dat.lis$row.meta
  } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
    if (is.null(assay.na)) {
      if (length(assays(data)) > 1) stop("Please specify which assay to use by assigning the assay name to 'assay.na'!") else if (length(assays(data)) == 1) assay.na <- 1
    }
    expr <- assays(data)[[assay.na]] 
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
      na <- names(par.list); if (!('method' %in% na)|is.null(par.list)) {
        par.list <- c(list(method='TMM'), par.list)
      }
      cat('Normalising:', norm.fun, '\n'); print(unlist(par.list))
      y <- DGEList(counts=expr); 
      y <- do.call(calcNormFactors, c(list(object=y), par.list))
      cat('Computing CPM ... \n')
      expr <- cpm(y, normalized.lib.sizes=TRUE, log=(log2.trans==TRUE))
    } else {
      pkg <- check_pkg('DESeq2'); if (is(pkg, 'character')) stop(pkg)
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=expr, colData=data.frame(col.dat=colnames(expr)), design=~1) # "design" does not affect "rlog" and "varianceStabilizingTransformation".
    }
    if (norm.fun=='ESF') {
      na <- names(par.list); if (!('type' %in% na)|is.null(par.list)) { par.list <- c(list(type='ratio'), par.list) }
      cat('Normalising:', norm.fun, '\n'); print(unlist(par.list))
      if (any(c('e', 'w') %in% check_pkg('DESeq2'))) stop('The package "DESeq2" is not detected!')
      dds <- do.call(DESeq2::estimateSizeFactors, c(list(object=dds), par.list))
      expr <- DESeq2::counts(dds, normalized=TRUE)
      if (log2.trans==TRUE) expr <- log2(expr+1)
    } else if (norm.fun=='VST') { 
      # Apply A Variance Stabilizing Transformation (VST) To The Count Data.
      if (is.null(par.list)) par.list <- list(fitType='parametric', blind=TRUE)
      na <- names(par.list); if (!('fitType' %in% na)) { par.list <- c(list(fitType='parametric'), par.list) }
      if (!('blind' %in% na)) { par.list <- c(list(blind=TRUE), par.list) }    
      cat('Normalising:', norm.fun, '\n'); print(unlist(par.list))
      if (any(c('e', 'w') %in% check_pkg('DESeq2'))) stop('The package "DESeq2" is not detected!')
      # Returns log2-scale data. 
      vsd <- do.call(DESeq2::varianceStabilizingTransformation, c(list(object=dds), par.list))
      expr <- assay(vsd)
      if (log2.trans==FALSE) expr <- 2^expr
    } else if (norm.fun=='rlog') { 
      if (is.null(par.list)) par.list <- list(fitType='parametric', blind=TRUE)
      na <- names(par.list); if (!('fitType' %in% na)) { par.list <- c(list(fitType='parametric'), par.list) }
      if (!('blind' %in% na)) { par.list <- c(list(blind=TRUE), par.list) }    
      cat('Normalising:', norm.fun, '\n'); print(unlist(par.list))
      if (any(c('e', 'w') %in% check_pkg('DESeq2'))) stop('The package "DESeq2" is not detected!')
      # Apply A 'Regularized Log' Transformation. 
      rld <- do.call(DESeq2::rlog, c(list(object=dds), par.list))
      expr <- assay(rld) 
      if (log2.trans==FALSE) expr <- 2^expr
    }
  } else if (min(expr)<0 | !all(round(expr)==expr)) stop('Nornalization only applies to data matrix of all non-negative integers! \n')
  if (log2.trans == FALSE) expr <- round(expr)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) { return(cbind(expr, ann)) } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) { assays(data)[[assay.na]] <- expr; return(data) }

} 



