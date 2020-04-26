#' Normalise Sequencing Count Matrix
#' 
#' This function normalizes sequencing count data. It accepts the count matrix and sample metadata in form of "SummarizedExperiment" object. In "colData" slot, at least replicates of "sample" and "condition" should be included.

#' @inheritParams filter_data   

#' @param norm.fun One of the normalisation functions: "CNF", "ESF", "VST", "rlog", and "none". Specifically, "CNF" stands for \code{\link[edgeR]{calcNormFactors}} from edgeR (McCarthy et al. 2012), and "EST", "VST", and "rlog" is equivalent to \code{\link[DESeq2]{estimateSizeFactors}}, \code{\link[DESeq2]{varianceStabilizingTransformation}}, and \code{\link[DESeq2]{rlog}} from DESeq2 respectively (Love, Huber, and Anders 2014). If "none", no normalisation is applied. Default is "CNF". The parameters of each normalisation function are specified through "parameter.list".

#' @param parameter.list A list of parameters for each normalisation function specified in "norm.fun". Default is NULL and it means list(method='TMM'), list(type='ratio'), list(fitType='parametric', blind=TRUE), list(fitType='parametric', blind=TRUE) is internally set for "CNF", "ESF", "VST", "rlog" respectively. Note the name of each element in the list is required. E.g. list(method='TMM') is expected while list('TMM') causes errors. \cr Complete parameters of CNF: https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors \cr Complete parameters of ESF: https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/estimateSizeFactors \cr Complete parameters of VST: https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation \cr Complete parameters of rlog: https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/rlog

#' @param data.trans One of "log2", "exp2", and "none", corresponding to transform the count matrix by log2, 2-based exponent, and no transformation respecitvely. Default is "none".

#' @return A "SummarizedExperiment" object containing normalised data matrix and metadata. The replicates of samples and conditions are concatenated by double underscore "__" to form "sample__condition" replicates. The column names of the output data matrix are replaced by "sample__condition" replicates. 

#' @seealso \code{\link[edgeR]{calcNormFactors}} in edgeR, and \code{\link[DESeq2]{estimateSizeFactors}}, \code{\link[DESeq2]{varianceStabilizingTransformation}}, \code{\link[DESeq2]{rlog}} in DESeq2.

#' @examples
#' # The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
#' library(ExpressionAtlas); library(SummarizedExperiment)
#' rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]; assay(rse.hum)[1:3, 1:3]
#'
#' # A targets file describing replicates of samples and conditions is required, which is made based on the "colData" slot in the downloaded "RangedSummarizedExperiment" and available in spatialHeatmap. See the "se" parameter for details. 
#' brain.pa <- system.file('extdata/shinyApp/example/target_brain.txt', package='spatialHeatmap')
#' target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
#' # The "organism_part" and "disease" column describes tissue and condition replicates respectively.  
#' target.hum[c(1:3, 41:42), 4:5]
#' # Place the targets file into "colData" slot as a DataFrame class. 
#' colData(rse.hum) <- DataFrame(target.hum)
#' 
#' # For users with little R expertise, if the gene expression matrix comes as a data frame, it should be placed into "SummarizedExperiment" before proceeding to next step. An example is shown below by borrowing a data matrix from the brain data.
#' # Borrow a data matrix.
#' df <- assay(rse.hum); df[1:2, 1:3]
#' # Place the data matrix and targets file (target.hum) into "SummarizedExperiment".
#' rse.hum <- SummarizedExperiment(assay=df, colData=target.hum, rowData=NULL)
#' 
#' # The count matrix is normalised with estimateSizeFactors (type=‘ratio’).
#' se.nor.hum <- norm_data(data=rse.hum, norm.fun='CNF', data.trans='log2')


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' Prudencio, Mercedes, Veronique V Belzil, Ranjan Batra, Christian A Ross, Tania F Gendron, Luc J Pregent, Melissa E Murray, et al. 2015. "Distinct Brain Transcriptome Profiles in C9orf72-Associated and Sporadic ALS." Nat. Neurosci. 18 (8): 1175–82
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97


#' @export norm_data
#' @importFrom SummarizedExperiment assay
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts varianceStabilizingTransformation rlog


norm_data <- function(data, norm.fun='CNF', parameter.list=NULL, data.trans='none') { 
  
  if (is(data, 'data.frame')|is(data, 'matrix')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- colnames(data) 
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    ann <- data[idx]; expr <- na[!idx]; colnames(expr) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) { expr <- SummarizedExperiment::assay(data) }
  
  if (is.null(norm.fun)) norm.fun <- 'no'
  if (!(norm.fun %in% c("CNF", "ESF", "VST", "rlog"))) { 

    if (data.trans=='log2') { 

      v.min <- min(expr); if (v.min<0) expr <- expr-v.min
      expr <- log2(expr+1) 

  } else if (data.trans=='exp2') expr <- 2^expr }
  
  if (min(expr)>=0 & all(round(expr)==expr)) {

    if (norm.fun=='CNF') {

      na <- names(parameter.list); if (!('method' %in% na)|is.null(parameter.list)) {

        parameter.list <- c(list(method='TMM'), parameter.list)

      }
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      y <- DGEList(counts=expr); 
      y <- do.call(calcNormFactors, c(list(object=y), parameter.list))
      expr <- cpm(y, normalized.lib.sizes=TRUE, log=(data.trans=='log2'))

    } else dds <- DESeqDataSetFromMatrix(countData=expr, colData=data.frame(col.dat=colnames(expr)), design=~1) # "design" does not affect "rlog" and "varianceStabilizingTransformation".

    if (norm.fun=='ESF') {

      na <- names(parameter.list); if (!('type' %in% na)|is.null(parameter.list)) { parameter.list <- c(list(type='ratio'), parameter.list) }
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      dds <- do.call(estimateSizeFactors, c(list(object=dds), parameter.list))
      expr <- counts(dds, normalized=TRUE)
      if (data.trans=='log2') expr <- log2(expr+1)

    } else if (norm.fun=='VST') { 
    
      # Apply A Variance Stabilizing Transformation (VST) To The Count Data.
      if (is.null(parameter.list)) parameter.list <- list(fitType='parametric', blind=TRUE)
      na <- names(parameter.list); if (!('fitType' %in% na)) { parameter.list <- c(list(fitType='parametric'), parameter.list) }
      if (!('blind' %in% na)) { parameter.list <- c(list(blind=TRUE), parameter.list) }    
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      # Returns log2-scale data. 
      vsd <- do.call(varianceStabilizingTransformation, c(list(object=dds), parameter.list))
      expr <- SummarizedExperiment::assay(vsd)
      if (data.trans=='exp2') expr <- 2^expr

    } else if (norm.fun=='rlog') { 
  
      if (is.null(parameter.list)) parameter.list <- list(fitType='parametric', blind=TRUE)
      na <- names(parameter.list); if (!('fitType' %in% na)) { parameter.list <- c(list(fitType='parametric'), parameter.list) }
      if (!('blind' %in% na)) { parameter.list <- c(list(blind=TRUE), parameter.list) }    
      cat('Normalising:', norm.fun, '\n'); print(unlist(parameter.list))
      # Apply A 'Regularized Log' Transformation. 
      rld <- do.call(rlog, c(list(object=dds), parameter.list))
      expr <- SummarizedExperiment::assay(rld) 
      if (data.trans=='exp2') expr <- 2^expr  

    }

  } else cat('Nornalisation only applies to data matrix with all non-negative values! \n')

  if (is(data, 'data.frame')|is(data, 'matrix')) { return(cbind(expr, ann)) } else if (is(data, 'SummarizedExperiment')) { SummarizedExperiment::assay(data) <- expr; return(data) }

} 



