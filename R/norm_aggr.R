#' Normalise Sequencing Count Matrix and Aggregate Sample Replicates

#' This function normalise count matrix from sequencing and averages sample replicates. It inputs the count matrix and sample metadata in form of "SummarizedExperiment". In "colData" slot, at least "sample" and "condition" information should be included.

#' @param se A "SummarizedExperiment" storing a gene expression matrix and metadata. The "assays" slot stores a expression matrix with row and column names being gene IDs and sample/conditions, respectively. The "rowData" can store a data frame of row (gene) anntation, but is optional. \cr The "colData" slot is required and contains a data frame with at least 2 columns corresponding sample and condition information respectively. It is important that only letters, digits, single underscore, dots, single space are allowed in the 2 columns. The 2 columns are ultimately concatenated by double underscore "__" to replace the original column names in the expression matrix. Thus the syntax of the column names in the output expression matrix is "sample__condition". E.g. "shoot_pGL2__hypoxia" (Mustroph et al. 2009), where "shoot_pGL2" is the sample and "hypoxia" is the condition. If the original column names in the expression matrix are already formatted as "sample__condition", then the "colData" slot is optional.  

#' @param method.norm A character of the normalisation methods. Options are "TMM", "TMMwsp", "RLE", "upperquartile" in \code{\link[edgeR]{calcNormFactors}} from edgeR (McCarthy et al. 2012), and "ratio" (i.e. "type" in \code{\link[DESeq2]{estimateSizeFactors}}), "iterate" (i.e. "type" in \code{\link[DESeq2]{estimateSizeFactors}}), "VST" (i.e. \code{\link[DESeq2]{varianceStabilizingTransformation}}), "\code{\link[DESeq2]{rlog}}" from DESeq2 (Love, Huber, and Anders 2014). If "none", no normalisation is applied. Default is "TMM".

#' @inheritParams DESeq2::estimateSizeFactors 

#' @inheritParams DESeq2::estimateDispersions

#' @param data.trans One of "log2", "exp2", and "none", corresponding to transform the count matrix by log2, 2-based exponent, and no transformation respecitvely. Default is "none".

#' @param sam.fctor A character. The column name corresponding to samples in the "colData" of "se" argument. Can be NULL if column names of expression matrix in "se" argument are already formatted as "sample__condition".

#' @param con.factor A character. The column name corresponding to conditions in the "colData" of "se" argument. Can be NULL if column names of expression matrix in "se" argument are already formatted as "sample__condition".

#' @param rep.aggr One of "mean", "median", and "none", methods of aggregating the replicates. Default is "mean". If "none", replicates are not aggregated.

#' @return A "SummarizedExperiment" object containing normalised data matrix and metadata. The column names of the expression matrix are formatted as "sample__condition". 

#' @seealso \code{\link[edgeR]{calcNormFactors}} in edgeR, and \code{\link[DESeq2]{estimateSizeFactors}}, \code{\link[DESeq2]{varianceStabilizingTransformation}}, \code{\link[DESeq2]{rlog}} in DESeq2.

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
#' se <- norm_aggr(se=se, method.norm='none', sam.factor='samples', con.factor='conditions', rep.aggr='mean')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.” Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. “Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation.” Nucleic Acids Research 40 (10): 4288–97
#' Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @export norm_aggr
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts varianceStabilizingTransformation rlog


norm_aggr <- function(se, method.norm='TMM', fitType='parametric', blind=TRUE, data.trans='none', sam.factor, con.factor, rep.aggr='mean') { 
  
  expr <- assay(se); if (is.null(method.norm)) method.norm <- 'no'
  if (!(method.norm %in% c('TMM', 'ratio', 'iterate', 'VST', 'rlog'))) { 

    if (data.trans=='log2') { 

      v.min <- min(expr); if (v.min<0) expr <- expr-v.min
      expr <- log2(expr+1) 

  } else if (data.trans=='exp2') expr <- 2^expr }
  
  if (min(expr)>=0 & all(round(expr)==expr)) {

    if (method.norm=='TMM') {

      cat('Normalising:', method.norm, '\n')
      y <- DGEList(counts=expr); y <- calcNormFactors(y, method=method.norm)
      expr <- cpm(y, normalized.lib.sizes=TRUE, log=(data.trans=='log2'))

    } else dds <- DESeqDataSetFromMatrix(countData=expr, colData=data.frame(col.dat=colnames(expr)), design=~1) # "design" does not affect "rlog" and "varianceStabilizingTransformation".

    if (method.norm %in% c('ratio', 'iterate')) {

      cat('Normalising:', method.norm, '\n') 
      # Estimates the size factors using the "median ratio method". 
      dds <- estimateSizeFactors(dds, type=method.norm)
      expr <- counts(dds, normalized=TRUE)
      if (data.trans=='log2') expr <- log2(expr+1)

    } else if (method.norm=='VST') { 
    
      # Apply A Variance Stabilizing Transformation (VST) To The Count Data.
      cat('Normalising:', method.norm, '\n')
      # Returns log2-scale data. 
      vsd <- varianceStabilizingTransformation(dds, blind=blind, fitType=fit.type); expr <- assay(vsd)
      if (data.trans=='exp2') expr <- 2^expr

    } else if (method.norm=='rlog') { 
  
      cat('Normalising:', method.norm, '\n') 
      # Apply A 'Regularized Log' Transformation. 
      rld <- rlog(dds, blind=blind, fitType=fit.type); expr <- assay(rld) 
      if (data.trans=='exp2') expr <- 2^expr  

    }

  } else cat('Nornalisation only applies to data matrix with all non-negative values! \n')

  if (!any(c('mean', 'median') %in% rep.aggr)) { assay(se) <- expr; return(se) }
  
  cat('Aggregating replicates...', '\n') 
  col.meta <- as.data.frame(colData(se))
  fct <- colnames(expr) <- paste0(make.names(col.meta[, sam.factor]), '__', make.names(col.meta[, con.factor]))
  # To keep colnames, "X" should be a character, not a factor.
  if (rep.aggr=='mean') mat <- sapply(X=unique(fct), function(x) rowMeans(expr[, fct==x, drop=FALSE]))
  if (rep.aggr=='median') { 
  
    mat <- sapply(X=unique(fct), function(x) Biobase::rowMedians(expr[, fct==x, drop=FALSE]))
    rownames(mat) <- rownames(se)

  }; col.meta <- col.meta[!duplicated(fct), ]; rownames(col.meta) <- NULL
  return(SummarizedExperiment(assays=list(expr=mat), rowData=rowData(se), colData=col.meta))
 
} 



