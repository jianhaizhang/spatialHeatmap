#' Normalise Sequencing Count Matrix
#' 
#' This function normalise count matrix from sequencing. It inputs the count matrix and sample metadata in form of "SummarizedExperiment". In "colData" slot, at least replicates of "sample" and "condition" should be included.

#' @inheritParams filter_data   

#' @param method.norm A character of the normalisation methods. Options are "TMM", "TMMwsp", "RLE", "upperquartile" in \code{\link[edgeR]{calcNormFactors}} from edgeR (McCarthy et al. 2012), and "ratio" (i.e. "type=ratio" in \code{\link[DESeq2]{estimateSizeFactors}}), "iterate" (i.e. "type=iterate" in \code{\link[DESeq2]{estimateSizeFactors}}), "VST" (i.e. \code{\link[DESeq2]{varianceStabilizingTransformation}}), \code{\link[DESeq2]{rlog}} from DESeq2 (Love, Huber, and Anders 2014). If "none", no normalisation is applied. Default is "TMM".

#' @inheritParams DESeq2::estimateSizeFactors 
#' @inheritParams DESeq2::estimateDispersions
#' @inheritParams DESeq2::varianceStabilizingTransformation

#' @param data.trans One of "log2", "exp2", and "none", corresponding to transform the count matrix by log2, 2-based exponent, and no transformation respecitvely. Default is "none".

#' @return A "SummarizedExperiment" object containing normalised data matrix and metadata. The replicates of samples and conditions are concatenated by double underscore "__" to form "sample__condition" replicates. The column names of the output data matrix are replaced by "sample__condition" replicates. 

#' @seealso \code{\link[edgeR]{calcNormFactors}} in edgeR, and \code{\link[DESeq2]{estimateSizeFactors}}, \code{\link[DESeq2]{varianceStabilizingTransformation}}, \code{\link[DESeq2]{rlog}} in DESeq2.

#' @examples
#' # The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
#' library(ExpressionAtlas); library(SummarizedExperiment)
#' rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]; assay(rse.hum)[1:3, 1:3]
#'
#' # A targets file describing replicates of samples and conditions is required, which is made based on the "colData" slot in the downloaded "RangedSummarizedExperiment" and available in spatialHeatmap. See the "se" parameter for details. 
#' brain.pa <- system.file('extdata/example_data/target_brain.txt', package='spatialHeatmap')
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
#' se.nor.hum <- norm_data(se=rse.hum, method.norm='ratio', data.trans='log2')


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' Prudencio, Mercedes, Veronique V Belzil, Ranjan Batra, Christian A Ross, Tania F Gendron, Luc J Pregent, Melissa E Murray, et al. 2015. "Distinct Brain Transcriptome Profiles in C9orf72-Associated and Sporadic ALS." Nat. Neurosci. 18 (8): 1175–82
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97

#' @export norm_data
#' @importFrom SummarizedExperiment assay
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts varianceStabilizingTransformation rlog


norm_data <- function(se, method.norm='TMM', fitType='parametric', blind=TRUE, data.trans='none') { 
  
  fit.type <- NULL
  expr <- SummarizedExperiment::assay(se); if (is.null(method.norm)) method.norm <- 'no'
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
      vsd <- varianceStabilizingTransformation(dds, blind=blind, fitType=fit.type); expr <- SummarizedExperiment::assay(vsd)
      if (data.trans=='exp2') expr <- 2^expr

    } else if (method.norm=='rlog') { 
  
      cat('Normalising:', method.norm, '\n') 
      # Apply A 'Regularized Log' Transformation. 
      rld <- rlog(dds, blind=blind, fitType=fit.type); expr <- SummarizedExperiment::assay(rld) 
      if (data.trans=='exp2') expr <- 2^expr  

    }

  } else cat('Nornalisation only applies to data matrix with all non-negative values! \n')
  SummarizedExperiment::assay(se) <- expr; return(se)

} 



