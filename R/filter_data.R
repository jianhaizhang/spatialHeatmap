#' Filter the Data Matrix
#' 
#' This function is designed to filter the numeric data in class of "data.frame" or "SummarizedExperiment". The filtering builds on two functions \code{\link[genefilter]{pOverA}} and \code{\link[genefilter]{cv}} from the package "genefilter" (Gentleman et al. 2018). 

#' @param data An object of \code{data.frame} or \code{SummarizedExperiment}. In either case, the columns and rows should be sample/conditions and assayed items (\emph{e.g.} genes, proteins, metabolites) respectively. If \code{data.frame}, the column names should follow the naming scheme "sample__condition". The "sample" is a general term and stands for cells, tissues, organs, \emph{etc}., where the values are measured. The "condition" is also a general term and refers to experiment treatments applied to "sample" such as drug dosage, temperature, time points, \emph{etc}. If certain samples are not expected to be colored in "spatial heatmaps" (see \code{\link{spatial_hm}}), they are not required to follow this naming scheme. In the downstream interactive network (see \code{\link{network}}), if users want to see node annotation by mousing over a node, a column of row item annotation could be optionally appended to the last column. \cr In the case of \code{SummarizedExperiment}, the \code{assays} slot stores the data matrix. Similarly, the \code{rowData} slot could optionally store a data frame of row item anntation, which is only relevant to the interactive network. The \code{colData} slot usually contains a data frame with one column of sample replicates and one column of condition replicates. It is crucial that replicate names of the same sample or condition must be identical. \emph{E.g.} If sampleA has 3 replicates, "sampleA", "sampleA", "sampleA" is expected while "sampleA1", "sampleA2", "sampleA3" is regarded as 3 different samples. If original column names in the \code{assay} slot already follow the "sample__condition" scheme, then the \code{colData} slot is not required at all. \cr In the function \code{\link{spatial_hm}}, this argument can also be a numeric vector. In this vector, every value should be named, and values expected to color the "spatial heatmaps" should follow the naming scheme "sample__condition". \cr In certain cases, there is no condition associated with data. Then in the naming scheme of \code{data frame} or \code{vector}, the "__condition" part could be discarded. In \code{SummarizedExperiment}, the "condition" column could be discarded in \code{colData} slot. \cr Note, regardless of data class the double underscore is a special string that is reserved for specific purposes in "spatialHeatmap", and thus should be avoided for naming feature/samples and conditions.

#' @param pOA It specifies parameters of the filter function \code{\link[genefilter]{pOverA}} from the package "genefilter" (Gentleman et al. 2018), where genes with expression values larger than "A" in at least the proportion of "P" samples are retained. The input is a vector of two numbers with the first being "P" and the second being "A". The default is c(0, 0), which means no filter is applied. \cr \emph{E.g.} c(0.1, 2) means genes with expression values over 2 in at least 10\% of all samples are kept. 

#' @param CV It specifies parameters of the filter function \code{\link[genefilter]{cv}} from the package "genefilter" (Gentleman et al. 2018), which filters genes according to the coefficient of variation (CV). The input is a vector of two numbers that specify the CV range. The default is c(-Inf, Inf), which means no filtering is applied. \cr \emph{E.g.} c(0.1, 5) means genes with CV between 0.1 and 5 are kept.

#' @param ann The column name of row item (gene, proteins, \emph{etc}.) annotation in the \code{rowData} slot of \code{SummarizedExperiment}. The default is NULL. In \code{\link{filter_data}}, this argument is only relevant if \code{dir} is specified, while in \code{\link{network}} it is only relevant if users want to see annotation when mousing over a node. 

#' @param sam.factor The column name corresponding to samples in the \code{colData} of \code{SummarizedExperiment}. If the original column names in the \code{assay} slot already follows the scheme "sample__condition", then the \code{colData} slot is not required and accordingly this argument could be NULL. 

#' @param con.factor The column name corresponding to conditions in the \code{colData} of \code{SummarizedExperiment}. Could be NULL if column names of in the \code{assay} slot already follows the scheme "sample__condition", or no condition is associated with the data.

#' @param dir The directory path where the filtered data matrix is saved as a TSV-format file "customData.txt", which is ready to upload to the Shiny app launched by \code{\link{shiny_all}}. In the "customData.txt", the rows are assayed items and column names are in the syntax "sample__condition". If gene annotation is provided to \code{ann}, it is appended to "customData.txt". The default is NULL and no file is saved. This argument is used only when the data is stored in \code{SummarizedExperiment} and need to be uploaded to the "customData" in the Shiny app.

#' @return The returned value is the same class with the input data, a \code{data.frame} or \code{SummarizedExperiment}. In either case, the column names of the data matrix follows the "sample__condition" scheme. If \code{dir} is specified, the filtered data matrix is saved in a TSV-format file "customData.txt". 

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
#' # A targets file describing samples and conditions is required for toy data2. It should be 
#' # made based on the experiment design, which is accessible through the accession number 
#' # "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in 
#' # this package and accessed below. 

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
#' # Filter out genes with low counts and low variance. Genes with counts over 5 (log2 unit) in 
#' # at least 1% samples (pOA), and coefficient of variance (CV) between 0.2 and 100 are retained.
#' # Filter toy data1.
#' df.fil.chk <- filter_data(data=df.chk, pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
#' # Filter toy data2.
#' se.fil.chk <- filter_data(data=se.chk, sam.factor='organism_part', con.factor='age', 
#' pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Gentleman, R, V Carey, W Huber, and F Hahne. 2018. "Genefilter: Methods for Filtering Genes from High-Throughput Experiments." http://bioconductor.uib.no/2.7/bioc/html/genefilter.html \cr Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9


#' @export filter_data
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment
#' @importFrom genefilter filterfun pOverA cv genefilter
#' @importFrom utils write.table

filter_data <- function(data, pOA=c(0, 0), CV=c(-Inf, Inf), ann=NULL, sam.factor, con.factor, dir=NULL) {

  options(stringsAsFactors=FALSE)
  # Process data.
  dat.lis <- check_data(data=data, sam.factor=sam.factor, con.factor=con.factor, usage='filter')
  expr <- dat.lis$dat; row.meta <- dat.lis$row.meta; col.meta <- dat.lis$col.meta
  ffun <- filterfun(pOverA(p=pOA[1], A=pOA[2]), cv(CV[1], CV[2]))
  filtered <- genefilter(expr, ffun); expr <- expr[filtered, , drop=FALSE] # Subset one row in a matrix, the result is a numeric vector not a matrix, so drop=FALSE.
  row.meta <- row.meta[filtered, , drop=FALSE]
  if (!is.null(dir)) { 

    dir <- normalizePath(dir, winslash="/", mustWork=FALSE)
    if (!dir.exists(dir)) stop(paste0(dir, ' does not exist!'))
    if (is(data, 'data.frame')|is(data, 'matrix')) {

      expr1 <- cbind.data.frame(expr, row.meta, stringsAsFactors=FALSE)

    }  else if (is(data, 'SummarizedExperiment')) {
      
      if (ncol(row.meta)>0 & !is.null(ann)) {

        expr1 <- cbind.data.frame(expr, row.meta[, ann], stringsAsFactors=FALSE)
        colnames(expr1)[ncol(expr1)] <- ann
      
      } else expr1 <- expr

    }
    write.table(expr1, paste0(dir, "/customData.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
      
  }

  if (is(data, 'data.frame')|is(data, 'matrix')) { return(cbind(expr, row.meta)) } else if (is(data, 'SummarizedExperiment')) {
  
    rownames(col.meta) <- NULL # If row names present in colData(data), if will become column names of assay(data).
    expr <- SummarizedExperiment(assays=list(expr=expr), rowData=row.meta, colData=col.meta); return(expr)

  }

}






