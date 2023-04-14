#' Filtering the Data Matrix
#' 
#' This function is designed to filter the numeric data in form of \code{data.frame} or \code{SummarizedExperiment}. The filtering builds on two functions \code{\link[genefilter]{pOverA}} and \code{\link[genefilter]{cv}} from the package \pkg{genefilter} (Gentleman et al. 2018). 

#' @param data 
#' \describe{
#'  \item{Terms}{
#'   spatial features: cells, tissues, organs, \emph{etc}; variables: experimental variables such as drug dosage, temperature, time points, \emph{etc}; biomolecules: genes, proteins, metabolites, \emph{etc}; spatial heatmap: SHM.
#'  }
#'  \item{`SummarizedExperiment`}{
#'   The \code{assays} slot stores the data matrix, where rows and columns are biomolecules and spatial featues respectively. Typically, at least two columns of spatial features and variables are stored in the \code{colData} slot respectively. When plotting SHMs, only identical spatial features between the data and aSVG will be colored according to the expression values of chosen biomolecules. Replicates of the same type in these two columns should be identical, \emph{e.g.} "tissueA", "tissueA" rather than "tissueA1", "tissueA2". If column names in the \code{assays} slot follow the "spatialFeature__variable" scheme, \emph{i.e.} spatial features and variables are concatenated by double underscore, then the \code{colData} slot is not required at all. If the data do not have experiment variables, the variable column in \code{colData} or the double underscore scheme is not required.   
#'  }
#'  \item{`data.frame`}{
#'   Rows and columns are biomolecules and spatial featues respectively. If there are experiment variables, the column names should follow the naming scheme "spatialFeature__variable". Otherwise, the column names should only include spatial features. The double underscore is a reserved string for specific purposes in \code{spatialHeatmap}, and thus should be avoided for naming spatial feature or variables. A column of biomolecule description can be included. This is only applicable in the interactive network graph (see \code{\link{network}}), where mousing over a node displays the corresponding description. 
#'  }
#'  \item{vector}{
#'   In the function \code{\link{shm}}, the data can be provided in a numeric \code{vector} for testing with a single gene. If so, the naming schme of the vector is the same with the \code{data.frame}.   
#'  }
#'  \item{Multiple variables}{
#'   For plotting SHMs, multiple variables contained in the data can be combined into a composite one, and the composite variable will be treated as a regular single variable. See the vignette for more details by running \code{browseVignettes('spatialHeatmap')} in R.  
#'  }
#' }

#' @param assay.na The name of target assay to use when \code{data} is \code{SummarizedExperiment}.  
#' @param pOA Parameters of the filter function \code{\link[genefilter]{pOverA}} from the package \pkg{genefilter} (Gentleman et al. 2018). Genes with expression values >= "A" at the proportion >= "P" of all samples are retained. It is a vector of two numbers, where the first and second is "P" and "A" respectively. The default is \code{c(0, 0)}, which means no filtering is applied.  
#' @param CV Parameters of the filter function \code{\link[genefilter]{cv}} from the package \pkg{genefilter} (Gentleman et al. 2018). Genes with coefficient of variation (CV) between CV1 and CV2 are retained. It is a vector of two numbers, where the first and second is CV1 and CV2 respectively. The default is \code{c(-Inf, Inf)}, which means no filtering is applied.  

#' @param top.CV The proportion (0-1) of top coefficient of variations (CVs). Only genes with CVs in this proportion are kept. \emph{E.g.} if \code{top.CV=0.7}, only genes with CVs ranked in the top 70\% are retained. Default is 1, which means all genes are retained. Note this argument takes precedence over \code{CV}.

#' @param desc A column name in the \code{rowData} slot of \code{SummarizedExperiment}. The default is NULL. In \code{\link{filter_data}}, this argument is only applicable if \code{file} is specified. 

#' @param sam.factor The column name corresponding to spatial features in \code{colData} of \code{SummarizedExperiment}. If the column names in the \code{assay} slot already follows the scheme "spatialFeature__variable", then the \code{colData} slot is not required and accordingly this argument could be NULL. 

#' @param con.factor The column name corresponding to experimental variables in \code{colData} of \code{SummarizedExperiment}. It can be NULL if column names of in the \code{assay} slot already follows the scheme "spatialFeature__variable", or no variable is associated with the data.

#' @param file The output file name for saving the filtered data matrix in TSV format, which is ready to upload in the Shiny App (see \code{\link{shiny_shm}}). In this file, the rows are biomoleclues and column names are in the scheme "spatialFeature__variable". If row annotation is provided to \code{desc}, it will be appended to the output file. The default is NULL and no file is saved. This argument is applicable only when the data is in \code{SummarizedExperiment} and need to be uploaded to the the Shiny App. 

#' @param verbose Logical. If TRUE (default), the summary of statistics is printed. 

#' @return An object of a \code{data.frame} or \code{SummarizedExperiment}, depending on the input data. 

#' @examples
#'
#' ## Two example data sets are showcased for the data formats of "data.frame" and 
#' ## "SummarizedExperiment" respectively. Both come from an RNA-seq analysis on 
#  ## development of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019).
#' ## For conveninece, they are included in this package. The complete raw count data are
#' ## downloaded using the R package ExpressionAtlas (Keays 2019) with the accession 
#' ## number "E-MTAB-6769". 
#'
#' # Access example data 1.
#' df.chk <- read.table(system.file('extdata/shinyApp/data/count_chicken_simple.txt', 
#' package='spatialHeatmap'), header=TRUE, row.names=1, sep='\t', check.names=FALSE)
#'
#' # Column names follow the naming scheme
#' # "spatialFeature__variable".  
#' df.chk[1:3, ]
#'
#' # A column of gene description can be optionally appended.
#' ann <- paste0('ann', seq_len(nrow(df.chk))); ann[1:3]
#' df.chk <- cbind(df.chk, ann=ann)
#' df.chk[1:3, ]
#'
#' # Access example data 2. 
#' count.chk <- read.table(system.file('extdata/shinyApp/data/count_chicken.txt', 
#' package='spatialHeatmap'), header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing spatial features and variables is required for example  
#' # data 2, which should be made based on the experiment design. 
#'
#' # Access the targets file. 
#' target.chk <- read.table(system.file('extdata/shinyApp/data/target_chicken.txt', 
#' package='spatialHeatmap'), header=TRUE, row.names=1, sep='\t')
#' # Every column in example data 2 corresponds with a row in the targets file. 
#' target.chk[1:5, ]
#' # Store example data 2 in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can optionally store a data frame of gene annotation.
#' rowData(se.chk) <- DataFrame(ann=ann)
#'
#' # Normalize data.
#' df.chk.nor <- norm_data(data=df.chk, norm.fun='CNF', log2.trans=TRUE)
#' se.chk.nor <- norm_data(data=se.chk, norm.fun='CNF', log2.trans=TRUE)
#'
#' # Aggregate replicates of "spatialFeature_variable", where spatial features are organs
#' # and variables are ages.
#' df.chk.aggr <- aggr_rep(data=df.chk.nor, aggr='mean')
#' df.chk.aggr[1:3, ]
#'
#' se.chk.aggr <- aggr_rep(data=se.chk.nor, sam.factor='organism_part', con.factor='age',
#' aggr='mean')
#' assay(se.chk.aggr)[1:3, 1:3]
#'
#' # Genes with experssion values >= 5 in at least 1% of all samples (pOA), and coefficient
#' # of variance (CV) between 0.2 and 100 are retained.
#' df.chk.fil <- filter_data(data=df.chk.aggr, pOA=c(0.01, 5), CV=c(0.2, 100))
#' se.chk.fil <- filter_data(data=se.chk.aggr, sam.factor='organism_part', con.factor='age', 
#' pOA=c(0.01, 5), CV=c(0.2, 100), file=NULL)
#'
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Gentleman, R, V Carey, W Huber, and F Hahne. 2018. "Genefilter: Methods for Filtering Genes from High-Throughput Experiments." http://bioconductor.uib.no/2.7/bioc/html/genefilter.html \cr Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' \cr Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x


#' @export filter_data
#' @importFrom SummarizedExperiment assays rowData colData SummarizedExperiment assayNames assayNames<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom genefilter filterfun pOverA cv genefilter
#' @importFrom utils write.table
#' @importFrom stats sd

filter_data <- function(data, assay.na=NULL, pOA=c(0, 0), CV=c(-Inf, Inf), top.CV=1, desc=NULL, sam.factor=NULL, con.factor=NULL, file=NULL, verbose=TRUE) {
  options(stringsAsFactors=FALSE)
  if (top.CV>1|top.CV<0) stop('"top.CV" should be between 0 and 1!')
  # Process data.
  dat.lis <- check_data(data=data, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, usage='filter')
  expr <- dat.lis$dat; row.meta <- dat.lis$row.meta; col.meta <- dat.lis$col.meta
  if (verbose==TRUE) { cat('All values before filtering:\n'); print(summary(unlist(as.data.frame(expr)))) }
  expr.t <- as.data.frame(t(expr)); cv.all <- sort(vapply(expr.t, sd, numeric(1))/vapply(expr.t, mean, numeric(1)), decreasing=TRUE)
  if (verbose==TRUE) { cat('All coefficient of variances (CVs) before filtering:\n'); print(summary(cv.all)) }
  if (top.CV<1) {
    cv.min <- cv.all[ceiling(length(cv.all)*top.CV)]
    CV <- c(cv.min, Inf)
  }
  ffun <- filterfun(pOverA(p=pOA[1], A=pOA[2]), cv(CV[1], CV[2]))
  filtered <- genefilter(expr, ffun); expr <- expr[filtered, , drop=FALSE] # Subset one row in a matrix, the result is a numeric vector not a matrix, so drop=FALSE.
  if (verbose==TRUE) { cat('All values after filtering:\n'); print(summary(unlist(as.data.frame(expr)))) }
  expr.t <- as.data.frame(t(expr))
  if (verbose==TRUE & nrow(expr)>0) {
    cv.all <- sort(vapply(expr.t, sd, numeric(1))/vapply(expr.t, mean, numeric(1)), decreasing=TRUE)
    cat('All coefficient of variances (CVs) after filtering:\n'); print(summary(cv.all))
  }
  row.meta <- row.meta[filtered, , drop=FALSE]
  rownames(row.meta) <- NULL # Some rownames could have been automatically mangled
                             # earlier to keep them unique so they don't necessarily
                             # match 'rownames(expr)' anymore. This will cause the
                             # call to 'SummarizedExperiment()' below to fail with
                             # SummarizedExperiment 1.23.2.
  if (!is.null(file)) {
    file <- normalizePath(file, winslash="/", mustWork=FALSE)
    if (file.exists(file)) stop(paste0(file, ' already exists!'))
    if (is(data, 'data.frame')|is(data, 'matrix')) {
      expr1 <- cbind.data.frame(expr, row.meta, stringsAsFactors=FALSE)
    }  else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
      
      if (ncol(row.meta)>0 & !is.null(desc)) {
        expr1 <- cbind.data.frame(expr, row.meta[, desc], stringsAsFactors=FALSE)
        colnames(expr1)[ncol(expr1)] <- desc
      } else expr1 <- expr

    }
    write.table(expr1, file, sep="\t", row.names=TRUE, col.names=TRUE)  
  }

  if (is(data, 'data.frame')|is(data, 'DFrame')|is(data, 'matrix')|is(data, 'dgCMatrix')) { return(cbind(expr, row.meta)) } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
    if (is.null(assay.na)) {
      if (length(assays(data)) > 1) stop("Please specify which assay to use by assigning the assay name to 'assay.na'!") else if (length(assays(data)) == 1) assay.na <- 1
    }
    rownames(col.meta) <- NULL # If row names present in colData(data), if will become column names of assay(data).
    # if (is(data, 'SingleCellExperiment')) { 
    # expr <- SingleCellExperiment(assays=list(expr=expr), rowData=row.meta, colData=col.meta)
     expr <- sce_sub(sce=data, mat=expr, cna=colnames(expr), assay.na=assay.na, row.idx=filtered, rdat=row.meta, cdat=col.meta)
     if (is.null(assayNames(data))) assayNames(expr)[1] <- 'expr'; return(expr)
#else {
#expr <- SummarizedExperiment(assays=list(expr=expr), rowData=row.meta, colData=col.meta) } 
#    }
  }

}
