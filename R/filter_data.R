#' Filter the Data Matrix
#' 
#' This function is designed to filter the data matrix, and the filtering is based on two functions \code{\link[genefilter]{pOverA}} and \code{\link[genefilter]{cv}} from the package "genefilter" (Gentleman et al. 2018). 

#' @param data An object of "data frame" or "SummarizedExperiment". If "data frame", the columns and rows should be sample/conditions and assayed items (e.g. genes, proteins, metabolites) respectively. The column names should follow the naming scheme "sample__condition". The "sample" is a general term and stands for cells, tissues, organs, etc. where the values are measured from. The "condition" is also a general term and refers to experiment treatments applied to "sample" such as drug dosage, temperature, time points, etc. If certain samples are not expected to be colored in the "spatial heatmap" (see \code{\link{spatial_hm}}), they are not required to follow this naming scheme. Since double underscore "__" is a reserved separator, it should not be used in any "sample" or "condition". In the downstream interactive network, if users want to see node annotation by mousing over a node, a column of row item annotation could be optionally appended to the last column. See \code{\link{network}} for details. \cr In the case of "SummarizedExperiment", the "assays" slot stores a data matrix with row and column names being assayed items and sample/conditions, respectively. Similarly, the "rowData" slot could optionally store a data frame of row item anntation. The "colData" slot contains a data frame and usually has a column of sample replicates and/or a column of condition replicates. Similarly, the reserved separator "__" should not be used in the sample and condition columns. It is crucial that replicate names of the same sample or condition must be identical. E.g. If sampleA has 3 replicates, "sampleA", "sampleA", "sampleA" is expected while "sampleA1", "sampleA2", "sampleA3" is regarded as 3 different samples. If original column names in the "assay" slot already follow the "sample__condition" scheme, then the "colData" slot is not required at all. \cr In the function \code{\link{spatial_hm}}, this argument can also be a numeric vector. In this vector, every value should be named, and values expected to color the "spatial heatmap" should follow the naming scheme "sample__condition". \cr In certain cases, there is no condition associated with data. Then in the naming scheme of "data frame" or "vector", the "__condition" part could be discarded. In "SummarizedExperiment", the "condition" column could be discarded in "colData" slot. Regardless of data class, "__" is still not allowed in "sample". 

#' @param pOA It specifies parameters of the filter function \code{\link[genefilter]{pOverA}} from the package "genefilter" (Gentleman et al. 2018), where genes with expression values larger than "A" in at least the proportion of "P" samples are retained. The input is a vector of two numbers with the first being "P" and the second being "A". The default is c(0, 0), which means no filter is applied. \cr E.g. c(0.1, 2) means genes with expression values over 2 in at least 10\% of all samples are kept. 

#' @param CV It specifies parameters of the filter function \code{\link[genefilter]{cv}} from the package "genefilter" (Gentleman et al. 2018), which filters genes according to the coefficient of variation (CV). The input is a vector of two numbers, specifying the CV range. The default is c(-Inf, Inf) and means no filtering is applied. \cr E.g. c(0.1, 5) means genes with CV between 0.1 and 5 are kept.

#' @param ann The column name corresponding to row item (gene, proteins, etc.) annotation in the "rowData" slot of "SummarizedExperiment". The default is NULL. In \code{\link{filter_data}}, this argument is only relevant if "dir" is specified, while in \code{\link{network}} it is only relevant if users want to see node annotation when mousing over a node. 

#' @param sam.factor The column name corresponding to samples in the "colData" of "SummarizedExperiment". If the original column names in the "assay" slot already follows the schme "sample__condition", then the "colData" slot is not required and accordingly this argument could be NULL. 

#' @param con.factor The column name corresponding to conditions in the "colData" of "SummarizedExperiment". Could be NULL if column names of in the "assay" slot already follows the schme "sample__condition", or no condition is associated with the data.

#' @param dir The directory path where the folder "local_mode_result/" is created automatically to save the filtered data matrix as tab-separated file "processed_data.txt", which is ready to upload to the Shiny App launched by \code{\link{shiny_all}}. In the "processed_data.txt", the rows are genes and column names are in the syntax "sample__condition". If gene annotation is provided to "ann", it is appended to the last column of "processed_data.txt". This parameter should be the same with that from the function \code{\link{adj_mod}} so that the files "adj.txt" and "mod.txt" generated by \code{\link{adj_mod}} are saved in the same directory. The default is NULL. 

#' @return The returned value is the same class with the input data, a "data frame" or "SummarizedExperiment". In either case, the column names of the data matrix follows the "sample__condition" scheme. If "dir" is specified and input data is "SummarizedExperiment", the filtered data matrix is saved in a "\\t" separated "txt" file (local_mode_result/processed_data.txt). 

#' @examples

#' ## In the following example, the 2 toy data come from an RNA-seq analysis on developments of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019). The complete raw count data are downloaded with the accession number "E-MTAB-6769" using the R package ExpressionAtlas (Keays 2019). Both toy data are generated by truncating the complete data. For conveninece, they are included in this package. Toy data1 is used as a "data frame" input to exemplify data with simple samples/conditions, while toy data2 as "SummarizedExperiment" to illustrate data involving complex samples/conditions.   
#'
#' ## Set up toy data.
#'
#' # Access the toy data1.
#' cnt.chk.simple <- system.file('extdata/shinyApp/example/count_chicken_simple.txt', package='spatialHeatmap')
#' df.chk <- read.table(cnt.chk.simple, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
#' # Note the naming scheme "sample__condition" in columns, where "sample" and "condition" stands for organs and time points respectively.
#' df.chk[1:3, ]

#' # A column of gene annotation can be appended to the data frame, but is not required.  
#' ann <- paste0('ann', seq_len(nrow(df.chk))); ann[1:3]
#' df.chk <- cbind(df.chk, ann=ann)
#' df.chk[1:3, ]
#'
#' # Access the toy data2. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]

#' # A targets file is required for toy data2. It should be made based on the experiment design, which is accessible through the accession number "E-MTAB-6769" in R package ExpressionAtlas. The completed targets file is included in this package. 

#' # Access the targets file. 
#' tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
#' target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
#' # Note every column in toy data2 corresponds with a row in targets file. 
#' target.chk[1:5, ]
#' # Store toy data2 in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene annotation, but not required.
#' rowData(se.chk) <- DataFrame(ann=ann)
#'
#' # Filter genes with low counts and low variance. Genes with counts over 5 (log2 unit) in at least 1% samples (pOA), and coefficient of variance (CV) between 0.2 and 100 are retained.
#' # Filter toy data1.
#' df.fil.chk <- filter_data(data=df.chk, pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
#' # Filter toy data2.
#' se.fil.chk <- filter_data(data=se.chk, sam.factor='organism_part', con.factor='age', pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Gentleman, R, V Carey, W Huber, and F Hahne. 2018. "Genefilter: Methods for Filtering Genes from High-Throughput Experiments." http://bioconductor.uib.no/2.7/bioc/html/genefilter.html \cr Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9


#' @export filter_data
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment
#' @importFrom genefilter filterfun pOverA cv genefilter
#' @importFrom utils write.table

filter_data <- function(data, pOA=c(0, 0), CV=c(-Inf, Inf), ann=NULL, sam.factor, con.factor,  dir=NULL) {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data))
    if (!identical(cna, colnames(data))) cat('Syntactically valid column names are made! \n')
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    row.meta <- data[idx]; expr <- na[!idx]; colnames(expr) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) {

    expr <- assay(data); col.meta <- as.data.frame(colData(data))
    row.meta <- as.data.frame(rowData(data), stringsAsFactors=FALSE)[, , drop=FALSE]
    # Factors teated by paste0/make.names are vectors.
    if (!is.null(sam.factor) & !is.null(con.factor)) { cna <- paste0(col.meta[, sam.factor], '__', col.meta[, con.factor]) } else if (!is.null(sam.factor) & is.null(con.factor)) { cna <- as.vector(col.meta[, sam.factor]) } else if (is.null(sam.factor) & !is.null(con.factor)) { cna <- as.vector(col.meta[, con.factor]) } else cna <- colnames(expr)
    colnames(expr) <- make.names(cna); if (!identical(cna, make.names(cna))) cat('Syntactically valid column names are made! \n')

  }

  ffun <- filterfun(pOverA(p=pOA[1], A=pOA[2]), cv(CV[1], CV[2]))
  filtered <- genefilter(expr, ffun); expr <- expr[filtered, ]
  row.meta <- row.meta[filtered, , drop=FALSE]

  if (!is.null(dir)) { 

    path <- paste0(dir, "/local_mode_result/"); if (!dir.exists(path)) dir.create(path)

    if (is(data, 'data.frame')|is(data, 'matrix')) {

      expr1 <- cbind.data.frame(expr, row.meta, stringsAsFactors=FALSE)

    }  else if (is(data, 'SummarizedExperiment')) {
      
      if (ncol(row.meta)>0 & !is.null(ann)) {
        
        expr1 <- cbind.data.frame(expr, row.meta[, ann], stringsAsFactors=FALSE)
        colnames(expr1)[ncol(expr1)] <- ann
      
      } else expr1 <- expr

    }
    write.table(expr1, paste0(path, "processed_data.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
      
  }

  if (is(data, 'data.frame')|is(data, 'matrix')) { return(cbind(expr, row.meta)) } else if (is(data, 'SummarizedExperiment')) {
  
    rownames(col.meta) <- NULL # If row names present in colData(data), if will become column names of assay(data).
    expr <- SummarizedExperiment(assays=list(expr=expr), rowData=row.meta, colData=col.meta); return(expr)

  }

}






