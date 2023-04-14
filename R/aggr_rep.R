#' Aggregate "Sample__Condition" Replicates in Data Matrix
#'
#' This function aggregates "sample__condition" (see \code{data} argument) replicates by mean or median. The input data is either a \code{data.frame} or \code{SummarizedExperiment}. 

#' @param aggr Aggregate "sample__condition" replicates by "mean" or "median". The default is "mean". If the \code{data} argument is a \code{SummarizedExperiment}, the "sample__condition" replicates are internally formed by connecting samples and conditions with "__" in \code{colData} slot, and are subsequently replace the original column names in \code{assay} slot. If no condition specified to \code{con.factor}, the data are aggregated by sample replicates. If "none", no aggregation is applied.  

#' @inheritParams filter_data

#' @return The returned value is the same class with the input data, a \code{data.frame} or \code{SummarizedExperiment}. In either case, the column names of the data matrix follows the "sample__condition" scheme.

#' @inherit filter_data examples

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Obenchain V, Hester J, Pagès H (2022). SummarizedExperiment: SummarizedExperiment container. R package version 1.28.0, <https://bioconductor.org/packages/SummarizedExperiment>.
#' R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x

#' @export 
#' @importFrom SummarizedExperiment assays assay rowData colData SummarizedExperiment assayNames<-
#' @importFrom SingleCellExperiment SingleCellExperiment


aggr_rep <- function(data, assay.na=NULL, sam.factor, con.factor=NULL, aggr='mean') {
  options(stringsAsFactors=FALSE)
  # Process data. "dgCMatrix" is converted to "matrix".
  dat.lis <- check_data(data=data, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, usage='aggr')
  mat <- assay(dat.lis$se); fct.cna <- dat.lis$fct.cna
  row.meta <- dat.lis$row.meta; col.meta <- dat.lis$col.meta
  rownames(data) <- rownames(row.meta)
  # To keep colnames, "X" should be a character, not a factor.
  if (aggr=='mean') mat <- vapply(unique(fct.cna), function(x) rowMeans(mat[, fct.cna==x, drop=FALSE]), numeric(nrow(mat)))
  if (aggr=='median') {
    mat <- vapply(unique(fct.cna), function(x) Biobase::rowMedians(mat[, fct.cna==x, drop=FALSE]), numeric(nrow(mat)))
    rownames(mat) <- rownames(data)
  }
  # sce is 'SummarizedExperiment', but not vice versa.
  # if (is(data, 'RangedSummarizedExperiment')) data <- as(data, 'SummarizedExperiment')
  if (is(data, 'data.frame') | is(data, 'matrix') | is(data, 'dgCMatrix')) { return(cbind(mat, row.meta)) } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {  
    col.meta <- col.meta[!duplicated(fct.cna), ]; rownames(col.meta) <- NULL
    # if (is(data, 'SingleCellExperiment')) {
      data <- sce_sub(sce=data, assay.na=assay.na, mat=mat, cna=colnames(mat), col.idx=!duplicated(fct.cna), cdat=col.meta)
      if (is.null(assayNames(data))) assayNames(data)[1] <- 'expr'; return(data)
      #data <- data[, !duplicated(fct.cna)]; colnames(data) <- colnames(mat)
      #assays(data)[[assay.na]] <- mat
      # Erase colnames in assay.
      #colData(data) <- as(col.meta, 'DataFrame'); colnames(data) <- colnames(mat); return(data)
    # } else { data <- SummarizedExperiment(assays=list(expr=mat), rowData=rowData(data), colData=col.meta); return(data) }
    }
}

