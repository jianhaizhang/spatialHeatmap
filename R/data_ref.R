#' Calculating relative expression values
#'
#' This function computes relative expression values for plotting spatial heatmaps (SHMs).  

#' @param se A `SummarizedExperiment` object containing non-log-transformed data. The `colData` slot contains a minimum of three columns: `spFeature` (spatial feature), `variable` (experiment variable), and `reference` (reference variable). The first two columns are explained in the `data` argument of the \code{\link{filter_data}} function. The `reference` indicates reference experimental variables that are comma-separated strings such as "variable1,variable2" (see the example below). Relative expression values will be computed based on the references. 
#' @param input.log Logical, `TRUE` or `FALSE`, indicating whether the input assay data are log-transformed or not, repectively. If `TRUE` this function will not compute relative expression values, since the log-transformed values significantly reduces real relative expression levels. Users are expected to ensure the input assay data are non-log-transformed and then select `FALSE`.  
#' @param output.log Logical. If `FALSE`, the output relative expression values are ratios such as treatment/control, while if `TRUE`, it would be log2 fold changes such as log2(treatment/control).  

#' @return A `SummarizedExperiment` object that contains relative expression values.   

#' @examples
#'
#' library(SummarizedExperiment)
#' # Access example data. 
#' se <- readRDS(system.file('extdata/shinyApp/data/mouse_organ.rds',
#' package='spatialHeatmap'))
#' colData(se)[1:4, 1:3]; assay(se)[1:3, 1:3]
#' # Data of relative expression values.
#' se.ref <- data_ref(se, input.log=FALSE)
#' colData(se.ref)[, 1:3]; assay(se.ref)[1:3, 1:3]

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} 

#' @references
#' Morgan M, Obenchain V, Hester J, Pagès H (2022). SummarizedExperiment: SummarizedExperiment container. R package version 1.28.0, <https://bioconductor.org/packages/SummarizedExperiment>.
#' R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export 
#' @importFrom SummarizedExperiment assay assay<- colData 

data_ref <- function(se, input.log, output.log=TRUE) {
  if (TRUE %in% input.log) {
    msg <- 'The input data should not be log-transformed, which reduces the real relative expression levels.'
      warning(msg); return(msg)
  }
  reference <- NULL
  pkg <- check_pkg('BiocGenerics'); if (is(pkg, 'character')) stop(pkg)
  colnames(se) <- cna <- paste0(se$spFeature, '__', se$variable)
  if (!'reference' %in% colnames(colData(se))) {
    msg <- 'The "reference" column is missing in "colData" slot!'
    warning(msg); return(msg)
  }
  se.ref <- subset(se, , !reference %in% 'none')
  cdat.ref <- colData(se.ref); se.ref.cal <- se[, 1]
  for (i in seq_len(ncol(se.ref))) {
    # Isolate the target and reference samples.
    ref0 <- cdat.ref$reference[i]
    ref0 <- unique(strsplit(ref0, ',| ')[[1]])
    if (!all(ref0 %in% se$variable)) {
      msg <- 'Ensure all reference names are valid!'
      warning(msg); return(msg)
    }
    ft.tar <- cdat.ref$spFeature[i]
    var.tar <- cdat.ref$variable[i]
    sams0 <- c(rownames(cdat.ref)[i], paste0(ft.tar, '__', ref0))
    se0 <- se[, colnames(se) %in% sams0]

    # The ratio of target/reference.
    asy <- assay(se0); asy <- asy[, 1]/asy
    if (TRUE %in% output.log) asy <- log2((asy[, 1]+1)/(asy+1))
    asy <- round(asy, 2); cdat0 <- colData(se0)
    se0$reference <- se0$variable
    # Target spFeature keeps the same, while the corresponding variables are different.
    se0$variable <- var.new <- paste0(var.tar, '_VS_', cdat0$variable)
    cna0 <- paste0(ft.tar, '__', var.new)
    colnames(asy) <- colnames(se0) <- cna0; assay(se0) <- asy
    se0 <- se0[, -1]; se.ref.cal <- BiocGenerics::cbind(se.ref.cal, se0)
  }; return(se.ref.cal[, -1])
}
