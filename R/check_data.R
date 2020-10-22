#' Check data validity and process dimension names
#'
#' @inheritParams spatial_hm
#' @param usage One of "shm", "aggr", "filter", "norm", "other". The default is "other". If the "data" is SummarizedExperiment, this argument is only relevant to "aggr_rep", "filter_data", "spatial_hm".
#' @return A list of 4 components: dat, fct.cna, col.meta, row.meta, con.na
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/

#' @importFrom SummarizedExperiment assay colData rowData

check_data <- function(data, sam.factor=NULL, con.factor=NULL, usage='other') {

  options(stringsAsFactors=FALSE)
  dat <- fct.cna <- col.meta <- row.meta <- con.na <- NULL
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data))
    if (!identical(cna, colnames(data))) cat('Syntactically valid column names are made! \n')
    # Only in replicate aggregation, data normalization, data filter, replicates are allowed. 
    if (!(usage %in% c('aggr', 'norm', 'filter'))) if (any(duplicated(cna))) stop('Please use function \'aggr_rep\' to aggregate replicates!') 
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    row.meta <- data[idx] # aggr_rep filter_data, submatrix 
    dat <- na[!idx]; colnames(dat) <- fct.cna <- cna[!idx]
    # spatial_hm
    if (usage=='shm') { form <- grepl("__", fct.cna); if (sum(form)==0) { colnames(dat) <- paste0(fct.cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE }

  } else if (is(data, 'SummarizedExperiment')) {

    dat <- assay(data); r.na <- rownames(dat); dat <- apply(dat, 2, as.numeric) # This step removes rownames.
    rownames(dat) <- r.na; cna <- colnames(dat) <- make.names(colnames(dat))
    if (!identical(cna, colnames(dat))) cat('Syntactically valid column names are made! \n')
    col.meta <- as.data.frame(colData(data))
    # filter_data
    row.meta <- as.data.frame(rowData(data))[, , drop=FALSE]
    # Factors teated by paste0/make.names are vectors.
    if (!is.null(sam.factor) & !is.null(con.factor)) { fct.cna <- colnames(dat) <- paste0(col.meta[, sam.factor], '__', col.meta[, con.factor]); con.na <- TRUE } 

    if (usage=='shm') {

      if (!is.null(sam.factor) & is.null(con.factor)) { sam.na <- as.vector(col.meta[, sam.factor]); fct.cna <- paste0(sam.na, "__", "con"); con.na <- FALSE } else if (is.null(sam.factor)) { form <- grepl("__", cna); if (sum(form)==0) { fct.cna <- paste0(cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE }
    
      if (!is.null(fct.cna)) { colnames(dat) <- make.names(fct.cna); if (!identical(fct.cna, make.names(fct.cna))) cat('Syntactically valid column names are made! \n') }
      if (any(duplicated(colnames(dat)))) stop('Please use function \'aggr_rep\' to aggregate \'sample__condition\' replicates!')
    
    } else if (usage %in% c('aggr', 'filter')) {
 
      if (!is.null(sam.factor) & is.null(con.factor)) {  fct.cna <- as.vector(col.meta[, sam.factor]) } else if (is.null(sam.factor) & !is.null(con.factor)) { fct.cna <- as.vector(col.meta[, con.factor]) } else fct.cna <- colnames(dat)
      if (!identical(fct.cna, make.names(fct.cna))) cat('Syntactically valid column names are made! \n')
      fct.cna <- colnames(dat) <- make.names(fct.cna) 

    }

  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", or "SummarizedExperiment", except that "spatial_hm" also accepts a "vector".')
 
  }; return(list(dat=dat, fct.cna=fct.cna, col.meta=col.meta, row.meta=row.meta, con.na=con.na)) 

}




