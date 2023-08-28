#' Converting gene ids using annotation databases
#'
#' This function is designed to convert one type of gene ids to another type, such as Ensembl, Entrez, UniProt.  

#' @param db A character of the annotation database name such as \code{'org.Hs.eg.db'}.
#' @param data A \code{data.frame} containing the gene ids in rownames to convert from.
#' @param from.id,to.id The type of ids to convert from (\code{'ENSEMBL'}, rownames in \code{data}) and to (\code{'SYMBOL'}) (UniProt) respectively. 
#' @param desc Logical. If \code{TRUE}, the description of each gene will be included.
#' @param other Other information to be added, such as `c('ENZYME', 'PATH')`. See `columns(org.Hs.eg.db)`.   

#' @return A \code{data.frame}.

#' @examples 

#' library(org.Hs.eg.db)
#' # Human Ensembl gene ids are rownames in the data frame.
#' data <- data.frame(tissue1=10:12, tissue2=20:22, row.names=c('ENSG00000006047', 'ENSG00000268433', 'ENSG00000268555'))
#' data <- cvt_id(db='org.Hs.eg.db', data=data, from.id='ENSEMBL', to.id='SYMBOL', desc=TRUE)
#' data

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Pagès H, Carlson M, Falcon S, Li N (2022). _AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor_. R package version 1.60.0, <https://bioconductor.org/packages/AnnotationDbi>.
#' Morgan M, Obenchain V, Hester J, Pagès H (2022). SummarizedExperiment: SummarizedExperiment container. R package version 1.28.  0, <https://bioconductor.org/packages/SummarizedExperiment>.

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<-

cvt_id <- function(db, data, from.id, to.id, desc=FALSE, other=NULL) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    msg <- 'Please install the "AnnotationDbi" package!'
    warning(msg); return(msg)
  }; if (TRUE %in% desc) desc <- 'GENENAME' else desc <- NULL
  if (!is(data, 'SummarizedExperiment')) if (is(as.data.frame(data), 'data.frame')) data <- SummarizedExperiment(assays=data)
  ann <- AnnotationDbi::select(get(db), keys=rownames(data), keytype=from.id, columns=c(from.id, to.id, desc, other))
  rdat <- rowData(data)
  if (ncol(rdat)>0) ann <- cbind(rdat[ann[, from.id], ], ann)
  # If 'to.id' is not available or duplicated, use 'from.id'.
  ids <- ann[, to.id]
  idx <- ids=='' | duplicated(ids) | is.na(ids)
  ann$to.id <- ann[, to.id]
  # Original ids are preserved, not filtered.
  ann$to.id[idx] <- ann[idx, from.id]
  # Useless, since all ids in from.id will be retained in ann even if not found in db.
  inter <- intersect(rownames(data), ann[, from.id])
  # ids from data: not available in the database.
  data.dif <- data[setdiff(rownames(data), inter), , drop=FALSE]
  # Convert ids.
  data <- data[inter, , drop=FALSE]
  ann <- subset(ann, get(from.id) %in% inter & (!duplicated(get(from.id))))
  data <- data[order(rownames(data)), , drop=FALSE]
  ann <- ann[order(ann[, from.id]), , drop=FALSE]
  rownames(data) <- ann$to.id
  if (is(data, 'SummarizedExperiment')|is(data, 'SingleCellExperiment')) {
    rowData(data) <- ann
    if (!is.null(desc)) rowData(data)$desc <- ann$GENENAME
    if (nrow(data.dif)>0) {
      if (!is.null(desc)) rowData(data.dif)$desc <- NA
      data <- rbind(data, data.dif)
    }; 
  } else {
    if (!is.null(desc)) data$desc <- ann$GENENAME 
    if (nrow(data.dif)>0) {
      if (!is.null(desc)) data.dif$desc <- NA
      data <- rbind(data, data.dif)
    }
  }; data
}


