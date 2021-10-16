#' A table of differentialy-expressed genes (DEGs) detected by different methods
#'
#' A table of up- and down-DEGs detected by different methods such as edgeR, limma, DEseq2.
#'
#' @docType data
#'
#' @usage data(deg.table)
#'
#' @format A table.
#'
#' @keywords datasets
#'
#' @references
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#'
#' @source 
#' \href{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6769/}{ExpressionAtlas E-MTAB-6769}
#'
#' @examples
#' data(deg.table)
#' deg.table[1:2, ]
"deg.table"
