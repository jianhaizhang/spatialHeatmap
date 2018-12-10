#' Filter the Data Matrix
#'
#' It designed to filter large data matrix (e.g.: > 10,000 rows) locally. The two main filter functions are "pOverA" and "cv", which are from the package "genefilter". The processed data matrix (processed_data.txt) is saved in the directory "local_mode_result", which should be uploaded to "Step 2: upload a gene expression file" in the Shiny app, which is launched by running "spatial.hm.all()".

#' @param data The path of the data matrix. In the example of gene expression matrix, the dimension names are gene IDs and sample/conditions. The sample/condition names MUST be fomatted this way: a sample name is followed by double underscore then the condition, such as "sample name__condition name". The meta data (e.g. gene annotation) can also be included in parallel with sample/condition. In the names of sample/condition and meta data, only letters, digits, single underscore, dots are allowed. \cr E.g.: system.file("extdata/example", "gene_expr_ann_row_gen.txt", package = "spatialHeatmap").

#' @param sep The seprator of the data matrix, e.g. ",", "\\t", ";".

#' @param isRowGene It specifies if the row names are used to display spatial heatmaps. The options are "TRUE" or "FALSE". For example, in a gene expression matrix genes are used to display heatmaps and the gene IDs are rows names, then the option is "TRUE".

#' @param pOA It specifies parameters of a filter function that filters according to the proportion of elements exceeding a threshold A. The input is a two-component vector, where the first one is the proportion and the second one is A, e.g.: c(0.1, 2). The default is c(0, 0), which means no filter is applied. Refer to "pOverA" from the package "genefilter". 

#' @param CV It specifies parameters of a filter function that filters according to the coefficient of variation (CV). The input is a two-component vector, where the first and second mean the lower and upper bound of CVs used to filter, e.g.: cv(0.1, 5). The default is cv(0, 10000), which tries to aviod filtering.  Refer to "cv" from the package "genefilter".

#' @param name The processed data matrix is a "\\t" separated "txt" file, and this argument is the name of the "txt file". E.g if name="processed_expression", the resulting file is "processed_expression.txt".

#' @return A data matrix filtered by "pOverA" and "cv" from the R package "genefilter", which is saved in the directory "local_mode_result". It is a "\\t" separated "txt" file and row names (e.g.: gene IDs) will be used to display the spatial heatmaps.

#' @examples

#' data.path <- system.file("extdata/example", "gene_expr_ann_row_gen.txt", package = "spatialHeatmap")
#' exp <- filter.data(data=data.path, sep="\t", isRowGene=TRUE, pOA=c(0, 0), CV=c(0.1, 10000), "processed_data")

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1 \cr Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table 

#' @export
#' @importFrom data.table fread
#' @importFrom genefilter filterfun pOverA cv genefilter
#' @importFrom utils write.table


filter.data <- function(data, sep, isRowGene, pOA=c(0, 0), CV=c(0, 10000), name) {

  # require(data.table); require(genefilter)
  if (!dir.exists("local_mode_result")) dir.create("local_mode_result")

  exp <- fread(data, sep=sep, header=T, fill=T); c.na <- colnames(exp)[-ncol(exp)]
  r.na <- as.data.frame(exp[, 1])[, 1];  df <- as.data.frame(exp, stringsAsFactors=F)[, -1] 
  rownames(df) <- r.na; colnames(df) <- c.na
  if(isRowGene==F) df <- t(df)
  idx <- grep("__", colnames(df)); idx1 <- setdiff(1:length(colnames(df)), idx)
  gene2 <- df[, idx, drop=F]; gene3 <- df[, idx1, drop=F]
  gene2 <- apply(gene2, 2, as.numeric) # This step removes rownames of gene2.

  ffun <- filterfun(pOverA(pOA[1], pOA[2]), cv(CV[1], CV[2]))
  filtered <- genefilter(gene2, ffun)
  gene2 <- gene2[filtered, ]; gene3 <- gene3[filtered, , drop=F]
  # gene2 have no rownames, after "cbind.data.frame" it has same names with gene3.
  df1 <- cbind.data.frame(gene2, gene3, stringsAsFactors=F) 
  write.table(df1, paste0("./local_mode_result/", name, ".txt"), sep="\t", row.names=T, col.names=T); return(df1)

}







