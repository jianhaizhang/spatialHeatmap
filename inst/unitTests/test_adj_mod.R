library(data.table); library(SummarizedExperiment)
data.path <- system.file("extdata/shinyApp/example", "root_expr_row_gen.txt", package = "spatialHeatmap")  

test_adj_mod <- function() {

  # Creat the "SummarizedExperiment" class. Refer to the R package "SummarizedExperiment" for more details. 
  ## The expression matrix, where the row and column names should be gene IDs and sample/conditions respectively. This data matrix is truncated from a GEO dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46205), which is already normalised.
  expr <- fread(data.path, sep='\t', header=TRUE, fill=TRUE)
  col.na <- colnames(expr)[-ncol(expr)]; row.na <- as.data.frame(expr[, 1])[, 1]
  expr <- as.matrix(as.data.frame(expr, stringsAsFactors=FALSE)[, -1])
  rownames(expr) <- row.na; colnames(expr) <- col.na
  col.met.path <- system.file("extdata/shinyApp/example", "col_metadata.txt", package = "spatialHeatmap") 
  ## Condition metadata is data frame. It has a column of tissues and a column of contidions, which correspond to columns of the data matrix "expr".
  col.metadata <- read.table(col.met.path, header=TRUE, row.names=NULL, sep='\t', stringsAsFactors=FALSE)
  row.met.path <- system.file("extdata/shinyApp/example", "row_metadata.txt", package = "spatialHeatmap")
  ## Row metadata is a data frame. It has a column of row (gene) annotations, which correspond to rows of the data matrix "expr".
  row.metadata <- read.table(row.met.path, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE)
  ## The expression matrix, row metadata, and column metadata are stored in a "SummarizedExperiment" object. The row metadata is optional while column metadata is mandatory. The column names in the expression matrix are not important, since they are ultimately renewed by column metadata.
  expr <- as.matrix(expr); se <- SummarizedExperiment(assays=list(expr=expr), rowData=row.metadata, colData=col.metadata)  
  exp <- filter_data(data=se, pOA=c(0, 0), CV=c(0.1, Inf), ann='ann', samples='sample', conditions='condition', dir=NULL) 

  adj_mod <- adj_mod(data=exp, type="signed", minSize=20)
  checkTrue(is(adj_mod, "list"))

}

