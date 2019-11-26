library(data.table); library(SummarizedExperiment)
data.path <- system.file("extdata/example", "root_expr_row_gen.txt", package = "spatialHeatmap")  

test_filter.data <- function() {
  ## The expression matrix, where the row and column names should be gene IDs and sample/conditions, respectively.
  expr <- fread(data.path, sep='\t', header=TRUE, fill=TRUE)
  col.na <- colnames(expr)[-ncol(expr)]; row.na <- as.data.frame(expr[, 1])[, 1]
  expr <- as.matrix(as.data.frame(expr, stringsAsFactors=FALSE)[, -1])
  rownames(expr) <- row.na; colnames(expr) <- col.na
  con.path <- system.file("extdata/example", "root_con.txt", package = "spatialHeatmap") 
  ## Condition is a single column data frame.
  con <- read.table(con.path, header=TRUE, row.names=NULL, sep='\t', stringsAsFactors=FALSE)
  ann.path <- system.file("extdata/example", "root_ann.txt", package = "spatialHeatmap")
  ## Gene annotation is a single column data frame.
  ann <- read.table(ann.path, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE)
  ## The expression matrix, gene annotation, and condition are stored in a "SummarizedExperiment" object. Gene annotation and condition are optional.
  expr <- SummarizedExperiment(assays=list(expr=expr), rowData=ann, colData=con)  
  exp <- filter.data(data=expr, pOA=c(0, 0), CV=c(0.1, 10000), dir=NULL) 
  checkTrue(is(exp, "SummarizedExperiment"))
}

test_adj.mod <- function() {
 ## The expression matrix, where the row and column names should be gene IDs and sample/conditions, respectively.
  expr <- fread(data.path, sep='\t', header=TRUE, fill=TRUE)
  col.na <- colnames(expr)[-ncol(expr)]; row.na <- as.data.frame(expr[, 1])[, 1]
  expr <- as.matrix(as.data.frame(expr, stringsAsFactors=FALSE)[, -1])
  rownames(expr) <- row.na; colnames(expr) <- col.na
  con.path <- system.file("extdata/example", "root_con.txt", package = "spatialHeatmap") 
  ## Condition is a single column data frame.
  con <- read.table(con.path, header=TRUE, row.names=NULL, sep='\t', stringsAsFactors=FALSE)
  ann.path <- system.file("extdata/example", "root_ann.txt", package = "spatialHeatmap")
  ## Gene annotation is a single column data frame.
  ann <- read.table(ann.path, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE)
  ## The expression matrix, gene annotation, and condition are stored in a "SummarizedExperiment" object. Gene annotation and condition are optional.
  expr <- SummarizedExperiment(assays=list(expr=expr), rowData=ann, colData=con)  
  exp <- filter.data(data=expr, pOA=c(0, 0), CV=c(0.1, 10000), dir=NULL) 
  adj_mod <- adj.mod(data=exp, type="signed", minSize=20)
  checkTrue(is(adj_mod, "list"))
}

