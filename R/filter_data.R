filter.data <- function(data, sep, isRowGene, pOA=c(0, 0), CV=c(0, 10000), name) {

  library(data.table); library(genefilter)
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
  write.table(df1, paste0("./local_mode_result/", name, ".txt"), sep="\t", row.names=T, 
  col.names=T); return(df1)

}







