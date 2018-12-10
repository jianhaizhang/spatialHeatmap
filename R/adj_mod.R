#' Compute Adjacency Matrix and Identify Modules
#'
#' It designed to compute adjacency matrix and define modules. The resulting adjacency matrix (adj.txt) and module definition (mod.txt) files are saved in the directory "local_mode_result", both of which should be uploaded to "Compute locally" in the Shiny app, which is launched by running "spatial.hm.all()".
#'
#' @param data The processed data matrix resulting from the function "filter.data()".

#' @param type There are two options in this argument, "signed" and "unsigned". "signed" means both positive and negative adjacency between genes are maintained in network module identification while "unsigned" takes the absolute values of negative adjacency. Refer to "WGCNA" for more details.

#' @param minSize The minimum module size in module identification. The default value is 15. Refer to "WGCNA" for more details.

#' @return The adjacency matrix file "adj.txt" and module assignment file "mod.txt", which are all automatically saved in the directory "local_mode_result."

#' @examples  
#' data.path <- system.file("extdata/example", "gene_expr_ann_row_gen.txt", package = "spatialHeatmap")
#' exp <- filter.data(data=data.path, sep="\t", isRowGen=TRUE, c(0, 0), c(0.1, 10000), "processed_data")
#' adj_mod <- adj.mod(data=exp, type="signed", minSize=15)

#' @section Details:
#' This function relies on "WGCNA". First, it computes the adjacency matrix based on the input data matrix, then converts the resulting adjacency matrix to a scale-free topology using a soft threshold (sft). If the type is "signed" sft=12, if the type is "unsigned" sft=6. Next, the modules are identified using dynamic tree cutting algorithm at four sensitivity levels 0, 1, 2, 3. From 0 to 3, more modules are identified but module sizes are decreasing. In order to enhance the performance of this Spatial Heatmap, modules are only identified at level 2 and 3, since the interactive network functionality of this app works better on smaller modules. Refer to "WGCNA" for more details.

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559 \cr Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. Journal of Statistical Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/ \cr Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. Journal of Statistical Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.

#' @export
#' @importFrom WGCNA adjacency TOMsimilarity 
#' @importFrom flashClust flashClust
#' @importFrom stats dist
#' @importFrom stats as.dist
#' @importFrom dynamicTreeCut cutreeHybrid


adj.mod <- function(data, type, minSize=15) {

  # require(WGCNA); require(flashClust)
  if (!dir.exists("local_mode_result")) dir.create("local_mode_result")

  if (type=="signed") sft <- 12; if (type=="unsigned") sft <- 6
  data <- t(data[, grep("__", colnames(data)), drop=F])
  adj=adjacency(data, power=sft, type=type); diag(adj)=0
  tom <- TOMsimilarity(adj, TOMType="signed")
  dissTOM=1-tom; tree.hclust=flashClust(as.dist(dissTOM), method="average")
  mcol <- NULL
  for (ds in 2:3) {

    tree <- cutreeHybrid(dendro=tree.hclust, pamStage=F, minClusterSize=(minSize-3*ds), cutHeight=0.99, deepSplit=ds, distM=dissTOM)
    mcol <- cbind(mcol, tree$labels)

  }; colnames(mcol) <- as.character(2:3); rownames(mcol) <- 1:nrow(mcol)
  write.table(adj, paste0("./local_mode_result/", "adj.txt"), sep="\t", row.names=T, col.names=T)
  write.table(mcol, paste0("./local_mode_result/", "mod.txt"), sep="\t", row.names=T, col.names=T)
  adj.mod <- list(adj=adj, mod=mcol); return(adj.mod)

}







