#' Cutting dendrograms
#'
#' This function is designed to cut dendrograms from heirarchical clustering at a certain height and returns the cluster containing the query biomolecule. 

#' @param dendro A dendrogram from heirarchical clustering.
#' @param h A numeric of height for cutting the dendrogram.
#' @param target The target biomoleclue. 

#' @return A vector of the cluster containing the query biomolecule.

#' @examples

#' blk.mus.pa <- system.file("extdata/shinyApp/data", "bulk_mouse_cocluster.rds", package="spatialHeatmap") 
#' blk.mus <- readRDS(blk.mus.pa)
#' res.hc <- matrix_hm(ID=c('Actr3b'), data=blk.mus, angleCol=60, angleRow=60, cexRow=0.8, cexCol=0.8, margin=c(10, 6), static=TRUE, arg.lis1=list(offsetRow=0.01, offsetCol=0.01))
#' cut_dendro(res.hc$rowDendrogram, h=1000, 'Actr3b')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Tal Galili (2015). dendextend: an R package for visualizing, adjusting, and comparing trees of hierarchical clustering. Bioinformatics. DOI: 10.1093/bioinformatics/btv428

#' @export

cut_dendro <- function(dendro, h, target) {
  pkg <- check_pkg('dendextend'); if (is(pkg, 'character')) {
    warning(pkg); return(pkg)
  }
  clus <- dendextend::cutree(dendro, h = h)[order.dendrogram(dendro)] 
  w <- which(names(clus)==target) 
  cl <- names(clus[clus==clus[w]]); sort(unique(cl)) 
} 

