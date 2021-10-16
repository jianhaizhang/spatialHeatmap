#' Plot Overlap of Spatially-Enriched Genes Across Methods 
#'
#' In \code{spatial_enrich}, the spatially-enriched genes are detected within each method (edgeR, limma, DESeq2, distinct). This function plot the overlap of these detected genes across methods in form of upset plot (Nils, 2019) and overlap matrix.  

#' @param lis.up.down The list of all up- and down-regulated genes organized by methods (edgeR, limma, DESeq2, distinct), which comes from the returned value by \code{spatial_enrich}.
#' @param type One of \code{up} (default) or \code{down}, which refers to up- or down-regulated genes.
#' @param plot One of \code{upset} (default) or \code{matrix}, which corresponds to upset plot or overlap matrix in the output plot. 
#' @inheritParams UpSetR::upset

#' @return An upset plot or matrix plot, which displays overlap of spatially-enriched genes across methods. 

#' @examples
#'
#' data(lis.deg.up.down)
#' # Overlap of up-regulated brain-specific genes across methods.
#' deg_ovl(lis.deg.up.down, type='up', plot='upset')
#' deg_ovl(lis.deg.up.down, type='up', plot='matrix')
#' # Overlap of down-regulated brain-specific genes across methods.
#' deg_ovl(lis.deg.up.down, type='down', plot='upset')
#' deg_ovl(lis.deg.up.down, type='down', plot='matrix')

#' # See detailed examples in the function spatial_enrich.


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' Nils Gehlenborg (2019). UpSetR: A More Scalable Alternative to Venn and Euler Diagrams for Visualizing Intersecting Sets. R package version 1.4.0. https://CRAN.R-project.org/package=UpSetR
 
#' @seealso \code{spatial_enrich}

#' @export deg_ovl
#' @importFrom UpSetR upset fromList

deg_ovl <- function(lis.up.down, type='up', plot='upset', order.by="degree", nintersects=40, point.size=3, line.size=1, mb.ratio=c(0.6, 0.4), text.scale=1.5) {
  if (type=='up') lis <- lis.up.down$up.lis else if (type=='down') lis <- lis.up.down$down.lis
  if (plot=='upset') {
    upset <- upset(fromList(lis), order.by=order.by, nintersects=nintersects, point.size=point.size, line.size=line.size, mb.ratio=mb.ratio, text.scale=text.scale)
    return(upset)
  } else if (plot=='matrix') {
    names(lis) <- sub('\\.up(\\.\\d+)*$|\\.down(\\.\\d+)*$', '', names(lis))
    g <- deg_ovl_mat(lis); return(g)
  }
}



#' Given a DEG list of different methods, plot the overlap matrix.
#'
#' @param deg.lis The list of all up- and down-regulated genes organized by methods, which comes from the returned value of \code{spatial_enrich}.

#' @return An image of ggplot.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' \cr Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber, Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2020). gplots: Various R Programming Tools for Plotting Data. R package version 3.1.1. https://CRAN.R-project.org/package=gplots 
#' \cr Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient theme_minimal theme element_text coord_fixed geom_text element_blank
 
deg_ovl_mat <- function(deg.lis) {
  Var1 <- Var2 <- value <- NULL
  mat <- vapply(names(deg.lis), function(x) vapply(names(deg.lis), function(y) length(intersect(deg.lis[[x]], deg.lis[[y]])), numeric(1)), numeric(length(deg.lis)))
  mel <- reshape2::melt(mat)
  g <- ggplot(data=mel, aes(x=Var1, y=Var2, fill=value))+geom_tile(colour="white")+scale_fill_gradient(low="lightcyan3", high="darkorange")+theme_minimal()+theme(axis.text=element_text(angle=45, vjust=1, size=10, hjust=1))+coord_fixed()+geom_text(aes(Var2, Var1, label=value), color="black", size=4)+theme(axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.ticks=element_blank()); return(g)
}
