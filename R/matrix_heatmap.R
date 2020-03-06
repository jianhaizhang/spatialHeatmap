#' Matrix Heatmap
#'
#' This function represents the input gene in the context of corresponding gene network module, where the rows and columns are sorted by hierarchical clustering dendrograms and the row of input gene is tagged by two lines. The matrix heatmap can be dispalyed in static or a web-browser based interactive mode. If the latter, users can zoom in and out by drawing a rectangle and by double clicking the image, respectively. Users can scale the expression values by gene or sample. \cr The network modules are identified at two alternative sensitivities levels (3, 2) by the function "adj.mod". From 3 to 2, the sensitivity decreases and results in less modules with larger sizes. The same module can also be displayed as an interactive form of a network through the "network" function.

#' @param geneID A gene ID from the expression matrix. 
#' @param data A "SummarizedExperiment" object containing the processed data matrix and metadata returned by the function \code{\link{filter_data}}. In the data matrix, rows are gene IDs and columns are samples/conditions.
#' @param adj.mod The list of "adjacency matrix" and "module" definition retured by the function "adj.mod".
#' @param ds The module identification sensitivity, either 2 or 3. See function "adj.mod" for details.
#' @param scale It specifies whether to scale the heatmap. There are three options: "row" means scale by row, "column" means scale by column, and "no" means no scale. 
#' @param col A vector of two colours. It is used for constructing the colour scale. The default is c('yellow', 'blue').
#' @param main The title of the matrix heatmap.
#' @param title.size A numeric, the size of the title font.
#' @param cexCol A numeric value, the size of column names. Default is 1.
#' @param cexRow A numeric value, the size of row names. Default is 1.
#' @param angleCol The angle of column names. The default is 45.
#' @param angleRow The angle of row names. The default is 45.
#' @param sepcolor The colour of the two lines labelling the target gene. The default is "black".
#' @param sep.width The width of two indicator lines of the input gene.
#' @param static "TRUE" gives the static mode and "FALSE" interactive mode.
#' @param margin A vector of two numbers, specifying bottom and right margins respectively. The default is c(10, 10).
#' @param arg.lis1 A list of additional arguments passed to the "heatmap.2" function from "gplots" package. E.g. list(xlab='sample', ylab='gene'). The default is an empty list.
#' @param arg.lis2 A list of additional arguments passed to the "ggplot" function from "gglot2" package. The default is an empty list. 
#' @return A static image or an interactive application lauched on the web browser. 
#' @examples
#' # Creat the "SummarizedExperiment" class. Refer to the R package "SummarizedExperiment" for more details.
#' data.path <- system.file("extdata/shinyApp/example", "root_expr_row_gen.txt", package = "spatialHeatmap")   
#' ## The expression matrix, where the row and column names should be gene IDs and sample/conditions respectively. This data matrix is truncated from a GEO dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46205), which is already normalised.
#' library(data.table); expr <- fread(data.path, sep='\t', header=TRUE, fill=TRUE)
#' col.na <- colnames(expr)[-ncol(expr)]; row.na <- as.data.frame(expr[, 1])[, 1]
#' expr <- as.matrix(as.data.frame(expr, stringsAsFactors=FALSE)[, -1])
#' rownames(expr) <- row.na; colnames(expr) <- col.na
#' col.met.path <- system.file("extdata/shinyApp/example", "col_metadata.txt", package = "spatialHeatmap") 
#' ## Condition metadata is data frame. It has a column of tissues and a column of contidions, which correspond to columns of the data matrix "expr".
#' col.metadata <- read.table(col.met.path, header=TRUE, row.names=NULL, sep='\t', stringsAsFactors=FALSE)
#' row.met.path <- system.file("extdata/shinyApp/example", "row_metadata.txt", package = "spatialHeatmap")
#' ## Row metadata is a data frame. It has a column of row (gene) annotations, which correspond to rows of the data matrix "expr".
#' row.metadata <- read.table(row.met.path, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE)
#' ## The expression matrix, row metadata, and column metadata are stored in a "SummarizedExperiment" object. The row metadata is optional while column metadata is mandatory. The column names in the expression matrix are not important, since they are ultimately renewed by column metadata.
#' library(SummarizedExperiment); expr <- as.matrix(expr); se <- SummarizedExperiment(assays=list(expr=expr), rowData=row.metadata, colData=col.metadata)  
#' exp <- filter_data(data=se, pOA=c(0, 0), CV=c(0.1, Inf), ann='ann', samples='sample', conditions='condition', dir=NULL) 

#' adj_mod <- adj_mod(data=exp, type="signed", minSize=15, dir=NULL)
#' # The gene "PSAC" is represented in the context of its gene module in the form of a static matrix heatmap.
#' matrix_heatmap(geneID="PSAC", data=exp, adj.mod=adj_mod, ds="2", scale="row", static=TRUE)
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr Andrie de Vries and Brian D. Ripley (2016). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1-20. https://CRAN.R-project.org/package=ggdendro \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Carson Sievert (2018) plotly for R. https://plotly-book.cpsievert.me \cr Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ \cr Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2019). gplots: Various R Programming Tools for Plotting Data. R package version 3.0.1.1.  https://CRAN.R-project.org/package=gplots \cr Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/ 

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom ggdendro dendro_data
#' @importFrom ggplot2 ggplot geom_segment geom_text position_dodge geom_rect theme theme_minimal geom_tile scale_fill_gradient geom_hline
#' @importFrom plotly layout subplot %>%
#' @importFrom stats hclust order.dendrogram as.dendrogram
#' @importFrom gplots heatmap.2
#' @importFrom reshape2 melt 
#' @importFrom graphics image mtext par plot title
#' @importFrom grDevices dev.off png

matrix_heatmap <- function(geneID, data, adj.mod, ds, scale, col=c('yellow', 'blue'), main=NULL, title.size=10, cexCol=1, cexRow=1, angleCol=45, angleRow=45, sepcolor="black", sep.width=0.02, static=TRUE, margin=c(10, 10), arg.lis1=list(), arg.lis2=list()) {

  mods <- adj.mod[["mod"]]; ds <- as.character(ds); gene <- assay(data)
  lab <- mods[, ds][rownames(gene)==geneID]
  if (lab=="0") stop("The input gene is not assigned to any module. Please input a different one.")
  mod <- as.matrix(gene[mods[, ds]==lab, ])
  
  if (static==TRUE) {

    tmp <- system.file("extdata/shinyApp/tmp", package="spatialHeatmap"); pa <- paste0(tmp, '/delete_hm.png')
    png(pa); hm <- heatmap.2(x=mod, scale=scale, main=main, trace="none"); dev.off()
    do.call(file.remove, list(pa))
    # Logical matrix with the same dimensions as module matrix.
    selection <- matrix(rep(FALSE, nrow(mod)*ncol(mod)), nrow=nrow(mod))
    # Select the row of target gene.  
    idx <- which(rev(colnames(hm$carpet) %in% geneID))
    lis1 <- c(arg.lis1, list(x=mod, scale=scale, main=main, margin=margin, col=colorRampPalette(col)(nrow(mod)*ncol(mod)), rowsep=c(idx-1, idx), cexCol=cexCol, cexRow=cexRow, srtRow=angleRow, srtCol=angleCol, dendrogram='both', sepcolor="black", sepwidth=c(sep.width, sep.width), key=TRUE, trace="none", density.info="none", Rowv=TRUE, Colv=TRUE))
    do.call(heatmap.2, lis1)

  } else if (static==FALSE) {

     x <- x1 <- x2 <- y <- y1 <- y2 <- xend <- yend <- value <- NULL 
     dd.gen <- as.dendrogram(hclust(dist(mod))); dd.sam <- as.dendrogram(hclust(dist(t(mod))))
     d.sam <- dendro_data(dd.sam); d.gen <- dendro_data(dd.gen)

     g.dengra <- function(df) {

       ggplot()+geom_segment(data=df, aes(x=x, y=y, xend=xend, yend=yend))+labs(x="", y="")+theme_minimal()+theme(axis.text= element_blank(), axis.ticks=element_blank(), panel.grid=element_blank())

     }

     p.gen <- g.dengra(d.gen$segments)+coord_flip(); p.sam <- g.dengra(d.sam$segments)
     gen.ord <- order.dendrogram(dd.gen); sam.ord <- order.dendrogram(dd.sam); mod.cl <- mod[gen.ord, sam.ord]
     if (scale=="column") mod.cl <- data.frame(scale(mod.cl), stringsAsFactors=FALSE); if (scale=="row") mod.cl <- data.frame(t(scale(t(mod.cl))), stringsAsFactors=FALSE)
     mod.cl$gene <- rownames(mod.cl)
     mod.m <- reshape2::melt(mod.cl, id.vars='gene', measure.vars=colnames(mod)); colnames(mod.m) <- c('gene', 'sample', 'value')
     # Use "factor" to re-order rows and columns as specified in dendrograms. 
     mod.m$gene <- factor(mod.m$gene, levels=rownames(mod.cl)); mod.m$sample <- factor(mod.m$sample, levels=colnames(mod.cl))
     # Plot the re-ordered heatmap.
     lis2 <- c(arg.lis2, list(data=mod.m, mapping=aes(x=sample, y=gene))) 
     g <- do.call(ggplot, lis2)+geom_tile(aes(fill=value), colour="white")+scale_fill_gradient(low=col[1], high=col[2])+theme(axis.text.x=element_text(size=cexRow*10, angle=angleCol), axis.text.y=element_text(size=cexCol*10, angle=angleRow))
     # Label target row/gene.
     g.idx <- which(rownames(mod.cl)==geneID)
     g <- g+geom_hline(yintercept=c(g.idx-0.5, g.idx+0.5), linetype="solid", color=sepcolor, size=sep.width*25)
     ft <- list(family = "sans serif", size=title.size, color='black')
     subplot(p.sam, ggplot(), g, p.gen, nrows=2, shareX=TRUE, shareY=TRUE, margin=0, heights=c(0.2, 0.8), widths=c(0.8, 0.2)) %>% plotly::layout(title=main, font=ft)

   }

}






