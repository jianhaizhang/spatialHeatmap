#' Hierarchical clustering combined with matrix heatmap
#'
#' Given a data matrix returned by \code{submatrix}, hierarchical clustering is performed on rows and columns respectively and the results are presented in a matrix heatmap, which supports static and interactive modes. In the matrix heatmap, rows and columns are sorted by hierarchical clustering dendrograms and rows of target biomolecules are tagged by black lines. In the interactive heatmap, users can zoom in and out by drawing a rectangle and by double clicking, respectively.

#' @param ID A vector of biomolecules of interest in the data matrix. 
#' @param data The subsetted data matrix returned by the function \code{\link{submatrix}}.
#' @param assay.na Applicable when \code{data} is `SummarizedExperiment` or `SingleCellExperiment`, where multiple assays could be stored. The name of target assay to use.
#' @param cut.h A numeric of the cutting height in the row dendrograms.
#' @param scale One of `row`, `column`, or `no` (default), corresponding to scale the heatmap by row, column, or no scaling respectively.
#' @param col A character vector of color ingredients for the color scale. The default is \code{c('yellow', 'orange', 'red')}.
#' @param col.n The number of colors in palette.
#' @param keysize A numeric value indicating the size of the color key.
#' @param main The title of the matrix heatmap.
#' @param title.size A numeric value of the title size.
#' @param cexCol A numeric value of column name size. Default is 1.
#' @param cexRow A numeric value of row name size. Default is 1.
#' @param angleCol The angle of column names. The default is 45.
#' @param angleRow The angle of row names. The default is 45.
#' @param sep.color The color of the two lines labeling the row of \code{ID}. The default is "black".
#' @param sep.width The width of two lines labeling the row of \code{ID}. The default is 0.02.
#' @param static Logical, `TRUE` and `FALSE` returns the static and interactive matrix heatmap respectively. 
#' @param margin A vector of two numbers, specifying bottom and right margins respectively. The default is c(10, 10).
#' @param arg.lis1 A list of additional arguments passed to the \code{\link[gplots]{heatmap.2}} function from "gplots" package. \emph{E.g.} `list(xlab='sample', ylab='gene')`.
#' @param arg.lis2 A list of additional arguments passed to the \code{\link[ggplot2]{ggplot}} function from "ggplot2" package.

#' @return A static or interactive matrix heatmap.

#' @inherit submatrix examples

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
#' Andrie de Vries and Brian D. Ripley (2016). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1-20. https://CRAN.R-project.org/package=ggdendro
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' Carson Sievert (2018) plotly for R. https://plotly-book.cpsievert.me
#' Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559
#' R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2019). gplots: Various R Programming Tools for Plotting Data. R package version 3.0.1.1.  https://CRAN.R-project.org/package=gplots
#' Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/ 
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' Matt Dowle and Arun Srinivasan (2019). data.table: Extension of `data.frame`. R package version 1.12.8. https://CRAN.R-project. org/package=data.table

#' @export
#' @importFrom SummarizedExperiment assays
#' @importFrom ggdendro dendro_data
#' @importFrom ggplot2 ggplot geom_segment geom_text position_dodge geom_rect theme theme_minimal geom_tile scale_fill_gradient geom_hline
#' @importFrom stats hclust order.dendrogram as.dendrogram
#' @importFrom gplots heatmap.2 
#' @importFrom graphics image mtext par plot title
#' @importFrom grDevices dev.off png

matrix_hm <- function(ID, data, assay.na=NULL, scale='row', col=c('yellow', 'red'), cut.h, col.n=200, keysize=1.8, main=NULL, title.size=10, cexCol=1, cexRow=1, angleCol=45, angleRow=45, sep.color="black", sep.width=0.02, static=TRUE, margin=c(10, 10), arg.lis1=list(), arg.lis2=list()) {
 # save(ID, data, assay.na, scale, col, col.n, keysize, main, title.size, cexCol, cexRow, angleCol, angleRow, sep.color, sep.width, static, margin, arg.lis1, arg.lis2, file='matrix.hm.arg')
  pkg <- check_pkg('dendextend'); if (is(pkg, 'character')) {
    warning(pkg); return(pkg)
  }
  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) {
    dat.lis <- check_data(data=data); gene <- dat.lis$dat
  } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
    if (is.null(assay.na)) {
      if (length(assays(data)) > 1) stop("Please specify which assay to use by assigning the assay name to 'assay.na'!") else if (length(assays(data)) == 1) assay.na <- 1
    }; gene <- assays(data)[[assay.na]] 
  } else { 
    stop('Accepted data classes are "data.frame", "matrix", "DFrame", "dgCMatrix", "SummarizedExperiment", or "SingleCellExperiment", except that "spatial_hm" also accepts a "vector".') }
  # The default hclust method is 'complete'.
  mod <- as.matrix(gene); dd.gen <- as.dendrogram(hclust(dist(mod))) 
  if (static==TRUE) {
    tmp <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE); pa <- paste0(tmp, '/delete_hm.png')
    png(pa); hm <- heatmap.2(x=mod, scale=scale, main=main, trace="none"); dev.off()
    do.call(file.remove, list(pa))
    # Select the row of target gene.  
    idx <- which(rev(colnames(hm$carpet) %in% ID))
    # If colour codes are more than 500, the colour key is blank.
    lis1 <- c(arg.lis1, list(x=mod, scale=scale, main=main, margin=margin, col=colorRampPalette(col)(col.n), keysize=keysize, rowsep=c(idx-1, idx), cexCol=cexCol, cexRow=cexRow, srtRow=angleRow, srtCol=angleCol, dendrogram='both', sepcolor=sep.color, sepwidth=c(sep.width, sep.width), key=TRUE, trace="none", density.info="none", Rowv=TRUE, Colv=TRUE))
    # The default hclust method is 'complete'.
    res <- do.call(heatmap.2, lis1)
    # Cutting row dendrograms.
    h.max <- ceiling(rev(dendextend::get_branches_heights(res$rowDendrogram))[1])
    if (missing(cut.h)) cut.h <- h.max * 0.2
    if (cut.h > h.max) cut.h <- h.max 
    clus <- cut_dendro(dd.gen, target=ID, h=cut.h)
    res$cluster <- clus; res$cut.h <- cut.h; invisible(res)
  } else if (static==FALSE) {
     pkg <- check_pkg('plotly'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
     x <- x1 <- x2 <- y <- y1 <- y2 <- xend <- yend <- value <- NULL 
     dd.sam <- as.dendrogram(hclust(dist(t(mod))))
     # Cutting row dendrograms.
     h.max <- ceiling(rev(dendextend::get_branches_heights(dd.gen))[1])
     if (missing(cut.h)) cut.h <- h.max * 0.2
     if (cut.h > h.max) cut.h <- h.max
     clus <- cut_dendro(dd.gen, target=ID, h=cut.h)
     
     d.sam <- dendro_data(dd.sam); d.gen <- dendro_data(dd.gen)
     g.dengra <- function(df) {
       ggplot()+geom_segment(data=df, aes(x=x, y=y, xend=xend, yend=yend))+labs(x="", y="")+theme_minimal()+theme(axis.text= element_blank(), axis.ticks=element_blank(), panel.grid=element_blank())
     }
     p.gen <- g.dengra(d.gen$segments) + geom_hline(yintercept=cut.h, linetype="solid", color = "red") + coord_flip()
     p.sam <- g.dengra(d.sam$segments)
     gen.ord <- order.dendrogram(dd.gen); sam.ord <- order.dendrogram(dd.sam); mod.cl <- mod[gen.ord, sam.ord]
     if (scale=="column") mod.cl <- scale(mod.cl); if (scale=="row") mod.cl <- t(scale(t(mod.cl)))
     mod.cl <- data.frame(mod.cl); mod.cl$gene <- rownames(mod.cl)
     mod.m <- reshape2::melt(mod.cl, id.vars='gene', measure.vars=colnames(mod)); colnames(mod.m) <- c('gene', 'sample', 'value')
     # Use "factor" to re-order rows and columns as specified in dendrograms. 
     mod.m$gene <- factor(mod.m$gene, levels=rownames(mod.cl)); mod.m$sample <- factor(mod.m$sample, levels=colnames(mod.cl))
     # Plot the re-ordered heatmap.
     lis2 <- c(arg.lis2, list(data=mod.m, mapping=aes(x=sample, y=gene))) 
     g <- do.call(ggplot, lis2)+geom_tile(aes(fill=value), colour="white")+scale_fill_gradient(low=col[1], high=col[2])+theme(axis.text.x=element_text(size=cexRow*10, angle=angleCol), axis.text.y=element_text(size=cexCol*10, angle=angleRow))
     # Label target row/gene.
     g.idx <- which(rownames(mod.cl) %in% ID)
     g <- g+geom_hline(yintercept=c(g.idx-0.5, g.idx+0.5), linetype="solid", color=sep.color, linewidth=sep.width*25)
     ft <- list(family = "sans serif", size=title.size, color='black')
     ply <- plotly::subplot(p.sam, ggplot(), g, p.gen, nrows=2, shareX=TRUE, shareY=TRUE, margin=0, heights=c(0.2, 0.8), widths=c(0.8, 0.2))
     ply <- plotly::layout(ply, title=main, font=ft)
     return(list(plot=ply, dendro.row=dd.gen, cluster=clus, cut.h=cut.h))
   }
}
