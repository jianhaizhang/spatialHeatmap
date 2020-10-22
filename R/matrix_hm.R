#' Matrix Heatmap
#'
#' This function visualizes the input assayed items (gene, protein, metabolite, \emph{etc}) in context of their nearest neighbors, which are subsetted by \code{submatrix}. The visualization is in form of  static or interactive matrix heatmap, where rows and columns are sorted by hierarchical clustering dendrograms and the row of target items are tagged by two lines. In the interactive heatmap, users can zoom in and out by drawing a rectangle and by double clicking the image, respectively.

#' @param ID A vector of target item identifiers in the data. 
#' @param data The subsetted data matrix returned by the function \code{\link{submatrix}}, where rows are assayed items and columns are samples/conditions.
#' @param scale One of "row", "column", or "no", corresponding to scale the heatmap by row, column, or no scale respectively. Default is "no".
#' @param col A character vector of color ingredients for constructing the color scale. The default is c('yellow', 'orange', 'red').
#' @param main The title of the matrix heatmap.
#' @param title.size A numeric value of the title size.
#' @param cexCol A numeric value of column name size. Default is 1.
#' @param cexRow A numeric value of row name size. Default is 1.
#' @param angleCol The angle of column names. The default is 45.
#' @param angleRow The angle of row names. The default is 45.
#' @param sep.color The color of the two lines labeling the row of \code{ID}. The default is "black".
#' @param sep.width The width of two lines labeling the row of \code{ID}. The default is 0.02.
#' @param static Logical, TRUE returns the static visualization and FALSE returns the interactive. 
#' @param margin A vector of two numbers, specifying bottom and right margins respectively. The default is c(10, 10).
#' @param arg.lis1 A list of additional arguments passed to the \code{\link[gplots]{heatmap.2}} function from "gplots" package. \emph{E.g.} list(xlab='sample', ylab='gene'). The default is an empty list.
#' @param arg.lis2 A list of additional arguments passed to the \code{\link[ggplot2]{ggplot}} function from "ggplot2" package. The default is an empty list. 
#' @return A static image or an interactive instance lauched on the web browser. 

#' @examples

#' ## In the following examples, the 2 toy data come from an RNA-seq analysis on development of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are
#' ## included in this package. The complete raw count data are downloaded using the R package
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as
#' ## a "data frame" input to exemplify data of simple samples/conditions, while toy data2 as
#' ## "SummarizedExperiment" to illustrate data involving complex samples/conditions.   

#' ## Set up toy data.
#' 
#' # Access toy data1.
#' cnt.chk.simple <- system.file('extdata/shinyApp/example/count_chicken_simple.txt', 
#' package='spatialHeatmap')
#' df.chk <- read.table(cnt.chk.simple, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
#' # Columns follow the namig scheme "sample__condition", where "sample" and "condition" stands
#' # for organs and time points respectively.
#' df.chk[1:3, ]
#'
#' # A column of gene annotation can be appended to the data frame, but is not required.  
#' ann <- paste0('ann', seq_len(nrow(df.chk))); ann[1:3]
#' df.chk <- cbind(df.chk, ann=ann)
#' df.chk[1:3, ]
#'
#' # Access toy data2. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing samples and conditions is required for toy data2. It should be 
#' # made based on the experiment design, which is accessible through the accession number
#' # "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in
#' # this package and accessed below. 

#' # Access the example targets file. 
#' tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
#' target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
#' # Every column in toy data2 corresponds with a row in targets file. 
#' target.chk[1:5, ]
#' # Store toy data2 in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene annotation, but not required.
#' rowData(se.chk) <- DataFrame(ann=ann)
#'
#' ## As conventions, raw sequencing count data should be normalized, aggregated, and filtered
#' ## to reduce noise.
#'
#' # Normalize count data.
#' # The normalizing function "calcNormFactors" (McCarthy et al. 2012) with default settings
#' # is used.  
#' df.nor.chk <- norm_data(data=df.chk, norm.fun='CNF', data.trans='log2')
#' se.nor.chk <- norm_data(data=se.chk, norm.fun='CNF', data.trans='log2')

#' # Aggregate count data.
#' # Aggregate "sample__condition" replicates in toy data1.
#' df.aggr.chk <- aggr_rep(data=df.nor.chk, aggr='mean')
#' df.aggr.chk[1:3, ]

#' # Aggregate "sample_condition" replicates in toy data2, where "sample" is "organism_part"
#' # and "condition" is "age". 
#' se.aggr.chk <- aggr_rep(data=se.nor.chk, sam.factor='organism_part', con.factor='age',
#' aggr='mean')
#' assay(se.aggr.chk)[1:3, 1:3]

#' # Filter out genes with low counts and low variance. Genes with counts over 5 (log2 unit) in
#' # at least 1% samples (pOA), and coefficient of variance (CV) between 0.2 and 100 are retained.
#' # Filter toy data1.
#' df.fil.chk <- filter_data(data=df.aggr.chk, pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
#' # Filter toy data2.
#' se.fil.chk <- filter_data(data=se.aggr.chk, sam.factor='organism_part', con.factor='age',
#' pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
#'
#' ## Select nearest neighbors for target genes 'ENSGALG00000019846' and 'ENSGALG00000000112',
#' ## which are usually genes visualized in spatial heatmaps.
#' # Toy data1.
#' df.sub.mat <- submatrix(data=df.fil.chk, ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), p=0.1)
#' # Toy data2.
#' se.sub.mat <- submatrix(data=se.fil.chk, ann='ann', ID=c('ENSGALG00000019846', 
#' 'ENSGALG00000000112'), p=0.1) 
#'
#' # In the following, "df.sub.mat" and "se.sub.mat" is used in the same way, so only 
#' # "se.sub.mat" illustrated.
#'
#' # The subsetted matrix is partially shown below.
#' se.sub.mat[c('ENSGALG00000019846', 'ENSGALG00000000112'), c(1:2, 63)]

#'
#' ## Matrix heatmap.
#' # Static matrix heatmap.
#' matrix_hm(ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), data=se.sub.mat, angleCol=80,
#' angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(8, 10), static=TRUE, 
#' arg.lis1=list(offsetRow=0.01, offsetCol=0.01))

#' # Interactive matrix heatmap.
#' \donttest{ matrix_hm(ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), data=se.sub.mat, 
#' angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(8, 10), static=FALSE, 
#' arg.lis1=list(offsetRow=0.01, offsetCol=0.01)) }

#' # In case the interactive heatmap is not automatically opened, run the following code snippet.
#' # It saves the heatmap as an HTML file according to the value assigned to the "file" argument.
#' \donttest{
#' mhm <- matrix_hm(ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), data=se.sub.mat, 
#' angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(8, 10), static=FALSE, 
#' arg.lis1=list(offsetRow=0.01, offsetCol=0.01))
#' htmlwidgets::saveWidget(widget=mhm, file='mhm.html', selfcontained=FALSE)
#' browseURL('mhm.html')
#' }


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr Andrie de Vries and Brian D. Ripley (2016). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1-20. https://CRAN.R-project.org/package=ggdendro \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Carson Sievert (2018) plotly for R. https://plotly-book.cpsievert.me \cr Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ \cr Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2019). gplots: Various R Programming Tools for Plotting Data. R package version 3.0.1.1.  https://CRAN.R-project.org/package=gplots \cr Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/ 
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export matrix_hm
#' @importFrom SummarizedExperiment assay
#' @importFrom ggdendro dendro_data
#' @importFrom ggplot2 ggplot geom_segment geom_text position_dodge geom_rect theme theme_minimal geom_tile scale_fill_gradient geom_hline
#' @importFrom plotly layout subplot %>%
#' @importFrom stats hclust order.dendrogram as.dendrogram
#' @importFrom gplots heatmap.2 
#' @importFrom graphics image mtext par plot title
#' @importFrom grDevices dev.off png

matrix_hm <- function(ID, data, scale='no', col=c('yellow', 'orange', 'red'), main=NULL, title.size=10, cexCol=1, cexRow=1, angleCol=45, angleRow=45, sep.color="black", sep.width=0.02, static=TRUE, margin=c(10, 10), arg.lis1=list(), arg.lis2=list()) {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) {
    dat.lis <- check_data(data=data); gene <- dat.lis$dat
  } else if (is(data, 'SummarizedExperiment')) { gene <- assay(data) } else { 
  stop('Accepted data classes are "data.frame", "matrix", "DFrame", or "SummarizedExperiment", except that "spatial_hm" also accepts a "vector".') }
  mod <- as.matrix(gene)
 
  if (static==TRUE) {

    tmp <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE); pa <- paste0(tmp, '/delete_hm.png')
    png(pa); hm <- heatmap.2(x=mod, scale=scale, main=main, trace="none"); dev.off()
    do.call(file.remove, list(pa))
    # Select the row of target gene.  
    idx <- which(rev(colnames(hm$carpet) %in% ID))
    # If colour codes are more than 500, the colour key is blank.
    lis1 <- c(arg.lis1, list(x=mod, scale=scale, main=main, margin=margin, col=colorRampPalette(col)(500), rowsep=c(idx-1, idx), cexCol=cexCol, cexRow=cexRow, srtRow=angleRow, srtCol=angleCol, dendrogram='both', sepcolor=sep.color, sepwidth=c(sep.width, sep.width), key=TRUE, trace="none", density.info="none", Rowv=TRUE, Colv=TRUE))
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
     g <- g+geom_hline(yintercept=c(g.idx-0.5, g.idx+0.5), linetype="solid", color=sep.color, size=sep.width*25)
     ft <- list(family = "sans serif", size=title.size, color='black')
     subplot(p.sam, ggplot(), g, p.gen, nrows=2, shareX=TRUE, shareY=TRUE, margin=0, heights=c(0.2, 0.8), widths=c(0.8, 0.2)) %>% plotly::layout(title=main, font=ft)

   }

}







