#' Matrix Heatmap
#'
#' This function displays the input gene in the context of corresponding gene network module, where the rows and columns are sorted by hierarchical clustering dendrograms and the row of input gene is tagged by two lines. The matrix heatmap can be dispalyed in static or a web-browser based interactive mode. If the latter, users can zoom in and out by drawing a rectangle and by double clicking the image, respectively. Users can scale the data matrix by gene or sample. The same module can also be displayed in form of a network through the \code{\link{network}} function.

#' @param geneID A gene ID from the expression data matrix. 
#' @param data A "SummarizedExperiment" object containing the processed data matrix and metadata returned by the function \code{\link{filter_data}}. In the data matrix, rows are gene IDs and columns are samples/conditions.
#' @param adj.mod The list of "adjacency matrix" and "module" definition retured by the function \code{\link{adj_mod}}.
#' @param ds The module identification sensitivity ds, either 2 or 3. See function \code{\link{adj_mod}} for details.
#' @param scale "row", "column", or "no", meaing scale the data matrix by row, column, or no scale. 
#' @param col A character vector of two colours. It is used for constructing the colour scale. The default is c('yellow', 'blue').
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
#' @param arg.lis1 A list of additional arguments passed to the \code{\link[gplots]{heatmap.2}} function from "gplots" package. E.g. list(xlab='sample', ylab='gene'). The default is an empty list.
#' @param arg.lis2 A list of additional arguments passed to the \code{\link[ggplot2]{ggplot}} function from "ggplot2" package. The default is an empty list. 
#' @return A static image or an interactive application lauched on the web browser. 
#' @examples

#' # The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
#' library(ExpressionAtlas); library(SummarizedExperiment)
#' rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]; assay(rse.hum)[1:3, 1:3]
#'
#' # A targets file describing replicates of samples and conditions is required, which is made based on the "colData" slot in the downloaded "RangedSummarizedExperiment" and available in spatialHeatmap. See the "se" parameter for details. 
#' brain.pa <- system.file('extdata/shinyApp/example/target_brain.txt', package='spatialHeatmap')
#' target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
#' # The "organism_part" and "disease" column describes tissue and condition replicates respectively.  
#' target.hum[c(1:3, 41:42), 4:5]
#' # Place the targets file into "colData" slot as a DataFrame class. 
#' colData(rse.hum) <- DataFrame(target.hum)
#' 
#' # For users with little R expertise, if the gene expression matrix comes as a data frame, it should be placed into "SummarizedExperiment" before proceeding to next step. An example is shown below by borrowing a data matrix from the brain data.
#' # Borrow a data matrix.
#' df <- assay(rse.hum); df[1:2, 1:3]
#' # Place the data matrix and targets file (target.hum) into "SummarizedExperiment".
#' rse.hum <- SummarizedExperiment(assay=df, colData=target.hum, rowData=NULL)
#' 
#' # The count matrix is normalised with estimateSizeFactors (type=‘ratio’).
#' se.nor.hum <- norm_data(data=rse.hum, method.norm='CNF', data.trans='log2')
#'
#' # Average replicates of concatenated sample__condition.
#' se.aggr.hum <- aggr_rep(data=se.nor.hum, sam.factor='organism_part', con.factor='disease', aggr='mean')
#' assay(se.aggr.hum)[49939:49942, ] # The concatenated tissue__conditions are the column names of the output data matrix.
#' 
#' # Genes with low expression level and low variantion are always filtered. 
#' se.fil.hum <- filter_data(data=se.aggr.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.01, 5), CV=c(0.3, 100), dir=NULL)

#' # Detect modules. 
#' adj.mod <- adj_mod(data=se.fil.hum, type="signed", minSize=15, dir=NULL)
#' # The first column is ds=2 while the second is ds=3. The numbers in each column are module labels with "0" meaning genes not assigned to any modules.
#' adj.mod[['mod']][1:3, ]

#' # Plot matrix heatmap on gene ENSG00000008196 with ds='3'. Set "static=TRUE" to launch the interactive mode. 
#' matrix_hm(geneID="ENSG00000008196", data=se.fil.hum, adj.mod=adj.mod, ds="3", scale="no", angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(10, 6), static=TRUE, arg.lis1=list(offsetRow=0.1, offsetCol=0.1))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr Andrie de Vries and Brian D. Ripley (2016). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1-20. https://CRAN.R-project.org/package=ggdendro \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Carson Sievert (2018) plotly for R. https://plotly-book.cpsievert.me \cr Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ \cr Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2019). gplots: Various R Programming Tools for Plotting Data. R package version 3.0.1.1.  https://CRAN.R-project.org/package=gplots \cr Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/ 
#' Prudencio, Mercedes, Veronique V Belzil, Ranjan Batra, Christian A Ross, Tania F Gendron, Luc J Pregent, Melissa E Murray, et al. 2015. "Distinct Brain Transcriptome Profiles in C9orf72-Associated and Sporadic ALS. Nat. Neurosci. 18 (8): 1175–82
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8


#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom ggdendro dendro_data
#' @importFrom ggplot2 ggplot geom_segment geom_text position_dodge geom_rect theme theme_minimal geom_tile scale_fill_gradient geom_hline
#' @importFrom plotly layout subplot %>%
#' @importFrom stats hclust order.dendrogram as.dendrogram
#' @importFrom gplots heatmap.2 
#' @importFrom graphics image mtext par plot title
#' @importFrom grDevices dev.off png

matrix_hm <- function(geneID, data, adj.mod, ds, scale, col=c('yellow', 'blue'), main=NULL, title.size=10, cexCol=1, cexRow=1, angleCol=45, angleRow=45, sepcolor="black", sep.width=0.02, static=TRUE, margin=c(10, 10), arg.lis1=list(), arg.lis2=list()) {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- colnames(data) 
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(se)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    gene <- na[!idx]; colnames(gene) <- cna[!idx]

  } else if (is(data, 'SummarizedExperiment')) { gene <- assay(data) }
  mods <- adj.mod[["mod"]]; ds <- as.character(ds)
  lab <- mods[, ds][rownames(gene)==geneID]
  if (lab=="0") stop("The input gene is not assigned to any module. Please input a different one.")
  mod <- as.matrix(gene[mods[, ds]==lab, ])
  
  if (static==TRUE) {

    tmp <- tempdir(); pa <- paste0(tmp, '/delete_hm.png')
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






