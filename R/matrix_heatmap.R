#' Matrix Heatmap
#'
#' This function represents the input gene in the context of corresponding gene network module as a web-browser based interactive matrix heatmap, where the rows and columns are sorted by hierarchical clustering dendrograms and the input gene is tagged by a red rectangle. To explore the results, users can zoom in and out by drawing a rectangle and by double clicking the image, respectively. Users can scale the expression values by gene or sample. \cr The network modules are identified at two alternative sensitivities levels (3, 2) by the function "adj.mod". From 3 to 2, the sensitivity decreases and results in less modules with larger sizes. The same module can also be displayed as an interactive network by the "network" function.

#' @param geneID A gene ID from the expression matrix. 
#' @param data The gene expression matrix, where rows are gene IDs and columns are samples/conditions. If annotation is included, it must be in the last column in parallel with samples/conditions.

#' @param adj.mod The "adjacency matrix" and "module" definition retured by the function "adj.mod".

#' @param ds The module identification sensitivity, either 2 or 3. See function "adj.mod" for details.

#' @param scale It specifies whether to scale the heatmap. There are three options: "row" means scale by row, "column" means scale by column, and "no" means no scale. 

#' @return An interactive matrix lauched on the web browser. 

#' @examples

#' data.path <- system.file("extdata/example", "root_expr_ann_row_gen.txt", package = "spatialHeatmap")
#' exp <- filter.data(data=data.path, sep="\t", isRowGen=TRUE, pOA=c(0, 0), CV=c(0.1, 10000), dir="./")
#' adj_mod <- adj.mod(data=exp, type="signed", minSize=15, dir="./")
#' # The gene "PSAC" is represented in the context of its network module in the form of an interactive matrix heatmap.
#' \donttest{ matrix.heatmap(geneID="PSAC", data=exp, adj.mod=adj_mod, ds="2", scale="row") }
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Andrie de Vries and Brian D. Ripley (2016). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1-20. https://CRAN.R-project.org/package=ggdendro \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Carson Sievert (2018) plotly for R. https://plotly-book.cpsievert.me \cr Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559 \cr Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. Journal of Statistical Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/ \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ \cr Carson Sievert (2018) plotly for R. https://plotly-book.cpsievert.me

#' @export
#' @importFrom ggdendro dendro_data
#' @importFrom ggplot2 ggplot geom_segment geom_text position_dodge geom_rect theme theme_minimal
#' @importFrom plotly plot_ly layout subplot %>%
#' @importFrom stats hclust order.dendrogram as.dendrogram


matrix.heatmap <- function(geneID, data, adj.mod, ds, scale) {

   x <- x1 <- x2 <- y <- y1 <- y2 <- xend <- yend <- NULL 
   adj <- adj.mod[["adj"]]; mods <- adj.mod[["mod"]]

   # The last column of data is expected to be annotation.
   col.len <- ncol(data)
   if (!is.numeric(data[, col.len])) gene <- data[, -col.len] else gene <- data
   lab <- mods[, ds][rownames(gene)==geneID]
    if (lab=="0") stop("The input gene is not assigned to any module. Please input a different one.")

      mod <- gene[mods[, ds]==lab, ]
      dd.gen <- as.dendrogram(hclust(dist(mod))); dd.sam <- as.dendrogram(hclust(dist(t(mod))))
      d.sam <- dendro_data(dd.sam); d.gen <- dendro_data(dd.gen)

      g.dengra <- function(df) {

        ggplot()+geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend))+labs(x = "", y = "") + theme_minimal()+ theme(axis.text = element_blank(), axis.ticks=element_blank(), panel.grid = element_blank())

      }

      p.gen <- g.dengra(d.gen$segments)+coord_flip(); p.sam <- g.dengra(d.sam$segments)
      df.gen <- data.frame(x=seq_len(length(labels(dd.gen))), y=0, lab=labels(dd.gen))
      gen.idx <- which(labels(dd.gen)==geneID)
      df.rec <- data.frame(x1=gen.idx-0.5, x2=gen.idx+0.5, y1=-1, y2=8)
      df.sam <- data.frame(x=seq_len(length(labels(dd.sam))), y=0, lab=labels(dd.sam)) 
      p.gen1 <- p.gen+geom_text(data=df.gen, aes(x=x, y=y, label=lab), position=position_dodge(0.9), vjust=0, hjust=-1, size=2, angle=0)+geom_rect(data=df.rec, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=NA, color="red", size=1, alpha=1)
      p.sam1 <- p.sam+geom_text(data=df.sam, aes(x=x, y=y, label=lab), 
      position=position_dodge(0.9), vjust=1, hjust=0, size=2, angle=90) 
      gen.ord <- order.dendrogram(dd.gen); sam.ord <- order.dendrogram(dd.sam)
      gene.clus <- rbind(Y=0, cbind(X=0, mod[gen.ord, sam.ord]))

      z <- as.matrix(gene.clus); if (scale=="column") z <- scale(z); if (scale=="row") z <- t(scale(t(z)))
      ply <- plot_ly(z=z, type="heatmap") %>% layout(yaxis=
      list(domain=c(0, 1), showticklabels=FALSE, showgrid=FALSE, ticks="", zeroline=FALSE), 
      xaxis=list(domain=c(0, 1), showticklabels=FALSE, showgrid=FALSE, ticks="", zeroline=FALSE))

      subplot(p.sam1, plot_ly(), ply, p.gen1, nrows=2, shareX=TRUE, shareY=TRUE, margin=0, heights=c(0.2, 0.8), widths=c(0.8, 0.2))

}


