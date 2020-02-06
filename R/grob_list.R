#' Draw each spatial heatmap and convert them to ggplot2 plot grob
#'

#' @param gene The gene expession matrix, where rows are genes and columns are tissue/conditions.
#' @param geneV The gene expression values used to construct the colour bar.
#' @param coord The coordidates extracted from the SVG file.
#' @param ID All gene ids selected after the App is launched.
#' @param cols All the colour codes used to construct the colour bar.
#' @param tis.path All the tissues/paths extracted from the SVG.
#' @param tis.trans The tissues selected by users that are expected to be transparent.
#' @param sub.title.size The subtitle font size of each individual spatial heatmap.
#' @param ... Other arguments passed to ggplot().

#' @return A list of spatial heatmaps in the form of ggplot2 plot grob.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @importFrom ggplot2 ggplot aes theme element_blank margin element_rect scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs 

grob_list <- function(gene, geneV, coord, ID, cols, tis.path, tis.trans=NULL, sub.title.size, ...) {

  x <- y <- tissue <- NULL
  # Map colours to samples according to expression level.
  g.list <- function(j) {

    g.col <- NULL; con.idx <- grep(paste0("^", j, "$"), con)
    tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]

    for (i in tis.path) {

      tis.idx <- which(tis.col1 %in% i); if (length(tis.idx)==1) { g.col <- c(g.col, scol1[tis.idx])
      } else if (length(tis.idx)==0) { g.col <- c(g.col, "white") }

    }
    names(g.col) <- tis.df <- unique(coord[, 'tissue']) # The colors might be internally re-ordered alphabetically during mapping, so give them names to fix the match with tissues. E.g. c('yellow', 'blue') can be re-ordered to c('blue', 'yellow'), which makes tissue mapping wrong. Correct: colours are not re-ordered. The 'tissue' in 'data=coord' are internally re-ordered according to a factor. Therfore, 'tissue' should be a factor with the right order. Otherwise, disordered mapping can happen.
    # Make selected tissues transparent by setting their colours as NA.
    if (!is.null(tis.trans)) for (i in tis.df) { if (sub('_\\d+$', '', i) %in% tis.trans) g.col[i] <- NA }
    g <- ggplot(...)+geom_polygon(data=coord, aes(x=x, y=y, fill=tissue), color="black")+scale_fill_manual(values=g.col, guide=FALSE)+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=sub.title.size))+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))+ggtitle(paste0(k, "_", j)); return(g)

  }
  cname <- colnames(gene); con <- gsub("(.*)(__)(.*)", "\\3", cname); con.uni <- unique(con)
  grob.na <- grob.lis <- NULL; for (k in ID) {

    scol <- NULL; for (i in gene[k, ]) { 
      ab <- abs(i-geneV); col.ind <- which(ab==min(ab))[1]; scol <- c(scol, cols[col.ind])
    }

    idx <- grep("__", cname); c.na <- cname[idx]
    tis.col <- gsub("(.*)(__)(.*)", "\\1", c.na); g.lis <- NULL
    grob.na0 <- paste0(k, "_", con.uni); g.lis <- lapply(con.uni, g.list)
    # Repress popups by saving it to a png file, then delete it.
    tmp <- system.file("extdata/shinyApp/tmp", package = "spatialHeatmap"); pa <- paste0(tmp, '/delete.png')
    png(pa); grob <- lapply(g.lis, ggplotGrob); dev.off(); if (file.exists(pa)) do.call(file.remove, list(pa))
    names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob) 

  }; return(grob.lis)

}
