#' Draw each spatial heatmap and convert them to ggplot2 plot grob
#'

#' @param gene The gene expession matrix, where rows are genes and columns are tissue/conditions.
#' @param geneV The gene expression values used to construct the colour bar.
#' @param coord The coordidates extracted from the SVG file.
#' @param ID All gene ids selected after the App is launched.
#' @param cols All the colour codes used to construct the colour bar.
#' @param tis.path All the tissues/paths extracted from the SVG.
#' @param tis.trans A character vector of tissue names. These tissues cover other tissues and should be set transparent. E.g c("epidermis", "cortex").
#' @param sub.title.size A numeric. The subtitle font size of each individual spatial heatmap. Default is 11.
#' @param sam.legend A character vector of tissue names in the SVG image. These tissues are shown in the legend. Default is "identical", meaning all the identical tissues between the data matrix and SVG image. Another value is 'all', meaning all tissues in the SVG image.
#' @param title.legend A character, the legend title. Default is NULL.
#' @param ncol.legend An integer, the column number of the items in the legend. Default is NULL.
#' @param nrow.legend An integer, the row number of the items in the legend. Default is NULL. 
#' @param pos.legend Legend position. One of "none", "top", "right", "bottom", "left". "None" meaning remove the legend. Default is "right".
#' @param legend.key.size A numeric (in "cm"). Default is 0.5. Size of the legend key.
#' @param legend.label.size A numeric. Default is 8. Size of the legend label.
#' @param legend.title.size A numeric. Default is 8. Size of the legend title.
#' @param line.size A numeric. The size of the polygon outline. Default is 0.2.
#' @param line.color A character. The color of polygon outline. Default is "grey70".
#' @param ... Other arguments passed to ggplot().



#' @return A list of spatial heatmaps in the form of ggplot2 plot grob.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @importFrom ggplot2 ggplot aes theme element_blank margin element_rect scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs guide_legend

grob_list <- function(gene, geneV, coord, ID, cols, tis.path, tis.trans=NULL, sub.title.size, sam.legend='identical', title.legend=NULL, ncol.legend=NULL, nrow.legend=NULL, pos.legend='right', legend.key.size=0.5, legend.label.size=8, legend.title.size=8, line.size=0.2, line.color='grey70', ...) {
  
  x <- y <- tissue <- NULL
  # Map colours to samples according to expression level.
  g.list <- function(j) {

    g.col <- NULL; con.idx <- grep(paste0("^", j, "$"), con)
    tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]

    for (i in tis.path) {

      tis.idx <- which(tis.col1 %in% i); if (length(tis.idx)==1) { g.col <- c(g.col, scol1[tis.idx])
      } else if (length(tis.idx)==0) { g.col <- c(g.col, NA) }

    }
    names(g.col) <- tis.df <- unique(coord[, 'tissue']) # The colors might be internally re-ordered alphabetically during mapping, so give them names to fix the match with tissues. E.g. c('yellow', 'blue') can be re-ordered to c('blue', 'yellow'), which makes tissue mapping wrong. Correct: colours are not re-ordered. The 'tissue' in 'data=coord' are internally re-ordered according to a factor. Therfore, 'tissue' should be a factor with the right order. Otherwise, disordered mapping can happen.
    # Make selected tissues transparent by setting their colours as NA.
    if (!is.null(tis.trans)) for (i in tis.df) { if (sub('_\\d+$', '', i) %in% tis.trans) g.col[i] <- NA }
    # Show selected or all samples in legend.
    if (length(sam.legend)==1) if (sam.legend=='identical') sam.legend <- unique(tis.path[!is.na(g.col)]) else if (sam.legend=='all') sam.legend <- unique(tis.path)
    leg.idx <- !duplicated(tis.path) & (tis.path %in% sam.legend)
    g <- ggplot()+geom_polygon(data=coord, aes(x=x, y=y, fill=tissue), color=line.color, size=line.size, linetype='solid')+scale_fill_manual(values=g.col, breaks=tis.df[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=title.legend, ncol=ncol.legend, nrow=nrow.legend))+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=sub.title.size), legend.position=pos.legend, legend.key.size=unit(legend.key.size, "cm"), legend.text=element_text(size=legend.label.size), legend.title=element_text(size=legend.title.size))+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))+ggtitle(paste0(k, "_", j)); return(g)


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
    tmp <- tempdir(); pa <- paste0(tmp, '/delete.png')
    png(pa); grob <- lapply(g.lis, ggplotGrob); dev.off(); if (file.exists(pa)) do.call(file.remove, list(pa))
    names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob) 

  }; return(grob.lis)

}
