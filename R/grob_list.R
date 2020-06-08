#' Draw Each Spatial Heatmap and Convert Them to Ggplot2 Plot Grob
#'
#' @param gene The gene expession matrix, where rows are genes and columns are tissue/conditions.
#' @param geneV The gene expression values used to construct the colour bar.
#' @param con.na Logical, TRUE or FALSE. Default is TRUE, meaning conditions are available.
#' @param coord The coordidates extracted from the SVG file.
#' @param ID All gene ids selected after the App is launched.
#' @param cols All the color codes used to construct the color bar.
#' @param tis.path All the tissues/paths extracted from the SVG.
#' @param tis.trans A character vector of tissue/spatial feature identifiers. These tissues may cover other tissues and should be set transparent. \emph{E.g} c("brain", "heart").
#' @param sub.title.size A numeric. The subtitle font size of each individual spatial heatmap. Default is 11.
#' @param sam.legend "identical", "all", or a character vector of tissue names from the aSVG image to show in the legend plot. Default is "identical", meaning all the identical/matching tissues between the data matrix and aSVG image. If "all", all tissues in the aSVG image are shown.
#' @param legend.col A character vector of colors for the legend keys. The lenght must be equal to the number of target samples shown in the legend. 
#' @param legend.title A character, the legend title. Default is NULL.
#' @param legend.ncol An integer, the total columns of items in the legend. Default is NULL.
#' @param legend.nrow An integer, the total rows of the items in the legend. Default is NULL. 
#' @param legend.key.size A numeric (in "cm"). Default is 0.5. Size of the legend key.
#' @param legend.label.size A numeric. Default is 8. Size of the legend label.
#' @param legend.title.size A numeric. Default is 8. Size of the legend title.
#' @param line.size A numeric. The size of the shape outlines. Default is 0.2.
#' @param line.color A character. The color of shape outlines. Default is "grey70".
#' @param ... Other arguments passed to \code{\link[ggplot2]{ggplot}}.
#' @inheritParams ggplot2::theme

#' @return A list of spatial heatmaps in the form of ggplot2 plot grob.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @importFrom ggplot2 ggplot aes theme element_blank margin element_rect scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs guide_legend alpha coord_fixed

grob_list <- function(gene, con.na=TRUE, geneV, coord, ID, cols, tis.path, tis.trans=NULL, sub.title.size, sam.legend='identical', legend.col, legend.title=NULL, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.5, legend.label.size=8, legend.title.size=8, line.size=0.2, line.color='grey70', ...) {

  g_list <- function(con, lgd=FALSE, ...) {

    x <- y <- tissue <- NULL; tis.df <- unique(coord[, 'tissue'])
    if (lgd==FALSE) { 

      g.col <- NULL; con.idx <- grep(paste0("^", con, "$"), cons)
      tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]

      for (i in tis.path) {

        tis.idx <- which(tis.col1 %in% i); if (length(tis.idx)==1) { g.col <- c(g.col, scol1[tis.idx])
        } else if (length(tis.idx)==0) { g.col <- c(g.col, NA) }

      }; names(g.col) <- tis.df
    # Make selected tissues transparent by setting their colours as NA.
    if (!is.null(tis.trans)) for (i in tis.df) { if (sub('_\\d+$', '', i) %in% tis.trans) g.col[i] <- NA }
    
    } else g.col <- rep(NA, length(tis.path))
    names(g.col) <- tis.df # The colors might be internally re-ordered alphabetically during mapping, so give them names to fix the match with tissues. E.g. c('yellow', 'blue') can be re-ordered to c('blue', 'yellow'), which makes tissue mapping wrong. Correct: colours are not re-ordered. The 'tissue' in 'data=coord' are internally re-ordered according to a factor. Therfore, 'tissue' should be a factor with the right order. Otherwise, disordered mapping can happen.

    if (lgd==FALSE) scl.fil <- scale_fill_manual(values=g.col, guide=FALSE) else { 

      # Show selected or all samples in legend.
      if (length(sam.legend)==1) if (sam.legend=='identical') sam.legend <- intersect(sam.uni, unique(tis.path)) else if (sam.legend=='all') sam.legend <- unique(tis.path)
      # Select legend key colours if identical samples between SVG and matrix have colors of "none".
      legend.col <- legend.col[sam.legend] 
      if (any(legend.col=='none')) {
       
         n <- sum(legend.col=='none'); col.all <- grDevices::colors()[grep('honeydew|aliceblue|white|gr(a|e)y', grDevices::colors(), invert=TRUE)]
         col.none <- col.all[seq(from=1, to=length(col.all), by=floor(length(col.all)/n))]
         legend.col[legend.col=='none'] <- col.none[seq_len(n)]

       }
       # Map legend colours to tissues.
       sam.legend <- setdiff(sam.legend, tis.trans) 
       leg.idx <- !duplicated(tis.path) & (tis.path %in% sam.legend)
       legend.col <- legend.col[sam.legend] # Exclude transparent tissues. 
       for (i in seq_along(g.col)) {

         g.col0 <- legend.col[sub('_\\d+', '', names(g.col)[i])]
         if (!is.na(g.col0)) g.col[i] <- g.col0

       }; scl.fil <- scale_fill_manual(values=g.col, breaks=tis.df[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=legend.title, ncol=legend.ncol, nrow=legend.nrow)) 

    }

    g <- ggplot(...)+geom_polygon(data=coord, aes(x=x, y=y, fill=tissue), color=line.color, size=line.size, linetype='solid')+scl.fil+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=sub.title.size))+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))
    if (con.na==FALSE) g.tit <- ggtitle(k) else g.tit <- ggtitle(paste0(k, "_", con)); g <- g+g.tit

    if (lgd==TRUE) {

      g <- g+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.01, 0.01, 0.2, 0, "npc"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=15, face="bold"), legend.position=legend.position, legend.direction=legend.direction, legend.background = element_rect(fill=alpha(NA, 0)), legend.key.size=unit(legend.key.size, "cm"), legend.text=element_text(size=legend.label.size), legend.title=element_text(size=legend.title.size), legend.margin=margin(l=0.1, r=0.1, unit='cm'))+ggtitle('Legend')

    }; return(g)

  }
  # Map colours to samples according to expression level.
  cname <- colnames(gene); form <- grep('__', cname) # Only take the column names with "__".
  cons <- gsub("(.*)(__)(.*)", "\\3", cname[form]); con.uni <- unique(cons)
  sam.uni <- unique(gsub("(.*)(__)(.*)", "\\1", cname)); tis.trans <- make.names(tis.trans)
  grob.na <- grob.lis <- NULL; for (k in ID) {

    scol <- NULL; for (i in gene[k, ]) { 
      ab <- abs(i-geneV); col.ind <- which(ab==min(ab))[1]; scol <- c(scol, cols[col.ind])
    }

    idx <- grep("__", cname); c.na <- cname[idx]
    tis.col <- gsub("(.*)(__)(.*)", "\\1", c.na); g.lis <- NULL
    grob.na0 <- paste0(k, "_", con.uni); g.lis <- lapply(con.uni, g_list)
    # Repress popups by saving it to a png file, then delete it.
    tmp <- tempfile()
    png(tmp); grob <- lapply(g.lis, ggplotGrob); dev.off(); if (file.exists(tmp)) do.call(file.remove, list(tmp))
    names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob) 

  }; g.lgd <- g_list(con=NULL, lgd=TRUE)
  return(list(grob.lis=grob.lis, g.lgd=g.lgd))

}

