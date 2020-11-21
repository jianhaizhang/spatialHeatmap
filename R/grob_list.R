#' Draw Each Spatial Heatmap and Convert Them to Ggplot2 Plot Grob
#'
#' @param gene The gene expession matrix, where rows are genes and columns are tissue/conditions.
#' @param geneV The gene expression values used to construct the colour bar.
#' @param con.na Logical, TRUE or FALSE. Default is TRUE, meaning conditions are available.
#' @param coord The coordidates extracted from the SVG file.
#' @param ID All gene ids selected after the App is launched.
#' @param cols All the color codes used to construct the color bar.
#' @param tis.path All the tissues/paths extracted from the SVG.
#' @param ft.trans A character vector of tissue/spatial feature identifiers that will be set transparent. \emph{E.g} c("brain", "heart"). This argument is used when target features are covered by  overlapping features and the latter should be transparent.
#' @param sub.title.size A numeric of the subtitle font size of each individual spatial heatmap. The default is 11.
#' @param ft.legend One of "identical", "all", or a character vector of tissue/spatial feature identifiers from the aSVG file. The default is "identical" and all the identical/matching tissues/spatial features between the data and aSVG file are indicated in the legend plot. If "all", all tissues/spatial features in the aSVG are shown. If a vector, only the tissues/spatial features in the vector are shown.
#' @param legend.col A character vector of colors for the keys in the legend plot. The lenght must be equal to the number of target samples shown in the legend. 
#' @param legend.ncol An integer of the total columns of keys in the legend plot. The default is NULL. If both \code{legend.ncol} and \code{legend.nrow} are used, the product of the two arguments should be equal or larger than the total number of shown spatial features.
#' @param legend.nrow An integer of the total rows of keys in the legend plot. The default is NULL. It is only applicable to the legend plot. If both \code{legend.ncol} and \code{legend.nrow} are used, the product of the two arguments should be equal or larger than the total number of matching spatial features.
#' @param legend.key.size A numeric of the legend key size ("npc"), applicable to the legend plot. The default is 0.02. 
#' @param legend.text.size A numeric of the legend label size, applicable to the legend plot. The default is 12.
#' @param line.size The thickness of each shape outline in the aSVG is maintained in spatial heatmaps, \emph{i.e.} the stroke widths in Inkscape. This argument is the extra thickness added to all outlines. Default is 0.2 in case stroke widths in the aSVG are 0.
#' @param line.color A character of the shape outline color. Default is "grey70".
#' @param mar.lb A two-component numeric vector. The first and second numeric is left/right and bottom/top margin (npc) respectively.
#' @param legend.plot.title The title of the legend plot. The default is NULL. 
#' @param legend.plot.title.size The title size of the legend plot. The default is 11.
#' @param ... Other arguments passed to \code{\link[ggplot2]{ggplot}}.
#' @inheritParams ggplot2::theme

#' @return A nested list of spatial heatmaps of ggplot2 plot grob, spatial heatmaps of ggplot, and legend plot of ggplot.
#' @keywords Internal
#' @noRd


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @importFrom ggplot2 ggplot aes theme element_blank margin element_rect scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs guide_legend alpha coord_fixed

grob_list <- function(gene, con.na=TRUE, geneV, coord, ID, cols, tis.path, ft.trans=NULL, sub.title.size, ft.legend='identical', legend.col, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, legend.plot.title=NULL, legend.plot.title.size=11, line.size=0.2, line.color='grey70', mar.lb=NULL, ...) {

  # save(gene, con.na, geneV, coord, ID, cols, tis.path, ft.trans, sub.title.size, ft.legend, legend.col, legend.ncol, legend.nrow, legend.position, legend.direction, legend.key.size, legend.text.size, legend.plot.title, legend.plot.title.size, line.size, line.color, mar.lb, file='all')
  
  # Main function to create SHMs and legend plot
  g_list <- function(con, lgd=FALSE, ...) {

    value <- feature <- x <- y <- tissue <- NULL; tis.df <- as.vector(unique(coord[, 'tissue']))
    # tis.path and tis.df have the same length by default.
    # Assign default colours to each path.
    g.col <- rep(NA, length(tis.path)); names(g.col) <- tis.df
    for (i in seq_along(g.col)) {
      g.col0 <- legend.col[sub('_\\d+$', '', names(g.col)[i])]
      if (g.col0=='none') next
      if (!is.na(g.col0)) g.col[i] <- g.col0
    }
    if (lgd==FALSE) {
      con.idx <- grep(paste0("^", con, "$"), cons)
      # Target tissues and colors in data columns.
      tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]
      for (i in unique(tis.path)) {
        # Map target colors to target tissues.
        tis.idx <- which(tis.col1 %in% i); if (length(tis.idx)==1) {
          pat <- paste0(paste0('^', i, '$'), '|', paste0('^', i, '_\\d+$'))
          g.col[grep(pat, tis.df)] <- scol1[tis.idx] # names(g.col) is tis.df 
        }
      }
    } # The colors might be internally re-ordered alphabetically during mapping, so give them names to fix the match with tissues. E.g. c('yellow', 'blue') can be re-ordered to c('blue', 'yellow'), which makes tissue mapping wrong. Correct: colours are not re-ordered. The 'tissue' in 'data=coord' are internally re-ordered according to a factor. Therfore, 'tissue' should be a factor with the right order. Otherwise, disordered mapping can happen.
    # Make selected tissues transparent by setting their colours as NA.
    if (!is.null(ft.trans)) for (i in tis.df) { if (sub('_\\d+$', '', i) %in% ft.trans) g.col[i] <- NA }
    # Show selected or all samples in legend.
    if (length(ft.legend)==1) if (ft.legend=='identical') ft.legend <- intersect(sam.uni, unique(tis.path)) else if (ft.legend=='all') ft.legend <- unique(tis.path)
    
    if (lgd==FALSE) { # Legend plot.
    
      ft.legend <- setdiff(ft.legend, ft.trans) 
      leg.idx <- !duplicated(tis.path) & (tis.path %in% ft.legend)
      # Bottom legends are set for each SHM and then removed in 'ggplotGrob', but a copy with legend is saved separately for later used in video.
      scl.fil <- scale_fill_manual(values=g.col, breaks=tis.df[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=NULL, ncol=legend.ncol, nrow=legend.nrow))
   
    } else { 

      # Assign legend key colours if identical samples between SVG and matrix have colors of "none".
      legend.col1 <- legend.col[ft.legend] 
      if (any(legend.col1=='none')) {
       
         n <- sum(legend.col1=='none'); col.all <- grDevices::colors()[grep('honeydew|aliceblue|white|gr(a|e)y', grDevices::colors(), invert=TRUE)]
         col.none <- col.all[seq(from=1, to=length(col.all), by=floor(length(col.all)/n))]
         legend.col1[legend.col1=='none'] <- col.none[seq_len(n)]

       }
       # Map legend colours to tissues.
       ft.legend <- setdiff(ft.legend, ft.trans) 
       leg.idx <- !duplicated(tis.path) & (tis.path %in% ft.legend)
       legend.col1 <- legend.col1[ft.legend] # Exclude transparent tissues. 
       # Copy colors across same numbered tissues.
       for (i in seq_along(g.col)) {
         if (!is.na(g.col[i])) next
         g.col0 <- legend.col1[sub('_\\d+$', '', names(g.col)[i])]
         if (!is.na(g.col0)) g.col[i] <- g.col0
       }; scl.fil <- scale_fill_manual(values=g.col, breaks=tis.df[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=NULL, ncol=legend.ncol, nrow=legend.nrow)) 

    }
    lgd.par <- theme(legend.position=legend.position, legend.direction=legend.direction, legend.background = element_rect(fill=alpha(NA, 0)), legend.key.size=unit(legend.key.size, "npc"), legend.text=element_text(size=legend.text.size), legend.margin=margin(l=0.1, r=0.1, unit='npc'))
    ## Add 'feature' and 'value' to coordinate data frame, since the resulting ggplot object is used in 'ggplotly'. Otherwise, the coordinate data frame is applied to 'ggplot' directly by skipping the following code.
    coord$gene <- k; coord$condition <- con; coord$value <- NA
    ft.pat <- paste0('(', paste0(unique(tis.path), collapse='|'), ')(_\\d+)$')
    coord$feature <- gsub(ft.pat, '\\1', coord$tissue)    
    # Assign values to each tissue.
    col.na <- paste0(coord$feature, '__', coord$condition)
    idx1 <- col.na %in% colnames(gene); df0 <- coord[idx1, ]
    df0$value <- unlist(gene[df0$gene[1], col.na[idx1]])
    coord[idx1, ] <- df0; coord <- line_size(coord, line.size)
    # If "data" is not in ggplot(), g$data slot is empty. x, y, and group should be in the same aes(). 
    g <- ggplot(data=coord, aes(x=x, y=y, value=value, group=tissue, text=paste0('feature: ', feature, '\n', 'value: ', value)), ...)+geom_polygon(aes(fill=tissue), color=line.color, size=coord$line.size, linetype='solid')+scl.fil+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), axis.title.x=element_text(size=16, face="bold"), plot.title=element_text(hjust=0.5, size=sub.title.size), legend.box.margin=margin(-20, 0, 2, 0, unit='pt'))+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))+lgd.par
    if (is.null(mar.lb)) g <- g+theme(plot.margin=margin(0.005, 0.005, 0.005, 0.005, "npc")) else g <- g+theme(plot.margin=margin(mar.lb[2], mar.lb[1], mar.lb[2], mar.lb[1], "npc"))
    if (con.na==FALSE) g.tit <- ggtitle(k) else g.tit <- ggtitle(paste0(k, "_", con)); g <- g+g.tit

    if (lgd==TRUE) {

      g <- g+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.005, 0.005, 0.2, 0, "npc"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=legend.plot.title.size))+lgd.par+ggtitle(legend.plot.title)

    }; return(g)

  }

  # Map colours to samples according to expression level.
  cname <- colnames(gene); form <- grep('__', cname) # Only take the column names with "__".
  cons <- gsub("(.*)(__)(.*)", "\\3", cname[form]); con.uni <- unique(cons)
  sam.uni <- unique(gsub("(.*)(__)(.*)", "\\1", cname)); ft.trans <- make.names(ft.trans)
  grob.na <- grob.lis <- g.lis.all <- NULL; for (k in ID) {

    scol <- NULL; for (i in gene[k, ]) { 
      ab <- abs(i-geneV); col.ind <- which(ab==min(ab))[1]; scol <- c(scol, cols[col.ind])
    }

    idx <- grep("__", cname); c.na <- cname[idx]
    tis.col <- gsub("(.*)(__)(.*)", "\\1", c.na) 
    tis.col.uni <- unique(tis.col); tis.path.uni <- unique(tis.path)
    tis.tar <- tis.col.uni[tis.col.uni %in% tis.path.uni]
    if (length(tis.tar)==0) return(NULL)
    c.na1 <- c.na[grepl(paste0('^(', paste0(tis.tar,collapse='|'), ')__'), c.na)]
    # Only conditions paired with valid tissues (have matching samples in data) are used. 
    con.vld <- gsub("(.*)(__)(.*)", "\\3", c.na1); con.vld.uni <- unique(con.vld)
    g.lis <- NULL; grob.na0 <- paste0(k, "_", con.vld.uni); g.lis <- lapply(con.vld.uni, g_list, ...)
    # Repress popups by saving it to a png file, then delete it.
    tmp <- normalizePath(tempfile(), winslash='/', mustWork=FALSE)
    png(tmp); grob <- lapply(g.lis, function(x) { x <- x+theme(legend.position="none"); ggplotGrob(x) })
    dev.off(); if (file.exists(tmp)) do.call(file.remove, list(tmp))
    names(g.lis) <- names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob); g.lis.all <- c(g.lis.all, g.lis)

  }; g.lgd <- g_list(con=NULL, lgd=TRUE, ...)
  return(list(grob.lis=grob.lis, g.lgd=g.lgd, g.lis.all=g.lis.all))

}


#' # Assign line size to each tissue.
#' @keywords Internal
#' @noRd
line_size <- function(coord, line.size) {
  coord$line.size <- 0; for (i in names(line.size)) {
    pat <- paste0('^', i, c('_\\d+', ''), '$', collapse='|')
    coord$line.size[grepl(pat, coord$tissue)] <- line.size[i]
  }; return(coord)
}




