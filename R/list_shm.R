#' Draw Each Spatial Heatmap and Legend Plot
#'
#' @param gene The gene expession matrix, where rows are genes and columns are tissue/conditions.
#' @param geneV The gene expression values used to construct the colour bar.
#' @param con.na Logical, TRUE or FALSE. Default is TRUE, meaning conditions are available.
#' @param coord The coordidates extracted from the SVG file.
#' @param ID All gene ids selected after the App is launched.
#' @param cols All the color codes used to construct the color bar.
#' @param tis.path All the tissues/paths extracted from the SVG.
#' @param lis.rematch A list for rematching features. In each slot, the slot name is an existing feature in the data, and the slot contains a vector of features in aSVG that will be rematched to the feature in the slot name. \emph{E.g.} \code{list(featureData1 = c('featureSVG1', 'featureSVG2'), featureData2 = c('featureSVG3'))}, where features \code{c('featureSVG1', 'featureSVG2')}, \code{c('featureSVG3')} in the aSVG are rematched to features \code{'featureData1'}, \code{'featureData2'} in data, respectively. 
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
#' @param aspect.ratio The ratio of width to height.
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
#' @importFrom parallel detectCores mclapply

gg_shm <- function(gene, con.na=TRUE, geneV, coord, ID, cols, tis.path, lis.rematch = NULL, ft.trans=NULL, sub.title.size, ft.legend='identical', legend.col, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, legend.plot.title=NULL, legend.plot.title.size=11, line.size=0.2, line.color='grey70', aspect.ratio = 1, ...) {

   # save(gene, con.na, geneV, coord, ID, cols, tis.path, lis.rematch, ft.trans, sub.title.size, ft.legend, legend.col, legend.ncol, legend.nrow, legend.position, legend.direction, legend.key.size, legend.text.size, legend.plot.title, legend.plot.title.size, line.size, line.color, aspect.ratio, file='all.shm')
  
  # Main function to create SHMs (by conditions) and legend plot.
  g_list <- function(con, lgd=FALSE, ...) {
    if (is.null(con)) cat('Legend plot ... \n') else cat(con, ' ')
    value <- feature <- x <- y <- tissue <- NULL; tis.df <- as.vector(unique(coord[, 'tissue']))
    legend.col[legend.col == 'none'] <- NA
    # tis.path and tis.df have the same length by default, but not entries, since tis.df is appended '__\\d+' at the end.
    # Assign default colours to each path.
    g.col <- rep(NA, length(tis.path)); names(g.col) <- tis.df
    if (lgd==FALSE) {
      # Keep global text/legend colors in the main SHM. No change on the local legend colors so they are all NA. The line widths of local legend will be set 0.
      g.col <- legend.col[grep('_globalLGD$', names(legend.col))][sub('__\\d+$', '', names(g.col))]

      names(g.col) <- tis.df # Resolves legend.col['tissue'] is NA by default. 
      # Un-related aSVG mapped to data.
      if (rematch.dif.svg) {
        # The data column index of data features that are assinged new aSVG features under a specific condition.
        tis.tar.idx <- cname %in% paste0(tis.tar, '__', con)
        # In data columns, feature__conditions and colors have the same index. 
        tis.tar.col <- color.dat[tis.tar.idx]; names(tis.tar.col) <- tis.tar
        # Copy colors of target features in data to corresponding aSVG features according to lis.rematch. 
        col.ft.rematch <- NULL; for (i in names(lis.rematch[tis.tar])) {
          ft.svg <- lis.rematch[tis.tar][[i]]
          col0 <- rep(tis.tar.col[i], length(ft.svg)); names(col0) <- ft.svg
          col.ft.rematch <- c(col.ft.rematch, col0)
        }
        for (i in names(col.ft.rematch)) g.col[tis.path %in% i] <- col.ft.rematch[i]
       } else {
        con.idx <- grep(paste0("^", con, "$"), cons)
        # Target tissues and colors in data columns.
        tis.col1 <- tis.col[con.idx]; color.dat1 <- color.dat[con.idx]
        for (i in unique(tis.path)) {
        # Map target colors to target tissues.
        tis.idx <- which(tis.col1 %in% i); if (length(tis.idx)==1) {
          # Account for single-shape tissue without '__\\d+$' and multi-shape tissue with '__\\d+$'.
          pat <- paste0(paste0('^', i, '$'), '|', paste0('^', i, '__\\d+$'))
          g.col[grep(pat, tis.df)] <- color.dat1[tis.idx] # names(g.col) is tis.df 
        }
      }
      # Rematch features between the same pair of data-aSVG.
      if (is.list(lis.rematch)) {
        for (i in seq_along(lis.rematch)) {
          lis0 <- lis.rematch[i]; if (!is.character(lis0[[1]])) next
          # Index of features will be rematched.
          idx.tis.rematch <- tis.path %in% lis0[[1]]
          # Index of color for rematching.
          idx.tis.rematch.color <- tis.path %in% names(lis0)
          # if (sum(idx.tis.rematch.color) == 0) stop(paste0('Feature "', names(lis0), '" is not detected in aSVG!'))
          if (sum(idx.tis.rematch) == 0 | sum(idx.tis.rematch.color) == 0) next
          g.col[idx.tis.rematch] <- unique(g.col[idx.tis.rematch.color])
         }
       }
     }
    } 
    # The colors might be internally re-ordered alphabetically during mapping, so give them names to fix the match with tissues. E.g. c('yellow', 'blue') can be re-ordered to c('blue', 'yellow'), which makes tissue mapping wrong. Correct: colours are not re-ordered. The 'tissue' in 'data=coord' are internally re-ordered according to a factor. Therfore, 'tissue' should be a factor with the right order. Otherwise, disordered mapping can happen. Alternatively, name the colors with corresponding tissue names.
    # aes() is passed to either ggplot() or specific layer. Aesthetics supplied to ggplot() are used as defaults for every layer. 
    # Show selected or all samples in legend.
    if (length(ft.legend)==1) if (ft.legend=='identical') {
      if (rematch.dif.svg) {
        ft.legend <- intersect(unlist(lis.rematch), unique(tis.path)) 
      } else ft.legend <- intersect(sam.uni, unique(tis.path)) 
    } else if (ft.legend=='all') ft.legend <- unique(tis.path)
    
    if (lgd==FALSE) { # Legend plot.
      # Make selected tissues transparent by setting their colours NA.
      if (!is.null(ft.trans)) g.col[sub('__\\d+$', '', tis.df) %in% ft.trans] <- NA # This step should not be merged with 'lgd=T'.
      ft.legend <- setdiff(ft.legend, ft.trans) 
      leg.idx <- !duplicated(tis.path) & (tis.path %in% ft.legend)
      # Bottom legends are set for each SHM and then removed in 'ggplotGrob', but a copy with legend is saved separately for later used in video.
      scl.fil <- scale_fill_manual(values=g.col, breaks=tis.df[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=NULL, ncol=legend.ncol, nrow=legend.nrow))
    } else { 
      # Assign legend key colours if identical samples between SVG and matrix have colors of "none".
      legend.col1 <- legend.col[ft.legend] # Only includes matching samples. 
      if (any(is.na(legend.col1))) {
         n <- sum(is.na(legend.col1)); col.all <- grDevices::colors()[grep('honeydew|aliceblue|white|gr(a|e)y', grDevices::colors(), invert=TRUE)]
         col.na <- col.all[seq(from=1, to=length(col.all), by=floor(length(col.all)/n))]
         legend.col1[is.na(legend.col1)] <- col.na[seq_len(n)]
       }
       # Map legend colours to tissues.
       # Exclude transparent tissues.
       ft.legend <- setdiff(ft.legend, ft.trans) 
       leg.idx <- !duplicated(tis.path) & (tis.path %in% ft.legend)
       legend.col1 <- legend.col1[ft.legend]
       # Keep all colors in the original SVG.
       g.col <- legend.col[sub('__\\d+$', '', names(g.col))]
       names(g.col) <- tis.df # Resolves "legend.col['tissue'] is NA" by default.
       # g.col[g.col=='none'] <- NA 
       # Make selected tissues transparent by setting their colours NA.
       if (!is.null(ft.trans)) g.col[sub('__\\d+$', '', tis.df) %in% ft.trans] <- NA
       # Copy colors across same numbered tissues. 
       g.col <- lapply(seq_along(g.col), function(x) {
         # In lapply each run must return sth.
         if (!is.na(g.col[x])) return(g.col[x])
         g.col0 <- legend.col1[sub('__\\d+$', '', names(g.col[x]))]
         if (!is.na(g.col0)) g.col[x] <- g.col0
         return(g.col[x])
         } 
       ); g.col <- unlist(g.col)
       # No matter the tissues in coordinate data frame are vector or factor, the coloring are decided by the named color vector (order of colors does not matter as long as names are right) in scale_fill_manual.
       scl.fil <- scale_fill_manual(values=g.col, breaks=tis.df[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=NULL, ncol=legend.ncol, nrow=legend.nrow)) 
    }
    lgd.par <- theme(legend.position=legend.position, legend.direction=legend.direction, legend.background = element_rect(fill=alpha(NA, 0)), legend.key.size=unit(legend.key.size, "npc"), legend.text=element_text(size=legend.text.size), legend.margin=margin(l=0.1, r=0.1, unit='npc'))
    ## Add 'feature' and 'value' to coordinate data frame, since the resulting ggplot object is used in 'ggplotly'. Otherwise, the coordinate data frame is applied to 'ggplot' directly by skipping the following code.
    coord$gene <- k; coord$condition <- con; coord$value <- NA
    coord$feature <- sub('__\\d+$', '', coord$tissue)
    # Assign values to each tissue.
    col.na <- paste0(coord$feature, '__', coord$condition)
    idx1 <- col.na %in% colnames(gene); df0 <- coord[idx1, ]
    # df0$value <- unlist(gene[df0$gene[1], col.na[idx1]]) # Slow in large data.
    col.na1 <- unique(col.na[idx1])
    tab0 <- table(col.na[idx1])[col.na1]
    lis.v <- lapply(col.na1, function(i) { rep(gene[df0$gene[1], i], tab0[i][[1]]) } )
    df0$value <- unlist(lis.v) 
    coord[idx1, ] <- df0; coord$line.size <- coord$line.size+line.size
    if (lgd==FALSE) { # Set line size as 0 for local legends. 
      idx.local <- grepl('_localLGD$|_localLGD__\\d+$', coord$tissue)
      coord$line.size[idx.local] <- 0
    }
    # If "data" is not in ggplot(), g$data slot is empty. x, y, and group should be in the same aes().
    g <- ggplot(data=coord, aes(x=x, y=y, value=value, group=tissue, text=paste0('feature: ', feature, '\n', 'value: ', value)), ...)+geom_polygon(aes(fill=tissue), color=line.color, size=coord$line.size, linetype='solid')+scl.fil+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.title=element_text(hjust=0.5, size=sub.title.size), legend.box.margin=margin(-20, 0, 2, 0, unit='pt'), plot.margin=margin(0.005, 0.005, 0.005, 0.005, "npc"), aspect.ratio = 1/aspect.ratio)+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))+lgd.par
    # The aspect ratio should not be calculated by margins that are inferred from original width/height (mar.lb <- (1-w.h/w.h.max*0.99)/2). It works on single plot, but will squeeze the subplots in arrangeGrob/grid.arrange. Rather the "aspect ratio" in theme should be used, which will not squeeze subplots.
    # After theme(aspect.ratio) is set, change in only top or left margin will not distort the image. Rather width/height will scale propotionally, since the aspect ratio is fixed.
    # if (is.null(mar.lb)) g <- g+theme(plot.margin=margin(0.005, 0.005, 0.005, 0.005, "npc")) else g <- g+theme(plot.margin=margin(mar.lb[2], mar.lb[1], mar.lb[2], mar.lb[1], "npc"))
    if (con.na==FALSE) g.tit <- ggtitle(k) else g.tit <- ggtitle(paste0(k, "_", con)); g <- g+g.tit
    if (lgd==TRUE) {
      g <- g+theme(plot.margin=margin(0.005, 0.005, 0.2, 0, "npc"), plot.title=element_text(hjust=0.5, size=legend.plot.title.size))+ggtitle(legend.plot.title)
    }; return(g)

  }

  # Map colours to samples according to expression level.
  # Column names without '__' are not excluded so as to keep one-to-one match with color.dat.
  cname <- colnames(gene)
  cons <- gsub("(.*)(__)(.*)", "\\3", cname)
  sam.uni <- unique(gsub("(.*)(__)(.*)", "\\1", cname)); ft.trans <- make.names(ft.trans)

  g.lis.all <- NULL; for (k in ID) {
    # Match color key values to selected genes.
    color.dat <- NULL; for (i in gene[k, ]) { 
      ab <- abs(i-geneV); col.ind <- which(ab==min(ab))[1]; color.dat <- c(color.dat, cols[col.ind])
    }

    # idx <- grep("__", cname); c.na <- cname[idx]
    # Column names without '__' are also included so as to keep one-to-one match with color.dat.
    tis.col <- gsub("(.*)(__)(.*)", "\\1", cname)
    # Valid/target tissues.
    tis.col.uni <- unique(tis.col); tis.path.uni <- unique(tis.path)
    tis.col.path.idx <- tis.col.uni %in% tis.path.uni
    tis.tar <- tis.col.uni[tis.col.path.idx]
    # Un-related aSVG mapped to data.
    rematch.dif.svg <- is.list(lis.rematch) & length(tis.tar)==0
    if (rematch.dif.svg) {
      # Take the features in data that are assinged new aSVG features.
      tis.tar <- unlist(lapply(names(lis.rematch), function(i) {
        vec0 <- lis.rematch[[i]]
        if (!is.null(vec0)) if (all(vec0 %in% tis.path.uni)) return(i)
      }))
    }; if (length(tis.tar)==0) return(NULL)
    cname1 <- cname[grepl(paste0('^(', paste0(tis.tar,collapse='|'), ')__'), cname)]
    # Only conditions paired with valid tissues (have matching samples in data) are used. 
    con.vld <- gsub("(.*)(__)(.*)", "\\3", cname1); con.vld.uni <- unique(con.vld)
    na0 <- paste0(k, "_", con.vld.uni); cat('ggplot: ', k, ', ', sep = ''); g.lis <- lapply(con.vld.uni, g_list, ...); cat('\n')
    names(g.lis) <- na0; g.lis.all <- c(g.lis.all, g.lis)
  }; g.lgd <- g_list(con=NULL, lgd=TRUE, ...)
  return(list(g.lgd = g.lgd, g.lis.all = g.lis.all))
}

#' Convert SHMs of Ggplot to Grobs
#'
#' @param gg.lis A list of SHMs of ggplot, returned by \code{\link{gg_shm}}.
#' @param cores The number of CPU cores for parallelization, relevant for aSVG files with size larger than 5M. The default is NA, and the number of used cores is 1 or 2 depending on the availability.
#' @return A list of grob SHMs.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @importFrom parallel detectCores mclapply

grob_shm <- function(gg.lis, cores=2) {
  tmp <- normalizePath(tempfile(), winslash='/', mustWork=FALSE)
  cat('Converting "ggplot" to "grob" ... \n')
  # mclapply does not give speed gain here.
  # Child jobs (on each core) are conducted in different orders, which can be reflected by "cat", but all final jobs are assembled in the original order.
  # Repress popups by saving it to a png file, then delete it.
  nas <- names(gg.lis)
  png(tmp); grob.lis <- mclapply(seq_along(gg.lis), function(x) {
    cat(nas[x], ' '); x <- gg.lis[[x]]
    # Remove legends in SHMs.
    x <- x+theme(legend.position="none"); ggplotGrob(x) 
  }, mc.cores=cores); dev.off(); cat('\n')
  names(grob.lis) <- nas; return(grob.lis)
}




