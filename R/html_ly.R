#' Build Interactive Spatial Heatmap in HTML
#'
#' @param gg A list of spatial heatmaps of ggplot.
#' @param cs.g The color key of ggplot.
#' @param sam.uni A vector of unique samples extracted from data matrix.
#' @inheritParams htmlwidgets::saveWidget
#' @inheritParams spatial_hm
#' @param ft.trans A character vector of tissue/spatial feature identifiers that will be set transparent. \emph{E.g} c("brain", "heart"). This argument is used when target features are covered by  overlapping features and the latter should be transparent.
#' @return HTML files of spatial heatmaps are saved in 'animaiton_shm'.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC Florida, 2020.
#' Ramnath Vaidyanathan, Yihui Xie, JJ Allaire, Joe Cheng and Kenton Russell (2019). htmlwidgets: HTML Widgets for R. R package version 1.5.1. https://CRAN.R-project.org/package=htmlwidgets

#' @importFrom plotly gg2list as_widget subplot 
#' @importFrom htmlwidgets saveWidget


html_ly <- function(gg, cs.g, ft.trans, sam.uni, anm.width, anm.height, selfcontained=FALSE, out.dir) {

  gg.na <- names(gg); cs.lis <- gg2list(cs.g, tooltip='color_scale')
  cs.lis$layout$title$text <- NULL 
  csly <- as_widget(cs.lis, tooltip='color_scale') 
  dir <- paste0(normalizePath(out.dir, winslash="/", mustWork=FALSE), '/html_shm')
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  rd1 <- '1. Double click the "html" files to display the interactive spatial heatmaps in a web browser.'
  rd2 <- '2. All files in the "lib" folder are required to display the spatial heatmaps, so the "html" files cannot work independently.'
  writeLines(text=c(rd1, rd2), con=paste0(dir, '/README.txt'))
  for (i in gg.na) {
    
    na.hl <- paste0(i, '.html')
    cat('Preparing', paste0("'", na.hl, "'"), '... \n')
    g <- gg[[i]]; lay.dat <- layer_data(g)
    dat <- g$data; g.col <- lay.dat$fill; names(g.col) <- dat$tissue
    g.col <- g.col[!duplicated(names(g.col))]; tis.path <- dat$feature
    ft.legend <- intersect(sam.uni, unique(tis.path))
    ft.legend <- setdiff(ft.legend, ft.trans) 
    leg.idx <- !duplicated(tis.path) & (tis.path %in% ft.legend)
    tis.show <- as.vector(dat$tissue)[leg.idx]
    tis.show1 <- tis.path[leg.idx]
    g2l <- gg2list(g, tooltip="text")
    cat('Preparing legend for', paste0("'", na.hl, "'"), '... \n')
    for (i in seq_along(g2l$data)) {
     
      idx.tis <- tis.show %in% g2l$data[[i]]$name
      if (any(idx.tis)) g2l$data[[i]]$name <- tis.show1[idx.tis] else g2l$data[[i]]$showlegend <- FALSE

    }; ggly <- as_widget(g2l)
    subly <- subplot(csly, ggly, nrows=1, shareX=FALSE, shareY=FALSE, margin=0, widths=c(0.05, 0.95))
    subly$width <- anm.width; subly$height <- anm.height
    saveWidget(subly, na.hl, selfcontained=selfcontained, libdir="lib")
    file.rename(na.hl, paste0(dir, '/', na.hl)) 

  }; unlink(paste0(dir, '/lib'), recursive=TRUE) 
  if (dir.exists('lib')) { lib1 <- paste0(dir, '/lib')
   if (dir.exists(lib1)) unlink(lib1, recursive=TRUE)
   file.rename('lib', lib1) 
  }

}

