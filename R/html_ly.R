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

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC Florida, 2020.
#' Ramnath Vaidyanathan, Yihui Xie, JJ Allaire, Joe Cheng and Kenton Russell (2019). htmlwidgets: HTML Widgets for R. R package version 1.5.1. https://CRAN.R-project.org/package=htmlwidgets

html_ly <- function(gg.all, cs.g, aspr=1, anm.scale=1, selfcontained=FALSE, out.dir) {
  # save(gg.all, cs.g, anm.scale, selfcontained, out.dir, file ='all.ly')
  pkg <- check_pkg('plotly'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  pkg <- check_pkg('htmlwidgets'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  
  dir <- paste0(normalizePath(out.dir, winslash="/", mustWork=FALSE), '/html_shm')
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  lib.fi <- file.path(dir, 'lib')
  if (dir.exists(lib.fi)) unlink(lib.fi, recursive=TRUE)

  nas <- names(gg.all)
  gg.na <- grep('^dim_', nas, invert=TRUE, value=TRUE)
  cs.g$theme$aspect.ratio <- NULL
  cs.lis <- plotly::gg2list(cs.g, tooltip='color_scale')
  cs.lis$layout$title$text <- NULL 
  csly <- plotly::as_widget(cs.lis, tooltip='color_scale') 

  rd1 <- '1. Double click the "html" files to display the interactive spatial heatmaps on a web browser.'
  rd2 <- '2. All files in the "lib" folder are required to display the spatial heatmaps, so the "html" files cannot work independently.'
  writeLines(text=c(rd1, rd2), con=paste0(dir, '/README.txt'))
  
  for (i in gg.na) {    
    na.hl <- paste0(i, '.html')
    cat('Preparing', paste0("'", na.hl, "'"), '... \n')
    g <- gg.all[[i]]; lay.dat <- layer_data(g)

    anm.height <- 550 * anm.scale; asp <- g$theme$aspect.ratio
    if (!is.null(asp)) asp <- 1; anm.width <- anm.height / asp

    dat <- g$data; tis.path <- dat$Feature
    # ft.legend <- intersect(sam.uni, unique(tis.path))
    # ft.legend <- setdiff(ft.legend, ft.trans) 
    # leg.idx <- !duplicated(tis.path) & (tis.path %in% ft.legend)
    # tis.show <- as.vector(dat$Feature)[leg.idx]
    # tis.show1 <- tis.path[leg.idx]
    g$theme$aspect.ratio <- NULL # Aspect ratio is not accepted.
    g2l <- plotly::gg2list(g, tooltip="text")
    cat('Preparing legend for', paste0("'", na.hl, "'"), '... \n')

    ft.show <- unique(dat$feature[grepl('^#', lay.dat$fill)])
    ft.all <- lapply(g2l$data, function(x) x$name)
    idx.show <- which(!duplicated(sub('__\\d+$', '', ft.all)) & ft.all %in% ft.show)
    for (j in seq_along(g2l$data)) { 
      na0 <- g2l$data[[j]]$name
      if (j %in% idx.show) g2l$data[[j]]$name <- sub('__\\d+$', '', na0) else g2l$data[[j]]$showlegend <- FALSE 
    }; ggly <- plotly::as_widget(g2l)
    na.dim <- paste0('dim_', i); if (na.dim %in% nas) {
      g.dim <- gg.all[[na.dim]]
      lay.dat.dim <- layer_data(g.dim); dat.dim <- g.dim$data
      cell.show <- unique(dat.dim$feature[grepl('^#', dat.dim$fill)])
      g.dim$theme$aspect.ratio <- NULL
      g2l.dim <- plotly::gg2list(g.dim, tooltip="text")
      for (k in seq_along(g2l.dim$data)) {
      # if (!g2l.dim$data[[k]]$legendgroup %in% cell.show) g2l.dim$data[[k]]$showlegend <- FALSE
      # If color cells independently, there are two many legends.
      g2l.dim$data[[k]]$showlegend <- FALSE
      }; ggly.dim <- plotly::as_widget(g2l.dim)
      ggly.dim <- plotly::layout(ggly.dim, legend = list(orientation = "h", xanchor = "center", x = 0.5))
    }
    if (na.dim %in% nas) { 
      subly <- plotly::subplot(csly, ggly.dim, ggly, nrows=1, shareX=TRUE, shareY=FALSE, margin=0, widths=c(0.05, 0.475, 0.475))
      subly$width <- anm.width*2*aspr; subly$height <- anm.height
    } else {
      subly <- plotly::subplot(csly, ggly, nrows=1, shareX=TRUE, shareY=FALSE, margin=0, widths=c(0.05, 0.95))
      subly$width <- anm.width*aspr; subly$height <- anm.height
    }
    htmlwidgets::saveWidget(subly, na.hl, selfcontained=selfcontained, libdir="lib")
    file.rename(na.hl, file.path(dir, na.hl)) 
  }; # unlink(paste0(dir, '/lib'), recursive=TRUE)
  if (dir.exists('lib')) file.rename('lib', lib.fi)
  return(list(width=subly$width, height=subly$height)) 
}
