#' The SVG class for storing annotated SVG (aSVG) instances
#'
#' The SVG class is designed to represent annotated SVG (aSVG) instances. 
#' 
#' @param coordinate A named \code{list} of x-y coordinates parsed from one or multile aSVG files respectively. Coordinates are represented in three columns \code{x}, \code{y}, and \code{feature} in form of \code{data.frame} or \code{tbl}, corresponding to x, y coordinates, and spatial features (cellular compartments, tissues, organs, etc.) in aSVGs respectively. The \code{list} name slots refer to aSVG instances respectively, e.g. list(SVGInstance1=coordinate1, SVGInstance2=coordinate2). 
#' @param attribute A named \code{list} of attributes of coordinates in \code{coordinate}. Attributes are represented in at least four columns \code{feature}, \code{id}, \code{fill} and \code{stroke} in form of \code{data.frame} or \code{tbl}, corresponding to \code{feature} in \code{coordinate}, ids of \code{feature}, fill colors of \code{feature}, and line widths of \code{feature} respectively. \code{id} can be the same as \code{feature} or ontology ids. The \code{list} name slots refer to aSVG instances respectively and must match those in \code{coordinate}, e.g. list(SVGInstance1=attribute1, SVGInstance2=attribute2).  
#' @param dimension A named \code{list} of width/height parsed from one or multile aSVG files respectively, which is calculated from \code{coordinate} automatically. Each pair of width/height is stored in a named vector. The \code{list} name slots refer to aSVG instances respectively and must match those in \code{coordinate}, e.g. list(SVGInstance1=c(with=100, height=8), SVGInstance2=c(with=20, height=15)).  
#' @param svg A named \code{list} of directory path(s) of one or multile aSVG files respectively. The \code{list} name slots refer to aSVG instances respectively and must match those in \code{coordinate}, e.g. list(SVGInstance1=svg.path1, SVGInstance2=svg.path2).  
#' @param raster A named \code{list} of directory path(s) of one or multile raster image files (jpg, png) respectively. This argument is relevant only when superimposing raster images with spatial heatmap plots that are created from aSVG images. The default is \code{NULL} for each aSVG instance. aSVG images are usually created by using these raster images as templates, otherwise spatial features between the two will not match. The \code{list} name slots refer to aSVG instances respectively and must match those in \code{coordinate}, e.g. list(SVGInstance1=raster.path1, SVGInstance2=raster.path2).  

#' @return A SVG object.
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @examples
#' 
#' # The first raste image used as a template to create an aSVG. 
#' raster.pa1 <- system.file('extdata/shinyApp/example/maize_leaf_shm1.png',
#' package='spatialHeatmap')
#' # The first aSVG created with the first template. 
#' svg.pa1 <- system.file('extdata/shinyApp/example/maize_leaf_shm1.svg',
#' package='spatialHeatmap')
#' # The second raster image used as a template to create an aSVG. 
#' raster.pa2 <- system.file('extdata/shinyApp/example/maize_leaf_shm2.png',
#' package='spatialHeatmap')
#' # The second aSVG created with the second template. 
#' svg.pa2 <- system.file('extdata/shinyApp/example/maize_leaf_shm2.svg',
#' package='spatialHeatmap')
#'
#' # Parse these two aSVGs without association with raster images.
#' svgs <- read_svg(svg.path=c(svg.pa1, svg.pa2), raster.path=NULL)
#'
#' # Parse these two aSVGs. The raster image paths are provide so as to 
#' # be associated with respective aSVGs, which will be used when 
#' # superimposing raster images with SHM plots.
#' svgs <- read_svg(svg.path=c(svg.pa1, svg.pa2), raster.path=c(raster.pa1, raster.pa2))
#'
#' # Two aSVG instances are stored in a "SVG" object of "svgs".
#' names(svgs)
#' # Access content of "svgs".
#' svgs[1, ] # The first aSVG instance
#' svgs[, 'coordinate'][1]; coordinate(svgs)[1] # The coordinates of the first aSVG instance
#' # Combine two "SVG" objects.
#' x <- svgs[1, ]; y <- svgs[2, ]; cmb(x, y)
#' # Extract slots from "svgs" and create a new "SVG" object.
#' lis <- list(cordn=coordinate(svgs), attrb=attribute(svgs), svg=svg(svgs), raster=raster(svgs))
#' new.svgs <- SVG(coordinate=lis$cordn, attribute=lis$attrb, svg=lis$svg, raster=lis$raster)
#' # Change aSVG instance names.
#' names(new.svgs) <- c('aSVG1', 'aSVG2')

#' @docType class
#' @export
#' @rdname coord

SVG <- function(coordinate=list(), attribute=list(), dimension=list(), svg=list(), raster=list()) {
  if (length(dimension)==0) {
    dimension <- lapply(coordinate, function(x) c(width=1, height=1))
  }
  if (length(raster)==0) {
    raster <- lapply(coordinate, function(x) return(NULL))
  }
  if (length(svg)==0) {
    svg <- lapply(coordinate, function(x) return(names(coordinate)[x]))
  }
  check_SVG(coord=coordinate, attr=attribute, wh=dimension, svg=svg, raster=raster)
  dimension <- lapply(seq_along(dimension), function(i, cordn) {
    cordn0 <- cordn[[i]] 
    w.h <- c(max(cordn0$x) - min(cordn0$x), max(cordn0$y) - min(cordn0$y))
    names(w.h) <- c('width', 'height'); w.h
  }, coordinate)
  names(dimension) <- names(coordinate)
  methods::new('SVG', coordinate=coordinate, attribute=attribute, dimension=dimension, svg=svg, raster=raster)
}

#' # Check validity of slots.
#' @keywords Internal
#' @noRd

check_SVG <- function(coord, attr, wh, svg, raster) {
  if (!is(coord, 'list') | !is(attr, 'list') | !is(wh, 'list')) stop('Each slot shoule be a named "list"!')
  len.inst <- unique(c(length(coord), length(attr), length(wh), length(svg), length(raster)))
  if (length(len.inst) > 1) stop('Each instance should have a value in each slot!')
  na.cor <- names(coord); na.att <- names(attr)
  na.wh <- names(wh); na.svg <- names(svg)
  na.raster <- names(raster)
  if (is.null(na.cor) | is.null(na.att) | is.null(na.wh) | is.null(na.svg) | is.null(na.raster)) stop('Each slot should be a "list" named by instance(s)!')
  for (i in seq_len(len.inst)) {
    inst <- table(c(na.cor[i], na.att[i], na.wh[i], na.svg[i], na.raster[i]))
    if (length(inst)!=1) stop('An instance should have the same name across slots!')
  }
  for (i in seq_along(coord)) {
    cordn <- coord[[i]]; attr0 <- attr[[i]]
    if (!is(cordn, 'data.frame') & !is(cordn, 'tbl')) {
      message(na.cor[i], ':')
      stop('The "coordinate" list should contain objects of "data.fram" or "tbl"!')
    }
    if (!is(attr0, 'data.frame') & !is(attr0, 'tbl')) {
      message(na.cor[i], ':')
      stop('The "attribute" list should contain objects of "data.fram" or "tbl"!')
    }
    if (!all(c('x', 'y', 'feature') %in% colnames(cordn))) stop('The table in "coordinate" should contain at least three columns: x, y, and feature!')
    if (!is.numeric(cordn$x)) stop('"x" in "coordinate" should be numeric!')
    if (!is.numeric(cordn$y)) stop('"y" in "coordinate" should be numeric!')
    if (!all(c('feature', 'id', 'fill', 'stroke') %in% colnames(attr0))) stop('The table in "attribute" should contain at least four columns: feature, id, fill, stroke!')
    if (!is.numeric(attr0$stroke)) stop('"stroke" in "attribute" should be numeric!')

    ft.cordn <- sort(unique(sub('__\\d+$', '', cordn$feature)))
    ft.attr <- sort(unique(attr0$feature))
    if (length(ft.cordn)!=length(ft.attr)) { 
      message(na.cor[i], ':')
      stop('Different features are detected between "coordinate" and "attribute"!')
    }
    if (!all(ft.cordn==ft.attr)) { 
      message(na.cor[i], ':')
      stop('Different features are detected between "coordinate" and "attribute"!')
    }
  }
}


setMethod("show", c(object="SVG"), function(object) {
  x <- object; nas <- names(x)
  if (length(nas) > 1) {
    na.pas <- paste(paste(nas[1:2], collapse=','), '...')
  } else na.pas <- nas
  message('Class "SVG": ', length(x), ' instance(s)', ' - ', na.pas)
  cat('Slots: coordinate, attribute, dimension, ... \n')
  message(nas[1], ':')
  cordn <- coordinate(x)[[1]]; attrb <- attribute(x)[[1]]
  if (nrow(cordn) > 3) r.cordn <- 1:3 else r.cordn <- seq_len(nrow(cordn))
  message('coordinate:'); print(cordn[r.cordn, ])
  if (nrow(attrb) > 3) r.attrb <- 1:3 else r.attrb <- seq_len(nrow(attrb))
  message('attribute:'); print(attrb[r.attrb, ])
  dim1 <- dimension(x)[[1]]
  cat('dimension: '); cat(names(dim1), '\n')
  cat('           ', dim1,'\n')
  raster <- raster(x)[[1]]
  if (!is.null(raster)) message('raster: ', raster)
})


