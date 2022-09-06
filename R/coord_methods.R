#' @title
#' Methods for S4 class \code{coord}
#'
#' @description
#' These are methods for subsetting, getting, setting, or combining \linkS4class{coord} objects.
#'
#' @section Main methods:
#' In the following code snippets, \code{cordn} is a \linkS4class{coord} object.
#' \describe{
#' \item{\code{cordn[i]}, \code{cordn[i, ]}}{
#'   Subsetting the ith aSVG instance.
#' }
#' \item{\code{cordn[i] <- cordn.new}}{
#'   Replacing the ith aSVG instance in \code{cordn} with a new \code{cordn} object \code{cordn.new}.
#' }
#' \item{\code{cordn[, j]}}{
#'   Subsetting the jth slot of all aSVG instances.
#' }
#' \item{\code{cordn[, 'coordinate']}, \code{coordinate(cordn)}}{
#'   Subsetting the \code{coordinate} slot that contains coordinates of all aSVG instances.
#' }
#' \item{\code{length(cordn)}}{
#'   Number of all aSVG instances.
#' }
#' \item{\code{names(cordn)}, \code{names(cordn)[1] <- 'newName'}}{
#'   Names of all aSVG instances, rename the first aSVG instance.
#' }
#' \item{\code{cbm(cordn1, cordn2)}}{
#'   Combining two aSVG instances.
#' }
#' }
#' 
#' @author
#' Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' 
#' @seealso
#' \code{\link{coord}}: creating \code{coord} objects.
#' 
#' @name coordMethods
#' @docType methods
#' @aliases coordinate attribute dimension svg raster cmb
#' @examples

#' # Create the first aSVG instance. 
#' svg.pa1 <- system.file('extdata/shinyApp/example/maize_leaf_shm1.svg',
#' package='spatialHeatmap')
#' svg1 <- parse_svg(svg.path=c(svg.pa1)); names(svg1); length(svg1); slotNames(svg1)
#' # Create the second aSVG instance. 
#' svg.pa2 <- system.file('extdata/shinyApp/example/maize_leaf_shm2.svg',
#' package='spatialHeatmap')
#' svg2 <- parse_svg(svg.path=c(svg.pa2)); names(svg2); length(svg2)
#' # Combine these two instances.
#' svg3 <- cmb(svg1, svg2); names(svg3); length(svg3)
#' # The first aSVG instance
#' svg3[1]
#' # Coordinates of the first aSVG instance 
#' svg3[, 'coordinate'][1]; coordinate(svg3)[1]
#' # Extract slots from "svg3" into a list and create a new "coord" object.
#' lis <- list(cordn=coordinate(svg3), attrb=attribute(svg3), svg=svg(svg3))
#' new.svgs <- coord(coordinate=lis$cordn, attribute=lis$attrb, svg=lis$svg)
#' # Change aSVG instance names.
#' names(new.svgs) <- c('aSVG1', 'aSVG2'); names(new.svgs)
#' # Replace an instance.
#' svg3[2] <- new.svgs[2]
#' # Replace a slot content.
#' coordinate(svg3)[[1]] <- coordinate(new.svgs)[[1]]
NULL

 
#' @rdname coordMethods
#' @export
setMethod("coordinate", "coord", function(x) { x@coordinate })
#' @rdname coordMethods
#' @export
setReplaceMethod("coordinate", "coord", function(x, value) {
    x@coordinate <- value; x
})

#' @rdname coordMethods
#' @export
setMethod("attribute", "coord", function(x) { x@attribute })

#' @rdname coordMethods
#' @export
setReplaceMethod("attribute", "coord", function(x, value) {
  x@attribute <- value; x
})

#' @rdname coordMethods
#' @export
setMethod("dimension", "coord", function(x) { x@dimension })
#' @rdname coordMethods
#' @export
setReplaceMethod("dimension", "coord", function(x, value) {
    x@dimension <- value; x
})


#' @rdname coordMethods
#' @export
setMethod("raster", "coord", function(x) { x@raster })
#' @rdname coordMethods
#' @export
setReplaceMethod("raster", "coord", function(x, value) {
    x@raster <- value; x
})

#' @rdname coordMethods
#' @export
setMethod("svg", "coord", function(x) { x@svg })
#' @rdname coordMethods
#' @export
setReplaceMethod("svg", "coord", function(x, value) {
    x@svg <- value; x
})


#' @rdname coordMethods
#' @export
setMethod("[", c("coord"), function(x, i, j) {
  if (!missing(i)) coord0 <- new('coord', coordinate=x@coordinate[i], attribute=x@attribute[i], dimension=x@dimension[i], svg=x@svg[i], raster=x@raster[i])
  if (!missing(j)) lis <- list(coordinate=x@coordinate, attribute=x@attribute, dimension=x@dimension, svg=x@svg, raster=x@raster)
  if (missing(j) & !missing(i)) return(coord0)
  if (missing(i) & !missing(j)) return(lis[[j]]) 
  if (!missing(i) & !missing(j)) {
    if (is(j, 'numeric')) return(slot(coord0, names(lis)[j]))
    if (is(j, 'character')) return(slot(coord0, lis[j]))
  }
  if (missing(i) & missing(j)) return(x)
})

#' @rdname coordMethods
#' @export
setMethod("[<-", c("coord"), function(x, i, value) {
  if (!is(value, 'list') & !is.null(value) & !is(value, 'coord')) stop('The replacement should be a coord, list or NULL!')
  if (is.null(value)) {
    slot(x, 'coordinate')[i] <- NULL
    slot(x, 'attribute')[i] <- NULL
    slot(x, 'dimension')[i] <- NULL
    slot(x, 'svg')[i] <- NULL
    slot(x, 'raster')[i] <- NULL
    return(x)
  } 
  if (is(value, 'coord')) value <- list(coordinate=value@coordinate, attribute=value@attribute, dimension=value@dimension, svg=value@svg, raster=value@raster)
  if (length(value$dimension)==0) {
    value$dimension <- lapply(value$coordinate, function(x) c(width=1, height=1))
  }
  if (length(value$svg)==0) {
    value$svg <- lapply(value$coordinate, function(x) return(names(value$coordinate)[x]))
  }
  if (length(value$raster)==0) {
    value$raster <- lapply(value$coordinate, function(x) return(NULL))
  }
  check_coord(coord=value$coordinate, attr=value$attribute, wh=value$dimension, svg=value$svg, raster=value$raster)
  
  lis <- list(coordinate=x@coordinate, attribute=x@attribute, dimension=x@dimension, svg=x@svg, raster=x@raster)
  for (k in seq_along(lis)) {
    na0 <- names(value[[k]])
    if (is.null(na0)) stop('The "list" provided to "value" should be named!')
    lis[[k]][i] <- value[[k]]; nas0 <- names(lis[[k]])
    if (is(i, 'numeric')) {
      if (na0 %in% nas0[-i]) {
        message(na0, ':')
        stop('Instance names should be unique!') 
      }
      names(lis[[k]])[i] <- na0
    } else if (is(i, 'character')) {
      if (is.null(nas0)) names(lis[[k]]) <- na0 else {
      names(lis[[k]])[nas0==i] <- na0
      }
    }
  }
  new('coord', coordinate=lis$coordinate, attribute=lis$attribute, dimension=lis$dimension, svg=lis$svg, raster=lis$raster)
})

#' @rdname coordMethods
#' @export
setMethod("length", "coord", function(x) {
  max(c(length(x@coordinate), length(x@attribute), length(x@dimension), length(x@svg), length(x@raster)))
})

#' @rdname coordMethods
#' @export
setMethod("names", "coord", function(x) {
  unique(c(names(x@coordinate), names(x@attribute), names(x@dimension), names(x@svg), names(x@raster)))
})


#' @rdname coordMethods
#' @export
setMethod("names<-", "coord", function(x, value) {
  # if (length(i)!=length(value)) stop('"i" and "value" should have the same size!')
  if (!is(x,  'coord')) stop('"x" should be an object of "coord"!')
  if (any(duplicated(value))) stop('Names should be unique!')
  lis <- list(coordinate=x@coordinate, attribute=x@attribute, dimension=x@dimension, svg=x@svg, raster=x@raster)
  lis <- lapply(lis, function(x) { names(x) <- value; x })
  new('coord', coordinate=lis$coordinate, attribute=lis$attribute, dimension=lis$dimension, svg=lis$svg, raster=lis$raster)
})

#' @rdname coordMethods
#' @export
setMethod("cmb", c(x="coord", y='coord'), function(x, y) {
  if (!is(x, 'coord') | !is(y, 'coord')) stop('The input should be coord classes!')
  if (length(intersect(names(x), names(y)))>0) stop('Instance names should be unique!')
  new('coord', coordinate=c(x@coordinate, y@coordinate), attribute=c(x@attribute, y@attribute), dimension=c(x@dimension, y@dimension), svg=c(x@svg, y@svg), raster=c(x@raster, y@raster))
})










