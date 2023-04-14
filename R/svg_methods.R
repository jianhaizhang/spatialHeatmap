#' @title
#' Methods for S4 class \code{SVG}
#'
#' @description
#' These are methods for subsetting, getting, setting, or combining \linkS4class{SVG} objects.
#'
#' @section Main methods:
#' In the following code snippets, \code{cordn} is a \linkS4class{SVG} object.
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
#' @return An object of \code{SVG}, \code{data.frame}, or \code{numeric}. 
#' @author
#' Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' 
#' @seealso
#' \code{\link{SVG}}: creating \code{SVG} objects.
#' 
#' @name SVGMethods
#' @rdname SVGMethods
#' @docType methods
#' @aliases coordinate coordinate<- attribute attribute<- dimension dimension<- svg_pa svg_pa<- raster_pa raster_pa<- cmb names names<-,SVG-method sub_sf angle angle<-
#' @examples
#'
#' # Create the first aSVG instance. 
#' svg.pa1 <- system.file('extdata/shinyApp/data/maize_leaf_shm1.svg',
#' package='spatialHeatmap')
#' svg1 <- read_svg(svg.path=c(svg.pa1)); names(svg1); length(svg1); slotNames(svg1)
#' # Create the second aSVG instance. 
#' svg.pa2 <- system.file('extdata/shinyApp/data/maize_leaf_shm2.svg',
#' package='spatialHeatmap')
#' svg2 <- read_svg(svg.path=c(svg.pa2)); names(svg2); length(svg2)
#' # Combine these two instances.
#' svg3 <- cmb(svg1, svg2); names(svg3); length(svg3)
#' # The first aSVG instance
#' svg3[1]
#' # Coordinates of the first aSVG instance 
#' svg3[, 'coordinate'][1]; coordinate(svg3)[1]
#' # Extract slots from "svg3" into a list and create a new "SVG" object.
#' lis <- list(cordn=coordinate(svg3), attrb=attribute(svg3), svg=svg_pa(svg3))
#' new.svgs <- SVG(coordinate=lis$cordn, attribute=lis$attrb, svg=lis$svg)
#' # Change aSVG instance names.
#' names(new.svgs) <- c('aSVG1', 'aSVG2'); names(new.svgs)
#' # Replace the second instance in "svg3".
#' svg3[2] <- new.svgs[2]
#' # Replace a slot content.
#' coordinate(svg3)[[1]] <- coordinate(new.svgs)[[1]]
NULL

 
#' @rdname SVGMethods
#' @export
setMethod("coordinate", "SVG", function(x) { x@coordinate })
#' @rdname SVGMethods

#' @export
#' @param value A value for replacement.
#' @importFrom methods slot slot<-

setReplaceMethod("coordinate", "SVG", function(x, value) {
  index <- NULL; x@coordinate <- value 
  for (i in seq_along(x)) { 
    svg0 <- x[i]; cordn0 <- coordinate(svg0)[[1]] 
    attr0 <- attribute(svg0)[[1]] 
    inter <- unique(intersect(cordn0$index, attr0$index)) 
    slot(x[i], 'coordinate')[[1]] <- subset(cordn0, index %in% inter)
    slot(x[i], 'attribute')[[1]] <- subset(attr0, index %in% inter) 
  }  
  check_SVG(coord=x@coordinate, attr=x@attribute, wh=x@dimension, svg=x@svg, raster=x@raster, angle=x@angle)
  x 
})

#' @rdname SVGMethods
#' @export
setMethod("attribute", "SVG", function(x) { x@attribute })

#' @rdname SVGMethods
#' @param value A value for replacement.
#' @references
#' Wickham H, François R, Henry L, Müller K (2022). _dplyr: A Grammar of Data Manipulation_. R package version 1.0.9, <https://    CRAN.R-project.org/package=dplyr>
#' @importFrom dplyr filter mutate %>%  
#' @export
#' @importFrom methods slot slot<-

setReplaceMethod("attribute", "SVG", function(x, value) {
  index <- feature <- NULL
  x@attribute <- value 
  for (k in seq_along(x)) { 
    svg0 <- x[k]; cordn0 <- coordinate(svg0)[[1]] 
    attr0 <- attribute(svg0)[[1]] 
    # Keep common subfeatures between coordinates and attributes. 
    inter <- unique(intersect(cordn0$index, attr0$index)) 
    cordn0 <- filter(cordn0, index %in% inter)
    attr0 <- filter(attr0, index %in% inter)

    # Identify regrouped subfeatures.
    cord0.uni <- filter(cordn0, !duplicated(index))
    df.sub.ft <- sub('__\\d+$' , '', cord0.uni$feature)
    sub.ft.attr <- attr0$feature
    w.regrp <- which(df.sub.ft!=sub.ft.attr)
    cord0.uni.idx <- cord0.uni$index
    # Change subfeatures identifiers in coordinates according to regrouping in attributes.
    if (sum(w.regrp)>0) for (i in w.regrp) {
      cordn0$feature <- as.vector(cordn0$feature)
      ft.attr0 <- sub.ft.attr[i]
      # Features change from factor to vector.
      cordn0 <- cordn0 %>% mutate(feature= replace(feature, index==cord0.uni.idx[i], ft.attr0))
      df.regrp <- filter(cordn0, grepl(paste0('^', ft.attr0, '$|^', ft.attr0, '__\\d+$'), feature))
      idx0 <- unique(df.regrp$index)
      # If the regrouped subfeatures have the same identifiers with other features, append '__\\d+$' to all these features.
      if (length(idx0)>1) {
        idx.sub <- seq_along(idx0)
        for (j in idx.sub) {
          df0 <- filter(df.regrp, index==idx0[j])
          # Features change from factor to vector.
          df0 <- df0 %>% mutate(df0, feature=paste0(sub('__\\d+', '', feature), '__', j))
          cordn0[cordn0$index %in% df0$index, ] <- df0
        }
      }
      # cordn0$feature <- factor(cordn0$feature, levels=unique(cordn0$feature))
    }
    slot(x[k], 'coordinate')[[1]] <- cordn0
    slot(x[k], 'attribute')[[1]] <- attr0 
    # Update dimension.
    w.h <- c(max(cordn0$x) - min(cordn0$x), max(cordn0$y) - min(cordn0$y))
    names(w.h) <- c('width', 'height')
    slot(x[k], 'dimension')[[1]] <- w.h 
  }  
  check_SVG(coord=x@coordinate, attr=x@attribute, wh=x@dimension, svg=x@svg, raster=x@raster, angle=x@angle)
  x 
})  


#' @rdname SVGMethods
#' @export
setMethod("dimension", "SVG", function(x) { x@dimension })

#' @rdname SVGMethods
#' @param value A value for replacement.
#' @export
setReplaceMethod("dimension", "SVG", function(x, value) {
    x@dimension <- value; x
})


#' @rdname SVGMethods
#' @export
setMethod("raster_pa", "SVG", function(x) { x@raster })

#' @rdname SVGMethods
#' @param value A value for replacement.
#' @export
setReplaceMethod("raster_pa", "SVG", function(x, value) {
    x@raster <- value; x
})

#' @rdname SVGMethods
#' @export
setMethod("svg_pa", "SVG", function(x) { x@svg })

#' @rdname SVGMethods
#' @param value A value for replacement.
#' @export
setReplaceMethod("svg_pa", "SVG", function(x, value) {
    x@svg <- value; x
})

#' @rdname SVGMethods
#' @export
setMethod("angle", "SVG", function(x) { x@angle })

#' @rdname SVGMethods
#' @param value A value for replacement.
#' @export
setReplaceMethod("angle", "SVG", function(x, value) {
    x@angle <- value; x
})

#' @rdname SVGMethods
#' @param i,j Two integers specifying an aSVG instance and a slot of the same aSVG respectively. 
#' @export
#' @importFrom methods new slot 

setMethod("[", c("SVG"), function(x, i, j) {
  if (!missing(i)) coord0 <- new('SVG', coordinate=x@coordinate[i], attribute=x@attribute[i], dimension=x@dimension[i], svg=x@svg[i], raster=x@raster[i], angle=x@angle[i])
  if (!missing(j)) lis <- list(coordinate=x@coordinate, attribute=x@attribute, dimension=x@dimension, svg=x@svg, raster=x@raster, angle=x@angle)
  if (missing(j) & !missing(i)) return(coord0)
  if (missing(i) & !missing(j)) return(lis[[j]]) 
  if (!missing(i) & !missing(j)) {
    if (is(j, 'numeric')) return(slot(coord0, names(lis)[j]))
    if (is(j, 'character')) return(slot(coord0, lis[j]))
  }
  if (missing(i) & missing(j)) return(x)
})

#' @rdname SVGMethods
#' @param value A value for replacement.
#' @export
#' @importFrom methods slot slot<- new

setMethod("[<-", c("SVG"), function(x, i, value) {
  if (!is(value, 'list') & !is.null(value) & !is(value, 'SVG')) stop('The replacement should be a coord, list or NULL!')
  if (is.null(value)) {
    slot(x, 'coordinate')[i] <- NULL
    slot(x, 'attribute')[i] <- NULL
    slot(x, 'dimension')[i] <- NULL
    slot(x, 'svg')[i] <- NULL
    slot(x, 'raster')[i] <- NULL
    slot(x, 'angle')[i] <- NULL
    return(x)
  } 
  if (is(value, 'SVG')) value <- list(coordinate=value@coordinate, attribute=value@attribute, dimension=value@dimension, svg=value@svg, raster=value@raster, angle=value@angle)
  if (length(value$dimension)==0) {
    value$dimension <- lapply(value$coordinate, function(x) c(width=1, height=1))
  }
  if (length(value$svg)==0) {
    value$svg <- lapply(value$coordinate, function(x) return(names(value$coordinate)[x]))
  }
  if (length(value$raster)==0) {
    value$raster <- lapply(value$coordinate, function(x) return(NULL))
  }
  if (length(value$angle)==0) {
    value$angle <- lapply(value$coordinate, function(x) return(NULL))
  }
  check_SVG(coord=value$coordinate, attr=value$attribute, wh=value$dimension, svg=value$svg, raster=value$raster, angle=value$angle)
  
  lis <- list(coordinate=x@coordinate, attribute=x@attribute, dimension=x@dimension, svg=x@svg, raster=x@raster, angle=x@angle)
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
  new('SVG', coordinate=lis$coordinate, attribute=lis$attribute, dimension=lis$dimension, svg=lis$svg, raster=lis$raster, angle=lis$angle)
})

#' @rdname SVGMethods
#' @export
setMethod("length", "SVG", function(x) {
  max(c(length(x@coordinate), length(x@attribute), length(x@dimension), length(x@svg), length(x@raster), length(x@angle)))
})

#' @rdname SVGMethods
#' @export
setMethod("names", "SVG", function(x) {
  unique(c(names(x@coordinate), names(x@attribute), names(x@dimension), names(x@svg), names(x@raster), names(x@angle)))
})

#' @rdname SVGMethods
#' @param value A value for replacement.
#' @export
#' @importFrom methods new
setReplaceMethod("names", "SVG", function(x, value) {
  # if (length(i)!=length(value)) stop('"i" and "value" should have the same size!')
  if (!is(x,  'SVG')) stop('"x" should be an object of "SVG"!')
  if (any(duplicated(value))) stop('Names should be unique!')
  lis <- list(coordinate=x@coordinate, attribute=x@attribute, dimension=x@dimension, svg=x@svg, raster=x@raster, angle=x@angle)
  lis <- lapply(lis, function(x) { names(x) <- value; x })
  new('SVG', coordinate=lis$coordinate, attribute=lis$attribute, dimension=lis$dimension, svg=lis$svg, raster=lis$raster, angle=lis$angle)
})

#' @rdname SVGMethods
#' @param x,y Two \code{SVG} objects.
#' @export
#' @importFrom methods new
setMethod("cmb", c(x="SVG", y='SVG'), function(x, y) {
  if (!is(x, 'SVG') | !is(y, 'SVG')) stop('The input should be SVG classes!')
  if (length(intersect(names(x), names(y)))>0) stop('Instance names should be unique!')
  new('SVG', coordinate=c(x@coordinate, y@coordinate), attribute=c(x@attribute, y@attribute), dimension=c(x@dimension, y@dimension), svg=c(x@svg, y@svg), raster=c(x@raster, y@raster), angle=c(x@angle, y@angle))
})


#' @rdname SVGMethods
#' @param svg An \code{SVG} object.
#' @param show,hide Two vectors of indexes in the \code{attribute} slot. aSVG features corresponding to these indexes will be shown or hidden in spatial heatmap plots respectively.

#' @references
#' Wickham H, François R, Henry L, Müller K (2022). _dplyr: A Grammar of Data Manipulation_. R package version 1.0.9, <https://CRAN.R-project.org/package=dplyr>

#' @export
#' @importFrom dplyr filter 

setMethod("sub_sf", c(svg="SVG"), function(svg, show=NULL, hide=NULL) {
  index <- NULL
  if (!is.null(show) & !is.null(hide)) stop('At least one of "show" and "hide" should be NULL!')
  if (is.null(show) & is.null(hide)) return(svg)
  svg0 <- svg[1]; df.attr <- attribute(svg0)[[1]]
  if (!is.null(show)) {
    if (any(!show %in% df.attr$index)) stop('Ensure all entries in "show" are from "index" of "attribute"!')
    df.attr <- filter(df.attr, index %in% show)
  }
  if (!is.null(hide)) {
    if (any(!hide %in% df.attr$index)) stop('Ensure all entries in "hide" are from "index" of "attribute"!')
    df.attr <- filter(df.attr, !index %in% hide)
  }
  attribute(svg0)[[1]] <- df.attr; svg0
})
