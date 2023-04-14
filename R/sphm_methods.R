#' @title
#' Methods for S4 class \code{SPHM}
#'
#' @description
#' These are methods for getting or setting \linkS4class{SPHM} objects.
#'
#' @section Main methods:
#' In the following code snippets, \code{x} is an \linkS4class{SPHM} object.
#' \describe{
#' \item{\code{x[i]}, svg(x), bulk(x), cell(x), match(x), output(x)}{
#'   Subsetting a slot of \code{x}.
#' }
#' \item{x[i] <- value, svg(x) <- value, bulk(x) <- value, cell(x) <- value, match(x) <- value}{
#'   Replacing a slot value in \code{x}.
#' }
#' \item{\code{name(x)}}{
#'   All slot names in \code{x}.
#' }
#' \item{\code{length(x)}}{
#'   Number of all slots in \code{x}.
#' }
#' }
#' 
#' @param x An `SPHM` object.

#' @return An object of \code{SVG}, \code{SummarizedExperiment}, \code{SingleCellExperiment}, \code{data.frame}, \code{numeric}, or \code{list}. 
#' @author
#' Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' 
#' @seealso
#' \code{\link{SPHM}}: creating \code{SPHM} objects.
#' 
#' @name SPHMMethods
#' @rdname SPHMMethods
#' @docType methods
#' @aliases svg svg<- bulk bulk<- cell cell<- match match<- output output<- name length
#' @examples
#'
#' # Import an aSVG image.
#' svg.hum.pa <- system.file("extdata/shinyApp/data", 'homo_sapiens.brain.svg', 
#' package="spatialHeatmap")
#' svg.hum <- read_svg(svg.hum.pa)
#' # aSVG features.
#' feature.hum <- attribute(svg.hum)[[1]]
#' set.seed(20)
#' unique(feature.hum$feature)[1:10] 
#' 
#' # Testing vector.
#' my_vec <- setNames(sample(1:100, 4), c('substantia.nigra', 'putamen', 'prefrontal.cortex', 
#' 'notMapped')); my_vec   
#' # An SPHM object for storing bulk data and aSVG.
#' dat.quick <- SPHM(svg=svg.hum, bulk=my_vec) 
#' # Getters and setters for SPHM objects.
#' svg(dat.quick); svg(dat.quick) <- svg.hum
#' dat.quick['bulk']; dat.quick['bulk'] <- my_vec

NULL

 
#' @rdname SPHMMethods
#' @export
setMethod("svg", "SPHM", function(x) { x@svg })

#' @rdname SPHMMethods
#' @param value A value for replacement.
#' @export
#' @importFrom methods slot<-

setReplaceMethod("svg", "SPHM", function(x, value) {
  check_SPHM(svg=value); slot(x, 'svg') <- value; x 
})

#' @rdname SPHMMethods
#' @export
setMethod("bulk", "SPHM", function(x) { x@bulk })

#' @rdname SPHMMethods
#' @param value A value for replacement.
#' @export
#' @importFrom methods slot<-

setReplaceMethod("bulk", "SPHM", function(x, value) {
  check_SPHM(bulk=value); slot(x, 'bulk') <- value; x 
})

#' @rdname SPHMMethods
#' @export
setMethod("cell", "SPHM", function(x) { x@cell })

#' @rdname SPHMMethods
#' @param value A value for replacement.
#' @export
#' @importFrom methods slot<-

setReplaceMethod("cell", "SPHM", function(x, value) {
  check_SPHM(cell=value); slot(x, 'cell') <- value; x 
})

#' @rdname SPHMMethods
#' @export
setMethod("match", "SPHM", function(x) { x@match })

#' @rdname SPHMMethods
#' @param value A value for replacement.
#' @export
#' @importFrom methods slot<-

setReplaceMethod("match", "SPHM", function(x, value) {
  check_SPHM(match=value); slot(x, 'match') <- value; x 
})

#' @rdname SPHMMethods
#' @export
setMethod("output", "SPHM", function(x) { x@output })

#' @rdname SPHMMethods
#' @param value A value for replacement.
#' @export
#' @importFrom methods slot<-

setReplaceMethod("output", "SPHM", function(x, value) {
  check_SPHM(output=value); slot(x, 'output') <- value; x 
})


#' @rdname SPHMMethods
#' @param i An integer specifying a slot of an \code{SPHM} object. 
#' @export
#' @importFrom methods slot 

setMethod("[", c("SPHM"), function(x, i) {
  nas <- name(x)
  if (is(i, 'numeric')) return(slot(x, nas[i]))
  if (is(i, 'character')) {
    if (! i %in% nas) stop('The provided slot is not found!')
    return(slot(x, i))
  }
})

#' @rdname SPHMMethods
#' @param value A value for replacement.
#' @export
#' @importFrom methods slot<-

setMethod("[<-", c("SPHM"), function(x, i, value) {
  nas <- name(x)
  if (is(i, 'numeric')) slot(x, nas[i]) <- value
  if (is(i, 'character')) {
    if (! i %in% nas) stop('The provided slot is not found!')
    slot(x, i) <- value
  }
  check_SPHM(svg=x@svg, bulk=x@bulk, cell=x@cell, match=x@match, output=x@output)
  x
})

#' @rdname SPHMMethods
#' @importFrom methods slotNames
#' @export
setMethod("length", "SPHM", function(x) { length(slotNames(x)) })

#' @rdname SPHMMethods
#' @export
#' @importFrom methods slotNames

setMethod("name", "SPHM", function(x) { slotNames(x) })




