########################################

# Methods for SPHM class.

#' @export
setGeneric("shm", function(data, ...) standardGeneric("shm"))
#' @export
setGeneric("covis", function(data, ...) standardGeneric("covis"))

#' @export
setGeneric("svg", function(x) standardGeneric("svg"))
#' @export
setGeneric("svg<-", function(x, value) standardGeneric("svg<-"))
#' @export
setGeneric("name", function(x) standardGeneric("name"))
#' @export
setGeneric("bulk", function(x) standardGeneric("bulk"))
#' @export
setGeneric("bulk<-", function(x, value) standardGeneric("bulk<-"))
#' @export
setGeneric("cell", function(x) standardGeneric("cell"))
#' @export
setGeneric("cell<-", function(x, value) standardGeneric("cell<-"))
#' @export
setGeneric("match", function(x) standardGeneric("match"))
#' @export
setGeneric("match<-", function(x, value) standardGeneric("match<-"))
#' @export
setGeneric("output", function(x) standardGeneric("output"))
#' @export
setGeneric("output<-", function(x, value) standardGeneric("output<-"))


# Methods for SVG class.


#' @export
setGeneric("spatial_hm", function(svg, ...) standardGeneric("spatial_hm"))

#' @export
setGeneric("cmb", function(x, y) standardGeneric("cmb"))

#' @export
setGeneric("coordinate", function(x) standardGeneric("coordinate"))

#' @export
setGeneric("coordinate<-", function(x, value) standardGeneric("coordinate<-"))

#' @export
setGeneric("attribute", function(x) standardGeneric("attribute"))

#' @export
setGeneric("attribute<-", function(x, value) standardGeneric("attribute<-"))

#' @export
setGeneric("dimension", function(x) standardGeneric("dimension"))

#' @export
setGeneric("dimension<-", function(x, value) standardGeneric("dimension<-"))

#' @export
setGeneric("svg_obj", function(x) standardGeneric("svg_obj"))
#' @export
setGeneric("svg_obj<-", function(x, value) standardGeneric("svg_obj<-"))
#' @export
setGeneric("raster_pa", function(x) standardGeneric("raster_pa"))
#' @export
setGeneric("raster_pa<-", function(x, value) standardGeneric("raster_pa<-"))

#' @export
setGeneric("angle", function(x) standardGeneric("angle"))
#' @export
setGeneric("angle<-", function(x, value) standardGeneric("angle<-"))

#' @export
setGeneric("sub_sf", function(svg, ...) standardGeneric("sub_sf"))


