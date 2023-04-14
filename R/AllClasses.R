#' @export
#' @rdname SVG
setClass("SVG", slots = c(coordinate='list', attribute='list', dimension='list', svg='list', raster='list', angle='list'))

#' @export
#' @rdname SPHM
setClass("SPHM", slots = c(svg='SVG', bulk='ANY', cell='ANY', match='list', output='list'))
