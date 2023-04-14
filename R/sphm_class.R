#' The SPHM class
#'
#' The SPHM class is designed to store numeric data and image objects for plotting spatial heatmaps.
#' 
#' @param svg An \code{SVG} object containing one or multiple aSVG instances (see \code{\link{SVG}} and \code{\link{read_svg}}). In the aSVGs, spatial features (tissues, organs, etc) having counterparts with the same identifiers in the `bulk` data will be colored accoording to expression profiles of chosen biomolecules (genes, proteins, etc). 
#' @param bulk The bulk data in form of numeric \code{vector}, \code{data.frame}, or \code{SummarizedExperiment}. See the `data` argument of the function \code{\link{filter_data}}.
#' @param cell The single-cell data in form of \code{SingleCellExperiment}. In the `colData` slot, a column that stores cell group labels is required.
#' @param match The \code{list} for re-matching in SHMs (only bulk data) or matching between single-cell and bulk data in co-visualization.
#' \describe{ 
#'  \item{SHMs}{ A named \code{list} for rematching spatial features between numeric data (ftA, ftB) and aSVGs (ftC, ftD, ftE). In each slot, the slot name is a spatial feature from the data and the corresponding element is one or multiple spatial features from the aSVG. \emph{E.g.} \code{list(ftA = c('ftC', 'ftD'), ftB = c('ftE'))}.
#'  } 
#'  \item{Co-visualization plots}{ Mapping cells to tissues: a named \code{list}, where cell group labels from \code{colData(sce.  dimred)[, 'cell.group']} are the name slots and aSVG features are the corresponding \code{list} elements. Mapping tissues to       cells: a named \code{list}, where tissues are the name slots and cells from \code{colData(sce.dimred)[, 'cell.group']} are the     corresponding \code{list} elements. Applicable when cell grouping methods are annodation labels, marker genes, clustering, or      manual assignments. 
#'  }
#' }
#' @param output A \code{list} of outputs, which is automatically generated. 

#' @return An SPHM object.
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @examples
#' 
#' library(SummarizedExperiment)
#' # Import single-cell data.
#' sce.pa <- system.file("extdata/shinyApp/data", "cell_mouse_brain.rds", 
#' package="spatialHeatmap")
#' sce <- readRDS(sce.pa)
#' # Pre-processing.
#' sce.dimred.quick <- process_cell_meta(sce, 
#' qc.metric=list(subsets=list(Mt=rowData(sce)$featureType=='mito'), threshold=1)) 
#' colData(sce.dimred.quick)[1:3, 1:2] 
#' sce.aggr.quick <- aggr_rep(sce.dimred.quick, assay.na='logcounts', sam.factor='label', 
#' aggr='mean')
#' # Import the aSVG image.   
#' svg.mus.brain.pa <- system.file("extdata/shinyApp/data", "mus_musculus.brain.svg", 
#' package="spatialHeatmap")
#' svg.mus.brain <- read_svg(svg.mus.brain.pa) 
#' # List for mapping single cells to bulk. 
#' lis.match.quick <- list(hypothalamus=c('hypothalamus'), cortex.S1=c('cerebral.cortex', 'nose'))
#'
#' # SPHM class for storing aSVG, bulk/sc data, and matching list. 
#' dat.quick <- SPHM(svg=svg.mus.brain, bulk=sce.aggr.quick, cell=sce.dimred.quick, 
#' match=lis.match.quick)
#'
#' # Co-visualization plot. 
#' # covis(data=dat.quick, ID=c('Apod'), dimred='PCA', cell.group='label', 
#' # tar.cell=names(lis.match.quick), assay.na='logcounts', bar.width=0.11, dim.lgd.nrow=1, 
#' # height=0.7, legend.r=1.5, legend.key.size=0.02, legend.text.size=12, legend.nrow=3) 

#' @docType class
#' @export
#' @rdname SPHM

SPHM <- function(svg=NULL, bulk=NULL, cell=NULL, match=list(), output=list()) {
  check_SPHM(svg=svg, bulk=bulk, cell=cell, match=match, output=output)
  methods::new('SPHM', svg=svg, bulk=bulk, cell=cell, match=match, output=output)
}

#' # Check validity of slots.
#' @keywords Internal
#' @noRd

check_SPHM <- function(svg=NULL, bulk=NULL, cell=NULL, match=NULL, output=NULL) {
  if (!is.null(svg)) if (!is(svg, 'SVG')) stop('"svg" should be an "SVG" object!')
  if (!is.null(bulk)) if (!is(bulk, 'numeric') & !is(bulk, 'data.frame') & !is(bulk, 'DFrame') & !is(bulk, 'SummarizedExperiment')) stop('"bulk" should be one of "numeric vector", "data.frame", "DFrame", or "SummarizedExperiment" object!')
  if (!is.null(cell)) if (!is(cell, 'SingleCellExperiment') & !is.null(cell) & !is(cell, 'Seurat')) stop('"cell" should be "NULL", "SingleCellExperiment", or "Seurat" object!')
  if (!is.null(match)) {
    if (!is(match, 'list')) stop('"match" should be a "list" object!') else {
      if (length(match) > 0 ) if (is.null(names(match)) | '' %in% names(match)) stop('"match" should be a named "list" object!')
    }
  }
  if (!is.null(output)) {
    if (!is(output, 'list')) stop('"output" should be a "list" object!') else {
      if (length(output) > 0 ) if (is.null(names(output)) | '' %in% names(output)) stop('"output" should be a named "list" object!')
    }
  }
}

setMethod("show", c(object="SPHM"), function(object) {
  message('Class "SPHM": slots - svg, bulk, cell, match, output')
  x <- object; svg <- x@svg
  message('svg:')
  if (is(svg, 'SVG')) {
    attrb <- attribute(svg)[[1]]
    if (nrow(attrb) > 3) r.attrb <- seq_len(3) else r.attrb <- seq_len(nrow(attrb))
    print(as.data.frame(attrb[r.attrb, c('feature', 'id', 'fill', 'stroke', 'sub.feature', 'index')]))
    message('...')
  }
  message('bulk:'); bulk <- x@bulk
  if (is(bulk, 'numeric')) {
    if (length(bulk) > 3) { print(bulk[seq_len(3)]); message('...') } else print(bulk)
  } else if (is(bulk, 'data.frame')) {
    if (nrow(bulk) > 3) { print(bulk[seq_len(3), , drop=FALSE]); message('...') } else print(bulk)
  } else if (is(bulk, 'SummarizedExperiment')) {
    print(bulk)
  }

  message('cell:'); cell <- x@cell
  if (is(cell, 'SingleCellExperiment') | is(cell, 'Seurat')) print(cell)
  message('match:'); match <- x@match
  if (is(match, 'list')) {
    if (length(match) > 3) print(match[seq_len(3)]) else if (length(match)>0) print(match)
  }
  message('output:'); output <- x@output
  if (is(output, 'list')) {
    if (length(output)>0) print(names(output))
  }
})



