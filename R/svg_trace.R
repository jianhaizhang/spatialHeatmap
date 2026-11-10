#' Trace the painted paths in an SVG file without Ghostscript
#'
#' A drop-in replacement for \code{rsvg_ps} followed by \code{grImport::PostScriptTrace}. The SVG is converted by librsvg/cairo via \code{rsvg_svg}, imported with \code{grImport2::readPicture}, and the flattened coordinates of the painted shapes are obtained from \code{grid::grobCoords}. The coordinates are written to an XML file in the RGML format of \code{PostScriptTrace}: one "path" node of type "fill" and/or "stroke" per painting operation, in painting order, with a "context" child followed by "move"/"line" children, and a final "summary" node with the path count. Coordinates use the same device space as \code{PostScriptTrace} (y axis upwards, 10 units per pixel).

#' @param svg.path The path of the input SVG file.
#' @param xml.path The path of the output XML file.
#' @return The path of the output XML file (invisibly).
#' @keywords Internal
#' @noRd

#' @references
#' Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg
#' Paul Murrell (2024). grImport2: Importing 'SVG' Graphics. R package version 0.3-3. https://CRAN.R-project.org/package=grImport2

#' @importFrom rsvg rsvg_svg
#' @importFrom grImport2 readPicture pictureGrob
#' @importFrom grid grobCoords
#' @importFrom grDevices pdf dev.off

svg_trace <- function(svg.path, xml.path) {
  cairo.svg <- tempfile(fileext='.svg'); on.exit(unlink(cairo.svg))
  rsvg_svg(svg.path, cairo.svg)
  # Newer cairo versions do not write the "surface" group checked by grImport2, which leads to a false warning.
  pic <- suppressWarnings(readPicture(cairo.svg))
  # Draw the picture on a device of w x h inches, so that grid coordinates in inches are pixels.
  w <- abs(diff(pic@summary@xscale)); h <- abs(diff(pic@summary@yscale))
  pdf(NULL, width=w, height=h); on.exit(dev.off(), add=TRUE)
  g <- pictureGrob(pic, x=0, y=0, width=w, height=h, default.units='in', just=c('left', 'bottom'), distort=TRUE, expansion=0)
  # Fills are closed paths and strokes are open polylines, which need different "closed" arguments.
  fills <- grobCoords(g, closed=TRUE); strokes <- grobCoords(g, closed=FALSE)
  out <- character(0)
  for (i in seq_along(fills)) for (j in seq_along(fills[[i]])) {
    # Each shape is a "picPath" (fill) followed by a "picPolyline" (stroke) if both are painted.
    stroke <- grepl('picPolyline', attr(fills[[i]][[j]], 'name'))
    part <- if (stroke) strokes[[i]][[j]] else fills[[i]][[j]]
    pts <- unlist(lapply(part, function(p) paste0('\t<', c('move', rep('line', length(p$x)-1)), " x='", signif(10*p$x, 7), "' y='", signif(10*p$y, 7), "'/>")))
    out <- c(out, paste0("<path type='", if (stroke) 'stroke' else 'fill', "' id='", length(out)+1, "'>\n\t<context/>\n", paste(pts, collapse='\n'), '\n</path>'))
  }
  cat("<?xml version='1.0' encoding='UTF-8' ?>\n<picture version='3'>\n", paste0(out, '\n'), "<summary count='", length(out), "'/>\n</picture>\n", file=xml.path, sep='')
  invisible(xml.path)
}
