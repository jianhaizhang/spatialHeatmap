#' Trace the painted paths in an SVG file without Ghostscript
#'
#' A drop-in replacement for \code{rsvg_ps} followed by \code{grImport::PostScriptTrace}. The SVG is normalised by librsvg/cairo via \code{rsvg_svg}, which resolves styles, groups, "use" elements and relative coordinates into plain absolute paths (only M/L/C/Z commands and an optional transform matrix). These paths are flattened in R and written to an XML file in the RGML format of \code{PostScriptTrace}: one "path" node of type "fill" and/or "stroke" per painting operation, in painting order, with a "context" child followed by "move"/"line" children, and a final "summary" node with the path count. Coordinates use the same device space as \code{PostScriptTrace} (y axis upwards, 10 units per pixel).

#' @param svg.path The path of the input SVG file.
#' @param xml.path The path of the output XML file.
#' @param tol The flattening tolerance of Bezier curves in pixels.
#' @return The path of the output XML file (invisibly).
#' @keywords Internal
#' @noRd

#' @references
#' Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2

#' @importFrom rsvg rsvg_svg
#' @importFrom xml2 read_xml xml_ns_strip xml_find_all xml_attr xml_parents xml_name

svg_trace <- function(svg.path, xml.path, tol=0.1) {
  cairo.svg <- tempfile(fileext='.svg'); on.exit(unlink(cairo.svg))
  rsvg_svg(svg.path, cairo.svg)
  doc <- read_xml(cairo.svg); xml_ns_strip(doc)
  # Map the user space of the cairo SVG to PostScript device space: y upwards, 10 units per pixel.
  w <- as.numeric(sub('[a-z]+$', '', xml_attr(doc, 'width')))
  h <- as.numeric(sub('[a-z]+$', '', xml_attr(doc, 'height')))
  vb <- as.numeric(strsplit(trimws(xml_attr(doc, 'viewBox')), '[ ,]+')[[1]])
  if (length(vb)!=4) vb <- c(0, 0, w, h)
  sx <- w/vb[3]; sy <- h/vb[4]
  dev <- c(10*sx, 0, 0, -10*sy, -10*sx*vb[1], 10*(h+sy*vb[2]))
  # Painted shapes only, i.e. not the clip paths, glyphs, patterns, etc. in "defs".
  nodes <- xml_find_all(doc, '//path[not(ancestor::defs)]|//rect[not(ancestor::defs)]')
  out <- character(0)
  for (node in nodes) {
    fill <- xml_attr(node, 'fill'); stroke <- xml_attr(node, 'stroke')
    # In SVG the default fill is black and the default stroke is none.
    types <- c(if (is.na(fill) || fill!='none') 'fill', if (!is.na(stroke) && stroke!='none') 'stroke')
    if (length(types)==0) next
    if (xml_name(node)=='rect') {
      x <- as.numeric(xml_attr(node, 'x')); y <- as.numeric(xml_attr(node, 'y'))
      x[is.na(x)] <- 0; y[is.na(y)] <- 0
      rw <- as.numeric(xml_attr(node, 'width')); rh <- as.numeric(xml_attr(node, 'height'))
      d <- sprintf('M %f %f L %f %f L %f %f L %f %f Z', x, y, x+rw, y, x+rw, y+rh, x, y+rh)
    } else d <- xml_attr(node, 'd')
    if (is.na(d)) next
    xy <- path_flatten(d, tol)
    if (nrow(xy)==0) next
    # Compose the transforms of the node and its ancestors.
    mat <- c(1, 0, 0, 1, 0, 0)
    for (nod in c(list(node), as.list(xml_parents(node)))) {
      tr <- xml_attr(nod, 'transform')
      if (!is.na(tr)) mat <- mat_mult(svg_transform(tr), mat)
    }
    mat <- mat_mult(dev, mat)
    x <- mat[1]*xy$x + mat[3]*xy$y + mat[5]
    y <- mat[2]*xy$x + mat[4]*xy$y + mat[6]
    pts <- paste0('\t<', ifelse(xy$move, 'move', 'line'), " x='", signif(x, 7), "' y='", signif(y, 7), "'/>", collapse='\n')
    for (type in types) {
      out <- c(out, paste0("<path type='", type, "' id='", length(out)+1, "'>\n\t<context/>\n", pts, '\n</path>'))
    }
  }
  # The bounding box in "summary" is not used downstream, so only the count is reported.
  cat("<?xml version='1.0' encoding='UTF-8' ?>\n<picture version='3'>\n", paste0(out, '\n'), "<summary count='", length(out), "'/>\n</picture>\n", file=xml.path, sep='')
  invisible(xml.path)
}

#' Flatten the path data of a cairo SVG path into line segments
#'
#' Cairo only writes absolute "M", "L", "C" and "Z" commands. Cubic Bezier curves are subdivided uniformly, where the number of segments is derived from the curvature (Wang's formula) so that the deviation is below \code{tol}.

#' @param d The "d" attribute of an SVG path.
#' @param tol The flattening tolerance.
#' @return A data frame of x/y coordinates, and a logical column indicating the start of subpaths.
#' @keywords Internal
#' @noRd

path_flatten <- function(d, tol=0.1) {
  tok <- regmatches(d, gregexpr('[A-Za-z]|[-+]?(\\d+\\.?\\d*|\\.\\d+)([eE][-+]?\\d+)?', d))[[1]]
  cmds <- grep('^[A-Za-z]$', tok)
  bad <- setdiff(unique(tok[cmds]), c('M', 'L', 'C', 'Z', 'z'))
  if (length(bad)>0) stop('Unsupported path command(s) in normalised SVG: ', paste0(bad, collapse=', '))
  x <- y <- numeric(0); move <- logical(0); cur <- start <- c(0, 0)
  for (i in seq_along(cmds)) {
    cmd <- tok[cmds[i]]; end <- if (i<length(cmds)) cmds[i+1]-1 else length(tok)
    v <- if (end>cmds[i]) as.numeric(tok[(cmds[i]+1):end]) else numeric(0)
    if (cmd=='M') {
      # Additional coordinate pairs after "M" are implicit "L" commands.
      pm <- matrix(v, nrow=2); cur <- start <- pm[, ncol(pm)]
      x <- c(x, pm[1, ]); y <- c(y, pm[2, ]); move <- c(move, TRUE, rep(FALSE, ncol(pm)-1))
    } else if (cmd=='L') {
      pm <- matrix(v, nrow=2); cur <- pm[, ncol(pm)]
      x <- c(x, pm[1, ]); y <- c(y, pm[2, ]); move <- c(move, rep(FALSE, ncol(pm)))
    } else if (cmd=='C') {
      pm <- matrix(v, nrow=6)
      for (j in seq_len(ncol(pm))) {
        p <- cbind(cur, pm[1:2, j], pm[3:4, j], pm[5:6, j])
        dd <- max(sqrt(colSums((p[, 1:2]-2*p[, 2:3]+p[, 3:4])^2)))
        n <- max(1, min(100, ceiling(sqrt(0.75*dd/tol))))
        t <- seq_len(n)/n
        b <- cbind((1-t)^3, 3*(1-t)^2*t, 3*(1-t)*t^2, t^3) %*% t(p)
        x <- c(x, b[, 1]); y <- c(y, b[, 2]); move <- c(move, rep(FALSE, n))
        cur <- pm[5:6, j]
      }
    } else if (cmd %in% c('Z', 'z')) {
      x <- c(x, start[1]); y <- c(y, start[2]); move <- c(move, FALSE); cur <- start
    }
  }; data.frame(x=x, y=y, move=move)
}

#' Parse an SVG transform attribute into an affine matrix c(a, b, c, d, e, f)
#' @param tr The "transform" attribute.
#' @keywords Internal
#' @noRd

svg_transform <- function(tr) {
  mat <- c(1, 0, 0, 1, 0, 0)
  ops <- regmatches(tr, gregexpr('[a-zA-Z]+\\s*\\([^)]*\\)', tr))[[1]]
  for (op in ops) {
    na <- trimws(sub('\\(.*', '', op))
    v <- as.numeric(strsplit(trimws(gsub('.*\\(|\\)', '', op)), '[ ,]+')[[1]])
    m <- switch(na,
      matrix=v,
      translate=c(1, 0, 0, 1, v[1], if (length(v)>1) v[2] else 0),
      scale=c(v[1], 0, 0, if (length(v)>1) v[2] else v[1], 0, 0),
      rotate={
        a <- v[1]*pi/180; r <- c(cos(a), sin(a), -sin(a), cos(a), 0, 0)
        if (length(v)==3) r <- mat_mult(mat_mult(c(1, 0, 0, 1, v[2], v[3]), r), c(1, 0, 0, 1, -v[2], -v[3]))
        r
      },
      skewX=c(1, 0, tan(v[1]*pi/180), 1, 0, 0),
      skewY=c(1, tan(v[1]*pi/180), 0, 1, 0, 0),
      c(1, 0, 0, 1, 0, 0))
    mat <- mat_mult(mat, m)
  }; mat
}

#' Multiply two affine matrices in the form of c(a, b, c, d, e, f)
#' @param m1,m2 The affine matrices. The result applies \code{m2} first.
#' @keywords Internal
#' @noRd

mat_mult <- function(m1, m2) {
  c(m1[1]*m2[1]+m1[3]*m2[2], m1[2]*m2[1]+m1[4]*m2[2],
    m1[1]*m2[3]+m1[3]*m2[4], m1[2]*m2[3]+m1[4]*m2[4],
    m1[1]*m2[5]+m1[3]*m2[6]+m1[5], m1[2]*m2[5]+m1[4]*m2[6]+m1[6])
}
