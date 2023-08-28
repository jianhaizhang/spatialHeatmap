#' Extract Attributes in an SVG File
#' 
#' This function extract basic attributes of nodes in an SVG file. The 'a' nodes are not deleted.

#' @param doc An object of class "xml_document", i.e. the imported SVG file in R.
#' @param feature A character vector of features/samples extracted from the data. If some of the input features are duplicated in SVG file, then a reminder message is returned.
#' @param br TRUE or FALSE. If TRUE, combined paths are broken apart and all styles are updated.
#' @return A 3-component list of attribute data frame, the outline node, and the sample node. 
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2
#' Müller K, Wickham H (2022). _tibble: Simple Data Frames_. R package version 3.1.7, <https://CRAN.R-project.org/package=tibble>  
#' Wickham H, François R, Henry L, Müller K (2022). _dplyr: A Grammar of Data Manipulation_. R package version 1.0.9, <https://    CRAN.R-project.org/package=dplyr> 

#' @importFrom xml2 xml_length xml_children xml_name xml_attr xml_remove xml_text
#' @importFrom tibble tibble 
#' @importFrom dplyr filter

svg_attr <- function(doc, feature, br=TRUE) {
  options(stringsAsFactors=FALSE); name <- element <- NULL
  len <- xml_length(doc); out <- xml_children(doc)[[len-1]]; ply <- xml_children(doc)[[len]]

  # If out is not a group, it is assigned an empty node. The "metadata" element should not be changed. "out" might be the "metadata" element, so before deleting elements in path_br make sure "out" has an empty element, so does "ply".
  if (xml_name(out)!='g') { xml_add_child(out, 'empty', .where=0); out1 <- xml_children(out)[[1]]; xml_remove(xml_children(out)[[1]], free=FALSE); out <- out1 }
  # If ply is not a group, it is assigned an empty node.
  if (xml_name(ply)!='g') { xml_add_child(ply, 'empty', .where=0); ply1 <- xml_children(ply)[[1]]; xml_remove(xml_children(ply)[[1]], free=FALSE); ply <- ply1 }
  # EBI SVG, if the outline shapes and tissue shapes are separate, they must be in two layers NOT two groups. Otherwise, 'fill' and 'stroke' in '.ps.xml' can be  messy.
  if (length(xml_children(out))>0 & length(xml_children(ply))>0) {
    g.mode <- xml_attr(out, 'groupmode')=='layer' & xml_attr(ply, 'groupmode')=='layer'
    msg <- 'If outline and regular spatial features are in two separate groups, the two groups must have the "groupmode" of "layer"!'
    if (is.na(g.mode)) return(msg) else if (g.mode==FALSE) return(msg)
  }
  # Break combined path to a group or siblings.
  if (br==TRUE) { 
    res1 <- path_br_all(out); if (is.character(res1)) return(res1)
    res2 <- path_br_all(ply); if (is.character(res2)) return(res2) 
  }
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  chdn.all <- c(chdn.out, chdn.ply)
  id.mat <- NULL; for (i in chdn.all) {
    idx <- grepl('^matrix',xml_attr(i, 'transform'))
     if (idx==TRUE) id.mat <- c(id.mat, xml_attr(i, 'id'))
  }
  if (grepl('^matrix',xml_attr(out, 'transform'))) id.mat <- c(id.mat, xml_attr(out, 'id'))
  if (grepl('^matrix',xml_attr(ply, 'transform'))) id.mat <- c(id.mat, xml_attr(ply, 'id'))
  if (!is.null(id.mat)) { cat('\n'); cat("Recommendation: please remove the 'transform' attribute with a 'matrix' value in the following groups by ungrouping and regrouping the respective groups in Inkscape. Otherwise, colors in spatial heatmap might be shifted! \n"); cat(id.mat, '\n\n') }
  
  ## Exrtact basic attributes into a data frame.
  idx <- seq_len(length(chdn.out)+length(chdn.ply))
  idx1 <- c(seq_len(length(chdn.out)), seq_len(length(chdn.ply)))
  parent <- make.names(c(rep(xml_attr(out, 'id'), length(chdn.out)), rep(xml_attr(ply, 'id'), length(chdn.ply))))
  nas <- c(xml_name(chdn.out), xml_name(chdn.ply))
  ids <- make.names(c(xml_attr(chdn.out, 'id'), xml_attr(chdn.ply, 'id')))
  if (any(duplicated(ids))) return(paste0('Duplicated node ids detected: ', paste0(ids[duplicated(ids)], collapse=' '), '!'))
  title <- make.names(c(vapply(chdn.out, tit_id, character(1)), vapply(chdn.ply, tit_id, character(1)))) # Use original names, no 'make.names'. '' after applied to 'make.names' becomes 'X'.
  w <- which(title=='X'|title==''); title[w] <- ids[w]
  # Duplicated titles.
  dup <- duplicated(title); if (any(dup)) { 
    tit.dup <- unique(title[dup]) 
    if (length(intersect(tit.dup, unique(make.names(feature))))>0) return(paste0('Duplicated title detected: ', paste0(tit.dup, collapse=' '), '!')) else {
      cat('Duplicated title text detected:', title[dup], '\n'); w <- title %in% tit.dup; title[w] <- paste0(title[w], seq_len(sum(w)))
    }
  }
  # Style inside groups are ignored. 
  sty <- c(xml_attr(chdn.out, 'style'), xml_attr(chdn.ply, 'style'))
  sty[!grepl(':', sty)] <- 'none'
  # Extract attribute values such as fill colors or stroke widths.
  sty_val <- function(style, attr) {

    w1 <- grepl(attr, style); st <- style[w1]; st <- strsplit(st, ';')
    st1 <- NULL; for (i in st) { st1 <- c(st1, i[grepl(attr, i)]) }; style[w1] <- st1
    # If only keep part of the string, the pattern should cover everything in the string, e.g. the '.*' on both ends.
    value <- gsub(paste0(".*(", attr, ")(.*).*"), '\\2', style); return(value)
  
  }; fil.cols <- sty_val(sty, 'fill:')
  fil.cols[!grepl('^#', fil.cols)] <- 'none'
  stro.w <- sty_val(sty, 'stroke-width:')
  # Convert real width to numeric.
  for (i in seq_along(stro.w)) {
    num <- tryCatch({ as.numeric(stro.w[i]) }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
    stro.w[i] <- ifelse(is.numeric(num), num, 0)
  }
  # index.all: counting groups of outlines and main shapes together.
  # index.sub: counting groups of outlines and main shapes independently.
  df.attr <- tibble(feature=title, id=ids, fill=fil.cols, stroke=as.numeric(stro.w), parent=parent, element=nas, index.all=idx, index.sub=idx1)
  df.attr <- filter(df.attr, element!='a')
  # 'fill' is not necessary. In Inkscape, resizing a "group" causes "matrix" in "transform" (relative positions) attribute, and this can lead to related polygons uncolored in the spatial heatmaps. Solution: ungroup and regroup to get rid of transforms and get absolute positions.
  # Change 'style' of all polygons. Since in SVG code, if no fill in style, no fill in ".ps.xml", so is the stroke.
  # "stroke" >= 0.51 px always introduces coordinates in .ps.xml, no matter "fill" is "none" or not. If "stroke" < 0.5 px, even though "fill" is not "none" there is no coordinates in ps.xml. E.g. irregular paths of dots.
  if (br==TRUE) {
    style <- 'fill:#46e8e8;fill-opacity:1;stroke:#000000;stroke-width:3;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1' 
    for (i in chdn.all) {
      xml_set_attr(i, 'style', style)
      if (xml_name(i)=='g') xml_set_attr(xml_children(i), 'style', style)
    }
  }; return(list(df.attr=df.attr, out=out, ply=ply))

}
