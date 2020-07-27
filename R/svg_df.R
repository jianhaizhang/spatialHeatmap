#' Extract Coordinates, Sample Names, and Colors from the SVG File
#'
#' @param svg.path The path of an SVG file.
#' @inheritParams svg_attr

#' @return A 3-length list, the first component is a data frame of the coordinates, the second is a vector of all sample/path names, and the third is the fill colors.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr    

#' @importFrom rsvg rsvg_ps 
#' @importFrom grImport PostScriptTrace 
#' @importFrom xml2 xml_length xml_children xml_name xml_attr xml_remove xml_text

svg_df <- function(svg.path, feature) {

  # Make sure the style is correct. If the stroke width is not the same across polygons such as '0.0002px', '0.216px', some stroke outlines cannot be recognised by 'PostScriptTrace'. Then some polygons are missing. Since the ggplot is based on 'stroke' not 'fill'.
  options(stringsAsFactors=FALSE)
  doc <- read_xml(svg.path); spa <- xml_attr(doc, 'space')
  if (!is.na(spa)) if (spa=='preserve') xml_set_attr(doc, 'xml:space', 'default')
  # w.h <- c(xml_attr(doc, 'width'), xml_attr(doc, 'height'))
  # w.h <- as.numeric(gsub("^(\\d+\\.\\d+|\\d+).*", "\\1", w.h))
  # names(w.h) <- c('width', 'height')
  svg.attr <- svg_attr(doc, feature=feature); if (is(svg.attr, 'character')) return(svg.attr)
  df.attr <- svg.attr[['df.attr']]; df.attr$title <- make.names(df.attr$title)
  out <- svg.attr[['out']]; ply <- svg.attr[['ply']]
  # Paths in 'a' node are recognised in .ps.xml file, so all 'a' nodes in or out groups are removed. 
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  chdn.all <- c(chdn.out, chdn.ply)
  for (i in chdn.all) {

    na <- xml_name(i)
    if (na=='a') xml_remove(i, free=FALSE) else if (na=='g') {

      chil <- xml_children(i); for (j in chil) {

        na1 <- xml_name(j); if (na1=='a') xml_remove(j, free=FALSE) else if (na1=='g') return(paste0('Nested group detected in ', xml_attr(i, 'id'), '!'))

      }

    }
 
  }
  # Renew the children after deletion of 'a' nodes.
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  chdn.all <- c(chdn.out, chdn.ply)

  # Get ids and titles for every path, including paths inside groups, except for 'a' nodes.
  tit <- id.all <- NULL; for (i in seq_along(chdn.all)) {

    if (df.attr[i, 'name']=='g') {

     na <- xml_name(xml_children(chdn.all[[i]]))
     tit0 <- rep(df.attr[i, 'title'], xml_length(chdn.all[[i]])-sum(na=='title')); tit0 <- paste0(tit0, '_', seq_along(tit0)); tit <- c(tit, tit0)
     id0 <- rep(df.attr[i, 'id'], xml_length(chdn.all[[i]])-sum(na=='title')); id.all <- c(id.all, id0)
     # If the styles in paths of a group are different with group style, they can lead to messy 'fill' and 'stroke' in '.ps.xml', so they are set NULL. This step is super important.
     xml_set_attr(xml_children(chdn.all[[i]]), 'style', NULL)

    } else if (df.attr[i, 'name']=='use') {

      ref <- paste0('#', df.attr[, 'id'])
      w <- which(ref %in% xml_attr(chdn.all[[i]], 'href'))
      # If reference is inside a group. A group contains no groups, so the use node has 1 shape.
      if (length(w)==0) { tit <- c(tit, df.attr[i, 'title']); id.all <- c(id.all, df.attr[i, 'id']) }
      # If reference is outside a group.
      if (length(w)>0) if (df.attr[w, 'name']=='g') {

        na <- xml_name(xml_children(chdn.all[[w]]))
        tit0 <- rep(df.attr[i, 'title'], xml_length(chdn.all[[w]])-sum(na=='title')); tit0 <- paste0(tit0, '_', seq_along(tit0)); tit <- c(tit, tit0)
        id0 <- rep(df.attr[i, 'id'], xml_length(chdn.all[[w]])-sum(na=='title')); id.all <- c(id.all, id0)

      } else { tit <- c(tit, df.attr[i, 'title']); id.all <- c(id.all, df.attr[i, 'id']) }

    } else { tit <- c(tit, df.attr[i, 'title']); id.all <- c(id.all, df.attr[i, 'id']) }

  }; tis.path <- gsub("_\\d+$", "", tit)  

  style <- 'fill:#46e8e8;fill-opacity:1;stroke:#000000;stroke-width:3;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1' # 'fill' is not necessary. In Inkscape, "group" or move an object adds transforms (relative positions), and this can lead to related polygons uncolored in the spatial heatmaps. Solution: ungroup and regroup to get rid of transforms and get absolute positions.
  # Change 'style' of all polygons.
  xml_set_attr(chdn.out, 'style', style); xml_set_attr(chdn.ply, 'style', style)  
  # xml_set_attr(out, 'style', style); xml_set_attr(ply, 'style', style)  
  # Export internal SVG.
  tmp <- tempdir(check=TRUE); svg.inter <- paste0(tmp, '/internal.svg')
  if (grepl("~", svg.inter)) svg.inter <- normalizePath(svg.inter)
  write_xml(doc, file=svg.inter)
  
  # SVG file conversion. 
  rsvg_ps(svg.inter, file=sub("svg$", "ps", svg.inter)) # Only the paths inside canvas of SVG are valid.
  p1 <- sub("svg$", "ps", svg.inter); p2 <- paste0(sub("svg$", "ps", svg.inter), ".xml"); PostScriptTrace(p1, p2) 
  chdn1 <- xml_children(read_xml(p2)) # Use internal svg to get coordinates.
     
  # Detect groups that use relative coordinates ("transform", "matrix" in Inkscape.), which leads to some plygons missing in ".ps.xml" file.
  # EBI SVG, if the outline shapes and tissue shapes are separate, they must be in two layers NOT two groups. Otherwise, 'fill' and 'stroke' in '.ps.xml' can be  messy.
  fil.stk <- xml_attr(chdn1[-length(chdn1)], 'type'); tab <- table(fil.stk)
  w <- which(fil.stk=='fill')%%2==0
  if (any(w) & tab['fill'] > tab['stroke']) { 
 
    # Index of wrong path.
    w1 <- which(w)[1]
    # Wrong path and related group.
    tis.wrg <- paste0(id.all[c(w1-1, w1)], collapse='; ')
    return(paste0("Error detected in '", tis.wrg, "' in SVG image. Please ungroup and regroup the respective group they belong to.")) 

  }
  
  # Get coordinates from '.ps.xml'.
  stroke <- chdn1[which(xml_attr(chdn1, 'type')=='stroke')]
  df <- NULL; for (i in seq_along(stroke)) {

    xy <- xml_children(stroke[[i]])[-1]
    x <- as.numeric(xml_attr(xy, 'x'))
    y <- as.numeric(xml_attr(xy, 'y'))
    df0 <- cbind(tissue=tit[i], data.frame(x=x, y=y), stringsAsFactors=TRUE) # The coordinates should not be factor.
    df <- rbind(df, df0)

  }; fil.cols <- df.attr$fil.cols; names(fil.cols) <- df.attr$title
  w.h <- c(max(abs(df$x)), max(abs(df$y))) 
  names(w.h) <- c('width', 'height')
  lis <- list(df=df, tis.path=sub('_\\d+$', '', tit), fil.cols=fil.cols, w.h=w.h); return(lis)

}


