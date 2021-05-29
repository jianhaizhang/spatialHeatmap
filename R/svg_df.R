#' Extract Coordinates, Sample Names, and Colors from the SVG File
#'
#' @param svg.path The path of an SVG file. 
#' @inheritParams spatial_hm 
#' @inheritParams svg_attr

#' @return A 3-length list, the first component is a data frame of the coordinates, the second is a vector of all sample/path names, and the third is the fill colors.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr    

#' @importFrom rsvg rsvg_ps 
#' @importFrom grImport PostScriptTrace 
#' @importFrom xml2 xml_length xml_children xml_name xml_attr xml_remove xml_text xml_attrs
#' @importFrom data.table setDF rbindlist
#' @importFrom parallel detectCores mclapply

svg_df <- function(svg.path, feature=NULL, cores) {
  # save(svg.path, feature, cores, file='sfc')
  # Make sure the style is correct. If the stroke width is not the same across polygons such as '0.0002px', '0.216px', some stroke outlines cannot be recognised by 'PostScriptTrace'. Then some polygons are missing. Since the ggplot is based on 'stroke' not 'fill'.
  options(stringsAsFactors=FALSE)
  doc <- read_xml(svg.path); spa <- xml_attr(doc, 'space')
  if (!is.na(spa)) if (spa=='preserve') xml_set_attr(doc, 'xml:space', 'default')
  # Even though 'out' and 'ply' are not returned by 'svg_attr', the paths in doc are broken accordingly, since the node in doc are pointed, any change on the node is actually changing the doc. 
  # Paths in 'a' node are recognised in .ps.xml file, so all 'a' nodes in ply and out groups are removed. 
  svg.attr <- svg_attr(doc, feature=feature); if (is(svg.attr, 'character')) return(svg.attr)
  df.attr <- svg.attr[['df.attr']]; df.attr$feature <- make.names(df.attr$feature)
  out <- svg.attr[['out']]; ply <- svg.attr[['ply']]
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  na.all <- c(xml_name(chdn.out), xml_name(chdn.ply))
  chdn.all <- c(chdn.out, chdn.ply)
  na.no <- na.all[!na.all %in% c('g', 'path', 'rect', 'ellipse', 'use', 'title')]
  if (length(na.no)>0) { cat('\n\n'); cat('Warning: accepted SVG elements are "g", "path", "rect", "ellipse", "use", and "title". Please remove these elements in Inkscape:', na.no, '\n\n') }
  # Get ids and titles for every path, including paths inside groups, except for 'a' nodes.
  tit <- id.all <- NULL; for (i in seq_along(chdn.all)) {

    if (df.attr[i, 'element']=='g') {

     na <- xml_name(xml_children(chdn.all[[i]]))
     len <- xml_length(chdn.all[[i]])-sum(na=='title')
     tit0 <- rep(df.attr[i, 'feature'], len)
     # Add the distinct pattern '__\\d+$' to each path in a group for easy recognition downstream.
     if (len>1) tit0 <- paste0(tit0, '__', seq_along(tit0))
     tit <- c(tit, tit0)
     id0 <- rep(df.attr[i, 'id'], len); id.all <- c(id.all, id0)
     # If the styles in paths of a group are different with group style, they can lead to messy 'fill' and 'stroke' in '.ps.xml', so they are set NULL. This step is super important.
     # xml_set_attr(xml_children(chdn.all[[i]]), 'style', NULL)

    } else if (df.attr[i, 'element']=='use') {

      ref <- paste0('#', df.attr[, 'id'])
      w <- which(ref %in% xml_attr(chdn.all[[i]], 'href'))
      # If reference is inside a group, since a group contains no nested groups, so the reference is a single path and the use node must has 1 shape.
      if (length(w)==0) { tit <- c(tit, df.attr[i, 'feature']); id.all <- c(id.all, df.attr[i, 'id']) }
      # If reference is outside a group.
      if (length(w)>0) if (df.attr[w, 'element']=='g') {

        na <- xml_name(xml_children(chdn.all[[w]]))
        # Length of the reference group (g).
        len.r <- xml_length(chdn.all[[w]])-sum(na %in% c('a', 'title'))
        tit0 <- rep(df.attr[i, 'feature'], len.r); tit0 <- paste0(tit0, '__', seq_along(tit0)); tit <- c(tit, tit0)
        id0 <- rep(df.attr[i, 'id'], len.r); id.all <- c(id.all, id0)

      } else { tit <- c(tit, df.attr[i, 'feature']); id.all <- c(id.all, df.attr[i, 'id']) }

    } else { tit <- c(tit, df.attr[i, 'feature']); id.all <- c(id.all, df.attr[i, 'id']) }

  }; tis.path <- gsub("__\\d+$", "", tit)
 # style <- 'fill:#46e8e8;fill-opacity:1;stroke:#000000;stroke-width:3;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1' # 'fill' is not necessary. In Inkscape, resizing a "group" causes "matrix" in "transform" (relative positions) attribute, and this can lead to related polygons uncolored in the spatial heatmaps. Solution: ungroup and regroup to get rid of transforms and get absolute positions.
  # Change 'style' of all polygons. Since in SVG code, if no fill in style, no fill in ".ps.xml", so is the stroke.
  # "stroke" >= 0.51 px always introduces coordinates in .ps.xml, no matter "fill" is "none" or not. If "stroke" < 0.5 px, even though "fill" is not "none" there is no coordinates in ps.xml. E.g. irregular paths of dots. 
  # xml_set_attr(chdn.out, 'style', style); xml_set_attr(chdn.ply, 'style', style)  
  # xml_set_attr(out, 'style', style); xml_set_attr(ply, 'style', style)  
  # Export internal SVG.
  tmp <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE);
  svg.inter <- paste0(tmp, '/internal.svg')
  if (grepl("~", svg.inter)) svg.inter <- normalizePath(svg.inter, winslash="/", mustWork=FALSE)
  write_xml(doc, file=svg.inter)
  
  # SVG file conversion.
  # In Inkscape, the path of letter "U" is combined by a long and short path: "M -145.74174, ..., 182.19339 Z M -141.22026,175.30862 Z". After broken apart, the short path "M -141.22026,175.30862 Z" still exists, but disappears after "rsvg_ps" and also "PostScriptTrace". As a result, the corresponding coordinates are gone. Since the short path introduces an id, the colors are shifted in SHM. Thus letters/text should be not be broken apart by "path_br".
  rsvg_ps(svg.inter, file=sub("svg$", "ps", svg.inter)) # Only the paths inside canvas of SVG are valid.
  p1 <- sub("svg$", "ps", svg.inter); p2 <- paste0(sub("svg$", "ps", svg.inter), ".xml"); PostScriptTrace(p1, p2) 
  chdn1 <- xml_children(read_xml(p2)) # Use internal svg to get coordinates.
  # Detect groups that use relative coordinates ("transform", "matrix" in Inkscape.), which leads to some plygons missing in ".ps.xml" file.
  # EBI SVG, if the outline shapes and tissue shapes are separate, they must be in two layers NOT two groups. Otherwise, 'fill' and 'stroke' in '.ps.xml' can be  messy.
  fil.stk <- xml_attr(chdn1[-length(chdn1)], 'type'); tab <- table(fil.stk)
  fil <- tab['fill']; stk <- tab['stroke']; tit.len <- length(tit)
  # 'eofill' is also accounted for.
  if (all(c('fill', 'stroke') %in% names(tab))) {
    if (fil==stk|stk==tit.len) nodeset <- chdn1[which(xml_attr(chdn1, 'type')=='stroke')] else if (fil==tit.len) nodeset <- chdn1[which(xml_attr(chdn1, 'type')=='fill')] else if (fil!=stk & ceiling(sum(tab)/2)==tit.len) {
      nodeset <- chdn1[seq(1, sum(tab), by=2)] 
    }
  } else if (sum(c('fill', 'stroke') %in% names(tab))==1) {
    nodeset <- chdn1[seq(1, sum(tab), by=2)]
  }
  w <- which(fil.stk=='fill')%%2==0
  # In Ubuntu 18.04.4, there are double 'fill', no 'stroke'.
  if (all(c('stroke', 'fill') %in% names(tab))) if (any(w) & tab['fill']!=tab['stroke']) { 
 
    w1 <- which(w)[1] # Index of wrong path.
    tis.wrg <- paste0(id.all[c(w1-1, w1)], collapse=';') # Wrong path and related group.
    cat('\n'); cat(paste0("Potential error detected in these elements: '", tis.wrg, "'! If they are groups, please remove the 'transform' attribute with a 'matrix' value by ungrouping and regrouping the respective groups in Inkscape. If individual paths, consider deleting them in Inkscape. Otherwise, colors in spatial heatmap might be shifted!"), '\n') 

  }
  # Assign stroke to every path including paths inside groups, since in cases of many shapes it is time-consuming to check every stroke against all tissues in the coordinate data frame to assign strokes after the coordinate data frame is done.
  stroke.w <- df.attr$stroke; names(stroke.w) <- df.attr$feature
  # Extract coordinates for each path independently. If many paths are included in an SVG, coordinates or fill/stroke order errors may arise if all coordinates are extracted as a whole. If extracted independently for each path, a little errors are raised, but the speed is slow. The coordinates (e.g. all y coord) of a shape may be slightly different (usually after decimal points) with extracted in a whole (i.e. use node), but no difference is observed on ggplot-plotted shapes extracted from the two contrasting methods.
  # Test if some paths/dots are missing: identical(as.vector(unique(df$tissue)), unique(tit))
  if (length(nodeset)==tit.len) df <- xy0(nodeset, tit, stroke.w, cores) else {
    cat('Extracting coordinates for each shape independently, which is slow ... \n')
    df.out <- cord_parent(svg.path, 'out', feature, stroke.w, cores)
    if (is(df.out, 'character')) return(df.out)
    df.ply <- cord_parent(svg.path, 'ply', feature, stroke.w, cores)
    if (is(df.ply, 'character')) return(df.ply)
    df <- rbind(df.out$df, df.ply$df); id.no <- c(df.out$ids, df.ply$ids)
    if (!is.null(id.no)) { cat('No coordinates were extracted for these element(s):', id.no, '!\n') }
  }
  # return("The 'transform' attribute with a 'matrix' value is not allowed in groups! Please remove them by ungrouping and regrouping the related groups in Inkscape if exist!") 
  # Get coordinates from '.ps.xml'.
  # nodeset <- chdn1[which(xml_attr(chdn1, 'type')=='stroke')]
  # In Ubuntu 18.04.4, there are double 'fill', no 'stroke', thus only use 'fill' of odd numbers.
  #if (xml_length(nodeset)[1]==0) {
  #  if (sum(tab) %% 2==0) nodeset <- chdn1[seq(1, sum(tab), by=2)] else return('Relative coordinates detected in aSVG file!')
  #}
  #if (length(tit)!=length(nodeset)) return('some shape(s) are missing!')

  # Move non-matching tissues on top of matching tissues in the data frame.
  idx.match <- sub('__\\d+$', '', df$tissue) %in% feature
  # In geom_polygon, the order to plot tissues is the factor level. If a tissue is the 1st according to factor level but is last in the coordinate data frame, it will be plotted first, and the 2nd tissue in the level can cover it if all tissues are colored.
  df <- rbind(df[!idx.match, ], df[idx.match, ])
  # Place some shapes on the top layer on purpose.
  idx.top <- grepl('_TOP$|_TOP__\\d+$', df$tissue)
  df <- rbind(df[!idx.top, ], df[idx.top, ])
  df$tissue <- factor(df$tissue, levels=unique(df$tissue))
  # Each entry in tis.path is represented by many x-y pairs in coordinate, and tissues in coord are tissues in tis.path appended '__\\d+$'.
  # Update tis.path.
  tis.path <- sub('__\\d+$', '', unique(df$tissue))
  fil.cols <- df.attr$color; names(fil.cols) <- df.attr$feature
  w.h <- c(max(df$x) - min(df$x), max(df$y) - min(df$y)); aspect.r <- w.h[1]/w.h[2]; names(aspect.r) <- NULL
  names(w.h) <- c('width', 'height')
  # tis.path=sub('_\\d+$', '', tit) introduces a potential bug, since the original single-path tissues can have '_\\d+$' pattern. Solution: in upstream append '__1', '__2', ... to the paths in a group.
  lis <- list(df=df, tis.path=tis.path, fil.cols=fil.cols, w.h = w.h, aspect.r = aspect.r, df.attr=df.attr); return(lis)

}

#' Extract children, id, element name from outline and tissue layer#' @param doc The document of SVG
#' @keywords Internal
#' @noRd

out_ply <- function(doc) {
  len <- xml_length(doc); out <- xml_children(doc)[[len-1]]
  ply <- xml_children(doc)[[len]]
  # If out is not a group, it is assigned an empty node.
  if (xml_name(out)!='g') { xml_add_child(out, 'empty', .where=0); out1 <- xml_children(out)[[1]]; xml_remove(xml_children(out)[[1]], free=FALSE); out <- out1 }
  # If ply is not a group, it is assigned an empty node.
  if (xml_name(ply)!='g') { xml_add_child(ply, 'empty', .where=0); ply1 <- xml_children(ply)[[1]]; xml_remove(xml_children(ply)[[1]], free=FALSE); ply <- ply1 }
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  id.all <- make.names(c(xml_attr(chdn.out, 'id'), xml_attr(chdn.ply, 'id')))
  na.all <- c(xml_name(chdn.out), xml_name(chdn.ply))
  return(list(out=out, ply=ply, chdn.out=chdn.out, chdn.ply=chdn.ply, id.all=id.all, na.all=na.all))
}

#' If the node is a single path or group, extract the tissue name from title or id
#' @param node A single node or group.
#' @keywords Internal
#' @noRd

tit_id <- function(node) {
  cld0 <- xml_children(node); na0 <- xml_name(cld0)
  idx <- grep('^title$', na0, ignore.case=TRUE)[1]
  if (!is.na(idx)) { tis <- xml_text(cld0[idx]); if (tis=='') tis <- xml_attr(node, 'id') } else tis <- xml_attr(node, 'id')
  return(tis)
}


#' Extract coordinates for a single node, or a group containing a single node, or a use node in the context of original SVG
#' @param doc The document of SVG. If use is FALSE, it contains empty outline and tissue layers.
#' @param parent The outline or tissue layer
#' @param node A single path or a group containing a single path
#' @param tis The title/id of the node. If the node is from a group, the group title/id is appended "_\\d+", and the "tis" is one of the appended title/ids.
#' @param use If TRUE, doc only contains the reference and use nodes. If the reference is a group, "tis" is title/id appended with "_\\d+". In the ps.xml file, only half of the second harf (use node) is extracted.
#' @param stroke.w A vector of all stroke widths extracted from the aSVG file, which is named by features in the aSVG file.
#' @param cores The number of CPU cores.
#' @keywords Internal
#' @noRd

xy <- function(doc, parent, node, tis, use=FALSE, stroke.w, cores) {        
  options(stringsAsFactors=FALSE)
  style <- 'fill:#46e8e8;fill-opacity:1;stroke:#000000;stroke-width:3;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1'
  # SVG file containing a single path.
  tmp <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE);
  svg.inter <- paste0(tmp, '/shm.svg')
  ps.inter <- paste0(tmp, '/shm.ps')
  xml.inter <- paste0(tmp, '/shm.ps.xml')
  # xml.inter <- paste0(tmp, '/', tis, '_shm.ps.xml')
  if (use==FALSE) {xml_set_attr(node, 'style', style); xml_add_child(parent, node) } 
  write_xml(doc, svg.inter)
  # Extract coordindates.
  rsvg_ps(svg.inter, ps.inter); # cat(xml.inter, '\n')
  PostScriptTrace(ps.inter, xml.inter)
  cld <- xml_children(read_xml(xml.inter))
  cnt <- as.numeric(xml_attr(cld[length(cld)], 'count'))
  if (use==TRUE) {
    # cat("Extracting coordinates for element 'use':", tis[1], '.. \n')
    # Reference and use nodes should generate the same coodinates.
    if (cnt %% 2 != 0) { cat(tis, ': problematic coordinates detected!\n'); return('no') }
    # The cooridnates at odd number.
    cld <- cld[seq(cnt/2 + 1, cnt, by=2)]; return(xy0(cld, tis, stroke.w, cores))
  } else {
    if (cnt==0) { cat('\n'); cat(tis, ': no coordinates detected!\n'); return('no') } else { cld <- cld[1]
    # return('yes') 
    }; return(xy0(cld, tis, stroke.w, cores))
  }
}

#' Extract coordinates for a nodeset
#' @param nodeset Node sets generated by xml_children.
#' @param tit.all The title/id corresponds to each node.
#' @param stroke.w A vector of line sizes corresponding to each tissue, which are named by tissue names extracted from the SVG file.
#' @param cores The number of CPU cores.
#' @keywords Internal
#' @noRd

xy0 <- function(nodeset, tit.all, stroke.w, cores) {
  # Cut node sets into chunks.
  idxs <- seq_along(tit.all); n <- ceiling(length(tit.all)/cores)
  chunk <- split(idxs, ceiling(idxs/n))
  # Extract coordinates of all tissues in each chunk and combine the all extracted coordinates.
  xy.mat <- setDF(rbindlist(mclapply(chunk,
    function(vec) {
      # Extract coordinates of all tissues in each chunk.
      # The class nodeset is like a list. xml_children(nodeset)[-1]: extract and combine the children of each node, [-1] removes the first child of each node.
      lis0 <- xml_attrs(xml_children(nodeset[vec])[-1], c('x','y'))
      mat0 <- as.data.frame(do.call("rbind", lis0)); return(mat0)
    }, mc.cores=cores)))
  xy.mat <- as.data.frame(apply(xy.mat, 2, as.numeric))
  # Vectorize tissue and line size and add them to the coordinate data frame.
  lens <- xml_length(nodeset)-1
  widths <- stroke.w[sub('__\\d+$', '', tit.all)]
  xy.mat$tissue <- rep(tit.all, lens)
  xy.mat$line.size <- rep(widths, lens); return(xy.mat)
}


#' Extract coordinates for each path in outline or tissue layer
#' @param doc The document of SVG containing all nodes.
#' @param out Outline layer.
#' @param ply Tissue layer.
#' @param parent The outline or tissue layer.
#' @param stroke.w A vector of all stroke widths extracted from the aSVG file, which is named by features in the aSVG file.
#' @param cores The number of CPU cores.
#' @keywords Internal
#' @noRd
#' @importFrom xml2 xml_new_root 

cord <- function(doc, out, ply, parent, stroke.w, cores) {
  if (xml_length(parent)==0) return(list(df=data.frame(), tits=NULL, ids=NULL))
  # The children of out cannot be inserted to ply to extract coordinates, vice versa, since coordinates may not be generated.
  doc0 <- xml_new_root(doc, .copy=TRUE); chdn0 <- xml_children(parent)
  # Each child of parent is inserted back to parent to extract coordinates.
  # tits: all tissues in original SVG. ids: all ids of paths not generating coordinates.
  cat('Extracting coordinates for:')
  df <- data.frame(); tits <- ids <- NULL; for (i in seq_along(chdn0)) { 
    nod0 <- chdn0[[i]]; na0 <- xml_name(nod0); id0 <- make.names(xml_attr(nod0, 'id'))
    tis <- make.names(tit_id(nod0)); tits <- c(tits, tis)
    if (na0 %in% c('a', 'title', 'text', 'flowRoot')) next
    cat(' ', tis)
    if (na0=='use') { # Use node
      # Make a copy of doc, since original doc might be emptied if the target node is not use.
      doc1 <- xml_new_root(doc0, .copy=TRUE)
      out.ply1 <- out_ply(doc1); chdn.all1 <- c(out.ply1$chdn.out, out.ply1$chdn.ply)
      id.all1 <- out.ply1$id.all; use.node <- chdn.all1[id.all1 %in% id0][[1]]
      len <- use(chdn.all1, id.all1, use.node) # Number of the paths in reference of a group.
      if (!is.null(len)) {
        if (len=='no') { ids <- c(ids, id0); next }
        if (is.numeric(len) & len > 0) tis <- paste0(tis, '__', seq_len(len))
      }
      df0 <- xy(doc=doc1, tis=tis, use=TRUE, stroke.w=stroke.w, cores=cores)
      if (!is(df0, 'data.frame')) { if (df0=='no') ids <- c(ids, id0) } else df <- rbind(df, df0)
    } else if (na0!='g'){ # If the child is not a group, inserted directly.
      xml_remove(xml_children(out)); xml_remove(xml_children(ply)) # Clean all parent children each time, since the children are accumulated otherwise.
      df0 <- xy(doc=doc, parent=parent, node=nod0, tis=tis, use=FALSE, stroke.w=stroke.w, cores=cores)
      if (!is(df0, 'data.frame')) { if (df0=='no') ids <- c(ids, id0) } else df <- rbind(df, df0)
    } else if (na0=='g') { # If the child is a group, each child is inserted back to the group and the group containing a single child is inserted back to parent.
      cld0 <- xml_children(nod0); nas0 <- xml_name(cld0)
      cld1 <- cld0[!nas0 %in% c('a', 'title', 'text', 'flowRoot', 'use')]
      if (length(cld1)>1) { # Each child of the group is inserted back to the group.
        tis.all <- paste0(tis, '__', seq_along(cld1))
        for (j in seq_along(cld1)) {
          xml_remove(xml_children(out)); xml_remove(xml_children(ply))
          xml_remove(xml_children(nod0)); xml_add_child(nod0, cld1[[j]])
          df0 <- xy(doc=doc, parent=parent, node=nod0, tis=tis.all[j], use=FALSE, stroke.w=stroke.w, cores=cores)
          if (!is(df0, 'data.frame')) { if (df0=='no') ids <- c(ids, xml_attr(cld1[[j]], 'id')) } else df <- rbind(df, df0)
        }
      } else if (length(cld1)==1) { # If the group contains only one child.
        xml_remove(xml_children(out)); xml_remove(xml_children(ply))
        xml_remove(xml_children(nod0))
        xml_add_child(nod0, cld1[[1]])
        df0 <- xy(doc=doc, parent=parent, node=nod0, tis=tis, use=FALSE, stroke.w=stroke.w, cores=cores)
        if (!is(df0, 'data.frame')) { if (df0=='no') ids <- c(ids, xml_attr(cld1[[j]], 'id')) } else df <- rbind(df, df0)
      } 
    }
  }; cat('\n'); return(list(df=df, tits=tits, ids=ids))
}

#' Keep only a reference and a use node in outline or tissue layer
#' @param chdn.all A list of all children nodes.
#' @param id.all A vector of all ids.
#' @param use.node The target use node.
#' @keywords Internal
#' @noRd

use <- function(chdn.all, id.all, use.node) {
  na.all <- NULL
  style <- 'fill:#46e8e8;fill-opacity:1;stroke:#000000;stroke-width:3;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1'
  href <- xml_attr(use.node, 'href')
  xml_set_attr(use.node, 'style', style)
  # Search for the reference node
  id <- xml_attr(use.node, 'id'); idx <- which(id.all==id) 
  w1 <- which(paste0('#', id.all) %in% href)
  if (length(w1)>0) { # The reference is not in a group.
    ref <- chdn.all[[w1[1]]]; xml_set_attr(ref, 'style', style)
    if (xml_name(ref)=='g') {
      xml_set_attr(xml_children(ref), 'style', style)
    }; for (i in chdn.all[-c(w1, idx)]) xml_remove(i) # Remove other nodes.
    nas0 <- xml_name(xml_children(ref))
    # Length of valid children in reference group.
    len <- sum(!nas0 %in% c('a', 'title', 'text', 'use', 'flowRoot'))
    if (len > 0) return(len)
  } else { # The reference is in a group.
    g <- chdn.all[na.all=='g']; if (length(g)==0) { cat('No reference element is detected for use', id, '\n'); return('no') }
    for (k in seq_along(g)) { # Search for the reference node in each group.
      g0 <- g[[k]]; cld.g <- xml_children(g0)
      w <- which(paste0('#', xml_attr(cld.g, 'id')) %in% href)
      if (sum(w)==0) next; xml_remove(cld.g[-w])
      w2 <- which(id.all==xml_attr(g0, 'id'))
      for (i in chdn.all[-c(w2, idx)]) xml_remove(i)
      ref <- chdn.all[[w2[1]]]; xml_set_attr(ref, 'style', style)
      xml_set_attr(xml_children(ref), 'style', style)
    }
  }
}


#' Extract coordinates for each path independently
#' @param svg.path The SVG file path.
#' @param parent The outline or tissue layer, where coordinates of each path will be extracted independently.
#' @param feature A character vector of features/samples extracted from the data. If some of the input features are duplicated in SVG file, then a reminder message is returned.
#' @param stroke.w A vector of all stroke widths extracted from the aSVG file, which is named by features in the aSVG file.
#' @param cores The number of CPU cores.
#' @keywords Internal
#' @noRd
#' @importFrom xml2 xml_text

cord_parent <- function(svg.path, parent, feature, stroke.w, cores) {
  options(stringsAsFactors=FALSE)
  doc <- read_xml(svg.path); spa <- xml_attr(doc, 'space')
  if (!is.na(spa)) if (spa=='preserve') xml_set_attr(doc, 'xml:space', 'default')
  svg.attr <- svg_attr(doc, feature=feature); if (is(svg.attr, 'character')) return(svg.attr)
  out.ply <- out_ply(doc); out <- out.ply$out; ply <- out.ply$ply
  chdn.out <- out.ply$chdn.out; chdn.ply <- out.ply$chdn.ply
  chdn.all <- c(chdn.out, chdn.ply)
  id.all <- out.ply$id.all; na.all <- out.ply$na.all

  # Append 1, 2, 3, ... to duplicated titles.
  tit.id <- c(vapply(chdn.out, tit_id, character(1)), vapply(chdn.ply, tit_id, character(1)))
  dup.tit <- tit.id[duplicated(tit.id)]
  if (length(dup.tit)>0) cat('Duplicated title text detected:', dup.tit, '\n')
  dup <- unique(tit.id[duplicated(tit.id)])
  for (i in dup) {
    w <- tit.id %in% dup; dup0 <- tit.id[w]; cld0 <- chdn.all[w] 
    dup1 <- paste0(dup0, seq_along(dup0))
    # Ids are not duplicated, so only update titles.
    for (j in seq_along(dup1)) {
      cld.all <- xml_children(cld0[[j]]); nas <- xml_name(cld.all)
      xml2::xml_text(cld.all[nas=='title']) <- dup1[j]
    }
  }
  if (parent=='out') df.par <- cord(doc, out, ply, out, stroke.w, cores) else if (parent=='ply') df.par <- cord(doc, out, ply, ply, stroke.w, cores)
  return(df.par)
}



