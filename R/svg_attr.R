#' Extract attributes in an SVG file
#' 
#' This function extract basic attributes of nodes in an SVG file. The 'a' nodes are not deleted.

#' @param doc An object of class "xml_document", i.e. the imported SVG file in R.
#' @param feature A character vector of features/samples extracted from the data. If some of the input features are duplicated in SVG file, then a reminder message is returned.

#' @return A 3-component list of attribute data frame, the outline node, and the tissue node. 

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2

#' @importFrom xml2 xml_length xml_children xml_name xml_attr xml_remove xml_text

svg_attr <- function(doc, feature) {

  options(stringsAsFactors=FALSE)
  len <- xml_length(doc); out <- xml_children(doc)[[len-1]]; ply <- xml_children(doc)[[len]]
  # Break combined path to a group or siblings.
  path_br_all(out); path_br_all(ply)


  # If out is not a group, it is assigned an empty node.
  if (xml_name(out)!='g') { xml_add_child(out, 'empty', .where=0); out1 <- xml_children(out)[[1]]; xml_remove(xml_children(out)[[1]], free=FALSE); out <- out1 }
  # If ply is not a group, it is assigned an empty node.
  if (xml_name(ply)!='g') { xml_add_child(ply, 'empty', .where=0); ply1 <- xml_children(ply)[[1]]; xml_remove(xml_children(ply)[[1]], free=FALSE); ply <- ply1 }
  chdn.out <- xml_children(out); chdn.ply <- xml_children(ply)
  
  ## Exrtact basic attributes into a data frame.
  idx <- seq_len(length(chdn.out)+length(chdn.ply))
  idx1 <- c(seq_len(length(chdn.out)), seq_len(length(chdn.ply)))
  parent <- c(rep(xml_attr(out, 'id'), length(chdn.out)), rep(xml_attr(ply, 'id'), length(chdn.ply)))
  nas <- c(xml_name(chdn.out), xml_name(chdn.ply))
  ids <- make.names(c(xml_attr(chdn.out, 'id'), xml_attr(chdn.ply, 'id')))
  if (any(duplicated(ids))) return(paste0('Duplicated node ids detected: ', paste0(ids[duplicated(ids)], collapse=' '), '!'))
  title <- c(xml_text(chdn.out), xml_text(chdn.ply)) # Use original names, no 'make.names'.
  w <- which(title=='X'|title==''); title[w] <- ids[w]
  dup <- duplicated(title); if (any(dup)) {
    
    if (length(intersect(title[dup], feature))>0) return(paste0('Duplicated title text detected: ', paste0(title[dup], collapse=' '), '!')) else {

      w <- title %in% title[dup]
      title[w] <- paste0(title[w], seq_len(sum(w)))

    }
  
  }
  # Style inside groups are ignored. 
  sty <- c(xml_attr(chdn.out, 'style'), xml_attr(chdn.ply, 'style'))
  sty[!grepl('fill:', sty)] <- 'none'
  w1 <- grepl(';', sty); st <- sty[w1]; st <- strsplit(st, ';')
  st1 <- NULL; for (i in st) { st1 <- c(st1, i[grepl('fill:', i)]) }; sty[w1] <- st1
  # If only keep part of the string, the pattern should cover everything in the string, e.g. the '.*' on both ends.
  fil.cols <- gsub('.*(fill:)(.*).*', '\\2', sty)
  df.attr <- data.frame(index=idx, index1=idx1, parent=parent, name=nas, id=ids, title=title, fil.cols=fil.cols)
  df.attr <- subset(df.attr, name!='a')    
  return(list(df.attr=df.attr, out=out, ply=ply))

}
