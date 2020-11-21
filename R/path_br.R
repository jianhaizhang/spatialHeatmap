#' Break a Combined Path into a Group or Siblings
#' 
#' The combined paths introduces connecting lines between the last point of one polygon and the first point of next polygon in the spatial heatmaps. Therefore, they should be broken apart into a group or siblings. This function checks if the input node is a combined path internally and breaks them apart if existing.
#' @param node An object of class "xml_node" without children nodes.
#' @param g Logical, TRUE or FALSE. Default is TRUE. If TRUE the combined path is broken into a group. Otherwise, as siblings.
#' @return Nothing is returned. The broken paths are updated in the root.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @noRd

#' @references
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2


#' @importFrom xml2 xml_attr xml_add_sibling xml_name xml_children xml_remove xml_add_child xml_set_attr

path_br <- function(node, g=TRUE) {

 na <- xml2::xml_name(node); if (na!='g') {
   id <- xml_attr(node, 'id')
   if (na=='a') { xml_remove(node); cat('Element "a" is removed:', id, '!\n'); return() }
   if (na %in% c('a', 'g', 'use', 'title', 'ellipse', 'rect')) return() 
   # Delete the node if tiny/dot path.
   dot <- rm_dot(node); if (dot=='yes') return() 

   d <- xml2::xml_attr(node, 'd') 
   if (grepl('m ', d)) return('Please use absolute coordinates for all paths!')
   if (grepl('Z M', d)) {
 
     z <- paste0(strsplit(d, 'Z')[[1]], 'Z')
     # Delete tiny/dot paths.
     for (i in seq_along(z)) {
       cords <- grep('\\d+', strsplit(z[i], ',|-| ')[[1]], value=TRUE)
       if (length(cords)<4) { z[i] <- NA; cat('Extracted tiny path from', id, 'is removed! \n') }
     }; z <- z[!is.na(z)]
       # If only one path is left, update the path and return.
       if (length(z)==1) { xml_attr(node, 'd'); return() }
       ids <- paste0(id, '_', seq_along(z))
       # Make node empty.
       xml2::xml_attr(node, 'd') <- NA
       # Break the combined path to a group.
      if (g==TRUE) {
        
        # Isolate 'title' node.
        na.chil <- xml_name(xml_children(node))
        w <- which(na.chil=='title')
        if (length(w)>0) { tit <- xml_children(node)[[w]]; xml_remove(xml_children(node)[w], free=FALSE) }

        # Add the empty node to itself as the first child.
        xml_add_child(node, node)
        # Copy the first child for length(z)-1 times.
        node1 <- xml_children(node)[[1]]
        for (j in seq_len(length(z)-1)) { xml_add_child(node, node1) }
        node.chl <- xml_children(node) # Function applies to 'nodeset' recusively. 
        # Set d and id for all childrend of node.
        xml_set_attr(node.chl, 'd', z)
        xml_set_attr(node.chl, 'id', ids)  
        # Name node 'g'.
        xml2::xml_name(node) <- 'g'; xml2::xml_attr(node, 'd') <- NULL
        if (length(w)>0) xml_add_child(node, tit, .where=0)

      } else {

        for (j in seq_along(z)[-1]) {

          # Copy node as its own siblings.
          xml_set_attr(node, 'd', z[j]); xml_set_attr(node, 'id', ids[j])
          xml_add_sibling(node, node, 'after')
          # Change 'd' in node at last.  
          xml_set_attr(node, 'd', z[1]); xml_set_attr(node, 'id', ids[1])

      }

    }

   }
    
  }

}

#' # Remove tiny/dot paths.
#' @keywords Internal
#' @noRd

rm_dot <- function(node) {
  d <- xml_attr(node, 'd'); if (!is.na(d)) cords <- grep('\\d+', strsplit(d, ',|-| ')[[1]], value=TRUE) else cords <- 0
  # Dots/tiny paths may only introduce "stoke" no "fill".
  if (length(cords)<4) { xml_remove(node); cat('Removing tiny paths and dots:', xml_attr(node, 'id'), '! \n'); return('yes') } else return('no')  
}
