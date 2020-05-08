#' Break a Combined Path into a Group or Siblings
#' 
#' The combined paths introduces connecting lines between the last point of one polygon and the first point of next polygon in the spatial heatmaps. Therefore, they should be broken apart into a group or siblings. This function checks if the input node is a combined path internally.
#' @param node An object of class "xml_node" without children nodes.
#' @param g Logical, TRUE or FALSE. Default is TRUE. If TRUE the combined path is broken into a group. Otherwise, as siblings.
#' @return Nothing is returned. The broken paths are updated in the root.

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2

#' @importFrom xml2 xml_name xml_attr xml_children xml_remove xml_add_child xml_set_attr

  path_br <- function(node, g=TRUE) {

    na <- xml_name(node); if (na!='g') {

      d <- xml_attr(node, 'd') 
      if (grepl('m ', d)) return('Please use absolute coordinates for all paths!')
      if (grepl('Z M', d)) {

        z <- paste0(strsplit(d, 'Z')[[1]], 'Z')
        ids <- paste0(xml_attr(node, 'id'), '_', seq_along(z))
        # Make node empty.
        xml_attr(node, 'd') <- NA
        
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
          xml_name(node) <- 'g'; xml_attr(node, 'd') <- NULL
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
