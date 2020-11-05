#' Break Children in a Node into a Group or Siblings
#' 
#' This function applies to the last two parent nodes in SVG file. Each child node in the input parent node is checked for combined path. The function \link{path_br} is used to break the combined paths.
#' @param node.parent An object of class "xml_node" with children, i.e. one of the last two nodes in SVG file, outline and tissue layer respectively.
#' @return Nothing is returned. The broken paths are updated in the root.

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @noRd

#' @references
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2

#' @importFrom xml2 xml_length xml_children xml_name xml_attr

path_br_all <- function(node.parent) {

  len <- xml_length(node.parent); chdn <- xml_children(node.parent)
  for (i in seq_len(len)) {

    nod0 <- chdn[[i]]; na <- xml_name(nod0); id <- xml_attr(nod0, 'id')
      # "xml_remove" does not affect "chdn" since it is "xml_nodeset", but affects "node.parent" and "doc" since they are "xml_node". There is no need to return "node.parent", since even in nested functions oprations on node are still effective on "doc" and "node.parent" such as "xml_remove" in "rm_dot".
      if (na %in% c('text', 'flowRoot')|(na=='g' & xml_length(nod0)==0)) { xml_remove(nod0); cat('Removing "text", "flowRoot" nodes, and empty "g" nodes:', id, '! \n'); next } # These nodes do not have coordinates after "PostScriptTrace" and thus lead to color shifting in SHM.
      rm_dot <- function(node) {

        d <- xml_attr(node, 'd'); cords <- grep('\\d+', strsplit(d, ',|-| ')[[1]], value=TRUE)
        # Dots/tiny paths may only introduce "stoke" no "fill".
        if (length(cords)<4) { xml_remove(node); cat('Removing tiny paths and dots:', xml_attr(node, 'id'), '! \n'); return('yes') } else return('no')
        
      } 
      if (na!='g') { dot <- rm_dot(nod0); if (dot=='yes') next }
      # 'a' node and 'text' path/group are ignored.
      if (na=='a'|grepl('^text', id, ignore.case=TRUE)) next
      if (na!='g') path_br(nod0, g=TRUE) else {

        nod0.chl <- xml_children(nod0); nas <- xml_name(nod0.chl)
        if ('g' %in% nas) return(paste0('Nested group detected in ', xml_attr(nod0, 'id'), '!'))
        if ('use' %in% nas) return(paste0('The "use" node is detected in ', xml_attr(nod0, 'id'), '!'))
        for (j in seq_along(nod0.chl)) {

          nod1 <- nod0.chl[[j]]; dot <- rm_dot(nod1)
          if (dot=='yes') next; if (xml_name(nod1)=='a') next
          path_br(nod1, g=FALSE)

        }

      }

  }

}


