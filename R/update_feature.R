#' Update SVG Images with Provided Features
#' 
#' Successful spatial heatmap plotting requires the feature identifiers of interest are identical between the data matrix and SVG image. This function is designed to replace existing features in SVG images with user-provided features. Note this function treats the first column in the feature data frame as user-provided features, so custom features must be the first column.

#' @param feature The data frame returned by \code{\link{return_feature}} with the user-provided features being the first column.
#' @param dir The directory where the SVG images are available. It should be same with "dir" in \code{\link{return_feature}}.
#' @return Nothing is returned. The SVG images of interest in "dir" are updated with new features, and are ready for plotting spatial heatmaps.
#' @examples
#' feature.df <- return_feature(feature='frontal cortex', species='homo sapiens', keywords.all=TRUE, desc=FALSE, return.all=FALSE, dir='.', remote=TRUE)
#' ft <- c('frontal.cortex', 'prefrontal.cortex', 'prefrontal.cortex', 'frontal.cortex', 'frontal.cortex', 'prefrontal.cortex')
#' # Add new features to the first column.
#' feature.df.new <- cbind(featureNew=ft, feature.df) 
#' update_feature(feature=ft2, dir='.')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2

#' @export update_feature
#' @importFrom xml2 read_xml xml_children xml_length xml_attr xml_set_attr xml_name xml_add_child xml_set_text write_xml

update_feature <- function(feature, dir) {

  dir.check <- !is.null(dir) 
  if (dir.check) dir.check <- !(is.na(dir)) else stop("\'dir\' is not valid!") 
  if (dir.check) { dir.check <- dir.exists(dir); if (!dir.check) stop("\'dir\' is not valid!") } else stop("\'dir\'is not valid!")
  feature[, 1] <- as.character(feature[, 1])
  feature$index <- as.numeric(feature$index)
  feature$index1 <- as.numeric(feature$index1)
  svgs.na <- unique(feature$SVG)
  for (i in svgs.na) {

    df0 <- subset(feature, SVG==i); dup <- duplicated(df0[, 1])
    if (any(dup)) stop(paste0("Duplicated feature \'", paste0(df0[, 1][dup], collapse=', '), "\' detected in ", i, "!"))
    path.in <- paste0(dir, '/', i); doc <- read_xml(path.in); len <- xml_length(doc)
    out <- xml_children(doc)[[len-1]]; ply <- xml_children(doc)[[len]]
    # Update features for two parent layers. 
    for (k in list(out, ply)) { 

      # Check if the parent needs to be updated. 
      df1 <- subset(df0, parent==xml_attr(k, 'id')); if (nrow(df1)==0) next
      cat(paste0('Setting titles \'', paste0(df1[, 1], collapse=', '), '\' in ', path.in, '\n'))
      # Update features by order.
      for (j in seq_len(nrow(df1))) {
    
        nod0 <- xml_children(k)[[df1$index1[j]]]
        chil0 <- xml_children(nod0); nas <- xml_name(chil0)
        # Add a 'title' node if the feature to update does not have it, and update 'chil0' and 'nas'.
        if (!('title' %in% nas)) { xml_add_child(nod0, 'title', .where=0); xml_set_attr(xml_children(nod0)[[1]], 'id', df1[j, 1]); chil0 <- xml_children(nod0); nas <- c('title', nas)
        }; xml_set_text(chil0[[which(nas=='title')]], df1[j, 1])

      }
  
    }; write_xml(doc, file=path.in)

  }

}

