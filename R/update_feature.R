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
#' @importFrom xml2 xml_children xml_length xml_attr xml_set_attr

update_feature <- function(feature, dir) {

  dir.check <- !is.null(dir) 
  if (dir.check) dir.check <- !(is.na(dir)) else stop("\'dir\' is not valid!") 
  if (dir.check) { dir.check <- dir.exists(dir); if (!dir.check) stop("\'dir\' is not valid!") } else stop("\'dir\'is not valid!")
  feature[, 1] <- as.character(feature[, 1])
  feature$index <- as.numeric(feature$index)
  svgs.na <- unique(feature$SVG)
  for (i in svgs.na) {

    df0 <- subset(feature, SVG==i); dup <- duplicated(df0[, 1])
    if (any(dup)) stop(paste0("Duplicated feature \'", paste0(df0[, 1][dup], collapse=', '), "\' detected in ", i, "!"))
    path.in <- paste0(dir, '/', i)  
    doc <- read_xml(path.in); chdn <- xml_children(doc)
    ply <- chdn[[xml_length(doc)]]; chdn1 <- xml_children(ply)
    cat(paste0('Setting \'', paste0(df0[, 1], collapse=', '), '\' in ', path.in, '\n'))
    xml_set_attr(chdn1[df0$index], 'id', df0[, 1])
    write_xml(doc, file=path.in)
  
  }

}

