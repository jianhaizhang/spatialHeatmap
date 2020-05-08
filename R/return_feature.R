#' Return SVG Images Containing Provided Features and Species
#' 
#' Successful spatial heatmap plotting requires the feature identifiers of interest are identical between the data matrix and SVG image. This function parses a collection of SVG images and returns existing features related to the provided keywords. If users want to use their custom feature identifiers, \code{\link{update_feature}} should be used. Otherwise the returned featrue identifiers should be used to replace the counterparts in the data matrix. Features denote tissues, cells, etc. 

#' @param feature A vector of keywords (case insentitive) of target feature(s), which is used to select SVG images from a collection. E.g. c("frontal", "cortex").
#' @param species A vector of keywords (case insentitive) of a target species, which is used to select SVG images from a collection. E.g. c("homo sapiens").
#' @param keywords.all Logical, TRUE or FALSE. Default is TRUE. If TRUE, every returned hit will contain all the provided keywords in "feature" and "species". Otherwise, every returned hit contain at least 1 keyword in "feature" and at least 1 keyword in "species".
#' @param remote Logical, FALSE or TRUE. Default is FALSE. If TRUE, the remote SVG repository "https://github.com/jianhaizhang/SVG_tutorial_file/tree/master/svg_repo" is used for query.
#' @param dir The directory where the SVG images are available. If "remote" is TRUE, the returned SVG images are saved in this directory. Note, in this case existing SVG images with identical names as returned ones are overwritten. If "remote" is FALSE, user-provided SVG images should be saved in this directory for query.
#' @param desc Logical, FALSE or TRUE. Default is FALSE. If TRUE, the feature descriptions from the package "rols" (Laurent Gatto 2019) are added. If too many features are returned, this process takes a long time.
#' @param match.only Logical, TRUE or FALSE. If TRUE (default), only features matching the feature keywords are returned. If FALSE, all features in the matching species are returned, and the matching features are listed on the top of the data frame. 
#' @param return.all Logical, FALSE or TRUE. Default is FALSE. If TRUE, all features together with all SVG images are returned, either remote or user-provided.

#' @return A feature data frame with columns corresponding to feature id, feature ontology id, feature is index in an SVG image, and/or feature description. 
#' @examples
#' feature.df <- return_feature(feature='frontal cortex', species='homo sapiens', keywords.all=TRUE, desc=FALSE, return.all=FALSE, dir='.', remote=TRUE)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Laurent Gatto (2019). rols: An R interface to the Ontology Lookup Service. R package version 2.14.0. http://lgatto.github.com/rols/
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2

#' @export return_feature
#' @importFrom xml2 read_xml
#' @importFrom rols term termDesc

return_feature <- function(feature, species, keywords.all=TRUE, remote=FALSE, dir=NULL, desc=FALSE, match.only=TRUE, return.all=FALSE) {

  options(stringsAsFactors=FALSE)
  dir.check <- !is.null(dir) 
  if (dir.check) dir.check <- !(is.na(dir)) else stop("\'dir\' is not valid!") 
  if (dir.check) { dir.check <- dir.exists(dir); if (!dir.check) stop("\'dir\' is not valid!") } else stop("\'dir\'is not valid!")

  # Parse and return features.
  ftr.return <- function(svgs, desc=desc) { 

    cat('Accessing features... \n')
    df <- NULL; for (i in svgs) {

      doc <- read_xml(i); df0 <- svg_attr(doc, feature=NULL)[['df.attr']]
      # Move ontology with NA or "NULL" to bottom.
      w.na <- which(is.na(df0$id)|df0$id=='NULL'|df0$id=='NA')
      if (length(w.na)>0) df0 <- rbind(df0[-w.na, ], df0[w.na, ])
      na <- strsplit(i, '/')[[1]]; na <- na[grep('.svg$', na)]
      cat(na); cat(', ')
      df0$SVG <- na
      df <- rbind(df, df0)

    }; cat('\n')
    cna <- colnames(df); colnames(df)[cna=='title'] <- 'feature'; df <- df[, c('feature', 'id', 'SVG', 'parent', 'index', 'index1')]

    if (desc==TRUE) {
  
      cat('Appending descriptions... \n')
      df$description <- NA; for (i in seq_len(nrow(df))) {

        ont <- df[i, 'id']; abbr <- tolower(sub('_.*', '', ont))
        trm <- tryCatch({ term(abbr, ont) }, error=function(e) { return(NA) })
        if (is(trm, 'Term')) { des <- termDesc(trm); if (!is.null(des)) df[i, 'description'] <- termDesc(trm) }
  
      }

    }; return(df)

  }

  if (remote==TRUE) {
  
    cat('Downloading SVG images... \n')
    tmp <- tempdir(check=TRUE); tmp1 <- paste0(tempdir(), '/git.zip')
    tmp2 <- paste0(tmp, '/git'); if (!dir.exists(tmp2)) dir.create(tmp2)
    download.file('https://github.com/ebi-gene-expression-group/anatomogram/archive/master.zip', tmp1); unzip(tmp1, exdir=tmp2)
    tmp3 <- paste0(tmp2, '/anatomogram-master/src/svg')
    svgs <- list.files(path=tmp3, pattern='.svg$', full.names=TRUE, recursive=TRUE)
    df <- ftr.return(svgs=svgs, desc=desc)
    svgs.na <- vapply(svgs, function(i) { str <- strsplit(i, '/')[[1]]; str[length(str)] }, FUN.VALUE=character(1))
    svgs1 <- list.files(path=dir, pattern='.svg$', full.names=TRUE)
    svgs.na1 <- list.files(path=dir, pattern='.svg$', full.names=FALSE)
    if (return.all==TRUE) { 

      svgs1.rm <- svgs1[svgs.na1 %in% svgs.na] 
      cat(paste0('Overwriting: ', svgs1.rm, '\n')) 
      # "file.copy" does not overwrite.
      sapply(svgs, function (i) file.copy(i, dir, overwrite=TRUE)); row.names(df) <- NULL;  return(df)

    }

  } else {

    svgs <- list.files(path=dir, pattern='.svg$', full.names=TRUE)
    df <- ftr.return(svgs=svgs, desc=desc)
    if (return.all==TRUE) { row.names(df) <- NULL; return(df) }

  }
 
  if (!is.null(species)) if (species[1]=='') species <- NULL
  if (!is.null(feature)) if (feature[1]=='') species <- NULL
  sp <- gsub(' |_|\\.|-|;|,|/', '|', make.names(species)); ft <- gsub(' |_|\\.|-|;|,|/', '|', make.names(feature))
  sp <- paste0(sp, collapse="|"); ft <- paste0(ft, collapse="|") 
  
  if (keywords.all==TRUE) {

    sp <- strsplit(sp, '\\|')[[1]]; ft <- strsplit(ft, '\\|')[[1]]
    df.idx <- vapply(sp, function(i) grepl(i, df$SVG, ignore.case=TRUE), FUN.VALUE=logical(nrow(df)))

    if (match.only==FALSE) {
      
      w1 <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w1)==0) return(data.frame())
      sp.df <- unique(df$SVG[w1])
      
      df.final <- NULL; for (i in sp.df) {

        df0 <- subset(df, SVG==i)
        df.idx1 <- vapply(ft, function(i) grepl(i, df0$feature, ignore.case=TRUE), FUN.VALUE=logical(nrow(df0)))
        w2 <- which(rowSums(df.idx1)==ncol(df.idx1))
        df.final <- rbind(df.final, df0[w2, ], df0[-w2, ])

      }

    } else {

      df.idx1 <- vapply(ft, function(i) grepl(i, df$feature, ignore.case=TRUE), FUN.VALUE=logical(nrow(df)))
      df.idx <- cbind(df.idx, df.idx1)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      df.final <- df[w, ]
    
    }

  } else {

    if (match.only==FALSE) {

      w1 <- grepl(sp, df$SVG, ignore.case=TRUE)
      if (sum(w1)==0) return(data.frame())
      df1 <- df[w1, ]; sp.df <- unique(df1$SVG)
      df.final <- NULL; for (i in sp.df) {

        df0 <- subset(df1, SVG==i)
        # Move queried features to top.
        w2 <- grep(ft, df0$feature, ignore.case=TRUE)
        # If w2 is empty, -w2 is also empty. Thus SVGs without input features are excluded. 
        df0 <- rbind(df0[w2, ], df0[-w2, ]); df.final <- rbind(df.final, df0)
      }

    } else {

      w <- grepl(sp, df$SVG, ignore.case=TRUE) & grepl(ft, df$feature, ignore.case=TRUE)
      df.final <- df[w, ]

    }

  }; rownames(df.final) <- NULL
  
  if (remote==TRUE) {

    svgs.cp <- svgs[svgs.na %in% df.final$SVG]
    svgs.rm <- svgs1[svgs.na1 %in% df.final$SVG]
    cat(paste0('Overwriting: ', svgs.rm, '\n'))
    sapply(svgs.cp, function (i) file.copy(i, dir, overwrite=TRUE))
    if (dir.exists(tmp)) unlink(tmp, recursive=TRUE)
　　
  }; return(df.final)

}



