#' Return aSVG Files Relevant to Target Features
#' 
#' This function parses a collection of aSVG files and returns those containing target features in a data frame. Successful spatial heatmap plotting requires the aSVG features of interest have matching samples (cells, tissues, \emph{etc}) in the data. To meet this requirement, the returned features could be used to replace target sample counterparts in the data. Alternatively, the target samples in the data could be used to replace matching features in the aSVG through function \code{\link{update_feature}}. Refer to function \code{\link{spatial_hm}} for more details on aSVG files.  

#' @param feature A vector of target feature keywords (case insentitive), which is used to select aSVG files from a collection. \emph{E.g.} c('heart', 'brain'). If NA or NULL, all features of all SVG files matching \code{species} are returned. 
#' @param species A vector of target species keywords (case insentitive), which is used to select aSVG files from a collection. \emph{E.g.} c('gallus'). If NA or NULL, all SVG files in \code{dir} are queried.
#' @param keywords.any Logical, TRUE or FALSE. Default is TRUE. The internal searching is case-insensitive. The space, dot, hypen, semicolon, comma, forward slash are treated as separators between words and not counted in searching. If TRUE, every returned hit contains at least one word in the \code{feature} vector and at least one word in the \code{species} vector, which means all the possible hits are returned. \emph{E.g.} "prefrontal cortex" in "homo_sapiens.brain.svg" would be returned if \code{feature=c('frontal')} and \code{species=c('homo')}. If FALSE, every returned hit contains at least one exact element in the \code{feature} vector and all exact elements in the \code{species} vector. \emph{E.g.} "frontal cortex" rather than "prefrontal cortex" in "homo_sapiens.brain.svg" would be returned if \code{feature=c('frontal cortex')} and \code{species=c('homo sapiens', 'brain')}.
#' @param remote Logical, FALSE or TRUE. If TRUE (default), the remote EBI aSVG repository \url{https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg} and spatialHeatmap aSVG Repository \url{https://github.com/jianhaizhang/SHM_SVG/tree/master/SHM_SVG_repo} developed in this project are queried.
#' @param dir The directory path of aSVG files. If \code{remote} is TRUE, the returned aSVG files are saved in this directory. Note existing aSVG files with same names as returned ones are overwritten. If \code{remote} is FALSE, user-provided (local) aSVG files should be saved in this directory for query. Default is NULL.
#' @param svg.path The path of a specific aSVG file. If the provided aSVG file exists, only features of this file are returned and there will be no querying process. Default is NULL.
#' @param desc Logical, FALSE or TRUE. Default is FALSE. If TRUE, the feature descriptions from the R package "rols" (Laurent Gatto 2019) are added. If too many features are returned, this process takes a long time.
#' @param match.only Logical, TRUE or FALSE. If TRUE (default), only target features are returned. If FALSE, all features in the matching aSVG files are returned, and the matching features are moved on the top of the data frame. 
#' @param return.all Logical, FALSE or TRUE. Default is FALSE. If TRUE, all features together with all respective aSVG files are returned, regardless of \code{feature} and \code{species}.

#' @return A data frame containing information on target features and aSVGs.

#' @examples

#' # This function is able to work on the EBI aSVG repository directly: https://github.com/
#' # ebi-gene-expression-group/anatomogram/tree/master/src/svg. The following shows how to
#' # download a chicken aSVG containing spatial features of 'brain' and 'heart'. An empty
#' # directory is recommended so as to avoid overwriting existing SVG files. 
#' # Here "~/test" is used. 
#'
#' \donttest{
#' # Make an empty directory "~/test" if not exist.
#' if (!dir.exists('~/test')) dir.create('~/test')
#' # Remote aSVG repos.
#' data(aSVG.remote.repo)
#' tmp.dir <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
#' tmp.dir.ebi <- paste0(tmp.dir, '/ebi.zip')
#' tmp.dir.shm <- paste0(tmp.dir, '/shm.zip')
#' # Download the remote aSVG repos as zip files. According to Bioconductor's 
#' # requirements, downloadings are not allowed inside functions, so the repos are 
#' # downloaded before calling "return_feature".  
#' download.file(aSVG.remote.repo$ebi, tmp.dir.ebi)
#' download.file(aSVG.remote.repo$shm, tmp.dir.shm)
#' remote <- list(tmp.dir.ebi, tmp.dir.shm)
#' # Query the remote aSVG repos.
#' feature.df <- return_feature(feature=c('heart', 'brain'), species=c('gallus'), dir='~/test',
#' match.only=FALSE, remote=remote)
#' feature.df
#' # The path of downloaded aSVG.
#' svg.chk <- '~/test/gallus_gallus.svg'
#' }
#'
#' # The spatialHeatmap package has a small aSVG collection and can be used to demonstrate the
#' # local query.
#' # Get the path of local aSVGs from the package.
#' svg.dir <- system.file("extdata/shinyApp/data", package="spatialHeatmap")
#' # Query the local aSVG repo. The "species" argument is set NULL on purpose so as to illustrate
#' # how to select the target aSVG among all matching aSVGs.
#' feature.df <- return_feature(feature=c('heart', 'brain'), species=NULL, dir=svg.dir, 
#' match.only=FALSE, remote=NULL)
#' # All matching aSVGs.
#' unique(feature.df$SVG)
#' # Select the target aSVG of chicken.
#' subset(feature.df, SVG=='gallus_gallus.svg')


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Laurent Gatto (2019). rols: An R interface to the Ontology Lookup Service. R package version 2.14.0. http://lgatto.github.com/rols/
#' \cr Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2
#' \cr R Core Team (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. "Gene Expression Across Mammalian Organ Development." Nature 571 (7766): 505-9

#' @export return_feature
#' @importFrom xml2 read_xml
#' @importFrom utils unzip

return_feature <- function(feature, species, keywords.any=TRUE, remote=NULL, dir=NULL, svg.path=NULL, desc=FALSE, match.only=TRUE, return.all=FALSE) {
  options(stringsAsFactors=FALSE)
  # Parse and return features.
  ftr_return <- function(svgs, desc=desc) {
    cat('Accessing features... \n'); df <- NULL; for (i in svgs) {
      doc <- read_xml(i); df0 <- svg_attr(doc, feature=NULL, TRUE)[['df.attr']]
      # Move ontology with NA or "NULL" to bottom.
      w.na <- which(is.na(df0$id)|df0$id=='NULL'|df0$id=='NA')
      if (length(w.na)>0) df0 <- rbind(df0[-w.na, ], df0[w.na, ])
      na <- strsplit(i, '/')[[1]]; na <- na[grep('\\.svg$', na)]
      cat(na); cat(', '); df0$SVG <- na; df <- rbind(df, df0)

    }; cat('\n')
    cna <- colnames(df); df <- df[, c('feature', 'stroke', 'fill', 'id', 'element', 'parent', 'index.sub', 'SVG')]
    colnames(df)[colnames(df)=='index.sub'] <- 'order'

    if (desc==TRUE) {
      pkg <- check_pkg('rols'); if (is(pkg, 'character')) stop(pkg)
      cat('Appending descriptions... \n')
      df$description <- NA; for (i in seq_len(nrow(df))) {
        ont <- df[i, 'id']; abbr <- tolower(sub('_.*', '', ont))
        trm <- tryCatch({ rols::Term(abbr, ont) }, error=function(e) { return(NA) })
        if (is(trm, 'Term')) { des <- rols::termDesc(trm); if (!is.null(des)) df[i, 'description'] <- rols::termDesc(trm) }
      }
    }; return(df)

  }
  # Parse the provided svg file.
  if (!is.null(svg.path)) {
    if (file.exists(svg.path)) return(ftr_return(svg.path, desc)) else stop('The provided aSVG file does not exists!')
  }
  # Query aSVG repo.
  SVG <- parent <- NULL
  dir <- normalizePath(dir, winslash = "/", mustWork=FALSE)
  if (!is.vector(feature) & is.vector(species)) stop('"feature" and "species" must be two vectors respectively!')
  if (length(species)==0) species <- NULL; if (length(feature)==0) feature <- NULL
  if (!is.null(species)) if (species[1]==''|is.na(species[1])) species <- NULL
  if (!is.null(feature)) if (feature[1]==''|is.na(feature[1])) feature <- NULL
  # Pattern for selecting relevant species.
  if (length(species)>0 & return.all==FALSE) { species1 <- gsub(' |_|\\.|-|;|,|/', '|', make.names(species))
  species1 <- paste0(species1, collapse="|") }
  if (is(remote, 'list')) { 
    tmp <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE) 
    # Extract paths of aSVG files in remote zip files.
    svgs <- svg_zpath(remote);
    # Only select relevant svgs to species.
    if (length(species)>0 & return.all==FALSE) {
      svgs.sp <- grep(species1, svgs, value=TRUE, ignore.case=TRUE)
      if (length(svgs.sp)>0) svgs <- svgs.sp
    }
    df <- ftr_return(svgs=svgs, desc=desc)
    svgs.na <- vapply(svgs, function(i) { str <- strsplit(i, '/')[[1]]; str[length(str)] }, FUN.VALUE=character(1))
    svgs1 <- list.files(path=dir, pattern='.svg$', full.names=TRUE)
    svgs.na1 <- list.files(path=dir, pattern='.svg$', full.names=FALSE)
    if (return.all==TRUE) { 

      svgs1.rm <- svgs1[svgs.na1 %in% svgs.na] 
      cat(paste0('Overwriting: ', svgs1.rm, '\n')) 
      vapply(svgs, function(i) file.copy(i, dir, overwrite=TRUE), logical(1)); row.names(df) <- NULL;  return(df)

    }

  } else {

    svgs <- list.files(path=dir, pattern='.svg$', full.names=TRUE)
    # Only select relevant svgs to species.
    if (length(species)>0 & return.all==FALSE) {
      svgs.sp <- grep(species1, svgs, value=TRUE, ignore.case=TRUE)
      if (length(svgs.sp)>0) svgs <- svgs.sp
    }
    df <- ftr_return(svgs=svgs, desc=desc)
    if (return.all==TRUE) { row.names(df) <- NULL; return(df) }

  }
  
  if (keywords.any==FALSE) {

    sp <- gsub(' |\\.|-|;|,|/', '_', make.names(species)); ft <- gsub(' |\\.|-|;|,|/', '_', make.names(feature))
    SVG <- paste0('_', gsub(' |\\.|-|;|,|/', '_', make.names(df$SVG)))
    df.idx <- vapply(sp, function(i) grepl(paste0('_', i, '_'), SVG, ignore.case=TRUE), FUN.VALUE=logical(nrow(df)))
    w1 <- which(rowSums(df.idx)==ncol(df.idx))
    if (length(w1)==0) return(data.frame())

    if (match.only==FALSE) {
      
      sp.df <- unique(df$SVG[w1])  
      df.final <- NULL; for (i in sp.df) {

        df0 <- subset(df, SVG==i)
        ft0 <- gsub(' |\\.|-|;|,|/', '_', make.names(df0$feature))
        df.idx1 <- vapply(ft, function(i) grepl(paste0('^', i, '$'), ft0, ignore.case=TRUE), FUN.VALUE=logical(nrow(df0)))
        w2 <- which(rowSums(df.idx1)>0)
        df.final <- rbind(df.final, df0[w2, ], df0[-w2, ])

      }

    } else {

      df0 <- df[w1, ]; ft0 <- gsub(' |\\.|-|;|,|/', '_', make.names(df0$feature))
      df.idx1 <- vapply(ft, function(i) grepl(paste0('^', i, '$'), ft0, ignore.case=TRUE), FUN.VALUE=logical(nrow(df0)))
      w2 <- which(rowSums(df.idx1)>0)
      df.final <- df0[w2, ]
    
    }

  } else {
  
    sp <- gsub(' |_|\\.|-|;|,|/', '|', make.names(species)); ft <- gsub(' |_|\\.|-|;|,|/', '|', make.names(feature))
    sp <- paste0(sp, collapse="|"); ft <- paste0(ft, collapse="|") 

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
  
  if (is(remote, 'list')) {
    svgs.cp <- svgs[svgs.na %in% unique(df.final$SVG)]
    svgs.rm <- svgs1[svgs.na1 %in% unique(df.final$SVG)]
    cat(paste0('Overwriting: ', svgs.rm, '\n'))
    vapply(svgs.cp, function (i) file.copy(i, dir, overwrite=TRUE), logical(1))
    if (dir.exists(tmp) & !is(remote, 'list')) unlink(tmp, recursive=TRUE)
  }; return(df.final)

}

#' Extract aSVG paths in remote zip files.
#'
#' @param zip.paths A list of paths of zipped aSVG files.
#' @return A vector of full paths of aSVG files.
#' @keywords Internal
#' @noRd

svg_zpath <- function(zip.path) { 
  if (!all(grepl('\\.zip$', zip.path))) stop('Please ensure the aSVG files are in a zip file with the extension ".zip".')
  tmp <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
  tmp.uz <- paste0(tmp, '/asvg'); tmp.zcopy <- paste0(tmp, '/zcopy')
  unlink(tmp.uz, recursive=TRUE); dir.create(tmp.uz)
  unlink(tmp.zcopy, recursive=TRUE); dir.create(tmp.zcopy)
 for (i in seq_along(zip.path)) {
   zip.pa <- zip.path[[i]]
   if (!file.exists(zip.pa)) stop(paste0(zip.pa, ': file does not exist!'))
   # Make a copy of provided SVG zip file.
   svg.zpa <- paste0(tmp.zcopy, '/svg.zip') 
   file.copy(zip.pa, svg.zpa); unzip(svg.zpa, exdir=tmp.uz)
 }
 svgs <- list.files(path=tmp.uz, pattern='.svg$', full.names=TRUE, recursive=TRUE)
 if (length(svgs)==0) stop('No aSVG files are found!')
 return(svgs)
}



