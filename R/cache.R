#' Save R Objects in Cache
#'
#' @param dir The directory path to save the cached data. Default is NULL and the cached data is stored in \code{~/.cache/shm}.
#' @param overwrite Logical, TRUE or FALSE. Default is TRUE and data in the cache with the same name of the object in \code{...} will be overwritten.
#' @param ... A single R object to be cached.

#' @return The directory path of the cache.

#' @examples
#' # Save the object "iris" in the default cache "~/.cache/shm".
#' cache.pa <- save_cache(dir=NULL, overwrite=TRUE, iris)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Lori Shepherd and Martin Morgan (2020). BiocFileCache: Manage Files Across Sessions. R package version 1.12.1.

#' @export save_cache
#' @importFrom BiocFileCache BiocFileCache bfcremove bfcquery bfcnew bfccache
#' @importFrom rappdirs user_cache_dir

save_cache <- function(dir=NULL, overwrite=TRUE, ...) {

  dir.con <- is.na(dir)|is.null(dir) 
  if (length(dir.con)==0|sum(dir.con)==1) {
      bfc <- BiocFileCache(user_cache_dir(appname="shm"), ask=FALSE)
  } else bfc <- BiocFileCache(dir, ask=FALSE)
  na <- deparse(substitute(...))
  if (overwrite==TRUE) bfcremove(bfc, bfcquery(bfc, na, exact=TRUE)$rid)
  path <- bfcnew(bfc, na); save(..., file=path)
  cat('Cache directory:', bfccache(bfc),  '\n')
  return(bfccache(bfc))

}


#' Read R Objects from Cache
#'
#' @param dir The directory path where cached data are located. It should be the path returned by \link{save_cache}.
#' @param name The name of the object to retrieve, which is one of the entries in the "rname" column returned by setting \code{info=TRUE}.
#' @param info Logical, TRUE or FALSE. If TRUE (default), the information of all tracked files in cache is returned in a table.

#' @return An R object retrieved from the cache.

#' @examples
#' # Save the object "iris" in the default cache "~/.cache/shm".
#' cache.pa <- save_cache(dir=NULL, overwrite=TRUE, iris)
#' # Retrieve "iris".
#' iris1 <- read_cache(cache.pa, 'iris')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Lori Shepherd and Martin Morgan (2020). BiocFileCache: Manage Files Across Sessions. R package version 1.12.1.

#' @export read_cache
#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcquery bfcinfo

read_cache <- function(dir, name, info=FALSE) {

    bfc <- tryCatch({ BiocFileCache(dir, ask=FALSE) }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
    if (!is(bfc, 'BiocFileCache')) { cat('No valid cache is detected in the provided "dir"! \n'); return() }
    if (info==TRUE) return(bfcinfo(bfc))
    rid <- bfcquery(bfc, name, exact=TRUE)$rid
    if (length(rid)==0) { cat('No valid record is detected for', name, '!\n'); return() }
    if (length(rid)>1) cat('Multiple files matched, the newest one is returned! \n\n')
    id <- rid[length(rid)]
    return(get(load(bfcrpath(bfc, rids=id))))

}




