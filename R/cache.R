#' Save R Objects in Cache
#'
#' @param dir The directory path to save the cached data. Default is NULL and the cached data is stored in \code{~/.cache/shm}.
#' @param overwrite Logical, TRUE or FALSE. Default is TRUE and data in the cache with the same name of the object in \code{...} will be overwritten.
#' @param obj,na A single R object (`obj`) to be cached into the file named `na`. 

#' @return The directory path of the cache.

#' @examples
#' # Save the object "iris" in the default cache "~/.cache/shm".
#' cache.pa <- save_cache(dir=NULL, overwrite=TRUE, iris)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Lori Shepherd and Martin Morgan (2020). BiocFileCache: Manage Files Across Sessions. R package version 1.12.1.

#' @export save_cache

save_cache <- function(dir=NULL, overwrite=TRUE, obj, na=NULL) {
  pkg <- check_pkg('BiocFileCache'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  pkg <- check_pkg('rappdirs'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  dir.con <- is.na(dir)|is.null(dir) 
  if (length(dir.con)==0|sum(dir.con)==1) {
      bfc <- BiocFileCache::BiocFileCache(rappdirs::user_cache_dir(appname="shm"), ask=FALSE)
  } else bfc <- BiocFileCache::BiocFileCache(dir, ask=FALSE)
  if (is.null(na)) na <- deparse(substitute(obj))
  if (overwrite==TRUE) BiocFileCache::bfcremove(bfc, BiocFileCache::bfcquery(bfc, na, exact=TRUE)$rid)
  path <- BiocFileCache::bfcnew(bfc, na); save(obj, file=path)
  pa <- BiocFileCache::bfccache(bfc)
  cat('Cache directory:', pa,  '\n'); return(pa)
}


#' Read R Objects from Cache
#'
#' @param dir The directory path where cached data are located. It should be the path returned by \link{save_cache}.
#' @param name The name of the object to retrieve, which is one of the entries in the "rname" column returned by setting \code{info=TRUE}.
#' @param info Logical, TRUE or FALSE. If TRUE (default), the information of all tracked files in cache is returned in a table.
#' @param day An integer representing cache age in days. If the cache is older than this age, NULL is returned.

#' @return An R object retrieved from the cache.

#' @examples
#' # Save the object "iris" in the default cache "~/.cache/shm".
#' cache.pa <- save_cache(dir=NULL, overwrite=TRUE, iris)
#' # Retrieve "iris".
#' iris1 <- read_cache(cache.pa, 'iris')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Lori Shepherd and Martin Morgan (2020). BiocFileCache: Manage Files Across Sessions. R package version 1.12.1.

#' @export read_cache

read_cache <- function(dir, name, day=180, info=FALSE) {
  pkg <- check_pkg('BiocFileCache'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  bfc <- tryCatch({ BiocFileCache::BiocFileCache(dir, ask=FALSE) }, error=function(e){ return('error') }, 
                  warning=function(w) { return('warning') } )
  if (!is(bfc, 'BiocFileCache')) { cat('No valid cache is detected in the provided "dir"! \n'); return() }
  tbl = BiocFileCache::bfcinfo(bfc)
  if (info==TRUE) return(tbl)
  rid <- BiocFileCache::bfcquery(bfc, name, exact=TRUE)$rid
  len = length(rid)
  if (length(rid)==0) { cat('No valid record is detected for', name, '!\n'); return() }
  if (len>1) {cat('Multiple files matched, the newest one is returned! \n\n')
    rid <- rid[len]
    file = tbl$rpath[tbl$rid==rid][len]
  } else if (len==1) file = tbl$rpath[tbl$rid==rid]
  old = file.info(file)$mtime < Sys.time() - as.difftime(day, units = "days")
  if (old) {
    return(NULL)
    # BiocFileCache::bfcremove(bfc, id); pa <- BiocFileCache::bfccache(bfc)
  } else {
    return(get(load(BiocFileCache::bfcrpath(bfc, rids=rid))))
  }
}

