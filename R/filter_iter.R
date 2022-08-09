#' Iteratively filter bulk and single data according to parameters in a data frame.
#'
#' In secondary filtering of coclustering optimization, iteratively filter normalized bulk and single data (in a list) according to parameter combinations in a data frame. 

#' @param bulk Normalized bulk data at log2-scale returned by \code{norm_multi}.
#' @param cell.lis Normalized single cell data at log2-scale in a named list returned by \code{norm_multi}.
#' @param df.par.fil A \code{data.frame} of filtering parameter settings that are passed to \code{filter_data} and \code{filter_cell} respectively. E.g. \code{df.par.fil <- data.frame(p=c(0.2, 0.3), A=rep(1, 4), cv1=c(0.2, 0.3), cv2=rep(100, 4), min.cnt=rep(1, 4), p.in.cell=c(0.25, 0.3), p.in.gen=c(0.05, 0.1))}.
#' @inheritParams filter_cell
#' @param wk.dir The work directory where filtered data are saved in ".rds" files by \code{saveRDS}.
#' @param norm.meth Methods used to normalize bulk and single cell data. One of \code{fct} and \code{cpm}, standing for \code{\link[scran]{computeSumFactors}} only and first \code{\link[scran]{computeSumFactors}} then further normalized by counts per million (cpm) respectively. No actual normalization is performed, only used in file names when saving filtered results.
#' @param verbose Logical. If \code{TRUE} (default), intermediate messages are printed.

#' @return Filtered data are save in \code{wk.dir}. 

#' @examples

#' # Example bulk data of Arabidopsis thaliana (Arabidopsis) root for coclustering optimization (Li et al 2016).
#' blk <- readRDS(system.file("extdata/cocluster/data", "bulk_cocluster.rds", package="spatialHeatmap"))
#'
#' # Example single cell data of Arabidopsis thaliana (Arabidopsis) root for coclustering optimization (Shahan et al 2020).
#' sc10 <- readRDS(system.file("extdata/cocluster/data", "sc10_cocluster.rds", package="spatialHeatmap"))
#'sc11 <- readRDS(system.file("extdata/cocluster/data", "sc11_cocluster.rds", package="spatialHeatmap"))
#'
#' \donttest{

#' # These example data are already pre-processed. To demonstrate the optimization process the pre-processing steps are perfomed again with few genes or cells removed.
#'
#' # Inital filtering before normalization. 
#' blk <- filter_data(data=blk, pOA=c(0.2, 15), CV=c(1.5, 100)); dim(blk)
#' 
#' fil.init <- filter_cell(lis=list(sc10=sc10, sc11=sc11), bulk=blk, gen.rm='^ATCG|^ATCG', min.cnt=1, p.in.cell=0.3, p.in.gen=0.1)
#' 
#' # Normalization.
#' # sum.factor.
#' norm.fct <- norm_multi(dat.lis=fil.init, cpm=FALSE)
#' # sum.factor + CPM.
#' norm.cpm <- norm_multi(dat.lis=fil.init, cpm=TRUE)
#'
#' # Secondary filtering.
#' # Filtering parameter sets.
#' df.par.fil <- data.frame(p=c(0.1, 0.2, 0.3, 0.4), A=rep(1, 4), cv1=c(0.1, 0.2, 0.3, 0.4), cv2=rep(100, 4), min.cnt=rep(1, 4), p.in.cell=c(0.1, 0.25, 0.3, 0.35), p.in.gen=c(0.01, 0.05, 0.1, 0.15))
#' df.par.fil
#'
#' # Filtered results are saved in "opt_res".
#' if (!dir.exists('opt_res')) dir.create('opt_res')
#' fct.fil.all <- filter_iter(bulk=norm.fct$bulk, cell.lis=list(sc10=norm.fct$sc10, sc11=norm.fct$sc11), df.par.fil=df.par.fil, gen.rm='^ATCG|^ATCG', wk.dir='opt_res')
#'
#' cpm.fil.all <- filter_iter(bulk=norm.cpm$bulk, cell.lis=list(sc10=norm.cpm$sc10, sc11=norm.cpm$sc11), df.par.fil=df.par.fil, gen.rm='^ATCG|^ATCG', wk.dir='opt_res')
#' }  

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}


#' @export filter_iter

filter_iter <- function(bulk, cell.lis, df.par.fil, gen.rm=NULL, wk.dir, norm.meth, verbose=TRUE) {
  fil.dir <- file.path(wk.dir, 'filter_res')
  if (!norm.meth %in% c('fct', 'cpm')) stop('"norm.meth" should be one of "fct" and "cpm"!')
  if (!dir.exists(fil.dir)) dir.create(fil.dir)
  lis <- NULL; for (i in seq_len(nrow(df.par.fil))) {
    df0 <- df.par.fil[i, ]
    blk <- filter_data(data=bulk, pOA=c(df0$p, df0$A), CV=c(df0$cv1, df0$cv2), verbose=verbose)
    if (nrow(blk)==0) {
      print(df0); stop('All genes are filtered out!')
    }
    dat.fil <- filter_cell(lis=cell.lis, bulk=blk, gen.rm=gen.rm, min.cnt=df0$min.cnt, p.in.cell=df0$p.in.cell, p.in.gen=df0$p.in.gen, verbose=verbose)
    lis <- c(lis, list(dat.fil))
    saveRDS(dat.fil, file=paste0(fil.dir, '/fil', i, '.', norm.meth, '.rds'))
  }; names(lis) <- paste0('fil', seq_len(nrow(df.par.fil)), '.', norm.meth)
  return(lis)
}
