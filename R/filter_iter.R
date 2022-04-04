#' Iteratively filter bulk and single data according to parameters in a data frame.
#'
#' In secondary filtering of coclustering optimization, iteratively filter normalized bulk and single data (in a list) according to parameter combinations in a data frame. 

#' @param bulk Normalized bulk data at log2-scale returned by \code{norm_multi}.
#' @param cell.lis Normalized single cell data at log2-scale in a named list returned by \code{norm_multi}.
#' @param df.par.fil A \code{data.frame} of filtering parameter settings that are passed to \code{filter_data} and \code{filter_cell} respectively. E.g. \code{df.par.fil <- data.frame(p=c(0.2, 0.3), A=rep(1, 4), cv1=c(0.2, 0.3), cv2=rep(100, 4), min.cnt=rep(1, 4), p.in.cell=c(0.25, 0.3), p.in.gen=c(0.05, 0.1))}.
#' @inheritParams filter_cell
#' @param norm.meth Methods used to normalize bulk and single cell data. One of \code{fct} and \code{cpm}, standing for \code{\link[scran]{computeSumFactors}} only and further normalized by counts per million (cpm) respectively. No actual normalization is performed, only used in file names when saving filtered results. 
#' @param wk.dir The work directory where filtered data are saved in ".rds" files by \code{saveRDS}.

#' @return Filtered data are save in \code{wk.dir}. 

#' @examples

#' See function "coclus_opt" by calling "?coclus_opt".   

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}


#' @export filter_iter

filter_iter <- function(bulk, cell.lis, df.par.fil, gen.rm=NULL, norm.meth, wk.dir) {
  fil.dir <- file.path(wk.dir, 'filter_res')
  if (!norm.meth %in% c('fct', 'cpm')) stop('"norm.meth" should be one of "fct" and "cpm"!')
  if (!dir.exists(fil.dir)) dir.create(fil.dir)
  lis <- NULL; for (i in seq_len(nrow(df.par.fil))) {
    df0 <- df.par.fil[i, ]
    blk <- filter_data(data=bulk, pOA=c(df0$p, df0$A), CV=c(df0$cv1, df0$cv2))
    dat.fil <- filter_cell(lis=cell.lis, bulk=blk, gen.rm=gen.rm, min.cnt=df0$min.cnt, p.in.cell=df0$p.in.cell, p.in.gen=df0$p.in.gen)
    lis <- c(lis, list(dat.fil))
    saveRDS(dat.fil, file=paste0(fil.dir, '/dat.', norm.meth, '.fil', i, '.rds'))
  }; names(lis) <- paste0('dat.', norm.meth, '.fil', seq_len(nrow(df.par.fil)))
    return(lis)
}
