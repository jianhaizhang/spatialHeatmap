#' Reduce sample replicates
#'
#' In an expression profile matrix such as RNA-seq count table, where columns and rows are samples and biological molecules respectively, reduce sample replicates according to sum of correlation coefficients (Pearson, Spearman, Kendall).

#' @param dat Abundance matrix in form of \code{data.frame} or \code{matrix}, where columns and rows are samples and biological molecules respectively. For example, gene exression matrix generated in RNA-seq.  
#' @param n An integer, the max number of replicates to keep per sample (e.g. tissue type). Within each sample, pairwise correlations are calculated among all replicates, and the correlations between one replicate and other relicates are summed. The replicates with top n largest sums are retained in each sample.
#' @param sim.meth One of \code{pearson} (default), \code{kendall}, or \code{spearman}, indicating which correlation coefficient method to use for calculating similarities between replicates.  

#' @return A \code{matrix}.

#' @examples
#' 
#' # Random abundance matrix.
#' dat <- matrix(rnorm(100), nrow=10)
#' # Two samples, each has 5 replicates.
#' colnames(dat) <- c(rep('sampleA', 5), rep('sampleB', 5))
#' rownames(dat) <- paste0('gene', seq_len(nrow(dat)))
#' reduce_rep(dat)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @export reduce_rep

reduce_rep <- function(dat, n=3, sim.meth='pearson') {
  dat <- as(dat, 'matrix')
  cna <- colnames(dat); cna.uni <- unique(cna)
  dat.red <- dat[, 1, drop=FALSE]; dat.red[, 1] <- NA
  for (i in cna.uni) {
    blk0 <- dat[, cna == i, drop = FALSE]
    # Select reps according to sum of correlations.
    idx <- rank(-colSums(cor(blk0, method=sim.meth))) <= n
    blk0 <- blk0[, idx, drop=FALSE]
    dat.red <- cbind(dat.red, blk0)
  }; dat.red <- dat.red[, -1]; return(dat.red)
}
