#' Generate random combinations of parameter settings
#'
#' In coclustering, generate random combinations of parameter settings for validating optimal settings. 
#' @param fil.set A character vector of filtering parameter set. E.g. \code{c('fil3', 'fil4')}.
#' @param norm A character vector of normalization methods. E.g. \code{c('cpm')}.
#' @param dimred A character vector of dimesionality reduction methods. E.g. \code{c('umap')}.
#' @param graph.meth A character vector of graph-building methods. E.g. \code{c('knn', 'snn')}.
#' @inheritParams coclus_meta
#' @param df.spd.opt A \code{data.frame} of optimized \code{spd.set} settings. These settings are avoided in the output random settings. 


#' @return A \code{data.frame}.

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @examples
#'
#' df.spd.opt <- data.frame(sim=c(0.2, 0.4, 0.3), sim.p=c(0.8, 0.6, 0.7), dim=c(12, 14, 13))
#' df.para.rdn <- random_para(fil.set=c('fil3', 'fil4'), norm='cpm', dimred='UMAP', graph.meth=c('knn', 'snn'), sim=round(seq(0.2, 0.8, by=0.1), 1), sim.p=round(seq(0.2, 0.8,by=0.1), 1), dim=seq(5, 40, by=1), df.spd.opt=df.spd.opt)

#' @export random_para
 
random_para <- function(fil.set, norm, dimred, graph.meth, sim=round(seq(0.2, 0.8, by=0.1), 1), sim.p=round(seq(0.2, 0.8,by=0.1), 1), dim=seq(5, 40, by=1), df.spd.opt) {
  spd.set <- NULL
  # All spd.set
  df.spd.all <- expand.grid(sim=sim, sim.p=sim.p, dim=dim, stringsAsFactors = FALSE)
  df.spd.all$spd.set <- paste0('s', df.spd.all$sim, 'p', df.spd.all$sim.p, 'd', df.spd.all$dim)
  if (!'spd.set' %in% colnames(df.spd.opt)) df.spd.opt$spd.set <- paste0('s', df.spd.opt$sim, 'p', df.spd.opt$sim.p, 'd', df.spd.opt$dim)
  df.spd.all <- subset(df.spd.all, !spd.set %in% df.spd.opt$spd.set)

  n <- nrow(df.spd.opt)
  # fil.set
  fil.set <- sample(fil.set, n, replace=TRUE)
  # Norm and dimred.
  norm <- sample(norm, n, replace=TRUE)
  dimred <- sample(dimred, n, replace=TRUE)

  df.spd.sel <- df.spd.all[sample(seq_len(nrow(df.spd.all)), n, replace=FALSE), ]
  # Clustering.
  graph.meth <- sample(graph.meth, n, replace=TRUE)
  df.par.rdn <- cbind(fil.set=fil.set, norm=norm, dimred=dimred, df.spd.sel[, c('sim', 'sim.p', 'dim', 'spd.set')], graph.meth=graph.meth)
  rownames(df.par.rdn) <- seq_len(n); return(df.par.rdn)
}


