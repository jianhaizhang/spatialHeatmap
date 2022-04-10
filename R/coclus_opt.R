#' Optimizing coclustering process with two levels of parallelizations available. 
#'
#' Optimizing coclustering process with two levels of parallelizations available.

#' @param wk.dir Working directory. Must be the same with that in \code{filter_iter}, since the filtered data will be used automatically.
#' @param parallel.info Logical. If \code{FALSE} (default), coclustering optimization is performed. If \code{TRUE}, parallelization guide is returned and no optimization. Users are advised to check the guide when editing the slurm template.  
#' @inheritParams cocluster
#' @param dimred Dimensionality redeuction methods to optimize: \code{c('PCA', 'UMAP')}.
#' @param graph.meth Graph-building methods in coclustering: \code{c('knn', 'snn')} (see \code{\link[scran]{buildKNNGraph}}), \code{\link[scran]{buildSNNGraph}}) respectively). The clusters are detected by first creating a nearest neighbor graph using \code{snn} or \code{knn} then partitioning the graph. 
#' @param sim,sim.p Used when refining cell clusters. Both are numeric scalars, ranging from 0 to 1. \code{sim} is a similarity (Spearman or Pearson correlation coefficient) cutoff between cells and \code{sim.p} is a proportion cutoff. In a certain cell cluster, cells having similarity >= \code{sim} with other cells in the same cluster at proportion >= \code{sim.p} would remain. Otherwise, they are discarded. The default of both is \code{seq(0.2, 0.8, by=0.1)} and can be customized.
#' @param dim Number of principle components (PCs, equivalent to genes) in combined bulk and single cell data. Used as the minimum number of PCs to retain in \code{\link[scran]{denoisePCA}} when coclustering bulk and single cells. The default is \code{seq(5, 40, by=1)}, and can be customized.
#' @inheritParams coclus_roc
#' @param batch.par The parameters for first-level parallelization. See \code{\link[BioParallel]{BatchtoolsParam}}. It works with the "slurm" scheduler, so "slurm" needs to be installed. If \code{NULL}, the first-level parallelization is skipped.
#' @param multi.core.par The parameters for second-level parallelization. See \code{\link[BioParallel]{MulticoreParam}}.
#' @param verbose Logical. If \code{TRUE} (default), intermediate messages are printed.

#' @return A list of normalized single cell and bulk data. 

#' @examples

#' # To obtain reproducible results, always start a new R session and set a fixed seed for Random Number Generator at the beginning, which is required only once in each R session.  
#' set.seed(10)
#' 
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
#' fct.fil.all <- filter_iter(bulk=norm.fct$bulk, cell.lis=list(sc10=norm.fct$sc10, sc11=norm.fct$sc11), df.par.fil=df.par.fil, gen.rm='^ATCG|^ATCG', wk.dir='opt_res', norm.meth='fct')
#'
#' cpm.fil.all <- filter_iter(bulk=norm.cpm$bulk, cell.lis=list(sc10=norm.cpm$sc10, sc11=norm.cpm$sc11), df.par.fil=df.par.fil, gen.rm='^ATCG|^ATCG', wk.dir='opt_res', norm.meth='cpm')
#'
#' # Matching table between bulk tissues and single cells.
#' match.pa <- system.file("extdata/cocluster/data", "match_arab_root_cocluster.txt", package="spatialHeatmap")
#' df.match.arab <- read.table(match.pa, header=TRUE, row.names=1, sep='\t')
#' df.match.arab[1:3, ]
#'
#' # Optimization. 
#' # Check parallelization guide.
#' coclus_opt(wk.dir='opt_res', parallel.info=TRUE, dimred=c('PCA', 'UMAP'), graph.meth=c('knn', 'snn'), sim=seq(0.2, 0.4, by=0.1), sim.p=seq(0.2, 0.4, by=0.1), dim=seq(5, 7, by=1))
#'
#' # The first-level parallel computing relies on the slurm scheduler (https://slurm.schedmd.com/documentation.html), so if it is available the whole optimization process could be parallelized at two levels. 
# Copy slurm template to current directory. Edit the template according to the parallelization guide and available computing resources.
#' file.copy(system.file("extdata/cocluster", "slurm.tmpl", package="spatialHeatmap"), './slurm.tmpl')
#' 
#' # The first- and second-level parallelizations are set 3 and 2 respectively.
#' library(BiocParallel)
#' opt <- coclus_opt(wk.dir='opt_res', dimred=c('PCA', 'UMAP'), graph.meth=c('knn', 'snn'), sim=seq(0.2, 0.4, by=0.1), sim.p=seq(0.2, 0.4, by=0.1), dim=seq(5, 7, by=1), df.match=df.match.arab, batch.par=BatchtoolsParam(workers=3, cluster="slurm", template='slurm.tmpl'), multi.core.par=MulticoreParam(workers=2))
#'
#' # If slurm is not available, parallelize the optimization only at the second-level through 2 workers. 
#' opt <- coclus_opt(wk.dir='opt_res', dimred=c('PCA', 'UMAP'), graph.meth=c('knn', 'snn'), sim=seq(0.2, 0.4, by=0.1), sim.p=seq(0.2, 0.4, by=0.1), dim=seq(5, 7, by=1), df.match=df.match.arab, batch.par=NULL, multi.core.par=MulticoreParam(workers=2))
#'
#' # The performaces of parameter settings are measured by AUC values in ROC curve. The following demonstrates how to visualize the AUCs and select optimal parameter settings.
# If one AUC is larger than 0.5, 0.6, 0.7, 0.8, or 0.9, and has total bulk assignments >= 500, true assignments >= 500, the corresponding parameter settings are extracted. 
#'
#' # Extract AUCs and other parameter settings for filtering parameter sets.
#' df.lis.fil <- auc_stat(wk.dir='opt_res', tar.par='filter', total.min=500, true.min=300, aucs=round(seq(0.5, 0.9, 0.1), 1))
#' df.lis.fil$df.auc.mean[1:3, ]
#' 
#' # Mean AUCs by each filtering settings and AUC cutoff.
#' mean_auc_bar(df.lis.fil[[1]], bar.width=0.07, title='Mean AUCs by filtering settings')
#' 
#' # All AUCs by each filtering settings and AUC cutoff. 
#' auc_violin(df.lis=df.lis.fil, xlab='Filtering settings')
#'
#' # Optimal filtering settings: fil1, fil2, fil3
#' df.par.fil[c(1, 2, 3), ]
#' 
#' # Extract AUCs and other parameter settings for normalization methods.
#' df.lis.norm <- auc_stat(wk.dir='opt_res', tar.par='norm', total.min=500, true.min=300, aucs=round(seq(0.5, 0.9, 0.1), 1))
#' df.lis.norm$df.auc.mean[1:3, ]
#'
#' # Mean AUCs by each normalization method and AUC cutoff.
#' mean_auc_bar(df.lis.norm[[1]], bar.width=0.07, title='Mean AUCs by normalization methods')
#' 
#' # All AUCs by each normalization method and AUC cutoff. 
#' auc_violin(df.lis=df.lis.norm, xlab='Normalization methods')
#' 
#' # Optimal normalization method: fct (computeSumFactors).
#' 
#' # Extract AUCs and other parameter settings for graph-building methods.
#' df.lis.graph <- auc_stat(wk.dir='opt_res', tar.par='graph', total.min=500, true.min=300, aucs=round(seq(0.5, 0.9, 0.1), 1))
#' df.lis.graph$df.auc.mean[1:3, ]
#' 
#' # Mean AUCs by each graph-building method and AUC cutoff.
#' mean_auc_bar(df.lis.graph[[1]], bar.width=0.07, title='Mean AUCs by graph-building methods')
#' 
#' # All AUCs by each graph-building method and AUC cutoff. 
#' auc_violin(df.lis=df.lis.graph, xlab='Graph-building methods')
#' 
#' # Optimal graph-building methods: knn (buildKNNGraph).
#' 
#' # Extract AUCs and other parameter settings for dimensionality reduction methods.
#' df.lis.dimred <- auc_stat(wk.dir='opt_res', tar.par='dimred', total.min=500, true.min=300, aucs=round(seq(0.5, 0.9, 0.1), 1))
#' df.lis.dimred$df.auc.mean[1:3, ]
#' 
#' # Mean AUCs by each dimensionality reduction method and AUC cutoff.
#' mean_auc_bar(df.lis.dimred[[1]], bar.width=0.07, title='Mean AUCs by dimensionality reduction methods')
#' 
#' # All AUCs by each dimensionality reduction method and AUC cutoff. 
#' auc_violin(df.lis=df.lis.dimred, xlab='Dimensionality reduction')
#' 
#' # Optimal dimensionality reduction method: pca (denoisePCA).
#' 
#' # Extract AUCs and other parameter settings for spd.sets.
#' df.lis.spd <- auc_stat(wk.dir='opt_res', tar.par='spd.set', total.min=500, true.min=300, aucs=round(seq(0.5, 0.9, 0.1), 1))
#' df.lis.spd$auc0.5$df.frq[1:3, ]
#'
#' # All AUCs of top spd.sets ranked by frequency. 
#' spd_auc_violin(df.lis=df.lis.spd, n=5, xlab='spd.sets', x.vjust=0.6)
#' 
#' # Optimal spd.sets.
#' df.spd.top <- rbind(df.lis.spd[[1]]$df.frq[1:5, ], df.lis.spd[[2]]$df.frq[1:5, ], df.lis.spd[[3]]$df.frq[1:5, ], df.lis.spd[[4]]$df.frq[1:5, ], df.lis.spd[[5]]$df.frq[1:5, ])
#'
#' df.spd.top$spd.set <- paste0('s', df.spd.top$sim, 'p', df.spd.top$sim.p, 'd', df.spd.top$dim)
#' df.spd.top <- subset(df.spd.top, !duplicated(spd.set))
#' df.spd.top[1:3, c('sim', 'sim.p', 'dim', 'spd.set')]

#' }

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references

#' Li, S., Yamada, M., Han, X., Ohler, U., and Benfey, P. N. (November, 2016) High-Resolution Expression Map of the Arabidopsis Root Reveals Alternative Splicing and lincRNA Regulation. Dev. Cell, 39(4), 508–522.
#' Shahan, R., Hsu, C.-W., Nolan, T. M., Cole, B. J., Taylor, I. W., Vlot, A. H. C., Benfey, P. N., and Ohler, U. (June, 2020) A single cell Arabidopsis root atlas reveals developmental trajectories in wild type and cell identity mutants.
#' Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2021). BiocParallel: Bioconductor facilities for parallel evaluation. R package version 1.28.3, https://github.com/Bioconductor/BiocParallel.

#' @export coclus_opt
#' @importFrom BiocParallel BatchtoolsParam MulticoreParam register bpRNGseed bpRNGseed<-

coclus_opt <- function(wk.dir, parallel.info=FALSE, sc.dim.min=10, max.dim=50, dimred=c('PCA', 'UMAP'), graph.meth=c('knn', 'snn'), sim=seq(0.2, 0.8, by=0.1), sim.p=seq(0.2, 0.8, by=0.1), dim=seq(5, 40, by=1), df.match, sim.meth='spearman', batch.par=BatchtoolsParam(workers=1, cluster="slurm", template='slurm.tmpl', RNGseed=100, stop.on.error = FALSE, log = TRUE, logdir=file.path(wk.dir, 'batch_log')), multi.core.par=MulticoreParam(workers=1, RNGseed=NULL, stop.on.error=FALSE, log=TRUE, logdir=file.path(wk.dir, 'multi_core_log')), verbose=TRUE) {
  fil.dir <- file.path(wk.dir, 'filter_res')
  bat.log.dir <- file.path(wk.dir, 'batch_log')
  mcore.log.dir <- file.path(wk.dir, 'multi_core_log')
  auc.dir <- file.path(wk.dir, 'auc_res')
  if (!dir.exists(fil.dir)) dir.create(fil.dir) 
  if (!dir.exists(bat.log.dir)) dir.create(bat.log.dir) 
  if (!dir.exists(mcore.log.dir)) dir.create(mcore.log.dir) 
  if (!dir.exists(auc.dir)) dir.create(auc.dir) 
  # Check random seed. 
  if (is(batch.par, 'BatchtoolsParam')) seed.bat <- bpRNGseed(batch.par) else seed.bat <- NULL
  if (is(multi.core.par, 'MulticoreParam')) seed.mul <- bpRNGseed(multi.core.par) else seed.mul <- NULL
  if (is.numeric(seed.bat) & !is.numeric(seed.mul)) { 
    bpRNGseed(multi.core.par) <- NULL
    message('"RNGseed" in MulticoreParam is set NULL, since it is already set in BatchtoolsParam.')
  }
  if (!is.numeric(seed.bat) & !is.numeric(seed.mul)) {
    message('"RNGseed" in BatchtoolsParam or MulticoreParam is NULL')
  }
  # All single cell data names.
  fil.nas <- list.files(fil.dir, '\\.fil\\d+\\.rds$')
  sc.na.all <- NULL; for (i in fil.nas) {
    pa <- file.path(fil.dir, i); dat.fil0 <- readRDS(pa)
    nas <- names(dat.fil0); sc.nas <- setdiff(nas, 'bulk')
    sc.na.all <- unique(c(sc.na.all, sc.nas))
  }

  df.par.com <- expand.grid(file=file.path(fil.dir, fil.nas), bulk='bulk', cell=sc.na.all, dimred=dimred, graph.meth=graph.meth, stringsAsFactors = FALSE)
  df.spd <- expand.grid(sim = sim, sim.p = sim.p, dim = dim, stringsAsFactors = FALSE)

  if (parallel.info==TRUE) {
    cat('\n Max first-level parallelizations across BatchtoolsParam workers: ', nrow(df.par.com), '\n', sep='')
    cat('\n Max second-level parallelizations across MulticoreParam workers (--cpus-per-task in slurm template) on each BatchtoolsParam worker: ', nrow(df.spd), '\n', sep='')
    return()
  }
  fun <- function(i, df.par.com, df.spd, df.match, sc.dim.min, max.dim, sim.meth, multi.core.par, auc.dir) {
    df0 <- df.par.com[i, ]
    df.para <- cbind(df.spd, df0, row.names=NULL)
    dat.fil <- readRDS(df0$file)
    blk <- dat.fil[['bulk']]; sc <- dat.fil[[df0$cell]]
    if (is.null(sc)) stop(paste0(df0$cell, ' is NULL in ', df0$file, '!'))  
    fil.na <- rev(strsplit(df.para[1, ]$file, '/')[[1]])[1]
    fil.na <- sub('\\.rds$', '', fil.na)
    out <- file.path(auc.dir, tolower(paste0('auc.', fil.na, '.', df0$cell, '.', df0$dimred, '.', df0$graph.meth)))
    # library(spatialHeatmap)
    df.para <- spatialHeatmap::cocluster(bulk=blk, cell=sc, df.match=df.match, df.para=df.para, sc.dim.min=sc.dim.min, max.dim=max.dim, sim.meth=sim.meth, return.all=FALSE, multi.core.par=multi.core.par, verbose=verbose, file=out); return('Done')
  }
  # Every row in df.par.com is combined with every row in df.spd for coclustering.
  if (is(batch.par, 'BatchtoolsParam')) {
    register(batch.par)
    res <- bplapply(seq_len(nrow(df.par.com)), fun, df.par.com=df.par.com, df.spd=df.spd, df.match=df.match, sc.dim.min=sc.dim.min, max.dim=max.dim, sim.meth=sim.meth, multi.core.par=multi.core.par,  auc.dir=auc.dir) 
  } else {
    res <- lapply(seq_len(nrow(df.par.com)), fun, df.par.com=df.par.com, df.spd=df.spd, df.match=df.match, sc.dim.min=sc.dim.min, max.dim=max.dim, sim.meth=sim.meth, multi.core.par=multi.core.par,  auc.dir=auc.dir)
  }
  return(res)
}



