#' Extract AUC statistics in coclustering optimization 
#'
#' Extract AUC statistics of target parameters in coclustering optimization.

#' @param wk.dir The working directory of coclustering.
#' @param tar.par The target parameter in optimization, one of \code{norm}, \code{filter}, \code{dimred}, \code{graph}, and \code{spd.set} corresponding to normalization methods, filtering parameter sets, dimensionality reduction methods, graph-building methods, and \code{sim}/\code{sim.p}/\code{dim} sets respectively.
#' @param total.min,true.min,aucs Cutoffs to extract AUCs. AUCs over a cutoff and having total bulk tissue assignments above \code{total} and true assignments above \code{true.min} are extracted. The default is \code{total.min=500}, \code{true.min=300}, \code{round(seq(0.5, 0.9, 0.1), 1)}. Each AUC cuttoff in \code{aucs} is used independently. 

#' @return An nested \code{list}.

#' @examples

#' # To obtain reproducible results, always start a new R session and set a fixed seed for Random Number Generator at the begaining, which is required only once in each R session.  
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
#' }

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}


#' @export auc_stat

auc_stat <- function(wk.dir, tar.par='norm', total.min=500, true.min=300, aucs=round(seq(0.5, 0.9, 0.1), 1)) { 
  dir <-file.path(wk.dir, 'auc_res')
  files <- list.files(file.path(wk.dir, 'auc_res'), '^auc.*\\.rds$')
  files <- sort(files)
  # All component names in ".rds" files.
  nas <- unique(unlist(strsplit(files, '\\.')))

  par.lis <- list(norm=c('fct', 'cpm'), dimred=c('pca', 'umap'), graph=c('knn', 'snn'), filter=grep('^fil\\d+', nas, value=TRUE))
  cell <- setdiff(nas, c('auc', 'dat', 'rds', unlist(par.lis)))
  if (tar.par[1]!='spd.set') tar.par <- par.lis[[tar.par]]

  # Frequency of "sim", "sim.p", "dim" at a single auc.thr.
  frq_spd <- function(files=files, auc.min=auc.min, total.min=total.min, true.min=true.min, dir) { 
    df.all <- NULL; for (i in files) {
      pa <- file.path(dir, i)
      if (!file.exists(pa)) {
        cat(pa, ' does not exist! \n', sep=''); next
      }; print(i)
      df0 <- sub_auc(pa, auc.min=auc.min, total.min=total.min, true.min=true.min)
      if (nrow(df0)==0) next
      df.all <- rbind(df.all, df0)
    }; print(dim(df.all))
    df.frq <- df.all
    # Combine "sim", "sim.p", "dim" and count freqeuncy.
    df.frq$spd.set <- paste0(df.frq$sim, df.frq$sim.p, df.frq$dim)
    tab <- table(df.frq$spd.set)
    # Add frequency column.
    df.frq$frequency <- tab[df.frq$spd.set]
    df.frq <- df.frq[order(df.frq$frequency, decreasing=TRUE), ]
    df.frq <- df.frq[!duplicated(df.frq$spd.set), ]
    df.frq$auc <- auc.min; df.frq$total <- total.min
    df.frq$true <- true.min
    df.frq <- df.frq[, setdiff(colnames(df.frq), c('thr', 'sens', 'spec', 'spd.set'))]
    cna <- colnames(df.frq)
    idx <- cna %in% c('auc', 'total', 'true') 
    cna[idx] <- paste0(cna[idx], '.thr')
    colnames(df.frq) <- cna
    return(list(df.frq=df.frq, df.all=df.all))
  }
 
  # The auc names are reserved and will be recognized in violin plots.
  names(aucs) <- paste0('auc', aucs)
  if (tar.par[1]=='spd.set') { # sim, sim.p, dim combination.
    df.lis <- lapply(aucs, function(i) frq_spd(files=files, auc.min=i, total.min=total.min, true.min=true.min, dir=dir))
    return(df.lis)
  }
  # auc stats of a single target.
  par_stat0 <- function(tar.par, files=files, auc.min=auc.min, total.min=total.min, true.min=true.min, dir) {
    # All target files.
    file.tar <- grep(paste0('\\.', tar.par, '\\.'), files, value=TRUE)
    df.all <- NULL
    for (i in file.tar) { pa <- file.path(dir, i)
      if (!file.exists(pa)) {
        cat(pa, ' does not exist! \n', sep=''); next
      }; print(paste0(i, ': auc.min ', auc.min))
      df0 <- sub_auc(pa, auc.min=auc.min, total.min=total.min, true.min=true.min)
      if (nrow(df0)==0) next
      df0$file <- i; df.all <- rbind(df.all, df0)
    }; print(dim(df.all)); return(df.all)
  }

  # The df names are reserved and will be recognized in violin plots.
  names(tar.par) <- paste0('df.', tar.par)
  # auc stats of a multiple targets at multiple auc.thrs.
  df.lis <- lapply(aucs, function(i) lapply(tar.par, function(j) par_stat0(tar.par=j, files=files, auc.min=i, total.min=total.min, true.min=true.min, dir=dir))
  )
  # Mean auc of each target-auc.thr pair.
  df.mean <- NULL; for (i in seq_along(df.lis)) {
    df.lis0 <- df.lis[[i]]
    for (j in seq_along(df.lis0)) { 
      # df.lis0 names have the same order of target.
      # df.lis names have the same order of aucs.
      df0 <- df.lis0[[j]] 
      df.mn <- data.frame(parameter=tar.par[j], auc.thr=aucs[i], mean=ifelse(is.null(df0), 0, mean(df0$auc)))
      df.mean <- rbind(df.mean, df.mn)
    }
  }; rownames(df.mean) <- seq_len(nrow(df.mean))
  df.lis <- c(list(df.auc.mean=df.mean), df.lis)
  return(df.lis)
}



#' Extract AUCs in a single file 
#'
#' Extract AUCs in a single ".rds" file generated in coclustering.

#' @param path The path of ".rds" file.
#' @param auc.min,total.min,true.min The cutoff of AUC, total bulk tissue assignments, and true assignments respectively.

#' @return A \code{data.frame}. 
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

sub_auc <- function(path, auc.min=0.1, total.min=300, true.min=100) {
  total <- true <- NULL
  df.par <- readRDS(path)
  df.par <- subset(df.par, auc >= auc.min & total >= total.min & true >= true.min)
  print(dim(df.par)); print(mean(df.par$auc))
  return(df.par)
}
