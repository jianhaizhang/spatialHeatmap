#' Optimization of co-clustering bulk and single cell data 
#'
#' This function is specialized in optimizing the co-clustering method that is able to automatically assign bulk tissues to single cells. A vignette is provide at \url{https://jianhaizhang.github.io/spatialHeatmap_supplement/cocluster_optimize.html}. 

#' @param dat.lis A two-level nested \code{list}. Each inner \code{list} consists of three slots of \code{bulk}, \code{cell}, and \code{df.match}, corresponding to bulk data, single cell data, and ground-truth matching between bulk and cells respectively. For example, list(dataset1=list(bulk=bulk.data1, cell=cell.data1, df.match=df.match1), dataset2=list(bulk=bulk.data2, cell=cell.data2, df.match=df.match2)). 
#' @param df.para A \code{data.frame} with each row corresponding to a combination of parameter settings in co-clustering. 
#' @param df.fil.set A \code{data.frame} of filtering settings. E.g. data.frame(p=c(0.1, 0.2), A=rep(1, 2), cv1=c(0.1, 0.2), cv2=rep(50, 2), cutoff=rep(1, 2), p.in.cell=c(0.15, 0.2), p.in.gen=c(0.05, 0.1), row.names=paste0('fil', seq_len(2))).
#' @param batch.par The parameters for first-level parallelization through a cluster scheduler such as SLURM, which is \code{\link[BiocParallel]{BatchtoolsParam}}. If \code{NULL} (default), the first-level parallelization is skipped.
#' @param multi.core.par The parameters for second-level parallelization, which is \code{\link[BiocParallel]{MulticoreParam}}.
#' @param wk.dir The working directory, where results will be saved.
#' @param verbose If \code{TRUE}, intermediate messages will be printed. 

#' @return A \code{data.frame}.

#' @examples
#' 
#' # Optimization includes many iterative runs of co-clustering. To reduce runtime, these runs 
#' # are parallelized with the package BiocParallel. 
#' library(BiocParallel)
#' # To obtain reproducible results, a fixed seed is set for generating random numbers.
#' set.seed(10)
#' 
#' # Read bulk (S. Li et al. 2016) and two single cell data sets (Shahan et al. 2020), all of
#' # which are from Arabidopsis root.
#' blk <- readRDS(system.file("extdata/cocluster/data", "bulk_cocluster.rds", 
#' package="spatialHeatmap")) # Bulk.
#' sc10 <- readRDS(system.file("extdata/cocluster/data", "sc10_cocluster.rds", 
#' package="spatialHeatmap")) # Single cell.
#' sc11 <- readRDS(system.file("extdata/cocluster/data", "sc11_cocluster.rds", 
#' package="spatialHeatmap")) # Single cell.
#' blk; sc10; sc11
#' 
#' # The ground-truth matching between bulk tissue and single cells needs to be defined in form 
#' # of a table so as to classify TRUE/FALSE assignments.
#' match.pa <- system.file("extdata/cocluster/data", "true_match_arab_root_cocluster.txt", 
#' package="spatialHeatmap")
#' df.match.arab <- read.table(match.pa, header=TRUE, row.names=1, sep='\t')
#' df.match.arab[1:3, ]
#' 
#' # Place the bulk, single cell data, and matching table in a list.
#' dat.lis <- list(
#'   dataset1=list(bulk=blk, cell=sc10, df.match=df.match.arab), 
#'   dataset1=list(bulk=blk, cell=sc11, df.match=df.match.arab) 
#' )
#' 
#' # Filtering settings. 
#' df.fil.set <- data.frame(p=c(0.1), A=rep(1, 1), cv1=c(0.1), cv2=rep(50, 1), cutoff=rep(1, 1),
#' p.in.cell=c(0.15), p.in.gen=c(0.05), row.names=paste0('fil', seq_len(1))) 
#' # Settings in pre-processing include normalization method (norm), filtering (fil). The 
#' # following optimization focuses on settings most relevant to co-clustering, including 
#' # dimension reduction methods (dimred), number of top dimensions for co-clustering (dims), 
#' # graph-building methods (graph), clustering methods (cluster). Explanations of these settings
#' # are provide in the help file of function "cocluster".  
#' norm <- c('FCT'); fil <- c('fil1'); dimred <- c('UMAP')
#' dims <- seq(5, 10, 1); graph <- c('knn', 'snn')
#' cluster <- c('wt', 'fg', 'le')
#' 
#' df.para <- expand.grid(dataset=names(dat.lis), norm=norm, fil=fil, dimred=dimred, dims=dims, 
#' graph=graph, cluster=cluster, stringsAsFactors = FALSE)
#' 
#' \donttest{
#' # Optimization is performed by calling "coclus_opt", and results to a temporary directory 
#' # "wk.dir".
#' wk.dir <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
#' df.res <- coclus_opt(dat.lis, df.para, df.fil.set, multi.core.par=MulticoreParam(workers=1, 
#' RNGseed=50), wk.dir=wk.dir, verbose=TRUE)
#' df.res[1:3, ]
#' }


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022). _BiocParallel: Bioconductor facilities for parallel evaluation_. R package version 1.30.3, <https://github.com/Bioconductor/BiocParallel>.
#' Li, Song, Masashi Yamada, Xinwei Han, Uwe Ohler, and Philip N Benfey. 2016. "High-Resolution Expression Map of the Arabidopsis Root Reveals Alternative Splicing and lincRNA Regulation." Dev. Cell 39 (4): 508–22
#' Shahan, Rachel, Che-Wei Hsu, Trevor M Nolan, Benjamin J Cole, Isaiah W Taylor, Anna Hendrika Cornelia Vlot, Philip N Benfey, and Uwe Ohler. 2020. "A Single Cell Arabidopsis Root Atlas Reveals Developmental Trajectories in Wild Type and Cell Identity Mutants." BioRxiv.

#' @export
#' @importFrom BiocParallel bpRNGseed bpRNGseed<- register bpnworkers bplapply MulticoreParam

coclus_opt <- function(dat.lis, df.para, df.fil.set, batch.par=NULL, multi.core.par=MulticoreParam(workers=1, RNGseed=50), wk.dir, verbose=TRUE) {
  e.w <- tryCatch( # Check if the cluster scheduler is installed.
    expr = { is(batch.par, 'BatchtoolsParam') }, 
    error = function(e){ return('e') }, warning = function(w){ return('w') } 
  ) 
  # Check random seed.
  seed.bat <- seed.mul <- NULL
  if (!e.w %in% c('e', 'w')) if (is(batch.par, 'BatchtoolsParam')) seed.bat <- bpRNGseed(batch.par)
  if (is(multi.core.par, 'MulticoreParam')) seed.mul <- bpRNGseed(multi.core.par)
  if (is.numeric(seed.bat) & is.numeric(seed.mul)) { 
    bpRNGseed(multi.core.par) <- NULL
    if (verbose==TRUE) message('"RNGseed" in MulticoreParam is set NULL, since it is already set in BatchtoolsParam.')
  }
  if (!is.numeric(seed.bat) & !is.numeric(seed.mul)) {
    if (verbose==TRUE) message('"RNGseed" in BatchtoolsParam and MulticoreParam is NULL')
  }

  nor.dir.bp <- opt_dir(wk.dir, sub.dir='norm_res', batch.par, multi.core.par)
  nor.dir <- nor.dir.bp$dir
  batch.par <- nor.dir.bp$batch.par
  multi.core.par <- nor.dir.bp$multi.core.par
  df.para$dims <- as.numeric(df.para$dims)
  df.para.nor <- df.para[, c('dataset', 'norm')]
  df.para.nor <- df.para.nor[!duplicated(df.para.nor), ]
  df.para.nor$file.norm <- file.path(nor.dir, tolower(paste0(df.para.nor$norm, '.', df.para.nor$dataset, '.rds')))

  if (is(batch.par, 'BatchtoolsParam')) {
    register(batch.par); workers <- bpnworkers(batch.par)
    rows <- seq_len(nrow(df.para.nor))
    split <- split(rows, ceiling(seq_along(rows) / ceiling(length(rows)/workers)))
    # Cut: is not able to deal with workers=1.
    # split <- split(rows, cut(seq_along(rows), workers, labels = FALSE))
    res.nor <- bplapply(seq_along(split), norm_fun, dat.lis=dat.lis, df.para=df.para.nor, split=split, multi.core.par=multi.core.par, com=FALSE, nor.dir=nor.dir)
  } else {
    split <- list('1'=seq_len(nrow(df.para.nor))) 
    res.nor <- norm_fun(1, dat.lis, df.para=df.para.nor, split, multi.core.par, com=FALSE, nor.dir=nor.dir)
  }

  if (all(res.nor=='done') & verbose==TRUE) message('Normalization is done!')

  fil.dir.bp <- opt_dir(wk.dir, sub.dir='filter_res', batch.par, multi.core.par)
  fil.dir <- fil.dir.bp$dir
  batch.par <- fil.dir.bp$batch.par
  multi.core.par <- fil.dir.bp$multi.core.par
  
  df.para.fil <- df.para[, c('dataset', 'norm', 'fil')]
  df.para.fil <- df.para.fil[!duplicated(df.para.fil), ]
  df.para.fil$file.norm <- tolower(paste0(df.para.fil$norm, '.', df.para.fil$dataset, '.rds'))
  df.para.fil$file.fil <- tolower(paste0(df.para.fil$fil, '.', df.para.fil$file.norm))
  df.para.fil$file.norm <- file.path(nor.dir, df.para.fil$file.norm)
  df.para.fil$file.fil <- file.path(fil.dir, df.para.fil$file.fil)
  if (is(batch.par, 'BatchtoolsParam')) { 
    register(batch.par); workers <- bpnworkers(batch.par)
    rows <- seq_len(nrow(df.para.fil))
    split <- split(rows, ceiling(seq_along(rows) / ceiling(length(rows)/workers)))
    # split <- split(rows, cut(seq_along(rows), workers, labels = FALSE))
    res.fil <- bplapply(seq_along(split), fil_fun, df.para=df.para.fil, split=split, bulk.aggr=TRUE, aggr='mean', df.fil.set, multi.core.par=multi.core.par, com=FALSE, fil.dir=fil.dir)
  } else {
    split <- list('1'=seq_len(nrow(df.para.fil)))
    res.fil <- fil_fun(1, df.para=df.para.fil, split=split, bulk.aggr=TRUE, aggr='mean', df.fil.set, multi.core.par=multi.core.par, com=FALSE, fil.dir=fil.dir)
  }
  if (all(res.fil=='done') & verbose==TRUE) message('Filtering is done!')
  
  coclus.dir.bp <- opt_dir(wk.dir, sub.dir='coclus_res', batch.par, multi.core.par)
  coclus.dir <- coclus.dir.bp$dir
  batch.par <- coclus.dir.bp$batch.par
  multi.core.par <- coclus.dir.bp$multi.core.par

  df.para.coclus <- df.para[, c('dataset', 'norm', 'fil', 'dimred', 'dims', 'graph', 'cluster')]
  df.para.coclus <- df.para.coclus[!duplicated(df.para.coclus), ]
  df.para.coclus$file.fil <- tolower(paste0(df.para.coclus$fil, '.', df.para.coclus$norm, '.', df.para.coclus$dataset, '.rds'))
  df.para.coclus$file.coclus <- tolower(paste0(df.para.coclus$cluster, '.', df.para.coclus$graph, '.', 'dims', df.para.coclus$dims, '.', df.para.coclus$dimred, '.', df.para.coclus$file.fil))
  df.para.coclus$file.fil <- file.path(fil.dir, df.para.coclus$file.fil)
  df.para.coclus$file.coclus <- file.path(coclus.dir, df.para.coclus$file.coclus)

  if (is(batch.par, 'BatchtoolsParam')) {
    register(batch.par); workers <- bpnworkers(batch.par)
    rows <- seq_len(nrow(df.para.coclus))
    split <- split(rows, ceiling(seq_along(rows) / ceiling(length(rows)/workers)))
    # split <- split(rows, cut(seq_along(rows), workers, labels = FALSE))
    res.coclus <- bplapply(seq_along(split), coclus_fun, dat.lis=dat.lis, df.para=df.para.coclus, split=split, multi.core.par=multi.core.par, coclus.dir=coclus.dir)
    df.res <- do.call('rbind', lapply(res.coclus, function(x) do.call("rbind", x)))
  } else {
    split <- list('1'=seq_len(nrow(df.para.coclus)))
    res.coclus <- coclus_fun(1, dat.lis=dat.lis, df.para=df.para.coclus, split=split, multi.core.par=multi.core.par, coclus.dir=coclus.dir) 
    df.res <- do.call('rbind', res.coclus)
  }
  if (verbose==TRUE) message('Co-clustering is done!')
  saveRDS(df.res, file=file.path(coclus.dir, 'df_res.rds'))
  df.res
}

#' Creating directories for each step in co-clustering optimization
#'
#' @param sce.coclus The coclustered bulk and single cell data in a \code{SingleCellExperiment}, where cocluster assignments are stored in the \code{cluster} column in \code{colData}.

#' @return A list of \code{roc} object and the data frame to create the \code{roc}.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#' Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022). _BiocParallel: Bioconductor facilities for parallel evaluation_. R package version 1.30.3, <https://github.com/Bioconductor/BiocParallel>.

#' @importFrom parallel detectCores
#' @importFrom BiocParallel bpnworkers bplog bplog<- bplogdir<-

opt_dir <- function(wk.dir, sub.dir, batch.par, multi.core.par) {
  dir <- file.path(wk.dir, sub.dir)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE) 
  if (is(multi.core.par, 'MulticoreParam')) {
    cpus <- detectCores(); workers <- bpnworkers(multi.core.par)
    if (workers <= 0) stop('The minimum worker(s) is 1 !')
    if (cpus - workers < 1) stop('The maximum worker(s) should be 1 less than all available workers (parallel::detectCores()) !')
    bplog(multi.core.par) <- TRUE 
    mcore.logdir <- file.path(dir, 'multi_core_log') 
    if (!dir.exists(mcore.logdir)) dir.create(mcore.logdir, recursive = TRUE) 
    bplogdir(multi.core.par) <- mcore.logdir
  }
  if (is(batch.par, 'BatchtoolsParam')) {
    if (bplog(batch.par)==FALSE) bplog(batch.par) <- TRUE
    bat.log.dir <- file.path(dir, 'batch_log')
    if (!dir.exists(bat.log.dir)) dir.create(bat.log.dir, recursive = TRUE)
    bplogdir(batch.par) <- bat.log.dir
  }; return(list(dir=dir, batch.par=batch.par, multi.core.par=multi.core.par))
}

#' Normalization function for co-clustering optimization
#'
#' @param sce.coclus The coclustered bulk and single cell data in a \code{SingleCellExperiment}, where cocluster assignments are stored in the \code{cluster} column in \code{colData}.

#' @return A list of \code{roc} object and the data frame to create the \code{roc}.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). "Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R." _Bioinformatics_, *33*, 1179-1186. doi:10.1093/bioinformatics/btw777 <https://doi.org/10.1093/bioinformatics/btw777>.
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). "Orchestrating single-cell analysis with Bioconductor." _Nature Methods_, *17*, 137-145. <https://www.nature.com/articles/s41592-019-0654-x>.
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/

#' @importFrom scuttle calculateCPM
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData colData<-

norm_opt <- function(cell, bulk, norm, com=FALSE) {
  bulkCell <- NULL
  if (norm %in% 'FCT') {
    # set.seed(10)
    res <- norm_cell(sce=cell, bulk=bulk, cpm=FALSE, com=com) 
  } else if (norm %in% c('CPM', 'TMM', 'VST')) {
    if (is(cell, 'data.frame')|is(cell, 'matrix')|is(cell, 'Matrix')) cell <- SingleCellExperiment(assays=list(counts=as.matrix(cell)))
    if (is(cell, 'SummarizedExperiment')) cell <- as(cell, 'SingleCellExperiment')  
    if (is(bulk, 'data.frame')|is(bulk, 'matrix')|is(bulk, 'Matrix')) bulk <- SingleCellExperiment(assays=list(counts=as.matrix(bulk)))
    if (is(bulk, 'SummarizedExperiment')) bulk <- as(bulk, 'SingleCellExperiment')  
    colData(cell) <- colData(bulk) <- NULL
    bulk$bulkCell <- 'bulk'; cell$bulkCell <- 'cell'
    bulk$sample <- colnames(bulk); cell$sample <- colnames(cell)
    int <- intersect(rownames(bulk), rownames(cell)) 
    sce <- cbind(bulk[int, ], cell[int, ]) 
    colnames(sce) <- seq_len(ncol(sce))  
    if ('CPM' %in% norm) { 
      cnt.cpm <- calculateCPM(sce); nor <- log2(cnt.cpm+1)
      nor <- SingleCellExperiment(assays=list(logcounts=as.matrix(nor)))
     colData(nor) <- colData(sce)
    } else if ('TMM' %in% norm) {
      nor <- norm_data(data=sce, norm.fun='CNF', parameter.list=list(method='TMM'))
    } else if ('VST' %in% norm) {
      nor <- norm_data(data=sce, norm.fun='VST')
    }
    colnames(nor) <- nor$sample
    assayNames(nor) <- 'logcounts'; res <- nor
    if (com==FALSE) {
      bulk <- subset(nor, , bulkCell=='bulk')
      cell <- subset(nor, , bulkCell=='cell')
      res <- list(bulk=bulk, cell=cell)
    }
  }; return(res)
}

#' Normalization function for co-clustering optimization
#'
#' @param sce.coclus The coclustered bulk and single cell data in a \code{SingleCellExperiment}, where cocluster assignments are stored in the \code{cluster} column in \code{colData}.

#' @return A list of \code{roc} object and the data frame to create the \code{roc}.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022). _BiocParallel: Bioconductor facilities for parallel evaluation_. R package version 1.30.3, <https://github.com/Bioconductor/BiocParallel>.

#' @importFrom BiocParallel bplapply

norm_fun <- function(i, dat.lis, df.para, split, multi.core.par, com=FALSE, nor.dir=NULL) {
  df.para0 <- df.para[split[[i]], ]
  lis <- bplapply(seq_len(nrow(df.para0)), BPPARAM=multi.core.par, FUN=function(j) {
  # lis <- bplapply(seq_len(nrow(df.para0)), FUN=function(j) {
  cell <- dat.lis[[df.para0$dataset[j]]]$cell
  bulk <- dat.lis[[df.para0$dataset[j]]]$bulk
  norm <- df.para0$norm[j]
  res <- norm_opt(cell=cell, bulk=bulk, norm, com=FALSE)
  if (!is.null(nor.dir)) { 
    if (!dir.exists(nor.dir)) dir.create(nor.dir, recursive = TRUE)
    saveRDS(res, file=df.para0$file.norm[j])
  }; return(res)
  }); return('done')
}


#' Filtering function for co-clustering optimization
#'
#' @param sce.coclus The coclustered bulk and single cell data in a \code{SingleCellExperiment}, where cocluster assignments are stored in the \code{cluster} column in \code{colData}.

#' @return A list of \code{roc} object and the data frame to create the \code{roc}.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022). _BiocParallel: Bioconductor facilities for parallel evaluation_. R package version 1.30.3, <https://github.com/Bioconductor/BiocParallel>.

#' @importFrom BiocParallel bplapply

fil_fun <- function(i, df.para, split, bulk.aggr=TRUE, aggr='mean', df.fil.set, multi.core.par, com=FALSE, fil.dir=NULL) { 
  df.para0 <- df.para[split[[i]], ]
  lis <- bplapply(seq_len(nrow(df.para0)), BPPARAM=multi.core.par, FUN=function(j) {
  # lis <- bplapply(seq_len(nrow(df.para0)), FUN=function(j) {
    blk.cell <- readRDS(df.para0$file.norm[j])
    if (bulk.aggr==TRUE) blk.aggr <- aggr_rep(data=blk.cell$bulk, assay.na='logcounts', sam.factor='sample', aggr=aggr)
    df.fil0 <- df.fil.set[df.para0$fil[j], ]
    # Filter bulk
    blk.fil <- filter_data(data=blk.aggr, pOA=c(df.fil0$p, df.fil0$A), CV=c(df.fil0$cv1, df.fil0$cv2), verbose=FALSE) 
    # Filter cell and subset bulk to genes in cell
    blk.sc.fil <- filter_cell(sce=blk.cell$cell, bulk=blk.fil, cutoff=df.fil0$cutoff, p.in.cell=df.fil0$p.in.cell, p.in.gen=df.fil0$p.in.gen, verbose=FALSE) 
    saveRDS(blk.sc.fil, file=df.para0$file.fil[j])
  }); return('done')
}


#' Co-clustering function for co-clustering optimization
#'
#' @param sce.coclus The coclustered bulk and single cell data in a \code{SingleCellExperiment}, where cocluster assignments are stored in the \code{cluster} column in \code{colData}.

#' @return A list of \code{roc} object and the data frame to create the \code{roc}.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022). _BiocParallel: Bioconductor facilities for parallel evaluation_. R package version 1.30.3, <https://github.com/Bioconductor/BiocParallel>.
#' Xavier Robin, Natacha Turck, Alexandre Hainard, Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics, 12, p. 77.  DOI: 10.1186/1471-2105-12-77 <http://www.biomedcentral.com/1471-2105/12/77/>

#' @importFrom BiocParallel bplapply

coclus_fun <- function(i, dat.lis=dat.lis, df.para, split, multi.core.par, coclus.dir=NULL) {
  df.para0 <- df.para[split[[i]], ]
  lis <- bplapply(seq_len(nrow(df.para0)), BPPARAM=multi.core.par, FUN=function(j) {
  # lis <- lapply(seq_len(nrow(df.para0)), FUN=function(j) {

    df0 <- df.para0[j, ]; print(df0); fil.file <- df0$file.fil
    df0$auc <- df0$total.assignment <- df0$true.assignment <- df0$threshold <- df0$sensitivity <- df0$specificity <- df0$accuracy <- 0
    df.match <- dat.lis[[gsub('.*(dataset\\d+).*', '\\1', fil.file)]]$df.match
    dat.fil <- readRDS(fil.file)
    df0$file.fil <- df0$file.coclus <- NULL
    bulk <- dat.fil$bulk; cell <- dat.fil$cell 
    if (nrow(bulk) <= 4) return(df0)
    # set.seed(10)
    res <- cocluster(bulk=bulk, cell=cell, min.dim=df0$dims, dimred=df0$dimred, graph.meth=df0$graph, cluster=df0$cluster, df.match=df.match)
    roc.obj <- res$roc.obj 
    if (!is.null(roc.obj)) {
      if (any(c('e', 'w') %in% check_pkg('pROC'))) stop('The package "pROC" is not detected!') 
      best <- round(pROC::coords(roc.obj, x='best', ret=c("threshold", "specificity", "sensitivity", "accuracy")), 3)
      if (nrow(best)>0) { 
        best <- best[1, ]  
        if (unique(max(abs(best$threshold)))!=Inf) {  
          cdat <- colData(res$sce.all)
          tab <- table(cdat$response)
          true <- tab['TRUE']; false <- tab['FALSE']
          if (is.na(true)) true <- 0
          if (is.na(false)) false <- 0
          df0$total.assignment <- true+false
          df0$true.assignment <- true
          df0$auc <- round(pROC::auc(roc.obj), 3)
          df0$threshold <- round(best$threshold, 3)  
          df0$sensitivity <- round(best$sensitivity, 3) 
          df0$specificity <- round(best$specificity, 3)  
          df0$accuracy <- round(best$accuracy, 3) 
        }   
      }    
    } 
    # df0 <- df0[, !grepl('^file\\.', colnames(df0))]
    return(df0)
  }); 
  df.res <- do.call('rbind', lis)
  saveRDS(df.res, file=file.path(coclus.dir, paste0('res', i, '.rds')))
  return(lis)
}


