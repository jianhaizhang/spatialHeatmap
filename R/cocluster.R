#' Coclustering bulk and single cell data in a single run
#'
#' Cluster single cells, refine single cell clusters, and cocluster bulk and single cell data. Accept multiple parameter settings combinations in form of \code{data.frame}, and allows for parallelization.

#' @param bulk Normalized and filtered bulk data at log2 scale in form of \code{matrix}, \code{data.frame}, \code{SingleCellExperiment}, or \code{SummarizedExperiment}.
#' @param cell Normalized and filtered bulk data at log2 scale in form of \code{matrix}, \code{data.frame}, \code{SingleCellExperiment}, or \code{SummarizedExperiment}.
#' @inheritParams coclus_roc

#' @param df.para A \code{data.frame} of paramerter settings in coclustering, where the parameter names are the column names and one row denotes one combination of settings. Missing parameters in the \code{data.frame} are replaced by their default settings internally. E.g. In \code{df.para = data.frame(sim = 0.2, sim.p = 0.8, dim = 12, graph.meth = 'knn'))}, the default settings of \code{sc.dim.min = 10; max.dim=50; dimred='PCA'; sim.meth='spearman'} are used internally. If multiple settings combinations are contained in multiple rows respectively, parallelization can be enabled through \code{multi.core.par}.
#' @param sc.dim.min Integer scalar specifying the minimum number of (principle components) PCs to retain in \code{\link[scran]{denoisePCA}} when clustering single cells without bulk data. The default is \code{10}.
#' @param max.dim Integer scalar specifying the maximum number of (principle components) PCs to retain in \code{\link[scran]{denoisePCA}} when clustering single cells without bulk data and coclustering single cells and bulk data. The default is \code{50}.
#' @inheritParams refine_cluster

#' @param dim Integer scalar specifying the minimum number of (principle components) PCs to retain in \code{\link[scran]{denoisePCA}} when coclustering single cells and bulk data. The default is \code{12}.
#' @param return.all Logical. If \code{TRUE}, single cell data after refining cluster, \code{roc} object and the \code{data.frame} to create the \code{roc} during coclustering are returned in a nested list. If \code{FALSE} (default), the \code{df.para}  table including coclustering statisctics are returned. \code{auc} denotes area under the curve. \code{true} and \code{total} indicate the total true assignments of bulk tissues and total assignments (true and false) respectively. \code{thr}, \code{spec}, and \code{sens} refer to the optimal similarity threshold between bulk and cells in coclusters, specificity, and sensitivity corresponding to \code{thr} respectively. 
#' @param multi.core.par The parallelization settings. Default is \code{MulticoreParam(workers=1, stop.on.error=FALSE, log=FALSE, logdir=NULL)}. See \code{\link[BiocParallel]{MulticoreParam}}.
#' @param file A file name without extension to save the table of parameter settings and coclustering statisctics if \code{return.all = FALSE}. The table is saved by \code{saveRDS} with extension \code{.rds}. Default is \code{NULL} and no file is saved.

#' @return A nested list or a table of coclustering results. 

#' @examples

#' # To obtain reproducible results, always start a new R session and set a fixed seed for Random Number Generator at the beginning, which is required only once in each R session.  
#' set.seed(10)
#' 
#' # Example bulk data of mouse brain for coclustering (Vacher et al 2021).
#' blk.mus.pa <- system.file("extdata/shinyApp/example", "bulk_mouse_cocluster.txt", package="spatialHeatmap") 
#'blk.mus <- as.matrix(read.table(blk.mus.pa, header=TRUE, row.names=1, sep='\t', check.names=FALSE))
#' blk.mus[1:3, 1:5]
#'
#' # Example single cell data for coclustering (Ortiz et al 2020).
#' sc.mus.pa <- system.file("extdata/shinyApp/example", "cell_mouse_cocluster.txt", package="spatialHeatmap") 
#'sc.mus <- as.matrix(read.table(sc.mus.pa, header=TRUE, row.names=1, sep='\t', check.names=FALSE))
#' sc.mus[1:3, 1:5]
#'
#' # Initial filtering. 
#' blk.mus <- filter_data(data=blk.mus, sam.factor=NULL, con.factor=NULL, pOA=c(0.1, 5), CV=c(0.2, 100), dir=NULL) 
#' dim(blk.mus)
#' mus.lis <- filter_cell(lis=list(sc.mus=sc.mus), bulk=blk.mus, gen.rm=NULL, min.cnt=1, p.in.cell=0.5, p.in.gen=0.1) 

#' \donttest{

#' # Normalization: bulk and single cell are combined and normalized, then separated.
#' mus.lis.nor <- norm_multi(dat.lis=mus.lis, cpm=FALSE)
#'
#' # Secondary filtering.
#' library(SingleCellExperiment)
#' blk.mus.fil <- filter_data(data=logcounts(mus.lis.nor$bulk), sam.factor=NULL, con.factor=NULL, pOA=c(0.1, 0.5), CV=c(0.2, 100)) 
#' dim(blk.mus.fil)
#' 
#' mus.lis.fil <- filter_cell(lis=list(sc.mus=logcounts(mus.lis.nor$sc.mus)), bulk=blk.mus.fil, gen.rm=NULL, min.cnt=1, p.in.cell=0.05, p.in.gen=0.02)
#'
#' # The aSVG file of mouse brain.
#' svg.mus <- system.file("extdata/shinyApp/example", "mus_musculus.brain.svg", package="spatialHeatmap")
#' # Spatial features.  
#' feature.df <- return_feature(svg.path=svg.mus) 
#'
#' # Matching table indicating true bulk tissues of each cell type and corresponding SVG bulk (spatial feature).
#' df.match.mus.pa <- system.file("extdata/shinyApp/example", "match_mouse_brain_cocluster.txt", package="spatialHeatmap")
#' df.match <- read.table(df.match.mus.pa, header=TRUE, row.names=1, sep='\t')
#' df.match
#'
#' # The SVG bulk tissues are in the aSVG file.  
#' df.match$SVGBulk %in% feature.df$feature
#'
#' # Cluster single cells.
#' clus.sc <- cluster_cell(data=mus.lis.fil$sc.mus, min.dim=10, max.dim=50, graph.meth='knn', dimred='PCA')
#' # Cluster labels are stored in the "cluster" column in "colData".
#'colData(clus.sc)[1:3, ]
#'
#' # Refine cell clusters.
#' cell.refined <- refine_cluster(clus.sc, sim=0.2, sim.p=0.8, sim.meth='spearman')
#'
#' # Include matching information in "colData".
#' cell.refined <- true_bulk(cell.refined, df.match)
#' colData(cell.refined)[1:3, ]
#'
#' # Cocluster bulk and single cells.
#' roc.lis <- coclus_roc(bulk=mus.lis.fil$bulk, cell.refined=cell.refined, df.match=df.match, min.dim=13, max.dim=50, graph.meth='knn', dimred='PCA', sim.meth='spearman') 
#'
#' # The colustering results. "predictor" is the similarity between bulk and cells within a co-cluster. "index" is the cell index in the "SingleCellExperiment" after cell clusters are refined. 
#' roc.lis$df.roc[1:3, ]
#' # ROC curve created according to "roc.lis$df.roc". 
#' plot(roc.lis$roc.obj, print.auc=TRUE)
#' # Incorporate "cell.refined" in "roc.lis" for downstream use in co-visualization.
#' res.lis <- c(list(cell.refined=cell.refined), roc.lis)
#' 
#' # The processes of clustering single cells, refining cell clusters, and coclustering bulk and single cells can be performed by a single function.
#' library(BiocParallel)
#' res.lis <- cocluster(bulk=mus.lis.fil$bulk, cell=mus.lis.fil$sc.mus, df.match=df.match, df.para=NULL, sim=0.2, sim.p=0.8, dim=13, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=TRUE)
#' res.lis <- res.lis[[1]]
#'
# Since "return.all" is "TRUE", the refined cell clusters, ROC object, and coclustering results are returned.  
#' names(res.lis)
#'
#' # "cocluster" accepts multiple combinations of parameter settings provided in a data frame, and coclustering on these combinations can be performed in parallel through "multi.core.par".
#' # Multiple combinations of parameter settings. If some parameters are not specified in this table such as "graph.meth", their default settings will be used. 
#' df.par <- data.frame(sim=c(0.2, 0.3), sim.p=c(0.8, 0.7), dim=c(12, 13))
#'
#' # The computation is parallelized on 2 cpu cores by "multi.core.par".
#' library(BiocParallel)
#' res.multi <- cocluster(bulk=mus.lis.fil$bulk, cell=mus.lis.fil$sc.mus, df.match=df.match, df.para=df.par, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=13, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=TRUE, multi.core.par=MulticoreParam(workers=2))
#'
#' # The results of auto-matching through coclustering can be tailored through "Lasso Select" on the convenience Shiny app (desired_bulk_shiny) or manually defining desired bulk. If the former, save "cell.refined" in an ".rds" file by "saveRDS(cell.refined, file='cell.refined.rds')" and upload "cell.refined.rds" to the Shiny app. 
#' df.desired.bulk <- NULL 
#' # Example of desired bulk downloaded from convenience Shiny app.
#' desired.blk.pa <- system.file("extdata/shinyApp/example", "selected_cells_with_desired_bulk.txt", package="spatialHeatmap")
#' df.desired.bulk <- read.table(desired.blk.pa, header=TRUE, row.names=1, sep='\t')
#' df.desired.bulk[1:3, ]
#'
#' # Desired bulk manually definded by x-y coordinates ranges.
#' plot_dim(res.lis$cell.refined, dim='PCA', color.by='cell', x.break=seq(-10, 10, 2), y.break=seq(-10, 10, 2))
#'
#' df.desired.bulk <- data.frame(x.min=c(2, 6), x.max=c(4, 10), y.min=c(6, 8), y.max=c(8, 10), desiredSVGBulk=c('cerebral.cortex', 'cerebral.cortex'), dimred='PCA')
#' df.desired.bulk
#'
#' # Extract cells with true bulk assignments. If "df.desired.bulk" is provided, the corrresponding assignments are incorporated and regarded true.
#' sce.lis <- sub_asg(res.lis=res.lis, df.desired.bulk=df.desired.bulk, df.match=df.match, true.only=TRUE)
#' 
#' # Aggregate the extracted cells and the aggregated abundance profiles are used to color matching bulk tissues in the aSVG image.
#' sce.aggr <- aggr_rep(data=sce.lis$cell.sub, assay.na='logcounts', sam.factor='SVGBulk', con.factor=NULL, aggr='mean')
#'
#' # Co-visualize bulk and single cells without abundance profiles.
#' shm.lis1 <- spatial_hm(svg.path=svg.mus, data=sce.aggr, ID=c('Adcy1'), legend.nrow=4, sce.dimred=sce.lis$cell.refined, dimred='PCA', tar.bulk=c('hippocampus'), profile=FALSE, dim.lgd.text.size=10, dim.lgd.nrow=1)
#'
#' # Co-visualize bulk and single cells with abundance profiles of gene "Adcy1".
#' shm.lis2 <- spatial_hm(svg.path=svg.mus, data=sce.aggr, ID=c('Adcy1'), legend.nrow=4, sce.dimred=sce.lis$cell.refined, dimred='PCA', assay.na='logcounts', tar.bulk=c('hippocampus'), profile=TRUE, dim.lgd.text.size=10, dim.lgd.nrow=1, bar.width=0.1)

#'}

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2021). BiocParallel: Bioconductor facilities for parallel evaluation. R package version 1.28.3, https://github.com/Bioconductor/BiocParallel.
#' R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Xavier Robin, Natacha Turck, Alexandre Hainard, Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics, 12, p. 77.  DOI: 10.1186/1471-2105-12-77 <http://www.biomedcentral.com/1471-2105/12/77/>
#' Vacher CM, Lacaille H, O'Reilly JJ, Salzbank J et al. Placental endocrine function shapes cerebellar development and social behavior. Nat Neurosci 2021 Oct;24(10):1392-1401. PMID: 34400844.
#' Ortiz C, Navarro JF, Jurek A, Märtin A et al. Molecular atlas of the adult mouse brain. Sci Adv 2020 Jun;6(26):eabb3446. PMID: 32637622
#' R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

#' @export cocluster
#' @importFrom BiocParallel bpnworkers bpworkers bpworkers<- bplapply
#' @importFrom parallel detectCores
#' @importFrom SummarizedExperiment assay 
#' @importFrom SingleCellExperiment logcounts 
#' @importFrom pROC auc coords 
#' @importFrom methods as 

cocluster <- function(bulk, cell, df.match, df.para=NULL, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=13, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=FALSE, multi.core.par=MulticoreParam(workers=1, stop.on.error=FALSE, log=FALSE), verbose=TRUE, file=NULL) {
  sc.par.com <- NULL
  # if (!dir.exists('./multi_core_log')) dir.create('./multi_core_log')
  cpus <- detectCores(); workers <- bpnworkers(multi.core.par)
  if (workers <= 0) stop('The minimum worker(s) is 1 !')
  if (cpus - workers < 1) stop('The maximum worker(s) should be 1 less than all available workers (parallel::detectCores()) !')
  if (is(bulk, 'SingleCellExperiment')) bulk <- logcounts(bulk)
  if (is(bulk, 'SummarizedExperiment')) bulk <- assay(bulk)
  if (is(cell, 'SingleCellExperiment')) cell <- logcounts(cell)
  if (is(cell, 'SummarizedExperiment')) cell <- assay(cell)
  if (all(round(bulk)==bulk)) stop('The bulk data should be in log2 scale!')
  if (all(round(cell)==cell)) stop('The single cell data should be in log2 scale!')
  inter <- intersect(rownames(bulk), rownames((cell)))
  bulk <- bulk[inter, , drop=FALSE]; cell <- cell[inter, , drop=FALSE]
  if (length(intersect(colnames(bulk), colnames(cell))) > 0) stop('Overlap between column names of bulk and single cell data is detected!')
  
  if (!is.null(df.para)) df.para <- as(df.para, 'data.frame')
  pars <- c(sc.dim.min=sc.dim.min, max.dim=max.dim, sim=sim, sim.p=sim.p, dim=dim, graph.meth=graph.meth, dimred=dimred, sim.meth=sim.meth)
  # Include missing parameters in the provided df.para.
  par.na <- setdiff(names(pars), colnames(df.para))
  if (length(par.na) > 0) {
    df0 <- data.frame(as.list(pars[par.na]))
    if (!is.null(df.para)) df.para <- cbind(df.para, df0) else df.para <- df0
  }
  par.num <- c('sc.dim.min', 'max.dim', 'sim', 'sim.p', 'dim')
  df.para$auc <- df.para$total <- df.para$true <- df.para$thr <- df.para$sens <- df.para$spec <- as.numeric(0)
  for (i in par.num) df.para[, i] <- as.numeric(df.para[, i])

  tmp.dir <- file.path(tempdir(check = TRUE), 'cocluster_para')
  if (!dir.exists(tmp.dir)) dir.create(tmp.dir)
  message('Temporary directory: ', tmp.dir)
  # Parallelization.
  # par.tar <- setdiff(colnames(df.para), c('sc.dim.min', 'max.dim'))

  # Divide parameters in df.para by same parameter combinations on clustering single cells. Multiple co-clusterings are under the same clustering on single cells instead of running the same single cell clustering for every downstream coclustering, which is time efficient. 
  df.para$sc.par.com <- paste0(df.para$sc.dim.min, df.para$max.dim,df.para$graph.meth, df.para$dimred)
  df.para$index <- seq_len(nrow(df.para))
  cna.val <- setdiff(colnames(df.para), 'sc.par.com')
  df.para.all <- data.frame(matrix(ncol=length(cna.val), nrow=0, dimnames=list(NULL, cna.val)))
  lis.all <- NULL
  for (i in unique(df.para$sc.par.com)) {
    df.para0 <- subset(df.para, sc.par.com == i)
    df.para0 <- df.para0[, cna.val, drop=FALSE]
    workers <- nrow(df.para0)
    clus.sc <- cluster_cell(data=cell, min.dim=df.para0$sc.dim.min[1], max.dim=df.para0$max.dim[1], pca=FALSE, graph.meth=df.para0$graph.meth[1], dimred=df.para0$dimred[1])
    if (workers < bpworkers(multi.core.par)) bpworkers(multi.core.par) <- workers
    # Split para combinations by # of workers so that the intermediate auc results could be saved more frequently. E.g. nrow(df.para0) is 100, workers are 2. If no splitting, the auc results are saved when every 50 cocusterings are done. If split, auc results are saved when every 2 cocusterings are done.  
    rows <- seq_len(nrow(df.para0)); row.gr <- split(rows, ceiling(seq_along(rows)/workers))
    for (k in seq_along(row.gr)) {
    df.para.gr <- df.para0[row.gr[[k]], , drop=FALSE]
    lis <- bplapply(seq_len(nrow(df.para.gr)), BPPARAM=multi.core.par, FUN=function(j)
    {
      df0 <- df.para.gr[j, , drop=FALSE]; if (verbose==TRUE) print(df0)
      # Cluster single cell only. min.dim/max.dim is independent from coclustering sc+bulk, which is time-efficient.
      #if (!is.null(sc.dim.min)) {
      if (is.null(clus.sc)) return(df0)
      #}
      # Cluster single cell only. min.dim/max.dim is the same with coclustering sc+bulk, which is time-inefficient.
      # if (is.null(sc.dim.min)) clus.sc <- cluster_cell(data=cell, min.dim=df0$dim, max.dim=max.dim, graph.meth=graph.meth, dimred=dimred)
      # if (is.null(clus.sc)) return(df0)
      # Filter cells by sim and proportion p.
      cell.refined <- refine_cluster(clus.sc, sim=df0$sim, sim.p=df0$sim.p, sim.meth=df0$sim.meth, verbose=verbose)
      cell.refined <- true_bulk(cell.refined, df.match)
      # Cocluster.
      roc.lis <- coclus_roc(bulk=bulk, cell.refined=cell.refined, df.match=df.match, min.dim=df0$dim, max.dim=df0$max.dim, graph.meth=df0$graph.meth, dimred=df0$dimred, sim.meth=df0$sim.meth) 
      if (return.all==TRUE) { 
        lis0 <- c(list(index=df0$index, cell.refined=cell.refined), roc.lis)
        return(lis0)
      }
      roc.obj <- roc.lis$roc.obj; if (is.null(roc.obj)) return(df0)
      df.roc <- roc.lis$df.roc
      df0$auc <- round(auc(roc.obj), 3); df0$total <- nrow(df.roc)
      df0$true <- sum(df.roc$response)
      best <- round(coords(roc.obj, x='best'), 3)
      if (nrow(best)==1) { 
        df0$thr <- round(best$threshold, 3)
        df0$sens <- round(best$sensitivity, 3)
        df0$spec <- round(best$specificity, 3)
      } else if (unique(max(abs(best$threshold)))==Inf) { 
      df0$thr <- df0$sens <- df0$spec <- Inf }
     if (verbose==TRUE) print(df0); 
      # pa <- paste0(tmp.dir, '/df.para_', j, '.rds')
      # if (!is.null(file)) saveRDS(df0, file=pa)
      return(df0)	
    })
    if (return.all==TRUE) {
      # Indexes in input df.para are preserved. 
      for (i in seq_along(lis)) {
        names(lis)[i] <- lis[[i]]$index; lis[[i]]$index <- NULL
      }
      lis.all <- c(lis.all, lis)
      lis.all <- lis.all[order(as.numeric(names(lis.all)))]
      if (!is.null(file)) saveRDS(lis.all, file=paste0(file, '.rds'))
    } else {
      # Indexes in input df.para are preserved. 
      df.para.all <- rbind(df.para.all, do.call('rbind', lis))
      df.para.all <- df.para.all[order(df.para.all$index), , drop=FALSE]
      if (!is.null(file)) saveRDS(df.para.all, file=paste0(file, '.rds'));
    }
  }

  }
  if (return.all==TRUE) {
    # The order is the same 
    return(lis.all) 
  } else { 
    # df.para <- do.call('rbind', lis.res)
    df.para.all <- df.para.all[, setdiff(colnames(df.para.all), 'index'), drop=FALSE]
    if (!is.null(file)) saveRDS(df.para.all, file=paste0(file, '.rds')); return(df.para.all) 
  }
}





