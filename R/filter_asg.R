#' Refine the bulk-cell assignments
#' 
#' Refine the bulk-cell assignments by including custom desired bulk tissues or subsetting the assignments according to a threshold, which is a similarity value between bulk and cells. 

#' @param res The coclustering results returned by \code{cocluster}.
#' @param min.sim The similarity cutoff for filterig bulk-cell assignments, which is a Pearson's or Spearman's correlation coefficient between bulk and cells. Only bulk-cell assignments with similarity values above the thresold would remain. The default is 0.

#' @return A \code{SingleCellExperiment} of remaining bulk-cell assignments. 

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
#' # Dimensionality reduction.
#' sce.dimred.sc <- reduce_dim(sce=mus.lis.fil$sc.mus, prop=0.1, min.dim=10, max.dim=50, de.pca=list(assay.type ="logcounts"))
#' # Cluster single cells.
#' clus.sc <- cluster_cell(sce=sce.dimred.sc, graph.meth='knn', dimred='PCA')
#' # Cluster labels are stored in the  "cluster" column in "colData".
#'colData(clus.sc)[1:3, ]
#'
#' # Refine cell clusters.
#' cell.refined <- refine_cluster(clus.sc, sim=0.2, sim.p=0.8, sim.meth='spearman')
#'
#' # Include matching information in "colData".
#' # cell.refined <- true_bulk(cell.refined, df.match)
#' # colData(cell.refined)[1:3, ]
#'
#' # Cocluster bulk and single cells.
#' roc.lis <- cocluster(bulk=mus.lis.fil$bulk, cell.refined=cell.refined, df.match=df.match, min.dim=12, max.dim=50, graph.meth='knn', dimred='PCA', sim.meth='spearman') 
#'
#' # The colustering results. "predictor" is the similarity between bulk and cells within a co-cluster. "index" is the cell index in the "SingleCellExperiment" after cell clusters are refined. 
#' roc.lis$df.asg[1:3, ]
#' # ROC curve created according to "roc.lis$df.asg". 
#' plot(roc.lis$roc.obj, print.auc=TRUE)
#' # Incorporate "cell.refined" in "roc.lis" for downstream use in co-visualization.
#' res.lis <- c(list(cell.refined=cell.refined), roc.lis)
#' 
#' # The processes of clustering single cells, refining cell clusters, and coclustering bulk and single cells can be performed by a single function.
#' library(BiocParallel)
#' res.lis <- coclus_meta(bulk=mus.lis.fil$bulk, cell=mus.lis.fil$sc.mus, df.match=df.match, df.para=NULL, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=TRUE)
#' res.lis <- res.lis[[1]]
#'
# Since "return.all" is "TRUE", the refined cell clusters, ROC object, and coclustering results are returned.  
#' names(res.lis)
#'
#' # "coclus_meta" accepts multiple combinations of parameter settings provided in a data frame, and coclustering on these combinations can be performed in parallel through "multi.core.par".
#' # Multiple combinations of parameter settings. If some parameters are not specified in this table such as "graph.meth", their default settings will be used. 
#' df.par <- data.frame(sim=c(0.2, 0.3), sim.p=c(0.8, 0.7), dim=c(12, 13))
#'
#' # The computation is parallelized on 2 cpu cores by "multi.core.par".
#' res.multi <- coclus_meta(bulk=mus.lis.fil$bulk, cell=mus.lis.fil$sc.mus, df.match=df.match, df.para=df.par, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=TRUE, multi.core.par=MulticoreParam(workers=2))
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
#' sce.lis <- refine_asg(res.lis=res.lis, df.desired.bulk=df.desired.bulk, df.match=df.match, true.only=TRUE)
#' }



#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.

#' @export
#' @importFrom SingleCellExperiment reducedDimNames 
#' @importFrom SummarizedExperiment colData colData<- 

filter_asg <- function(res, min.sim=0) {  
  if (is(res, 'SingleCellExperiment')) sce.all <- res else if (is(res, 'list')) sce.all <- res$sce.all
  # Set 'none' to assignments filtered out.
  cdat <- colData(sce.all)
  cdat0 <-  subset(cdat, bulkCell=='cell' & similarity!='none')
  cdat0$similarity <- as.numeric(cdat0$similarity)
  idx <- cdat0$similarity < min.sim
  cdat0$assignedBulk[idx] <- cdat0$similarity[idx] <- 'none'
  # Optimization.
  if (is(res, 'list')) cdat0$response[idx] <- cdat0$trueBulk[idx] <- 'none'
  cdat[cdat0$index, ] <- cdat0; colData(sce.all) <- cdat

  if (is(res, 'SingleCellExperiment')) res <- sce.all else if (is(res, 'list')) {
    cdat <- colData(sce.all)
    cdat0 <- subset(cdat, similarity!='none')
    roc.obj <- roc(cdat0$response, as.numeric(cdat0$similarity), smoothed = TRUE, ci=TRUE, ci.alpha=0.9, stratified=FALSE, plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=TRUE, show.thres=TRUE, direction='<', levels=c('FALSE', 'TRUE'))
    res$sce.all <- sce.all; res$roc.obj <- roc.obj
  }
  return(res)
}
