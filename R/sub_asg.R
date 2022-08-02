#' Subset the bulk-cell assignments
#' 
#' Subset the bulk-cell assignments according to a threshold, which is a similarity value between bulk and cells. 

#' @param res.lis The result list of coclustering, which is the output of \code{tests} and comprises three slots \code{sce}, \code{roc.obj}, \code{df.roc}.
#' @param thr The threshold for subsetting bulk-cell assignments, which is a similarity value (Pearson's or Spearman's correlation coefficient) between bulk and cells. Only bulk-cell assignments with similarity values above the thresold would remain. The default is 0.
#' @param df.desired.bulk A "data.frame" of desired bulk for some cells. The cells could be specified by providing x-y axis ranges in an embedding plot ("UMAP", "PCA", "TSNE") returned by \code{plot_dim}. E.g. \code{df.desired.bulk <- data.frame(x.min=c(4, -6), x.max=c(5, -5), y.min=c(-2.5, 2), y.max=c(-2, 2.5), desiredSVGBulk=c('CORT', 'STELE'), dimred='UMAP')}, where columns \code{x.min}, \code{x.max}, \code{y.min}, \code{y.max}, \code{desiredSVGBulk}, \code{dimred} are required. In this example, cells located in 4 <= x <= 5 and -2.5 <= y <= -2 in the "UMAP" plot are assigned "STELE", and cells located in -6 <= x <= -5 and 2 <= y <= 2.5 in the "UMAP" plot are assigned "CORT". \cr Alternatively, the "data.frame" could be downloaded from the Shiny app launched by \code{desired_bulk_shiny}. \cr \code{df.desired.bulk} and \code{df.match} together are used to tailor the co-clustering results. That is to say additional true bulk-cell assignments are created and included in the final assignments. If these assignments conflict with the co-colustering results the latter would be overwritten.
#' @param df.match The ground-truth matching between cells and bulk. See the example of \code{data(df.match)}.
#' @param true.only Logical. If \code{TRUE}, only the true assignments are returned in the subset \code{SingleCellExperiment}. This argument affects the values for plottting SHMs.

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
#' # Cluster single cells.
#' clus.sc <- cluster_cell(data=mus.lis.fil$sc.mus, min.dim=10, max.dim=50, graph.meth='knn', dimred='PCA')
#' # Cluster labels are stored in "label" column in "colData".
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
#' roc.lis <- coclus_roc(bulk=mus.lis.fil$bulk, cell.refined=cell.refined, df.match=df.match, min.dim=12, max.dim=50, graph.meth='knn', dimred='PCA', sim.meth='spearman') 
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
#' res.lis <- cocluster(bulk=mus.lis.fil$bulk, cell=mus.lis.fil$sc.mus, df.match=df.match, df.para=NULL, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=TRUE)
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
#' res.multi <- cocluster(bulk=mus.lis.fil$bulk, cell=mus.lis.fil$sc.mus, df.match=df.match, df.para=df.par, sc.dim.min=10, max.dim=50, sim=0.2, sim.p=0.8, dim=12, graph.meth='knn', dimred='PCA', sim.meth='spearman', return.all=TRUE, multi.core.par=MulticoreParam(workers=2))
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
#' }



#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.

#' @export sub_asg
#' @importFrom SingleCellExperiment colLabels colLabels<- reducedDimNames 
#' @importFrom SummarizedExperiment colData 

sub_asg <- function(res.lis, thr=0, df.desired.bulk=NULL, df.match=NULL, true.only=TRUE) {
  # save(res.lis, thr, df.desired.bulk, df.match, true.only, file='sub.asg.arg')
  x <- y <- desiredSVGBulk <- SVGBulk <- index <- predictor <- total <- true <- NULL
  # Validate labels.
  sce <- res.lis$cell.refined; labs <- colLabels(sce)
  er.wa <- check(as.vector(labs), as.numeric)
  if (is.numeric(er.wa) & length(er.wa) > 0) colLabels(sce) <- as.character(labs)
  # If desired bulk is provided with x-y coordinate range, convert "df.desired.bulk" to the form downloaded from Shiny app.
  xy.ran.cna <- c('x.min', 'x.max', 'y.min', 'y.max') 
  if (any(xy.ran.cna %in% colnames(df.desired.bulk))) {
    if (!all(c(xy.ran.cna, 'desiredSVGBulk', 'dimred') %in% colnames(df.desired.bulk))) stop(paste0('Make sure the "df.desired.bulk" has these columns: ', paste0(xy.ran.cna, collapse=', '), '!'))
    x.min <- df.desired.bulk$x.min
    x.max <- df.desired.bulk$x.max
    y.min <- df.desired.bulk$y.min
    y.max <- df.desired.bulk$y.max
    dimred <- unique(df.desired.bulk$dimred)
    if (!(all(x.min <= x.max) & all(y.min <= y.max))) stop('Make sure x.min <= x.max and y.min <= y.max !')
    if (! dimred %in% reducedDimNames(sce)) stop(paste0(dimred, ' is not found!'))
    gg <- plot_dim(sce, dim=dimred, color.by=colnames(colData(sce))[1])
    dat <- gg$data; df.desire <- NULL
    blk <- df.desired.bulk$desiredSVGBulk
    # Assign each bulk to cells in corresponding x-y range.
    for (i in seq_len(nrow(df.desired.bulk))) { 
      df.desired.bulk0 <- df.desired.bulk[i, ]
      x.min0 <- df.desired.bulk0$x.min
      x.max0 <- df.desired.bulk0$x.max
      y.min0 <- df.desired.bulk0$y.min
      y.max0 <- df.desired.bulk0$y.max
      df0 <- subset(dat, x >= x.min0 & x <= x.max0 & y >= y.min0 & y <= y.max0)
      if (nrow(df0)==0) {
        cat('No cells selected for row: ', i, '!\n', sep=''); next
      }
      # Same form as downloaded from Shiny app.
      df.desire0 <- data.frame(x=df0$x, y=df0$y, key=df0$key, desiredSVGBulk=blk[i], dimred=dimred) 
      df.desire <- rbind(df.desire, df.desire0)   
    }; df.desired.bulk <- df.desire
  }

  df.roc <- res.lis$df.roc
  # Include desired bulk to "df.roc".
  if (!is.null(df.desired.bulk)) {
    xy.cna <- c('x', 'y', 'key', 'desiredSVGBulk', 'dimred') 
    if (!all(xy.cna %in% colnames(df.desired.bulk))) stop(paste0('Make sure "df.desired.bulk" has these columns: ', paste0(xy.cna, collapse=', '), '!'))
    key <- as.numeric(df.desired.bulk$key)
    # "none" indicates no change.
    desire.blk <- setdiff(unique(df.desired.bulk$desiredSVGBulk), 'none')
    if (length(desire.blk)>0) {
      for (i in desire.blk) {
      df.desired.bulk0 <- subset(df.desired.bulk, desiredSVGBulk==i)
      df.match0 <- subset(df.match, SVGBulk==i)[1, ]
      # Use "df.match" and "desire.blk" as gold standard to construct true coclustering results "df.roc0" for each "desire.blk", then append "df.roc0" to "df.roc".
      df.roc0 <- data.frame(matrix(ncol=ncol(df.roc), nrow=nrow(df.desired.bulk0), dimnames=list(NULL, colnames(df.roc))))
      df.roc0[, colnames(df.match0)] <- df.match0 
      df.roc0$assignedBulk <- df.roc0$trueBulk
      df.roc0$response <- TRUE; df.roc0$predictor <- 1
      df.roc0$index <- df.desired.bulk0$key
      df.roc <- rbind(subset(df.roc, !index %in% df.roc0$index), df.roc0)
      }; res.lis$df.roc <- df.roc
    }
  }
 
  # Subset bulk-cell assignment according to thr.
  # False assignments are aslo included and will be aggregated for plotting SHMs.
  df.roc <- subset(df.roc, predictor >= thr)
  # Cells without bulk assigments have original labels.
  cdat <- colData(sce)
  cdat$assignedBulk <- cdat$response <- 'none'
  cdat$index <- seq_len(ncol(sce))
  cdat$assignedBulk[df.roc$index] <- df.roc$assignedBulk
  cdat$response[df.roc$index] <- df.roc$response
  cdat$index[df.roc$index] <- df.roc$index
  # In consideration of rematching in coclustering.
  cdat$cell[df.roc$index] <- df.roc$cell
  cdat$SVGBulk[df.roc$index] <- df.roc$SVGBulk
  cdat$trueBulk[df.roc$index] <- df.roc$trueBulk

  cdat <- cdat[, c('label', 'cell', 'assignedBulk', 'response', 'SVGBulk', 'trueBulk', 'index')]
  colData(sce) <- cdat
  sce.sub <- sce[, df.roc$index]
  if (true.only==TRUE) sce.sub <- sce.sub[, colData(sce.sub)$response==TRUE]
  # Indicate the dim plot through which desired bulk is assigned, so that in SHMs the same dim plot is generated.
  if (!is.null(df.desired.bulk)) colData(sce.sub)$dimred <- unique(df.desired.bulk$dimred)
  return(list(cell.refined=sce, cell.sub=sce.sub, df.roc=df.roc))
}
