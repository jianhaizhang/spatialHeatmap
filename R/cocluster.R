#' Co-clustering bulk and single cell data 
#'
#' Automatically assigns bulk tissues to single cells through co-clustering.

#' @param bulk The normalized bulk data (log2-scale) in form of \code{SingleCellExperiment}.
#' @param cell The normalized single cell data in form of \code{SingleCellExperiment}. 
#' @param df.match A \code{data.frame} specifying ground-truth matching between cells and bulk, applicable in co-clustering optimization. 
#' @param sim.meth Method to calculate similarities between bulk and cells in each cocluster when assigning bulk to cells. \code{spearman} (default) or \code{pearson}.
#' @inheritParams reduce_dim
#' @inheritParams cluster_cell

#' @return A list of coclustering results in \code{SingleCellExperiment} and an \code{roc} object (relevant in optimization).

#' @examples

#' # To obtain reproducible results, a fixed seed is set for generating random numbers.
#' set.seed(10); library(SummarizedExperiment)
#' # Example bulk data of mouse brain for coclustering (Vacher et al 2021).
#' blk.mus.pa <- system.file("extdata/shinyApp/data", "bulk_mouse_cocluster.rds", 
#' package="spatialHeatmap") 
#' blk.mus <- readRDS(blk.mus.pa)
#' assay(blk.mus)[1:3, 1:5]
#'
#' # Example single cell data for coclustering (Ortiz et al 2020).
#' sc.mus.pa <- system.file("extdata/shinyApp/data", "cell_mouse_cocluster.rds", 
#' package="spatialHeatmap") 
#' sc.mus <- readRDS(sc.mus.pa)
#' colData(sc.mus)[1:3, , drop=FALSE]
#'
#' \donttest{
#' # Normalization: bulk and single cell are combined and normalized, then separated.
#' mus.lis.nor <- norm_cell(sce=sc.mus, bulk=blk.mus, com=FALSE)
#'
#' # Aggregate bulk replicates.
#' blk.mus.aggr <- aggr_rep(data=mus.lis.nor$bulk, assay.na='logcounts', sam.factor='sample', 
#' aggr='mean')
#' # Filter bulk
#' blk.mus.fil <- filter_data(data=blk.mus.aggr, pOA=c(0.1, 1), CV=c(0.1, 50), verbose=FALSE) 
#' # Filter cell and subset bulk to genes in cell
#' blk.sc.mus.fil <- filter_cell(sce=mus.lis.nor$cell, bulk=blk.mus.fil, cutoff=1, p.in.cell=0.1,
#' p.in.gen=0.01, verbose=FALSE) 
#' # Co-cluster bulk and single cells.
#' coclus.mus <- cocluster(bulk=blk.sc.mus.fil$bulk, cell=blk.sc.mus.fil$cell, min.dim=12, 
#' dimred='PCA', graph.meth='knn', cluster='wt')
#' # Co-clustering results. The 'cluster' indicates cluster labels, the 'bulkCell' indicates bulk
#' # tissues or single cells, the 'sample' suggests original labels of bulk and cells, the 
#' # 'assignedBulk' refers to bulk tissues assigned to cells with none suggesting un-assigned, 
#' # and the 'similarity' refers to Spearman's correlation coefficients for assignments between 
#' # bulk and cells, which is a measure of assignment strigency.
#' colData(coclus.mus)
#' 
#' # Filter bulk-cell assignments according a similarity cutoff (min.sim).
#' coclus.mus <- filter_asg(coclus.mus, min.sim=0.1)
#' 
#' # Tailor bulk-cell assignments in R.
#' plot_dim(coclus.mus, dim='UMAP', color.by='sample', x.break=seq(-10, 10, 1), 
#' y.break=seq(-10, 10, 1), panel.grid=TRUE)
#' # Define desired bulk tissues for selected cells.
#' df.desired.bulk <- data.frame(x.min=c(-8), x.max=c(-3.5), y.min=c(-2.5), y.max=c(0.5), 
#' desiredBulk=c('hippocampus'), dimred='UMAP')
#' df.desired.bulk
#' # Tailor bulk-cell assignments.
#' coclus.mus.tailor <- refine_asg(sce.all=coclus.mus, df.desired.bulk=df.desired.bulk)
#' 
#' # Define desired bulk tissues for selected cells on a Shiny app.
#' # Save "coclus.mus" using "saveRDS" then upload the saved ".rds" file to the Shiny app.
#' saveRDS(coclus.mus, file='coclus.mus.rds')
#' 
#' # Start the Shiny app.
#' desired_bulk_shiny()
#' }


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Vacher, Claire-Marie, Helene Lacaille, Jiaqi J. O’Reilly, Jacquelyn Salzbank, Dana Bakalar, Sonia Sebaoui, Philippe Liere, et al. 2021. “Placental Endocrine Function Shapes Cerebellar Development and Social Behavior.” Nature Neuroscience 24 (10): 1392–1401.
#' Ortiz, Cantin, Jose Fernandez Navarro, Aleksandra Jurek, Antje Märtin, Joakim Lundeberg, and Konstantinos Meletis. 2020. “Molecular Atlas of the Adult Mouse Brain.” Science Advances 6 (26): eabb3446.
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and        environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.

#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment logcounts

cocluster <- function(bulk, cell, df.match=NULL, min.dim=11, max.dim=50, dimred='PCA', graph.meth='knn', knn.gr=list(), snn.gr=list(), cluster='fg', wt.arg=list(steps = 4), fg.arg=list(), sl.arg=list(spins = 25), le.arg=list(), eb.arg=list(), sim.meth='spearman') {
  # save(bulk, cell, df.match, min.dim, max.dim, dimred, graph.meth, knn.gr, snn.gr, cluster, wt.arg, fg.arg, sl.arg, le.arg, eb.arg, sim.meth, file='cocluster.arg')
  # Two or more cells are needed.
  if (length(unique(colnames(cell)))<=1) {
    cat('At least two cells are needed!\n'); return()
  }
  if (is(bulk, 'SingleCellExperiment')) {
    bulk.log <- logcounts(bulk)
    if (all(round(bulk.log)==bulk.log)) { 
      message('The bulk data should be in log2 scale!')
      return()
    }
  }
  if(length(intersect(colnames(bulk), colnames(cell)))>0) stop('Common identifiers are detected between bulk and cell data!')
  # Combine bulk and single cell data.
  inter <- intersect(rownames(bulk), rownames(cell))
  blk.kp <- bulk[inter, ]; sc.kp <- cell[inter, ]
  pkg <- check_pkg('BiocGenerics'); if (is(pkg, 'character')) stop(pkg)
  com.kp <- BiocGenerics::cbind(blk.kp, sc.kp)
  com.kp$index <- seq_len(ncol(com.kp))
  com.kp$sample <- colnames(com.kp)
  # Cluster cells+bulk.
  message('Dimension reduction ...')
  sce.dimred <- reduce_dim(sce=com.kp, min.dim=min.dim, max.dim=max.dim)
  if (is.null(sce.dimred)) return()
  sce.tsne.com <- cluster_cell(sce=sce.dimred, graph.meth=graph.meth, knn.gr=knn.gr, snn.gr=snn.gr, dimred=dimred, cluster=cluster, wt.arg=wt.arg, fg.arg=fg.arg, sl.arg=sl.arg, le.arg=le.arg, eb.arg=eb.arg)
  if (is.null(sce.tsne.com)) return()
  # sce.tsne.com <- clus.com$sce.tsne; sce.com.nor <- clus.com$sce.nor
  # data.com: logcounts(sce.com.nor), reducedDim(sce.tsne.com, 'PCA')
  # roc <- com_roc(dat.com=logcounts(sce.com.nor), cdat.com, dat.blk=bulk)
  res <- com_roc(sce.coclus=sce.tsne.com, dimred=dimred, dat.blk=bulk, df.match=df.match, sim.meth=sim.meth) 
  # Include assignment info to "colData".  
  # cdat <- colData(cell)
  # if (!'dataBulk' %in% colnames(cdat)) cdat$dataBulk <- 'none'
  # if (!'SVGBulk' %in% colnames(cdat)) cdat$SVGBulk <- 'none'
  # colData(cell) <- cdat

  # Isolate bulk data.
  # sce.bulk <- sce.tsne.com[, !colnames(sce.tsne.com) %in% cdat$cell]
  # cdat.blk <- colData(sce.bulk)
  # cdat.blk$cluster <- NULL
  # names(cdat.blk)[names(cdat.blk)=='cell'] <- 'dataBulk'
  # svg.blk <- df.match$SVGBulk
  # data.blk <- df.match$dataBulk
  # names(svg.blk) <- data.blk
  # cdat.blk$SVGBulk <- svg.blk[cdat.blk$dataBulk]
  # colData(sce.bulk) <- cdat.blk

  # sce.lis <- refine_asg(res.lis=c(list(cell=cell, sce.bulk.cell=sce.tsne.com, sce.bulk=sce.bulk), res.lis), thr=-Inf, df.desired.bulk=NULL, df.match=df.match)
  return(res)
}

#' Calculate ROC/AUC for the combined bulk and single cell data
#'
#' @param sce.coclus The coclustered bulk and single cell data in a \code{SingleCellExperiment}, where cocluster assignments are stored in the \code{cluster} column in \code{colData}.
#' @param dimred The reduced dimensionality to use for calculating similarities between bulk and cells in each cocluster. \code{PCA} or \code{UMAP}.
#' @param dat.blk The bulk data at log2 scale. 
#' @param df.match The matching data frame between cells and true bulk.
#' @param sim.meth Similarity method between bulk and cells in each cocluster. \code{spearman} (default) or \code{pearson}.

#' @return A list of \code{roc} object and the data frame to create the \code{roc}.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.  0, https://bioconductor.org/packages/SummarizedExperiment. 
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Xavier Robin, Natacha Turck, Alexandre Hainard, Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics, 12, p. 77.  DOI: 10.1186/1471-2105-12-77 <http://www.biomedcentral.com/1471-2105/12/77/>

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim

com_roc <- function(sce.coclus, dimred, dat.blk, df.match=NULL, sim.meth='spearman') {
  pkg <- check_pkg('pROC'); if (is(pkg, 'character')) stop(pkg)
  # save(sce.coclus, dimred, dat.blk, df.match, sim.meth, file='com.roc.arg')
  if (!is.null(df.match)) {
  # if (!'SVGBulk' %in% colnames(df.match)) df.match$SVGBulk <- 'none' 
    df.match <- df.match[!duplicated(df.match), , drop=FALSE]
  }
  response <- NULL; dat.com <- t(reducedDim(sce.coclus, dimred))
  cdat.com <- colData(sce.coclus) 
  if (nrow(dat.com)<5) stop('Dimensions should be >= 5!')
  dat.com <- as(dat.com, 'matrix')
  # Co-cluster labels.
  lab.com <- cdat.com$cluster
  labs.com <- as.character(sort(lab.com))
  lab.com.uni <- unique(labs.com)
  # Co-clusters in a list.
  coclus <- lapply(lab.com.uni, function(x) sort(table(rownames(cdat.com)[lab.com==x])))
  names(coclus) <- lab.com.uni; bulk.na <- unique(colnames(dat.blk))
  # Index all bulk and cells for convenience to select cells with a bulk tissue assignment.
  cna.all <- colnames(dat.com)
  names(cna.all) <- seq_along(cna.all)

  # In each co-cluster, compute PCC between cells and bulk.
  df.asg <- NULL; for (i in lab.com.uni) {
    # Cell and bulk names in a co-cluster.
    co.na <- names(coclus[[i]])
    blk.na0 <- unique(bulk.na[bulk.na %in% co.na])
    if (length(blk.na0)==0) next
    cat('Bulk in cocluster ', i, ':\n', sep='')
    cat(blk.na0, '\n', sep=' ')

    # Subset a cocluster in the combined data frame.
    clus0 <- dat.com[, lab.com==i, drop=FALSE]
    cna <- colnames(clus0); cna.all0 <- cna.all[lab.com==i]
    # Aggregate bulk in the subsetted co-cluster.
    mat <- vapply(blk.na0, function(x) rowMeans(clus0[, cna==x, drop=FALSE]), numeric(nrow(clus0)))

    # Combine aggregated bulk back with cells. In practice, cell ids are not known.
    sc0 <- clus0[, !colnames(clus0) %in% blk.na0, drop=FALSE]
    cna.sc0 <- cna.all0[!colnames(clus0) %in% blk.na0] 
    sc0.na.uni <- unique(colnames(sc0))
    blk.ag.sc <- cbind(mat, sc0)
    # Correlation between aggregated bulk and cells.
    sim <- cor(blk.ag.sc, method=sim.meth)
    sim.sub <- sim[rownames(sim) %in% blk.na0, !colnames(sim) %in% blk.na0, drop=FALSE]
    sc0.na <- colnames(sim.sub)
    # In each cocluster, compare each cell with each bulk according to sims, TRUE/FALSE is decided by the bulk with the largest sim with the cell.
   for (j in seq_len(ncol(sim.sub))) {
     # All sims of one cell between each bulk. Select the largest.
     blk.sel <- sort(sim.sub[, j], decreasing=TRUE)[1]
     if (length(blk.na0)==1) names(blk.sel) <- blk.na0
     blk.asg <- names(blk.sel) 
     # Matching between the cells and bulk in the co-cluster.
     # match.lis0 <- match.lis[sc0.na[j]]
     if (is.null(df.match)) {
       #asg.idx <- unlist(lapply(df.match$dataBulk, function(x) {
       # blk.asg %in% strsplit(x, ',|;| ')[[1]] }))
       # if (sum(asg.idx)==0) { 
         # message(paste0('Assigned bulk for cell "', sc0.na[j], '" is not found in "df.match": ', blk.asg)); 
       #  next
       # }
       # index is based on bulk+cell.
       df0 <- data.frame(cell=sc0.na[j], assignedBulk=blk.asg,  similarity=blk.sel, index=names(cna.sc0[j]), row.names=NULL)  
     } else if (is(df.match, 'data.frame')) {    
       true.bulk.pas <- df.match$trueBulk[df.match$cell==sc0.na[j]]
       if (length(true.bulk.pas)==0) { # No matching bulk.
         # message(paste0('Assigned bulk for cell "', sc0.na[j], '" is not found in "df.match": ', blk.asg));
         next
       }
       true.bulk <- strsplit(true.bulk.pas, ',|;| ')[[1]]
       true.bulk <- true.bulk[true.bulk!='']
       # svg.bulk <- df.match$SVGBulk[df.match$cell==sc0.na[j]]
       # Check if the cell matches with the true bulk.
       # index is based on bulk+cell.
       df0 <- data.frame(cell=sc0.na[j], assignedBulk=blk.asg, trueBulk=true.bulk.pas, response=blk.asg %in% true.bulk, similarity=blk.sel, index=names(cna.sc0[j]), row.names=NULL)
       }; df.asg <- rbind(df.asg, df0)
   }
  }
  if (!is.null(df.asg)) { 
    df.asg$index <- as.numeric(df.asg$index)
    df.asg$similarity <- round(df.asg$similarity, 3)
    sce.coclus$assignedBulk <- sce.coclus$similarity <- 'none'
    idx <- df.asg$index
    sce.coclus$assignedBulk[idx] <- df.asg$assignedBulk 
    sce.coclus$similarity[idx] <- df.asg$similarity 
  }
  roc.obj <- NULL; if (is.null(df.asg)) {
    message('No bulk tissues assigned!')
    return(list(sce=sce.coclus, roc.obj=roc.obj))
  }
  if ('cell' %in% colnames(df.match)) {
    # If a cell is not assigned bulk, its true bulk is also 'none'.
    sce.coclus$trueBulk <- sce.coclus$response <- 'none'
    idx <- df.asg$index
    sce.coclus$trueBulk[idx] <- df.asg$trueBulk 
    sce.coclus$response[idx] <- df.asg$response

    df.tr <- subset(df.asg, response==TRUE)
    cat('TRUE and FALSE:\n'); print(dim(df.asg))
    cat('TRUE:\n'); print(dim(df.tr))
    if (length(unique(df.asg$response))!=2) return()
    # index is based on cells only.
    # index.sc.asg <- seq_along(dat.com)[as.numeric(df.asg$index)]
    # index.all <- c(rep(0, ncol(dat.blk)), seq_along(colnames(dat.com)[!colnames(dat.com) %in% bulk.na]))
    # df.asg$index <- index.all[index.sc.asg]
    roc.obj <- pROC::roc(df.asg$response, df.asg$similarity, smoothed = TRUE, ci=TRUE, ci.alpha=0.9, stratified=FALSE, plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=TRUE, show.thres=TRUE, direction='<', levels=c('FALSE', 'TRUE'))
  }
  cdat <- colData(sce.coclus); cdat.na <- colnames(cdat)
  sel.na <- c('cluster', 'bulkCell', 'sample', 'assignedBulk', 'similarity', 'index')
  if ('cell' %in% colnames(df.match)) sel.na <- c('cluster', 'bulkCell', 'sample', 'assignedBulk', 'trueBulk', 'response', 'similarity', 'index')
  cdat <- cdat[, c(sel.na, setdiff(cdat.na, sel.na))]
  colData(sce.coclus) <- cdat
  if ('cell' %in% colnames(df.match)) res <- list(sce.all=sce.coclus, roc.obj=roc.obj) else res <- sce.coclus
  return(res)
  # df.asg$index <- as.numeric(df.asg$index) - ncol(dat.blk)
  # return(list(df.asg=df.asg, roc.obj=roc.obj))
}
