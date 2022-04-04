#' Co-cluster bulk and single cell data Calculate ROC/AUC for the combined bulk and single cell data
#'
#' Co-cluster bulk and refined single cell data and assign bulk to single cells. Since the identities of bulk tissues and single cells are labeled, ROC/AUC are calculated to evaluate the co-clustering performance.

#' @param bulk The normalized and filtered bulk data in form of \code{SingleCellExperiment}, \code{matrix} (log2-scale), or \code{data.frame} (log2-scale).
#' @param cell.refined The refined cell data in form of \code{SingleCellExperiment}, which is returned by \code{refine_cluster}. 
#' @param df.match The \code{data.frame} specifying matching between cells and true bulk.
#' @inheritParams cluster_cell
#' @param sim.meth Method to calculate similarities between bulk and cells in each cocluster when assigning bulk to cells. \code{spearman} (default) or \code{pearson}.
#' @inheritParams cluster_cell

#' @return A list of \code{roc} object and the data frame to create the \code{roc}.

#' @examples

#' See function "cocluster" by running "?cocluster".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods,    17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.

#' @export coclus_roc 
#' @importFrom SingleCellExperiment logcounts

coclus_roc <- function(bulk, cell.refined, df.match, min.dim=10, max.dim=50, graph.meth='snn', dimred='PCA', sim.meth='spearman', seed=1000) {
  # Combine bulk and single cell data.
  dat.kp.sc <- logcounts(cell.refined)
  cat('Single cell data: \n'); print(dim(dat.kp.sc))
  # Two or more cells are needed.
  if (length(unique(colnames(dat.kp.sc)))<=1) {
    cat('At least two cells are needed!\n'); return()
  }
  if (is(bulk, 'SingleCellExperiment')) bulk <- logcounts(bulk)
  if (all(round(bulk)==bulk)) stop('The bulk data should be in log2 scale!')
  inter <- intersect(rownames(bulk), rownames(dat.kp.sc))
  blk.kp <- as.matrix(bulk[inter, ]); sc.kp <- dat.kp.sc[inter, ]
  
  com.kp <- as.matrix(cbind(blk.kp, sc.kp)); dim(com.kp)
  table(colnames(blk.kp)); table(colnames(sc.kp))

  # Cluster cells+bulk.
  sce.tsne.com <- cluster_cell(data=com.kp, min.dim=min.dim, max.dim=max.dim, graph.meth=graph.meth, dimred=dimred, seed=seed)
  if (is.null(sce.tsne.com)) return()
  # sce.tsne.com <- clus.com$sce.tsne; sce.com.nor <- clus.com$sce.nor
  # data.com: logcounts(sce.com.nor), reducedDim(sce.tsne.com, 'PCA')
  # roc <- com_roc(dat.com=logcounts(sce.com.nor), cdat.com, dat.blk=bulk)
  roc.lis <- com_roc(sce.coclus=sce.tsne.com, dimred=dimred, dat.blk=bulk, df.match=df.match, sim.meth=sim.meth)
  return(roc.lis)
}



#' Calculate ROC/AUC for the combined bulk and single cell data
#'
#' @param sce.coclus The coclustered bulk and single cell data in a \code{SingleCellExperiment}, where cocluster assignments are stored in the \code{label} column in \code{colData}.
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
#' @importFrom pROC roc

com_roc <- function(sce.coclus, dimred, dat.blk, df.match, sim.meth='spearman') { 
  # save(sce.coclus, dimred, dat.blk, df.match, sim.meth, file='com.roc.arg')
  dat.com <- t(reducedDim(sce.coclus, dimred))
  cdat.com <- colData(sce.coclus) 
  if (nrow(dat.com)<5) stop('Dimensions should be >= 5!')
  dat.com <- as(dat.com, 'matrix')
  # Co-cluster labels.
  lab.com <- cdat.com$label
  labs.com <- as.character(sort(lab.com))
  lab.com.uni <- unique(labs.com)
  # Co-clusters in a list.
  coclus <- lapply(lab.com.uni, function(x) sort(table(rownames(cdat.com)[lab.com==x])))
  names(coclus) <- lab.com.uni; bulk.na <- unique(colnames(dat.blk))
  # Index all bulk and cells for convenience to selected cells with a bulk tissue assignment.
  cna.all <- colnames(dat.com)
  names(cna.all) <- seq_along(cna.all)

  # In each co-cluster, compute PCC between cells and bulk.
  df.roc <- NULL; for (i in lab.com.uni) {
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

     # Matching between the cells and bulk in the co-cluster.
     # match.lis0 <- match.lis[sc0.na[j]]
     true.bulk.pas <- df.match$trueBulk[df.match$cell==sc0.na[j]]
     if (length(true.bulk.pas)==0) next # No matching bulk.
     true.bulk <- strsplit(true.bulk.pas, ',|;| ')[[1]]
     true.bulk <- true.bulk[true.bulk!='']
     svg.bulk <- df.match$SVGBulk[df.match$cell==sc0.na[j]]
     # Check if the cell matches with the true bulk.
     # index is based on bulk+cell.
     df0 <- data.frame(assignedBulk=names(blk.sel), cell=sc0.na[j], response=names(blk.sel) %in% true.bulk, predictor=blk.sel, index=names(cna.sc0[j]), trueBulk=true.bulk.pas, SVGBulk=svg.bulk)
     df.roc <- rbind(df.roc, df0)
   }

  }
  if (is.null(df.roc)) return()
  df.tr <- subset(df.roc, response==TRUE)
  cat('TRUE and FALSE:\n'); print(dim(df.roc))
  cat('TRUE:\n'); print(dim(df.tr))
  if (length(unique(df.roc$response))!=2) return()
  # index is based on cells only.
  # index.sc.asg <- seq_along(dat.com)[as.numeric(df.roc$index)]
  # index.all <- c(rep(0, ncol(dat.blk)), seq_along(colnames(dat.com)[!colnames(dat.com) %in% bulk.na]))
  # df.roc$index <- index.all[index.sc.asg]
  df.roc$index <- as.numeric(df.roc$index) - ncol(dat.blk)

  roc.obj <- roc(df.roc$response, df.roc$predictor, smoothed = TRUE, ci=TRUE, ci.alpha=0.9, stratified=FALSE, plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=TRUE, show.thres=TRUE, direction='<', levels=c('FALSE', 'TRUE'))
  return(list(roc.obj=roc.obj, df.roc=df.roc))
}
