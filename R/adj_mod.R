#' Computing adjacency matrix and identifying network modules
#'
#' This function is designed to identify network modules in a data matrix, e.g. returned by \code{\link{submatrix}}, which builds on the \pkg{WGCNA} algorithm (Langfelder and Horvath 2008; Ravasz et al. 2002). Each module consists of \emph{i.e.} groups of biomolecules showing most similar coexpression profiles.

#' @inheritParams matrix_hm
#' @param type The network type, one of `unsigned`, `signed`, `signed hybrid`, or `distance`. Correlations and distances are transformed as follows: for \code{type="unsigned"}, \code{adjacency=|cor|^power}; for \code{type="signed"}, \code{adjacency=(0.5 * (1+cor) )^power}; for \code{type="signed hybrid"}, if cor>0 \code{adjacency=cor^power}, otherwise adjacency=0; and for \code{type="distance"}, \code{adjacency=(1-(dist/max(dist))^2)^power}. Refer to \pkg{WGCNA} (Langfelder and Horvath 2008) for more details. 
#' @param power A numeric of soft thresholding power for generating the adjacency matrix. The default is 1 for \code{type=='distance'} and 6 for other network types.
#' @param arg.adj A list of additional arguments passed to \code{\link[WGCNA]{adjacency}}, \emph{e.g.} \code{list(corFnc='cor')}. The default is an empty list \code{list()}.
#' @inheritParams WGCNA::TOMsimilarity
#' @param arg.tom A list of additional arguments passed to \code{\link[WGCNA]{TOMsimilarity}}, \emph{e.g.} \code{list(verbose=1)}. The default is an empty list \code{list()}.
#' @inheritParams flashClust::flashClust
#' @param ds One of 0, 1, 2, or 3. The sensitivity of module identification. The smaller, the less modules identified with larger sizes. 
#' @param minSize The expected minimum module size. The default is 15. Refer to \pkg{WGCNA} for more details.
#' @param arg.cut A list of additional arguments passed to \code{\link[dynamicTreeCut]{cutreeHybrid}}, \emph{e.g.} \code{list(verbose=2)}. The default is an empty list \code{list()}.

#' @param dir The directory to save the results, where the adjacency matrix and module assignments are saved as TSV-format files "adj.txt" and "mod.txt" respectively. 

#' @return A list containing the adjacency matrix and module assignments.

#' @section Details:
#' To identify modules, first a similarity matrix including similarities between any pair of biomolecules is computed using distance or correlation-based similarity metrics. Second, the similarity matrix is transformed into an adjacency matrix. Third, the adjacency matrix is used to calculate a topological overlap matrix (TOM) where shared neighborhood information among biomolecules is used to preserve robust connections, while removing spurious connections. Fourth, the distance transformed TOM is used for hierarchical clustering with the "flashClust" package (Langfelder and Horvath 2012). Fifth, network modules are identified with the "dynamicTreeCut" package (Langfelder, Zhang, and Steve Horvath 2016). Its \code{ds} (deepSplit) argument can be assigned integer values from 0 to 3, where higher values increase the stringency of the module identification process. Since this is a coexpression analysis, total samples should be at least 5. Otherwise, identified modules are not reliable.
#' 
#' These procedures are wrapped in \code{adj_mod} for convenience. The result is a list containing the adjacency matrix and the final module assignments stored in a \code{data.frame}. Since the interactive network feature (see \code{network}) used in the downstream visualization performs best on smaller modules, only modules obtained with stringent ds settings (here \code{ds=2} and \code{ds=3}) are returned. 

#' @inherit submatrix examples

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559
#' Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. Journal of Statistical Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/ 
#' R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ 
#' Peter Langfelder, Bin Zhang and with contributions from Steve Horvath (2016). dynamicTreeCut: Methods for Detection of Clusters in Hierarchical Clustering Dendrograms. R package version 1.63-1. https://CRAN.R-project.org/package=dynamicTreeCut
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' Ravasz, E, A L Somera, D A Mongru, Z N Oltvai, and A L Barabási. 2002. “Hierarchical Organization of Modularity in Metabolic Networks.” Science 297 (5586): 1551–5. 

#' @export
#' @importFrom stats quantile dist as.dist
#' @importFrom SummarizedExperiment SummarizedExperiment


adj_mod <- function(data, assay.na=NULL, type='signed', power=if (type=='distance') 1 else 6, arg.adj=list(), TOMType='unsigned', arg.tom=list(), method='complete', ds=0:3, minSize=15, arg.cut=list(), dir=NULL) {
  pkg <- check_pkg('WGCNA'); if (is(pkg, 'character')) stop(pkg)
  pkg <- check_pkg('flashClust'); if (is(pkg, 'character')) stop(pkg)
  pkg <- check_pkg('dynamicTreeCut'); if (is(pkg, 'character')) stop(pkg)
  options(stringsAsFactors=FALSE)
  # Get data matrix.
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) {
    dat.lis <- check_data(data=data); data <- t(dat.lis$dat)
  } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
    if (is.null(assay.na)) {
      if (length(assays(data)) > 1) stop("Please specify which assay to use by assigning the assay name to 'assay.na'!") else if (length(assays(data)) == 1) assay.na <- 1
    }
    data <- t(assays(data)[[assay.na]])
    if (any(duplicated(rownames(data)))) stop('Please use function \'aggr_rep\' to aggregate replicates!')

  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", "dgCMatrix", "SummarizedExperiment", or "SingleCellExperiment", except that "spatial_hm" also accepts a "vector".') }
  if (nrow(data)<5) cat('Warning: variables of sample/condition are less than 5! \n')
  if (ncol(data)>10000) cat('More than 10,000 rows are detected in data. Computation may take a long time! \n')
  
  # Compute adjacency matrix.
  arg.adj <- c(list(datExpr=data, power=power, type=type), arg.adj)
  adj <- do.call(WGCNA::adjacency, arg.adj)
  # Compute TOM and hierarchical clustering.
  arg.tom <- c(list(adjMat=adj, TOMType=TOMType), arg.tom)
  tom <- do.call(WGCNA::TOMsimilarity, arg.tom)
  dissTOM <- 1-tom; tree.hclust <- flashClust::flashClust(d=as.dist(dissTOM), method=method)
  # Cut the tree to get modules.
  cutHeight <- quantile(tree.hclust[['height']], probs=seq(0, 1, 0.05))[19]
  arg.cut1 <- list(dendro=tree.hclust, pamStage=FALSE, cutHeight=cutHeight, distM=dissTOM)
  na.cut1 <- names(arg.cut1); na.cut <- names(arg.cut)
  if ('deepSplit' %in% na.cut) arg.cut <- arg.cut[!(na.cut %in% 'deepSplit')]
  w <- na.cut1 %in% na.cut
  arg.cut1 <- c(arg.cut1[!w], arg.cut)
  mcol <- NULL; for (d in ds) {       
    min <- as.numeric(minSize)-3*d; if (min < 5) min <- 5
    arg.cut <- c(list(minClusterSize=min, deepSplit=d), arg.cut1)
    tree <- do.call(dynamicTreeCut::cutreeHybrid, arg.cut)
    mcol <- cbind(mcol, tree$labels); arg.cut <- list()

  }; colnames(mcol) <- as.character(ds); rownames(mcol) <- colnames(adj) <- rownames(adj) <- colnames(data)
  if (!is.null(dir)) {  
    dir <- normalizePath(dir, winslash="/", mustWork=FALSE)
    if (!dir.exists(dir)) dir.create(dir)
    write.table(adj, file.path(dir, "adj.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    write.table(mcol, file.path(dir, "mod.txt"), sep="\t", row.names=TRUE, col.names=TRUE) 
  }; return(list(adj=adj, mod=mcol))

}






