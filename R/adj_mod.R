#' Compute Adjacency Matrix and Identify Modules
#'
#' The objective is to explore target items (gene, protein, metabolite, \emph{etc}) in context of their neighbors sharing highly similar abundance profiles in a more advanced approach than \code{\link{matrix_hm}}. This advanced approach is the WGCNA algorithm (Langfelder and Horvath 2008; Ravasz et al. 2002). It takes the assay matrix subsetted by \code{\link{submatrix}} as input and splits the items into network modules, \emph{i.e.} groups of items showing most similar coexpression profiles.

#' @inheritParams matrix_hm
#' @param type The network type, one of "unsigned", "signed", "signed hybrid", "distance". Correlation and distance are transformed as follows: for \code{type="unsigned"}, \code{adjacency=|cor|^power}; for \code{type="signed"}, \code{adjacency=(0.5 * (1+cor) )^power}; for \code{type="signed hybrid"}, if cor>0 \code{adjacency=cor^power}, otherwise adjacency=0; and for \code{type="distance"}, \code{adjacency=(1-(dist/max(dist))^2)^power}. Refer to "WGCNA" (Langfelder and Horvath 2008) for more details. 
#' @param power A numeric of soft thresholding power for generating the adjacency matrix. The default is 1 for \code{type=='distance'} and 6 for other network types.
#' @param arg.adj A list of additional arguments passed to \code{\link[WGCNA]{adjacency}}, \emph{e.g.} \code{list(corFnc='cor')}. The default is an empty list \code{list()}.
#' @inheritParams WGCNA::TOMsimilarity
#' @param arg.tom A list of additional arguments passed to \code{\link[WGCNA]{TOMsimilarity}}, \emph{e.g.} \code{list(verbose=1)}. The default is an empty list \code{list()}.
#' @inheritParams flashClust::flashClust
#' @param minSize The expected minimum module size. The default is 15. Refer to "WGCNA" for more details.
#' @param arg.cut A list of additional arguments passed to \code{\link[dynamicTreeCut]{cutreeHybrid}}, \emph{e.g.} \code{list(verbose=2)}. The default is an empty list \code{list()}.

#' @param dir The directory to save the results. In this directory, a folder "customComputedData" is created automatically, where the adjacency matrix and module assignments are saved as TSV-format files "adj.txt" and "mod.txt" respectively. This argument should be the same with the \code{dir} in \code{\link{submatrix}} so that the "sub_matrix.txt" generated in \code{\link{submatrix}} is saved in the same folder. This argument is designed since the computation is intensive for large data matrix (\emph{e.g.} > 10,000 genes). Therefore, to avoid system crash when using the Shiny app (see \code{\link{shiny_all}}), "adj.txt" and "mod.txt" can be computed in advance and then uploaded to the app. In addition, the saved files can be used repetitively and therefore avoid repetitive computation. The default is NULL and no file is saved. This argument is used only when the "customComputedData" is chosen in the Shiny app. \cr The large matrix issue could be resolved by increasing the subsetting strigency to get smaller matrix in \code{\link{submatrix}} in most cases. Only in rare cases users cannot avoid very large subsetted matrix, this argument is recommended. 

#' @return A list containing the adjacency matrix and module assignment, which should be provided to \code{\link{network}}. The module assignment is a data frame. The first column is \code{ds=2} while the second is \code{ds=3} (see the "Details" section). The numbers in each column are module labels, where "0" means items not assigned to any modules. If \code{dir} is specified, both adjacency matrix and module assignment are automatically saved in the folder "customComputedData" as "adj.txt" and "mod.txt" respectively, which can be uploaded under "customComputedData" in the Shiny app (see \code{\link{shiny_all}}).

#' @section Details:
#' To identify modules, first a correlation matrix is computed using distance or correlation-based similarity metrics. Second, the obtained matrix is transformed into an adjacency matrix defining the connections among items. Third, the adjacency matrix is used to calculate a topological overlap matrix (TOM) where shared neighborhood information among items is used to preserve robust connections, while removing spurious connections. Fourth, the distance transformed TOM is used for hierarchical clustering. To maximize time performance, the hierarchical clustering is performed with the flashClust package (Langfelder and Horvath 2012). Fifth, network modules are identified with the dynamicTreeCut package (Langfelder, Zhang, and Steve Horvath 2016). Its \code{ds} (deepSplit) argument can be assigned integer values from 0 to 3, where higher values increase the stringency of the module identification process. Since this is a coexpression analysis, variables of sample/condition should be at least 5. Otherwise, identified modules are not reliable.
#' These procedures are wrapped in \code{adj_mod} for convenience. The result is a list containing the adjacency matrix and the final module assignments stored in a \code{data.frame}. Since the interactive network feature (see \code{network}) used in the downstream visualization performs best on smaller modules, only modules obtained with stringent ds settings (here \code{ds=2} and \code{ds=3}) are returned. 

#' @examples

#' ## In the following examples, the 2 toy data come from an RNA-seq analysis on development of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are 
#' ## included in this package. The complete raw count data are downloaded using the R package 
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as a 
#' ## "data frame" input to exemplify data of simple samples/conditions, while toy data2 as 
#' ## "SummarizedExperiment" to illustrate data involving complex samples/conditions.   

#' ## Set up toy data.
#' 
#' # Access toy data1.
#' cnt.chk.simple <- system.file('extdata/shinyApp/example/count_chicken_simple.txt', 
#' package='spatialHeatmap')
#' df.chk <- read.table(cnt.chk.simple, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
#' # Columns follow the namig scheme "sample__condition", where "sample" and "condition" stands 
#' # for organs and time points respectively.
#' df.chk[1:3, ]
#'
#' # A column of gene annotation can be appended to the data frame, but is not required.  
#' ann <- paste0('ann', seq_len(nrow(df.chk))); ann[1:3]
#' df.chk <- cbind(df.chk, ann=ann)
#' df.chk[1:3, ]
#'
#' # Access toy data2. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing samples and conditions is required for toy data2. It should be 
#' # made based on the experiment design, which is accessible through the accession number 
#' # "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in this
#' # package and accessed below. 

#' # Access the example targets file. 
#' tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
#' target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
#' # Every column in toy data2 corresponds with a row in targets file. 
#' target.chk[1:5, ]
#' # Store toy data2 in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene annotation, but not required.
#' rowData(se.chk) <- DataFrame(ann=ann)
#'
#' ## As conventions, raw sequencing count data should be normalized, aggregated, and filtered to
#' ## reduce noise.
#'
#' # Normalize count data.
#' # The normalizing function "calcNormFactors" (McCarthy et al. 2012) with default settings 
#' # is used.  
#' df.nor.chk <- norm_data(data=df.chk, norm.fun='CNF', data.trans='log2')
#' se.nor.chk <- norm_data(data=se.chk, norm.fun='CNF', data.trans='log2')

#' # Aggregate count data.
#' # Aggregate "sample__condition" replicates in toy data1.
#' df.aggr.chk <- aggr_rep(data=df.nor.chk, aggr='mean')
#' df.aggr.chk[1:3, ]

#' # Aggregate "sample_condition" replicates in toy data2, where "sample" is "organism_part" and
#' # "condition" is "age". 
#' se.aggr.chk <- aggr_rep(data=se.nor.chk, sam.factor='organism_part', con.factor='age', 
#' aggr='mean')
#' assay(se.aggr.chk)[1:3, 1:3]

#' # Filter out genes with low counts and low variance. Genes with counts over 5 (log2 unit) in 
#' # at least 1% samples (pOA), and coefficient of variance (CV) between 0.2 and 100 are retained.
#' # Filter toy data1.
#' df.fil.chk <- filter_data(data=df.aggr.chk, pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
#' # Filter toy data2.
#' se.fil.chk <- filter_data(data=se.aggr.chk, sam.factor='organism_part', con.factor='age', 
#' pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
#'
#' ## Select nearest neighbors for target genes 'ENSGALG00000019846' and 'ENSGALG00000000112', 
#' ## which are usually genes visualized in spatial heatmaps.
#' # Toy data1.
#' df.sub.mat <- submatrix(data=df.fil.chk, ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), p=0.1)
#' # Toy data2.
#' se.sub.mat <- submatrix(data=se.fil.chk, ann='ann', ID=c('ENSGALG00000019846', 
#' 'ENSGALG00000000112'), p=0.1) 
#'
#' # In the following, "df.sub.mat" and "se.sub.mat" is used in the same way, so only 
#' # "se.sub.mat" illustrated.
#'
#' # The subsetted matrix is partially shown below.
#' se.sub.mat[c('ENSGALG00000019846', 'ENSGALG00000000112'), c(1:2, 63)]

#' ## Adjacency matrix and module identification

#' # The modules are identified by "adj_mod". It returns a list containing an adjacency matrix and
#' # a data frame of module assignment. 
#' adj.mod <- adj_mod(data=se.sub.mat)

#' # The adjacency matrix is a measure of co-expression similarity between genes, where larger 
#' # value denotes higher similarity.
#' adj.mod[['adj']][1:3, 1:3]

#' # The modules are identified at two alternative sensitivity levels (ds=2 or 3). From 2 to 3, 
#' # more modules are identified but module sizes are smaller. The two sets of module assignment
#' # are returned in a data frame. The first column is ds=2 while the second is ds=3. The numbers
#' # in each column are module labels, where "0" means genes not assigned to any module.
#' adj.mod[['mod']][1:3, ]

#' @author Jianhai Zhang \email{zhang.jianhai@@hotmail.com; jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559 \cr Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. Journal of Statistical Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/ \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ \cr Peter Langfelder, Bin Zhang and with contributions from Steve Horvath (2016). dynamicTreeCut: Methods for Detection of Clusters in Hierarchical Clustering Dendrograms. R package version 1.63-1. https://CRAN.R-project.org/package=dynamicTreeCut \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' \cr Ravasz, E, A L Somera, D A Mongru, Z N Oltvai, and A L Barabási. 2002. “Hierarchical Organization of Modularity in Metabolic Networks.” Science 297 (5586): 1551–5. 

#' @export adj_mod
#' @importFrom WGCNA adjacency TOMsimilarity 
#' @importFrom flashClust flashClust
#' @importFrom stats quantile dist as.dist
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom dynamicTreeCut cutreeHybrid


adj_mod <- function(data, type='signed', power=if (type=='distance') 1 else 6, arg.adj=list(), TOMType='unsigned', arg.tom=list(), method='complete', minSize=15, arg.cut=list(), dir=NULL) {

  options(stringsAsFactors=FALSE)
  # Get data matrix.
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) {
    dat.lis <- check_data(data=data); data <- t(dat.lis$dat)
  } else if (is(data, 'SummarizedExperiment')) { data <- t(assay(data)) 

    if (any(duplicated(rownames(data)))) stop('Please use function \'aggr_rep\' to aggregate replicates!')

  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", or "SummarizedExperiment", except that "spatial_hm" also accepts a "vector".') }
  if (nrow(data)<5) cat('Warning: variables of sample/condition are less than 5! \n')
  if (ncol(data)>10000) cat('More than 10,000 rows are detected in data. Computation may take a long time! \n')
  
  # Compute adjacency matrix.
  arg.adj <- c(list(datExpr=data, power=power, type=type), arg.adj)
  adj <- do.call(adjacency, arg.adj)
  # Compute TOM and hierarchical clustering.
  arg.tom <- c(list(adjMat=adj, TOMType=TOMType), arg.tom)
  tom <- do.call(TOMsimilarity, arg.tom)
  dissTOM <- 1-tom; tree.hclust <- flashClust(d=as.dist(dissTOM), method=method)
  # Cut the tree to get modules.
  cutHeight <- quantile(tree.hclust[['height']], probs=seq(0, 1, 0.05))[19]
  arg.cut1 <- list(dendro=tree.hclust, pamStage=FALSE, cutHeight=cutHeight, distM=dissTOM)
  na.cut1 <- names(arg.cut1); na.cut <- names(arg.cut)
  if ('deepSplit' %in% na.cut) arg.cut <- arg.cut[!(na.cut %in% 'deepSplit')]
  w <- na.cut1 %in% na.cut
  arg.cut1 <- c(arg.cut1[!w], arg.cut)
  mcol <- NULL; for (ds in 2:3) {
         
    min <- as.numeric(minSize)-3*ds; if (min < 5) min <- 5
    arg.cut <- c(list(minClusterSize=min, deepSplit=ds), arg.cut1)
    tree <- do.call(cutreeHybrid, arg.cut)
    mcol <- cbind(mcol, tree$labels); arg.cut <- list()

  }; colnames(mcol) <- as.character(2:3); rownames(mcol) <- colnames(adj) <- rownames(adj) <- colnames(data)
  if (!is.null(dir)) { 
  
    dir <- normalizePath(dir, winslash="/", mustWork=FALSE)
    if (!dir.exists(dir)) stop(paste0(dir, ' does not exist!'))
    path <- paste0(dir, "/customComputedData/"); if (!dir.exists(path)) dir.create(path)
    write.table(adj, paste0(path, "adj.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    write.table(mcol, paste0(path, "mod.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
    
  }; return(list(adj=adj, mod=mcol))

}






