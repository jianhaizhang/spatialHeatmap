#' Subsetting data matrix 
#'
#' Given one or multiple biomolecules (gene, protein, metabolite, \emph{etc}) from a data matrix, this function subsets other biomolecules that have similar expression profiles with each of the given biomolecule independently. The subset data matrix of each biomolecule is combined in a row-wise manner and returned.   

#' @param data A `data.frame`, `SummarizedExperiment`, or `SingleCellExperiment` object, where the columns and rows of the data matrix are samples and biomolecules respectively. Since this function builds on co-expression analysis, samples should be at least 5, otherwise, the results are not reliable.
#' @param assay.na Applicable when \code{data} is `SummarizedExperiment` or `SingleCellExperiment`, where multiple assays could be stored. The name of assay to use.
#' @param ID A vector of biomolecules of interest.
#' @param p The proportion of top biomolecules with most similar expression profiles with a given biomolecule. Only biomolecules within this proportion are returned. It applies to each biomolecule independently and selected biomolecules of each given biomolecule are returned together.
#' @param n An integer of top biomolecules with most similar expression profiles with a given biomolecule. Only biomolecules within this number are returned. It applies to each biomolecule independently and selected biomolecules of each given biomolecule are returned together.
#' @param v A cutoff of correlation coefficient (CC, -1 to 1) or distance (>=0) for subsetting biomolecules sharing the most similar expression profiles with a given biomolecule. If \code{fun='cor'}, only biomolecules with CC larger than \code{v} are returned. If \code{fun='dist'}, only biomolecules with distance less than \code{v} are returned. It applies to each given biomolecule independently and selected biomolecules of each given biomolecule are returned together.
#' @param fun The function to calculate similarity/distance measures, `cor` (default) or `dist`, corresponding to \code{\link[stats]{cor}} or \code{\link[stats]{dist}} from the "stats" package respectively.
#' @param cor.absolute Logical. If `TRUE`, absolute correlation coefficients (CCs) are used. Only applies to \code{fun='cor'}. Default is `FALSE`, meaning the CCs preserve the negative sign when subsetting biomolecules. 
#' @param arg.cor A list of arguments passed to \code{\link[stats]{cor}} in the "stats" package. Default is \code{list(method="pearson")}.
#' @param arg.dist A list of arguments passed to \code{\link[stats]{dist}} in the "stats" package. Default is \code{list(method="euclidean")}.
#' @param file The file name to save subset data matrix. 

#' @return The subset data matrix in form of `data.frame`, `SummarizedExperiment`, or `SingleCellExperiment`.

#' @examples
#'
#' ## The example data included in this package come from an RNA-seq analysis on 
#' ## development of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019). 
#' ## The complete raw count data are downloaded using the R package ExpressionAtlas
#' ## (Keays 2019) with the accession number "E-MTAB-6769". 
#'
#' # Access example count data. 
#' count.chk <- read.table(system.file('extdata/shinyApp/data/count_chicken.txt', 
#' package='spatialHeatmap'), header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing spatial features and variables is made based on the 
#' # experiment design.
#' target.chk <- read.table(system.file('extdata/shinyApp/data/target_chicken.txt', 
#' package='spatialHeatmap'), header=TRUE, row.names=1, sep='\t')
#' # Every column in example data 2 corresponds with a row in the targets file. 
#' target.chk[1:5, ]
#' # Store example data in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#'
#' # Normalize data.
#' se.chk.nor <- norm_data(data=se.chk, norm.fun='CNF', log2.trans=TRUE)
#'
#' # Aggregate replicates of "spatialFeature_variable", where spatial features are organs
#' # and variables are ages.
#' se.chk.aggr <- aggr_rep(data=se.chk.nor, sam.factor='organism_part', con.factor='age',
#' aggr='mean')
#' assay(se.chk.aggr)[1:3, 1:3]
#'
#' # Genes with experssion values >= 5 in at least 1% of all samples (pOA), and coefficient
#' # of variance (CV) between 0.2 and 100 are retained.
#' se.chk.fil <- filter_data(data=se.chk.aggr, sam.factor='organism_part', con.factor='age', 
#' pOA=c(0.01, 5), CV=c(0.2, 100), file=NULL)
#' 
#' ## Subset the data matrix for gene 'ENSGALG00000019846' and 'ENSGALG00000000112'.
#' se.sub.mat <- submatrix(data=se.chk.fil, ID=c('ENSGALG00000019846', 
#' 'ENSGALG00000000112'), p=0.1) 
#' 
#' ## Hierarchical clustering. 
#' library(dendextend)
#' # Static matrix heatmap.
#' mhm.res <- matrix_hm(ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), data=se.sub.mat, 
#' angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(8, 10), static=TRUE, 
#' arg.lis1=list(offsetRow=0.01, offsetCol=0.01))
#' # Clusters containing "ENSGALG00000019846".
#' cut_dendro(mhm.res$rowDendrogram, h=15, 'ENSGALG00000019846')
#'
#' # Interactive matrix heatmap.
#' \donttest{ matrix_hm(ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), data=se.sub.mat, 
#' angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(8, 10), static=FALSE, 
#' arg.lis1=list(offsetRow=0.01, offsetCol=0.01)) 
#' }
#'
#' # In case the interactive heatmap is not automatically opened, run the following code snippet.
#' # It saves the heatmap as an HTML file that is assigned to the "file" argument.
#' \donttest{
#' mhm <- matrix_hm(ID=c('ENSGALG00000019846', 'ENSGALG00000000112'), data=se.sub.mat, 
#' angleCol=80, angleRow=35, cexRow=0.8, cexCol=0.8, margin=c(8, 10), static=FALSE, 
#' arg.lis1=list(offsetRow=0.01, offsetCol=0.01))
#' htmlwidgets::saveWidget(widget=mhm, file='mhm.html', selfcontained=FALSE)
#' browseURL('mhm.html')
#' }
#'
#' ## Adjacency matrix and module identification 
#' adj.mod <- adj_mod(data=se.sub.mat)
#'
#' # The adjacency is a measure of co-expression similarity between genes, where larger
#' # value denotes higher similarity.
#' adj.mod[['adj']][1:3, 1:3]
#'
#' # The modules are identified at four sensitivity levels (ds=0, 1, 2, or 3). From 0 to 3, 
#' # more modules are identified but module sizes are smaller. The 4 sets of module 
#' # assignments are returned in a data frame, where column names are sensitivity levels. 
#' # The numbers in each column are module labels, where "0" means genes not assigned to 
#' # any module.
#' adj.mod[['mod']][1:3, ]
#'
#' # Static network graph. Nodes are genes and edges are adjacencies between genes. 
#' # The thicker edge denotes higher adjacency (co-expression similarity) while larger node
#' # indicates higher gene connectivity (sum of a gene's adjacencies with all its direct 
#' # neighbors). The target gene is labeled by "_target".
#' network(ID="ENSGALG00000019846", data=se.sub.mat, adj.mod=adj.mod, adj.min=0, 
#' vertex.label.cex=1.5, vertex.cex=4, static=TRUE)
#'
#' # Interactive network. The target gene ID is appended "_target".  
#' \donttest{ network(ID="ENSGALG00000019846", data=se.sub.mat, adj.mod=adj.mod, static=FALSE) }
#'
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

#' @export
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom stats cor dist

submatrix <- function(data, assay.na=NULL, ID, p=0.3, n=NULL, v=NULL, fun='cor', cor.absolute=FALSE, arg.cor=list(method="pearson"), arg.dist=list(method="euclidean"), file=NULL) {
  options(stringsAsFactors=FALSE)
  # Process data.
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) {
    data <- s4se(data); if (is(data, 'character')) { warning(data); return(data) }
    # dat.lis <- check_data(data=data); data <- dat.lis$dat; ann <- dat.lis$row.meta
  }; dat.ipt <- data
  if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
    if (is.null(assay.na)) {
      if (length(assays(data)) > 1) stop("Please specify which assay to use by assigning the assay name to 'assay.na'!") else if (length(assays(data)) == 1) assay.na <- 1
    }
    data <- as.data.frame(assays(data)[[assay.na]]) 
    if (any(duplicated(rownames(data)))) stop('Please use function \'aggr_rep\' to aggregate replicates!')
  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", "dgCMatrix", "SummarizedExperiment", or "SingleCellExperiment" except that "spatial_hm" also accepts a "vector".') }

  if (nrow(data)<5) cat('Warning: variables of sample/condition are less than 5! \n')
  na <- NULL; len <- nrow(data)
  if (len>=50000) cat('More than 50,000 rows are detected in data. Computation may take a long time! \n')
  if (sum(c(!is.null(p), !is.null(n), !is.null(v)))>1) return('Please only use one of \'p\', \'n\', \'v\' as the threshold!')

  if (fun=='cor') {
    dat.t <- t(data)
    m <- do.call(cor, c(x=list(dat.t[, ID]), y=list(dat.t), arg.cor))
    if (length(ID)==1) rownames(m) <- ID; m <- t(m)
    if (cor.absolute==TRUE) { m1 <- m; m <- abs(m) }
    na <- sub_na(mat=m, ID=ID, p=p, n=n, v=v)
    if (cor.absolute==TRUE) m <- m1
    m <- t(m)
  } else if (fun=='dist') { 
    m <- data.frame(); for (i in ID) {
      dis <- NULL; for (j in rownames(data)) {
        # Distances are transformed to negative.
        dis0 <- -do.call(dist, c(x=list(data[c(i, j), ]), arg.dist))[1]
        dis <- c(dis, dis0)
      }; m <- rbind(m, dis)
    }; row.names(m) <- ID; colnames(m) <- rownames(data)
    if (!is.null(v)) v <- -v
    na <- sub_na(mat=t(m), ID=ID, p=p, n=n, v=v); m <- -m
  }
  if (!is.null(file)) {
    file <- normalizePath(file, winslash="/", mustWork=FALSE)
    if (file.exists(file)) stop(paste0(file, ' already exists!'))
    write.table(data[na, ], file, col.names=TRUE, row.names=TRUE, sep='\t')
  }; return(dat.ipt[na, ])
  # return(cbind(data[na, ], ann[na, , drop=FALSE]))
}







