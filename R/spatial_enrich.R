#' Identifying spatially enriched or depleted biomolecules 
#'
#' @description
#' The spatial enrichment (SpEn) is designed to detect spatially enriched or depleted biomolecules (genes, proteins, etc) for chosen spatial features (cellular compartments, tissues, organs, \emph{etc}). It compares each feature with all other reference features. The biomolecules significantly up- or down-regulated in one feature relative to reference features are denoted spatially enriched or depleted respectively. The underlying differential expression analysis methods include edgeR (Robinson et al, 2010), limma (Ritchie et al, 2015), DESeq2 (Love et al, 2014), and distinct (Tiberi et al, 2020). By querying a feature of interest from the enrichment results, the enriched or depleted biomolecules will be returned. \cr In addition, the SpEn is also able to identify biomolecules enriched or depleted in experiment vairables in a similar manner.   
#'
#' `sf_var()` subsets data according to given spatial features and variables.
#' 
#' `spatial_enrich()` detects enriched or depleted biomolecules for each given spatial feature.
#' 
#' `query_enrich()` queries enriched or depleted biomolecules in the enrichment results returned by \code{spatial_enrich} for a chosen spatial feature.
#'
#' `ovl_enrich()` plots overlap of enrichment results across spatial features in form an upset plot, overlap matrix, or Venn diagram.
#'
#' `graph_line()` plots expression values of chosen biomolecules in a line graph.

#' @param data
#' \describe{
#'  \item{\code{sf_var}}{
#'   A \code{SummarizedExperiment} object. The \code{colData} slot is required to contain at least two columns of spatial features and experiment variables respectively.   
#'  }
#'  \item{\code{spatial_enrich}}{
#'  A \code{SummarizedExperiment} object returned by \code{sf_var}.
#'  } 
#'  \item{\code{graph_line}}{
#'  A \code{data.frame}, where rows are biomolecules and columns are spatial features.
#'  } 
#' }

#' @return 
#' \describe{
#'  \item{`sf_var`}{
#'     A \code{SummarizedExperiment} object.
#'    }
#'  \item{`spatial_enrich`}{
#'     A \code{list} object.
#'    }
#'  \item{`query_enrich`}{
#'     A \code{SummarizedExperiment} object.
#'    }
#'  \item{`ovl_enrich`}{
#'     An UpSet plot, overlap matrix plot, or Venn diagram.
#'    }
#'  \item{`graph_line`}{
#'     A ggplot.
#'    }
#' }

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @name SpatialEnrichment
#' @rdname SpatialEnrichment
#' @aliases sf_var spatial_enrich query_enrich ovl_enrich graph_line

#' @references
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
#' Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#' Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
#' Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' Simone Tiberi and Mark D. Robinson. (2020). distinct: distinct: a method for differential analyses via hierarchical permutation tests. R package version 1.2.0. https://github.com/SimoneTiberi/distinct
#' Nils Gehlenborg (2019). UpSetR: A More Scalable Alternative to Venn and Euler Diagrams for Visualizing Intersecting Sets. R package version 1.4.0. https://CRAN.R-project.org/package=UpSetR
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

#' @examples
#'
#' ## In the following examples, the toy data come from an RNA-seq analysis on development of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, it is
#' ## included in this package. The complete raw count data are downloaded using the R package
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769".   
#'
#' library(SummarizedExperiment) 
#' # Access the count table. 
#' cnt.chk <- read.table(system.file('extdata/shinyApp/data/count_chicken.txt', package='spatialHeatmap'), header=TRUE, row.names=1,sep='\t')
#' cnt.chk[1:3, 1:5]
#' # A targets file describing spatial features and conditions is required for toy data. It should be made
#' # based on the experiment design, which is accessible through the accession number 
#' # "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in this
#' # package and accessed below. 
#'
#' # Access the example targets file. 
#' tar.chk <- read.table(system.file('extdata/shinyApp/data/target_chicken.txt', package='spatialHeatmap'), header=TRUE, row.names=1, sep='\t') 
#' # Every column in count table corresponds with a row in targets file. 
#' tar.chk[1:5, ]
#' # Store count data and targets file in "SummarizedExperiment".
#' se.chk <- SummarizedExperiment(assay=cnt.chk, colData=tar.chk)
#' # The "rowData" slot can store a data frame of gene metadata, but not required. Only the 
#' # column named "metadata" will be recognized. 
#' # Pseudo row metadata.
#' metadata <- paste0('meta', seq_len(nrow(cnt.chk))); metadata[1:3]
#' rowData(se.chk) <- DataFrame(metadata=metadata)
#'
#' # Subset the count data by features (brain, heart, kidney) and variables (day10, day12).
#' # By setting com.by='ft', the subsequent spatial enrichment will be performed across 
#' # features with the variables as replicates. 
#' data.sub <- sf_var(data=se.chk, feature='organism_part', ft.sel=c('brain', 'kidney',
#'  'heart', 'liver'), variable='age', var.sel=c('day10', 'day35'), com.by='ft')
#' 
#' ## As conventions, raw sequencing count data should be normalized and filtered to
#' ## reduce noise. Since normalization will be performed in spatial enrichment, only filtering
#' ## is required.  
#'
#' # Filter out genes with low counts and low variance. Genes with counts over 5 in
#' # at least 10% samples (pOA), and coefficient of variance (CV) between 3.5 and 100 are 
#' # retained.
#' data.sub.fil <- filter_data(data=data.sub, sam.factor='organism_part', con.factor='age',
#' pOA=c(0.1, 5), CV=c(0.7, 100))
#'
#' # Spatial enrichment for every spatial feature with 1 outlier allowed.  
#' enr.res <- spatial_enrich(data.sub.fil, method=c('edgeR'), norm='TMM', log2.fc=1, fdr=0.05, outliers=1)
#' # Overlaps of enriched genes across features.
#' ovl_enrich(enr.res, type='up', plot='upset')
#' # Query the results for brain.
#' en.brain <- query_enrich(enr.res, 'brain')
#' rowData(en.brain)[1:3, c('type', 'total', 'method')] 
#'
#' # Read aSVG image into an "SVG" object.
#' svg.chk <- system.file("extdata/shinyApp/data", "gallus_gallus.svg", 
#' package="spatialHeatmap")
#' svg.chk <- read_svg(svg.chk)
#' # Plot an enrichment SHM.
#' dat.enrich <- SPHM(svg=svg.chk, bulk=en.brain)
#' shm(data=dat.enrich, ID=rownames(en.brain)[1], legend.r=1, legend.nrow=7, sub.title.size=10, ncol=2, bar.width=0.09, lay.shm='gene')
#' # Line graph of gene expression profile.
#' graph_line(assay(en.brain[1, , drop=FALSE]), lgd.pos='bottom')
NULL

#' @rdname SpatialEnrichment
#' @param feature The column name in the \code{colData} slot of \code{SummarizedExperiment} that contains spatial features.
#' @param ft.sel A vector of spatial features to choose. 
#' @param variable The column name in the \code{colData} slot of \code{SummarizedExperiment} that contains experiment variables.
#' @param var.sel A vector of variables to choose. 
#' @param com.by One of \code{ft}, \code{var}, or \code{ft.var}. If \code{ft}, the enrichment is performed for each spatial feature and the variables are treated as replicates. If \code{var} the enrichment is performed for each variable and spatial features are treated as replicates. If \code{ft.var}, spatial features (tissue1, tissue2) and variables (var1, var2) are combined such as tissue1__var1, tissue1_var2, tissue2__var1, tissue2_var2. The enrichment is performed for each combination.
 
#' @export
#' @importFrom SummarizedExperiment colData colData<- assayNames assayNames<-

sf_var <- function(data, feature, ft.sel=NULL, variable=NULL, var.sel=NULL, com.by='ft') {
  # save(data, feature, ft.sel, variable, var.sel, com.by, file='sf.var.arg')
  fts <- vars <- NULL; cdat <- colData(data)
  if (!is.null(feature)) fts <- cdat[, feature] <- make.names(cdat[, feature])
  if (!is.null(variable)) vars <- cdat[, variable] <- make.names(cdat[, variable])
  
  if (is.null(ft.sel)) ft.sel <- unique(fts)[seq_len(2)] else if (ft.sel[1]=='all') ft.sel <- unique(fts)
  lgc.var <- !is.null(variable) & !is.null(var.sel)
  if (lgc.var) if (var.sel[1]=='all') var.sel <- unique(cdat[, variable])  
  # Subset data according to selected features and variables.
  if (lgc.var) idx <- cdat[, feature] %in% ft.sel & cdat[, variable] %in% var.sel else idx <- cdat[, feature] %in% ft.sel

  if (sum(idx)==0) { 
    msg <- 'Ensure "ft.sel" and "var.sel" are correct!'
    warning(msg); return(msg)
  }
  colData(data) <- cdat; data <- data[, idx]
  cdat <- colData(data); cdat.na <- colnames(cdat)
  # Re-order colData slot.
  if (! com.by %in% c('ft', 'var', 'ft.var')) stop('"com.by" is one of "ft", "var", or "ft.var"!')
  if (com.by=='ft') {
    cdat$com.by <- cdat[, feature]
    cdat <- cdat[, c('com.by', feature, variable, setdiff(cdat.na, c(feature, variable)))]
  } else if (com.by=='var') {
    cdat$com.by <- cdat[, variable]
    cdat <- cdat[, c('com.by', variable, feature, setdiff(cdat.na, c(feature, variable)))]
  } else if (com.by=='ft.var') { # Combine features and variables.
    ft.fct <- paste0(cdat[, feature], '__', cdat[, variable])
    cdat$com.by <- ft.fct
    cdat <- cdat[, c('com.by', feature, variable, setdiff(cdat.na, c(feature, variable)))]
  }
  colData(data) <- cdat
  if (!is.null(variable)) colnames(data) <- paste0(cdat[, feature], '__', cdat[, variable]) else colnames(data) <- cdat[, feature]
  # Name the assay: required in distinct.
  if (is.null(assayNames(data))) assayNames(data) <- 'count'
  return(data)
}

#' @rdname SpatialEnrichment
#' @param method One of \code{edgeR}, \code{limma}, \code{DESeq2}, \code{distinct}. 
#' @param norm The normalization method (\code{TMM}, \code{RLE}, \code{upperquartile}, \code{none}) in edgeR. The default is \code{TMM}. Details: https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors. 

#' @param log2.trans.dis Logical, only applicable when \code{method='distinct'}. If \code{TRUE} the count data is transformed to log-2 scale.

#' @param aggr One of \code{mean} (default) or \code{median}. The method to aggregated replicates in the assay data.  
#' @param log2.trans Logical. If \code{TRUE} (default), the aggregated data (see \code{aggr}) is transformed to log2-scale and will be further used for plotting SHMs. 
#' @param p.adjust The method (\code{holm}, \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH}, \code{BY}, \code{fdr}, or \code{none}) for adjusting p values in multiple hypothesis testing. The default is \code{BH}.
#' @param log2.fc The log2-fold change cutoff. The default is 1.
#' @param fdr The FDR cutoff. The default is 0.05.
#' @param outliers The number of outliers allowed in the references. If there are too many references, there might be no enriched/depleted biomolecules in the query feature. To avoid this, set a certain number of outliers.  

#' @export
#' @importFrom SummarizedExperiment colData
 
spatial_enrich <- function(data, method=c('edgeR'), norm='TMM', log2.trans.dis=TRUE, log2.fc=1, p.adjust='BH', fdr=0.05, outliers=0, aggr='mean', log2.trans=TRUE) {
  #save(data, method, norm, log2.trans.dis, log2.fc, p.adjust, fdr, aggr, log2.trans, file='spatial.enrich.arg')
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'dgCMatrix')|is(data, 'DataFrame')) {                 
    data <- SummarizedExperiment(assays=list(data=data))                                                           
  }
  edg <- dsq <- lim <- dis <- NULL
  if ('edgeR' %in% method) { cat('edgeR ... \n')
    edg <- edgeR(data, method.norm=norm, com.factor='com.by', method.adjust=p.adjust, return.all=FALSE, log2.fc=log2.fc, fdr=fdr, outliers=outliers)
    cat('Done! \n')
  }
  if ('limma' %in% method) { cat('limma ... \n')
    lim <- limma(data, m.array=FALSE, method.norm=norm, com.factor='com.by', method.adjust=p.adjust, return.all=FALSE, log2.fc=log2.fc, fdr=fdr, outliers=outliers)
    cat('Done! \n') 
  }
  if ('DESeq2' %in% method) { cat('DESeq2 ... \n')
    dsq <- deseq2(data, com.factor='com.by', method.adjust=p.adjust, return.all=FALSE, log2.fc=log2.fc, fdr=fdr, outliers=outliers)
    cat('Done! \n')
  }
  if ('distinct' %in% method) { cat('distinct ... \n')
    dis <- distt(data, norm.fun='CNF', par.list=list(method=norm), log2.trans=log2.trans.dis, com.factor='com.by', return.all=FALSE, log2.fc=log2.fc, fdr=fdr, outliers=outliers)
    cat('Done! \n')
  }

  lis <- list(edgeR=edg, limma=lim, DESeq2=dsq, distinct=dis)[c('edgeR', 'limma', 'DESeq2', 'distinct') %in% method]
  names(lis) <- 'result'
  dat.nor <- norm_data(data, norm.fun='CNF', par.list=list(method=norm), log2.trans=log2.trans)
  dat.aggr <- aggr_rep(dat.nor, sam.factor=NULL, con.factor=NULL, aggr=aggr)
  lis$data <- dat.aggr; return(lis)
}


#' @rdname SpatialEnrichment
#' @param res Enrichment results returned by \code{spatial_enrich}.
#' @param query A spatial feature for query.
 
#' @export
#' @importFrom SummarizedExperiment rowData rowData<-

query_enrich <- function(res, query) {
  up <- res$result[[query]]$up; dn <- res$result[[query]]$down
  if (length(intersect(rownames(up), rownames(dn)))>0) {
    msg <- 'Same biomolecules detected as enriched and depleted! Please reduce the outliers.'
    warning(msg); return(msg)
  }
  df.deg <- DataFrame(rbind(up, dn))
  if (nrow(df.deg)==0) { 
    msg <- 'No enriched/depleted results are detected!'
    warning(msg); return(msg)
  }
  df.deg$total <- as.numeric(df.deg$total)
  na.deg <- rownames(df.deg)
  data <- res$data; data <- data[na.deg, ]
  rdat <- rowData(data)
  if (nrow(rdat) > 0) rowData(data) <- cbind(df.deg, rdat)
  if (nrow(rdat) == 0) rowData(data) <- df.deg
  cat('Done! \n')
  # The "feature__factor" columns in df.deg will be extacted in 'spatial_hm'.
  return(data)
}


#' @rdname SpatialEnrichment
#' @param type One of \code{up} (default) or \code{down}, which refers to up- or down-regulated biomolecules.
#' @param plot One of \code{upset} (default), \code{matrix}, or \code{venn}, corresponding to upset plot, overlap matrix, or Venn diagram respectively.
#' @inheritParams UpSetR::upset
#' @param upset.arg A \code{list} of additional arguments passed to \code{\link[UpSetR]{upset}}.
#' @inheritParams gplots::venn
#' @param venn.arg A \code{list} of additional arguments passed to \code{\link[gplots]{venn}}.
#' @param axis.agl The angle of axis text in overlap matrix.
#' @param font.size The font size of all text in overlap matrix.
#' @param line.size The line thickness in overlap matrix. 
#' @param cols A vector of two colors indicating low and high values in the overlap matrix respectively. The default is \code{c("lightcyan3", "darkorange")}.

 
#' @export
#' @importFrom UpSetR upset fromList

ovl_enrich <- function(res, type='up', plot='upset', order.by="freq", nintersects=40, point.size=3, line.size=1, mb.ratio=c(0.6, 0.4), text.scale=1.5, upset.arg=list(), show.plot=TRUE, venn.arg=list(), axis.agl=45, font.size=5, cols=c("lightcyan3", "darkorange")) {
  sams <- names(res$result)
  if (type=='up') lis <- lapply(sams, function(x) {rownames(res$result[[x]]$up)} )
  if (type=='down') lis <- lapply(sams, function(x) {rownames(res$result[[x]]$down)} )
  names(lis) <- sams
  if (plot=='upset') {
    upset.arg <- c(list(data=fromList(lis), order.by=order.by, nintersects=nintersects, point.size=point.size, line.size=line.size, mb.ratio=mb.ratio, text.scale=text.scale), upset.arg)
    ups <- do.call('upset', upset.arg)
    return(ups)
  } else if (plot=='matrix') {
    g <- deg_ovl_mat(lis, axis.agl=axis.agl, font.size=font.size, cols=cols); return(g)
  } else if (plot=='venn') {
    venn.arg <- c(list(data=lis, show.plot=show.plot), venn.arg)
    inter <- do.call('venn', venn.arg)
    invisible(inter)
  }
}

#' Plot the overlap matrix for enrichment results across spatial features
#'
#' @param deg.lis The list of all up- and down-regulated biomolecules organized by spatial features returned by \code{spatial_enrich}. 
#' @param axis.agl The angle of axis text.
#' @param font.size The font size of all text in overlap matrix.
#' @param cols A vector of two colors indicating low and high values in the overlap matrix respectively. The default is \code{c("lightcyan3", "darkorange")}.

#' @return An image of ggplot.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' \cr Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber, Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2020). gplots: Various R Programming Tools for Plotting Data. R package version 3.1.1. https://CRAN.R-project.org/package=gplots 
#' \cr Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient theme_minimal theme element_text coord_fixed geom_text element_blank

deg_ovl_mat <- function(deg.lis, axis.agl, font.size, cols) {
  Var1 <- Var2 <- value <- NULL
  mat <- vapply(names(deg.lis), function(x) vapply(names(deg.lis), function(y) length(intersect(deg.lis[[x]], deg.lis[[y]])), numeric(1)), numeric(length(deg.lis)))
  mel <- reshape2::melt(mat) 
  g <- ggplot(data=mel, aes(x=Var1, y=Var2, fill=value))+geom_tile(colour="white")+scale_fill_gradient(low=cols[1], high=cols[2])+ theme_minimal()+theme(axis.text=element_text(angle=axis.agl, vjust=1, size=font.size+12, hjust=1))+coord_fixed()+geom_text(aes(Var2, Var1, label=value), color="black", size=font.size)+theme(axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.ticks=element_blank(), text = element_text(size=font.size+17)); return(g)
}


#' @rdname SpatialEnrichment
#' @param scale The method to scale the data. If \code{none} (default), no scaling. If \code{row}, each row is scaled independently. If \code{all}, all rows are scaled as a whole.
#' @param x.title,y.title The title of X-axis and Y-axis respectively.
#' @param text.size The font size of all text.
#' @param text.angle The angle of axis text.
#' @param lgd.pos The position of legend. The default is \code{right}.
#' @param lgd.guide The \code{\link[ggplot2]{guides}} function in \code{ggplot2} for customizing legends. 

#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_line theme labs element_text element_rect element_line

graph_line <- function(data, scale='none', x.title='Samples/conditions', y.title='Assay values', text.size=15, text.angle=60, lgd.pos='right', lgd.guide=guides(color=guide_legend(nrow=1, byrow=TRUE, title=NULL))) {
  Samples <- Value <- Genes <- NULL 
  if (all(c('type', 'total', 'method') %in% colnames(data))) { # Data frame of spatial enrichment.
    data <- data[, !colnames(data) %in% c('type', 'total', 'metadata', 'method'), drop=FALSE] }
  if (scale=='row') data <- t(scale(t(data))) else if (scale=='all') data <- scale_all(data)
  # convert to long format. 
  df.long <- reshape2::melt(as.matrix(data))
  colnames(df.long) <- c('Genes', 'Samples', 'Value')
  # The levels in melted data frame has the same order (left to right) with row order (top to bottom) in original data frame before melted. 
  # Colours map to variables in original data frame before melted.
  # Possible: the colour order (left to right) matches with the row order (top to bottom) in original data frame before melted, but the coloured lined is plotted in the order of levels (left to right) in melted data frame. 
  # if (length(cols)<nrow(data)) cols <- diff_cols(nrow(data)) 
  # Custom colours: scale_color_manual(values=cols) 
  g <- ggplot(data=df.long, aes(x=Samples, y=Value, colour=Genes, group=Genes))+geom_line()+labs(title="", x=x.title, y=y.title)+theme(legend.position=lgd.pos, axis.text=element_text(size=text.size), axis.title=element_text(size=text.size, face="bold"), axis.text.x=element_text(angle=text.angle, hjust=1), legend.title=element_text(size=text.size-3), legend.text=element_text(size=text.size-3), panel.background = element_rect(fill = "gray95", colour = "gray95", linewidth = 0.5, linetype = "solid"), panel.grid.major = element_line(linewidth= 0.5, linetype = 'solid', colour = "white"), panel.grid.minor = element_line(linewidth = 0.5, linetype = 'solid', colour ="white")) + lgd.guide 
  return(g)
}
