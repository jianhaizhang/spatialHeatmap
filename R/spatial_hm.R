#' Plot Spatial Heatmaps
#'
#' The input are a pair of aSVG image and formatted data (\code{vector}, \code{data frame}, \code{SummarizedExperiment}). The former is an annotated SVG file (aSVG) where spatial features are represented by shapes and assigned unique identifiers, while the latter are numeric values measured from these spatial features and organized in specific formats. In biological cases, aSVGs are anatomical or cell structures, and data are measurements of genes, proteins, metabolites, \emph{etc}. in different samples (\emph{e.g.} cells, tissues). Data are mapped to the aSVG according to identifiers of assay samples and aSVG features. Only the data from samples having matching counterparts in aSVG features are mapped. The mapped features are filled with colors translated from the data, and the resulting images are termed spatial heatmaps. Note, "sample" and "feature" are two equivalent terms referring to cells, tissues, organs etc. where numeric values are measured. Matching means a target sample in data and a target spatial feature in aSVG image have the same identifier. \cr This function is designed as much flexible as to achieve optimal visualization. For example, subplots of spatial heatmaps can be organized by gene or condition for easy comparison, in multi-layer anotomical structures selected tissues can be set transparent to expose burried features, color scale is customizable to highlight difference among features. This function also works with many other types of spatial data, such as population data plotted to geographic maps.

#' @param svg.path The path of the aSVG image. \emph{E.g.}: system.file("extdata/shinyApp/example", "gallus_gallus.svg", package="spatialHeatmap"). See \code{\link{return_feature}} for details on how to download aSVGs from the EBI aSVG repository (https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg) directly.

#' @inheritParams filter_data
#' @inheritParams grob_list
#' @inheritParams col_bar

#' @param ID A character vector of assyed items (\emph{e.g.} genes, proteins) whose values are used to color the aSVG.
#' @param col.com A character vector of the color components used to build the color scale. The default is c('purple', 'yellow', 'blue').
#' @param col.bar "selected" or "all", the former uses values of input assayed items to build the color scale while the latter uses all values from the data. The default is "selected".
#' @param data.trans "log2", "exp2", or NULL. If colors across samples cannot distinguish due to low variance or outliers, transform the data by "log2" or "2-base expoent" (exp2). Default is NULL (data will not be transformed).
#' @param bar.width The width of color bar. Default if 0.08.
#' @param width A numeric of each subplot width. The default is 1.
#' @param height A numeric of each subplot height. The default is 1.
#' @param legend.r A numeric to adjust the dimension of the legend plot. Default is 1. The larger, the higher ratio of width to height.
#' @param lay.shm "gene" or "con", the former organizes spatial heatmaps by genes, proteins, metabolites, \emph{etc}. while the latter by conditions/treatments applied to experiments.
#' @param ncol An integer, number of columns to display the spatial heatmaps, not including the legend plot.
#' @param verbose Logical, FALSE or TRUE. If TRUE the samples in data not colored in spatial heatmaps are printed to R console. Default is TRUE.

#' @return An image of spatial heatmap(s), a two-component list of the spatial heatmap(s) in \code{ggplot} format and a data frame of mapping between assayed samples and aSVG features.

#' @section Details:
#' See the package vignette (\code{browseVignettes('spatialHeatmap')}).  

#' @examples

#' ## In the following examples, the 2 toy data come from an RNA-seq analysis on developments of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are included in this package. The complete raw count data are downloaded using the R package ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as a "data frame" input to exemplify data with simple samples/conditions, while toy data2 as "SummarizedExperiment" to illustrate data involving complex samples/conditions.     
#'
#' ## Set up toy data.
#'
#' # Access toy data1.
#' cnt.chk.simple <- system.file('extdata/shinyApp/example/count_chicken_simple.txt', package='spatialHeatmap')
#' df.chk <- read.table(cnt.chk.simple, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
#' # Note the naming scheme "sample__condition" in columns, where "sample" and "condition" stands for organs and time points respectively.
#' df.chk[1:3, ]

#' # A column of gene annotation can be appended to the data frame, but is not required.  
#' ann <- paste0('ann', seq_len(nrow(df.chk))); ann[1:3]
#' df.chk <- cbind(df.chk, ann=ann)
#' df.chk[1:3, ]
#'
#' # Access toy data2. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]

#' # A targets file describing samples and conditions is required for toy data2. It should be made based on the experiment design, which is accessible through the accession number "E-MTAB-6769" in R package ExpressionAtlas. An example targets file is included in this package.  

#' # Access the targets file. 
#' tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
#' target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
#' # Note every column in toy data2 corresponds with a row in targets file. 
#' target.chk[1:5, ]
#' # Store toy data2 in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene annotation, but not required.
#' rowData(se.chk) <- DataFrame(ann=ann)
#'
#' ## As conventions, raw sequencing count data should be normalized, aggregated, and filtered fo reduce noise.
#'
#' # Normalize count data.
#' # The normalizing function "calcNormFactors" (McCarthy et al. 2012) with default settings is used.  
#' df.nor.chk <- norm_data(data=df.chk, norm.fun='CNF', data.trans='log2')
#' se.nor.chk <- norm_data(data=se.chk, norm.fun='CNF', data.trans='log2')

#' # Aggregate count data.
#' # Aggregate "sample__condition" replicates in toy data1.
#' df.aggr.chk <- aggr_rep(data=df.nor.chk, aggr='mean')
#' df.aggr.chk[1:3, ]

#' # Aggregate "sample_condition" replicates in toy data2, where "sample" is "organism_part" and "condition" is "age". 
#' se.aggr.chk <- aggr_rep(data=se.nor.chk, sam.factor='organism_part', con.factor='age', aggr='mean')
#' assay(se.aggr.chk)[1:3, 1:3]

#' # Filter genes with low counts and low variance. Genes with counts over 5 (log2 unit) in at least 1% samples (pOA), and coefficient of variance (CV) between 0.2 and 100 are retained.
#' # Filter toy data1.
#' df.fil.chk <- filter_data(data=df.aggr.chk, pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
#' # Filter toy data2.
#' se.fil.chk <- filter_data(data=se.aggr.chk, sam.factor='organism_part', con.factor='age', pOA=c(0.01, 5), CV=c(0.2, 100), dir=NULL)
#'
#' ## Spatial heatmaps.
#'
#' # The target chicken aSVG is downloaded from the EBI aSVG repository (https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg) directly with function "return_feature". It is included in this package and accessed as below. Details on how this aSVG is selected are documented in function "return_feature".
#' svg.chk <- system.file("extdata/shinyApp/example", "gallus_gallus.svg", package="spatialHeatmap")

#' # Plot spatial heatmaps on gene "ENSGALG00000019846".
#' # Toy data1. 
#' spatial_hm(svg.path=svg.chk, data=df.fil.chk, ID='ENSGALG00000019846', height=0.4, legend.r=1.7, sub.title.size=9, ncol=3, legend.nrow=2)
#' # Toy data2.
#' spatial_hm(svg.path=svg.chk, data=se.fil.chk, ID='ENSGALG00000019846', legend.r=1.5, sub.title.size=9, ncol=3, legend.nrow=2)
#'
#' # The data can also come as as a simple named vector. The following gives an example on a vector of 3 random values. 
#' # Random values.
#' vec <- sample(1:100, 3)
#' # Name the vector. The last name is assumed as a random sample without a matching feature in aSVG.
#' names(vec) <- c('brain', 'heart', 'notMapped')
#' vec
#' # Plot.
#' spatial_hm(svg.path=svg.chk, data=vec, ID='geneX', height=0.7, legend.r=1.5, sub.title.size=9, ncol=1, legend.nrow=2)


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1 \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg 
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' Guangchuang Yu (2020). ggplotify: Convert Plot to 'grob' or 'ggplot' Object. R package version 0.0.5. https://CRAN.R-project.org/package=ggplotify
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggplot geom_bar aes theme element_blank margin element_rect coord_flip scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs
#' @importFrom rsvg rsvg_ps 
#' @importFrom grImport PostScriptTrace 
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid grobTree unit
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' @importFrom ggplotify as.ggplot

spatial_hm <- function(svg.path, data, sam.factor=NULL, con.factor=NULL, ID, col.com=c('purple', 'yellow', 'blue'), col.bar="selected", bar.width=0.08, bar.title.size=10, data.trans=NULL, tis.trans=NULL, width=1, height=1, legend.r=1, sub.title.size=11, lay.shm="gene", ncol=2, sam.legend='identical', legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.5, legend.label.size=8, line.size=0.2, line.color='grey70', verbose=TRUE, ...) {

  x <- y <- color_scale <- tissue <- line.type <- NULL  
  # Extract and filter data.
  options(stringsAsFactors=FALSE) 
  if (is.vector(data)) {

    vec.na <- make.names(names(data)); if (is.null(vec.na)) stop("Please provide names for the input data!")
    if (!identical(vec.na, names(data))) cat('Syntactically valid column names are made! \n')
    if (any(duplicated(vec.na))) stop('Please make sure data names are unique!')
    form <- grepl("__", vec.na); if (sum(form)==0) { vec.na <- paste0(vec.na, '__', 'con'); con.na <- FALSE } else con.na <- TRUE
    data <- tryCatch({ as.numeric(data) }, warning=function(w) { stop("Please make sure input data are numeric!") }, error=function(e) { stop("Please make sure input data are numeric!") })
    if (is.null(ID)) stop('Please provide a name for the data!')
    gene <- as.data.frame(matrix(data, nrow=1, dimnames=list(ID, vec.na)))

  } else if (is(data, 'data.frame')|is(data, 'matrix')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- make.names(colnames(data))
    if (!identical(cna, colnames(data))) cat('Syntactically valid column names are made! \n')
    if (any(duplicated(cna))) stop('Please make sure column names are unique!')
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    gene <- na[!idx]; colnames(gene) <- cna <- cna[!idx]
    form <- grepl("__", cna); if (sum(form)==0) { colnames(gene) <- paste0(cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE

  } else if (is(data, 'SummarizedExperiment')) {

    gene <- assay(data); r.na <- rownames(gene); gene <- apply(gene, 2, as.numeric) # This step removes rownames of gene2.
    rownames(gene) <- r.na; cna <- colnames(gene) <- make.names(colnames(gene))
    if (!identical(cna, colnames(gene))) cat('Syntactically valid column names are made! \n')
    col.meta <- as.data.frame(colData(data), stringsAsFactors=FALSE)
    # Factors teated by paste0/make.names are vecters.
    if (!is.null(sam.factor) & !is.null(con.factor)) { cna.tar <- paste0(col.meta[, sam.factor], '__', col.meta[, con.factor]); con.na <- TRUE } else if (!is.null(sam.factor) & is.null(con.factor)) { sam.na <- as.vector(col.meta[, sam.factor]); cna.tar <- paste0(sam.na, "__", "con"); con.na <- FALSE } else if (is.null(sam.factor)) { form <- grepl("__", cna); if (sum(form)==0) { cna.tar <- paste0(cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE }
    if (exists('cna.tar')) { colnames(gene) <- make.names(cna.tar); if (!identical(cna.tar, make.names(cna.tar))) cat('Syntactically valid column names are made! \n') }
    if (any(duplicated(colnames(gene)))) stop('Please use function \'aggr_rep\' to aggregate \'sample__condition\' replicates!')

  }; gene <- as.data.frame(gene)
  if (!is.null(data.trans)) if (data.trans=='log2') { 
          
      g.min <- min(gene) 
      if (g.min<0) gene <- gene-g.min+1; if (g.min==0) gene <- gene+1; gene <- log2(gene)  

    } else if (data.trans=='exp2') gene <- 2^gene
 
    # Color bar.
    bar.len=1000
    if (col.bar=="all") geneV <- seq(min(gene), max(gene), len=bar.len) else if (col.bar=="selected") geneV <- seq(min(gene[ID, , drop=FALSE]), max(gene[ID, , drop=FALSE]), len=bar.len)
    col <- colorRampPalette(col.com)(length(geneV))
    cs.g <- col_bar(geneV=geneV, cols=col, width=1, bar.title.size=bar.title.size, mar=c(3, 0.1, 3, 0.1)); cs.grob <- ggplotGrob(cs.g)    

    # Only take the column names with "__".
    cname <- colnames(gene); form <- grepl('__', cname)
    con <- gsub("(.*)(__)(.*)", "\\3", cname[form]); con.uni <- unique(con)
    sam.uni <- unique(gsub("(.*)(__)(.*)", "\\1", cname))

    df_tis <- svg_df(svg.path=svg.path, feature=sam.uni)
    if (is.character(df_tis)) stop(df_tis)
    g.df <- df_tis[['df']]; tis.path <- df_tis[['tis.path']]

    not.map <- setdiff(sam.uni, unique(tis.path)); if (verbose==TRUE & length(not.map)>0) cat('Enrties not mapped:', paste0(not.map, collapse=', '), '\n')
    sam.com <- intersect(unique(tis.path), sam.uni) 

    idx.com <- vapply(cname, function(i) grepl(paste0(sam.com, '__', collapse='|'), i), FUN.VALUE=logical(1))
    map.gene <- gene[idx.com]; cna <- colnames(map.gene)
    map.sum <- data.frame(); for (i in ID) {

      df0 <- as.data.frame(t(map.gene[i, ])); colnames(df0) <- 'value'
      featureSVG <- gsub("(.*)(__)(.*)", "\\1", cna)
      rowID <- i; df1 <- data.frame(rowID=rowID, featureSVG=featureSVG)
      if (con.na==TRUE) { condition <- gsub("(.*)(__)(.*)", "\\3", cna); df1 <- cbind(df1, condition=condition) }
      df1 <- cbind(df1, df0); map.sum <- rbind(map.sum, df1)

    }; row.names(map.sum) <- NULL    
    grob.lis <- grob_list(gene=gene, con.na=con.na, geneV=geneV, coord=g.df, ID=ID, cols=col, legend.col=df_tis[['fil.cols']], tis.path=tis.path, tis.trans=tis.trans, sub.title.size=sub.title.size, sam.legend=sam.legend, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.label.size=legend.label.size, line.size=line.size, line.color=line.color, line.type=line.type, ...)
    g.arr <- lay_shm(lay.shm=lay.shm, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis[['grob.lis']], width=width, height=height, shiny=FALSE)
    cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc")) # "mm" is fixed, "npc" is scalable.
    g.lgd <- grob.lis[['g.lgd']]; grob.lgd <- ggplotGrob(g.lgd)
    # Layout matrix of legend.
    if (lay.shm=='gene') lay.lgd <- matrix(seq_len(ceiling(length(con.uni)/ncol)), byrow=FALSE)
    if (lay.shm=='con') lay.lgd <- matrix(seq_len(ceiling(length(ID)/ncol)), byrow=FALSE)
    lgd.arr <- arrangeGrob(grobs=list(grobTree(grob.lgd)), layout_matrix=lay.lgd, widths=unit(width, "npc"), heights=unit(width/legend.r, "npc"))
    shm <- grid.arrange(cs.arr, g.arr, lgd.arr, ncol=3, widths=unit(c(bar.width-0.005, (1-bar.width)/(ncol+1)*ncol, (1-bar.width)/(ncol+1)), 'npc')); shm
    lis <- list(spatial_heatmap=as.ggplot(shm), mapped_feature=map.sum); invisible(lis)

}
