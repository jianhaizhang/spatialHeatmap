#' Create Spatial Heatmaps
#'
#' The input are a pair of annotated SVG (aSVG) file and formatted data (\code{vector}, \code{data.frame}, \code{SummarizedExperiment}). In the former, spatial features are represented by shapes and assigned unique identifiers, while the latter are numeric values measured from these spatial features and organized in specific formats. In biological cases, aSVGs are anatomical or cell structures, and data are measurements of genes, proteins, metabolites, \emph{etc}. in different samples (\emph{e.g.} cells, tissues). Data are mapped to the aSVG according to identifiers of assay samples and aSVG features. Only the data from samples having matching counterparts in aSVG features are mapped. The mapped features are filled with colors translated from the data, and the resulting images are termed spatial heatmaps. Note, "sample" and "feature" are two equivalent terms referring to cells, tissues, organs \emph{etc.} where numeric values are measured. Matching means a target sample in data and a target spatial feature in aSVG have the same identifier. \cr This function is designed as much flexible as to achieve optimal visualization. For example, subplots of spatial heatmaps can be organized by gene or condition for easy comparison, in multi-layer anotomical structures selected tissues can be set transparent to expose burried features, color scale is customizable to highlight difference among features. This function also works with many other types of spatial data, such as population data plotted to geographic maps.

#' @param svg.path The path of aSVG file(s). \emph{E.g.}: system.file("extdata/shinyApp/example", "gallus_gallus.svg", package="spatialHeatmap"). Multiple aSVGs are also accepted, such as aSVGs depicting organs development across mutiple times. In this case, the aSVGs should be indexed with suffixes "_shm1", "_shm2", ..., such as "arabidopsis.thaliana_organ_shm1.svg", "arabidopsis.thaliana_organ_shm2.svg", and the paths of these aSVGs be provided in a character vector.  \cr See \code{\link{return_feature}} for details on how to directly download aSVGs from the EBI aSVG repository \url{https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg} and spatialHeatmap aSVG Repository \url{https://github.com/jianhaizhang/spatialHeatmap_aSVG_Repository} developed in this project.

#' @param bar.title.size A numeric of color bar title size. The default is 0. 
#' @param bar.value.size A numeric of value size in the color bar y-axis. The default is 10.

#' @inheritParams filter_data
#' @inheritParams htmlwidgets::saveWidget
#' @inheritParams ggplot2::theme

#' @param position.2nd The position of the secondary legend. One of "top", "right", "bottom", "left", or a two-component numeric vector. The default is "bottom". Applies to the static image and video.
#' @param legend.nrow.2nd An integer of rows of the secondary legend keys. Applies to the static image and video.
#' @param legend.ncol.2nd An integer of columns of the secondary legend keys. Applies to the static image and video.
#' @param legend.key.size.2nd A numeric of legend key size. The default is 0.03. Applies to the static image and video.
#' @param add.feature.2nd Logical TRUE or FALSE. Add feature identifiers to the secondary legend or not. The default is FALSE. Applies to the static image.
#' @param legend.text.size.2nd A numeric of the secondary legend text size. The default is 10. Applies to the static image and video.
#' @param angle.text.key.2nd A value of angle of key text in the secondary legend. Default is 0. Applies to the static image and video.
#' @param position.text.key.2nd The position of key text in the secondary legend, one of "top", "right", "bottom", "left". Default is "right". Applies to the static image and video.


#' @param angle.text.key A value of key text angle in legend plot. The default is NULL, equivalent to 0.
#' @param position.text.key The position of key text in legend plot, one of "top", "right", "bottom", "left". Default is NULL, equivalent to "right".
#' @param legend.value.vdo Logical TRUE or FALSE. If TRUE, the numeric values of matching spatial features are added to video legend. The default is NULL.
#' @param label Logical. If TRUE, spatial features having matching samples are labeled by feature identifiers. The default is FALSE. It is useful when spatial features are labeled by similar colors. 
#' @param label.size The size of spatial feature labels in legend plot. The default is 4.
#' @param label.angle The angle of spatial feature labels in legend plot. Default is 0.
#' @param hjust The value to horizontally adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param vjust The value to vertically adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param opacity The transparency of colored spatial features in legend plot. Default is 1. If 0, features are totally transparent.
#' @param key Logical. The default is TRUE and keys are added in legend plot. If \code{label} is TRUE, the keys could be removed. 


#' @param ft.trans A character vector of tissue/spatial feature identifiers that will be set transparent. \emph{E.g} c("brain", "heart"). This argument is used when target features are covered by  overlapping features and the latter should be transparent.
#' @param lis.rematch A list for rematching features. In each slot, the slot name is an existing feature in the data, and the slot contains a vector of features in aSVG that will be rematched to the feature in the slot name. \emph{E.g.} \code{list(featureData1 = c('featureSVG1', 'featureSVG2'), featureData2 = c('featureSVG3'))}, where features \code{c('featureSVG1', 'featureSVG2')}, \code{c('featureSVG3')} in the aSVG are rematched to features \code{'featureData1'}, \code{'featureData2'} in data, respectively.
#' @param tis.trans This argument is deprecated and replaced by \code{ft.trans}. 
#' @param sub.title.size A numeric of the subtitle font size of each individual spatial heatmap. The default is 11.
#' @param ft.legend One of "identical", "all", or a character vector of tissue/spatial feature identifiers from the aSVG file. The default is "identical" and all the identical/matching tissues/spatial features between the data and aSVG file are colored in the legend plot. If "all", all tissues/spatial features in the aSVG are shown. If a vector, only the tissues/spatial features in the vector are shown.
#' @param legend.ncol An integer of the total columns of keys in the legend plot. The default is NULL. If both \code{legend.ncol} and \code{legend.nrow} are used, the product of the two arguments should be equal or larger than the total number of shown spatial features.
#' @param legend.nrow An integer of the total rows of keys in the legend plot. The default is NULL. It is only applicable to the legend plot. If both \code{legend.ncol} and \code{legend.nrow} are used, the product of the two arguments should be equal or larger than the total number of matching spatial features.
#' @param legend.key.size A numeric of the legend key size ("npc"), applicable to the legend plot. The default is 0.02. 
#' @param legend.text.size A numeric of the legend label size, applicable to the legend plot. The default is 12.
#' @param line.size The thickness of each shape outline in the aSVG is maintained in spatial heatmaps, \emph{i.e.} the stroke widths in Inkscape. This argument is the extra thickness added to all outlines. Default is 0.2 in case stroke widths in the aSVG are 0.
#' @param line.color A character of the shape outline color. Default is "grey70".
#' @param legend.plot.title The title of the legend plot. The default is 'Legend'. 
#' @param legend.plot.title.size The title size of the legend plot. The default is 11.


#' @param ID A character vector of assyed items (\emph{e.g.} genes, proteins) whose abudance values are used to color the aSVG.
#' @param col.com A character vector of the color components used to build the color scale. The default is c('yellow', 'orange', 'red').
#' @param col.bar One of "selected" or "all", the former uses values of \code{ID} to build the color scale while the latter uses all values from the data. The default is "selected".
#' @param cores The number of CPU cores for parallelization, relevant for aSVG files with size larger than 5M. The default is NA, and the number of used cores is 1 or 2 depending on the availability. 
#' @param trans.scale One of "log2", "exp2", "row", "column", or NULL, which means transform the data by "log2" or "2-base expoent", scale by "row" or "column", or no manipuation respectively. This argument should be used if colors across samples cannot be distinguished due to low variance or outliers. 
#' @param bar.width The width of color bar that ranges from 0 to 1. The default is 0.08.

#' @param legend.plot A vector of suffix(es) of aSVG file name(s) such as \code{c('shm1', 'shm2')}. Only aSVG(s) whose suffix(es) are assigned to this arugment will have a legend plot on the right. The default is \code{all} and each aSVG will have a legend plot. If NULL, no legend plot is shown.
#' @param legend.r A numeric (between -1 and 1) to adjust the legend plot size. The default is 0.2.
#' @param legend.2nd Logical, TRUE or FALSE. If TRUE, the secondary legend is added to each spatial heatmap, which are the numeric values of each matching spatial features. The default its FALSE. Only applies to the static image.
#' @param lay.shm One of "gene", "con", or "none". If "gene", spatial heatmaps are organized by genes proteins, or metabolites, \emph{etc.} and conditions are sorted whithin each gene. If "con", spatial heatmaps are organized by the conditions/treatments applied to experiments, and genes are sorted winthin each condition. If "none", spaital heatmaps are organized by the gene order in \code{ID} and conditions follow the order they appear in \code{data}. 
#' @param ncol An integer of the number of columns to display the spatial heatmaps, which does not include the legend plot.

#' @param verbose Logical, FALSE or TRUE. If TRUE the samples in data not colored in spatial heatmaps are printed to R console. Default is TRUE.
#' @param out.dir The directory to save interactive spatial heatmaps as independent HTML files and videos. Default is NULL, and the HTML files and videos are not saved.
#' @param animation.scale A numeric to scale the spatial heatmap size in the HTML files. The default is 1, and the height is 550px and the width is calculated according to the original aspect ratio in the aSVG file.
#' @param video.dim A single character of the dimension of video frame in form of 'widthxheight', such as '1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480' (default). The aspect ratio of spatial heatmaps are decided by \code{width} and \code{height}. 
#' @param interval The time interval (seconds) between spatial heatmap frames in the video. Default is 1.
#' @param framerate An integer of video framerate in frames per seconds. Default is 1. Larger values make the video smoother. 
#' @param res Resolution of the video in dpi.
#' @return An image of spatial heatmap(s), a three-component list of the spatial heatmap(s) in \code{ggplot} format, a data frame of mapping between assayed samples and aSVG features, and a data frame of feature attributes.

#' @section Details:
#' See the package vignette (\code{browseVignettes('spatialHeatmap')}).  

#' @examples

#' ## In the following examples, the 2 toy data come from an RNA-seq analysis on development of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are
#' ## included in this package. The complete raw count data are downloaded using the R package
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as
#' ## a "data frame" input to exemplify data of simple samples/conditions, while toy data2 as
#' ## "SummarizedExperiment" to illustrate data involving complex samples/conditions.   
#'
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
#' # A targets file describing samples and conditions is required for toy data2. It should be made
#' # based on the experiment design, which is accessible through the accession number 
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
#' df.nor.chk <- norm_data(data=df.chk, norm.fun='CNF', log2.trans=TRUE)
#' se.nor.chk <- norm_data(data=se.chk, norm.fun='CNF', log2.trans=TRUE)

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
#' ## Spatial heatmaps.
#'
#' # The target chicken aSVG is downloaded from the EBI aSVG repository
#' # (https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg) directly with
#' # function "return_feature". It is included in this package and accessed as below. Details on
#' # how this aSVG is selected are documented in function "return_feature".
#' svg.chk <- system.file("extdata/shinyApp/example", "gallus_gallus.svg",
#' package="spatialHeatmap")

#' # Plot spatial heatmaps on gene "ENSGALG00000019846".
#' # Toy data1. 
#' spatial_hm(svg.path=svg.chk, data=df.fil.chk, ID='ENSGALG00000019846', height=0.4,
#' legend.r=1.9, sub.title.size=7, ncol=3)
#' # Save spaital heatmaps as HTML and video files by assigning "out.dir" "~/test". 
#' \donttest{
#' if (!dir.exists('~/test')) dir.create('~/test')
#' spatial_hm(svg.path=svg.chk, data=df.fil.chk, ID='ENSGALG00000019846', height=0.4,
#' legend.r=1.9, sub.title.size=7, ncol=3, out.dir='~/test')
#' }
#' # Toy data2.
#' spatial_hm(svg.path=svg.chk, data=se.fil.chk, ID='ENSGALG00000019846', legend.r=1.9,
#' legend.nrow=2, sub.title.size=7, ncol=3)
#'
#' # The data can also come as as a simple named vector. The following gives an example on a
#' # vector of 3 random values. 
#' # Random values.
#' vec <- sample(1:100, 3)
#' # Name the vector. The last name is assumed as a random sample without a matching feature
#' # in aSVG.
#' names(vec) <- c('brain', 'heart', 'notMapped')
#' vec
#' # Plot.
#' spatial_hm(svg.path=svg.chk, data=vec, ID='geneX', height=0.6, legend.r=1.5, ncol=1)
#'
#' # Plot spatial heatmaps on aSVGs of two Arabidopsis thaliana development stages.
#' 
#' # Make up a random numeric data frame.
#' df.test <- data.frame(matrix(sample(x=1:100, size=50, replace=TRUE), nrow=10))
#' colnames(df.test) <- c('shoot_totalA__condition1', 'shoot_totalA__condition2', 
#' 'shoot_totalB__condition1', 'shoot_totalB__condition2', 'notMapped')
#' rownames(df.test) <- paste0('gene', 1:10) # Assign row names 
#' df.test[1:3, ]

#' # aSVG of development stage 1.
#' svg1 <- system.file("extdata/shinyApp/example", "arabidopsis.thaliana_organ_shm1.svg",
#' package="spatialHeatmap")
#' # aSVG of development stage 2.
#' svg2 <- system.file("extdata/shinyApp/example", "arabidopsis.thaliana_organ_shm2.svg",
#' package="spatialHeatmap")
#' # Spatial heatmaps. 
#' spatial_hm(svg.path=c(svg1, svg2), data=df.test, ID=c('gene1'), height=0.8, legend.r=1.6,
#' preserve.scale=TRUE) 
#'
#' # Multiple development stages can also be arranged in a single aSVG image, but the 
#' # samples, stages, and conditions should be formatted in different ways. See the vignette
#' # for details by running "browseVignette('spatialHeatmap')" in R. 


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1 \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg 
#' \cr Yu, G., 2020. ggplotify:  Convert Plot to ’grob’ or ’ggplot’ Object. R package version 0.0.5.URLhttps://CRAN.R-project.org/package=ggplotify30
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' \cr Guangchuang Yu (2020). ggplotify: Convert Plot to 'grob' or 'ggplot' Object. R package version 0.0.5. https://CRAN.R-project.org/package=ggplotify
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

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

spatial_hm <- function(svg.path, data, sam.factor=NULL, con.factor=NULL, ID, lay.shm="gene", ncol=2, col.com=c('yellow', 'orange', 'red'), col.bar='selected', cores=NA, bar.width=0.08, bar.title.size=0, trans.scale=NULL, ft.trans=NULL, tis.trans=ft.trans, lis.rematch = NULL, legend.r=0.2, sub.title.size=11, legend.plot='all', ft.legend='identical', bar.value.size=10, legend.plot.title='Legend', legend.plot.title.size=11, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, angle.text.key=NULL, position.text.key=NULL, legend.2nd=FALSE, position.2nd='bottom', legend.nrow.2nd=NULL, legend.ncol.2nd=NULL, legend.key.size.2nd=0.03, legend.text.size.2nd=10, angle.text.key.2nd=0, position.text.key.2nd='right', add.feature.2nd=FALSE, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, line.size=0.2, line.color='grey70', relative.scale = NULL, verbose=TRUE, out.dir=NULL, animation.scale = 1, selfcontained=FALSE, video.dim='640x480', res=500, interval=1, framerate=1, legend.value.vdo=NULL, ...) {

  calls <- names(vapply(match.call(), deparse, character(1))[-1])
  if("tis.trans" %in% calls) { ft.trans <- tis.trans; warning('"tis.trans" is deprecated and replaced by "ft.trans"! \n') }
  x <- y <- color_scale <- tissue <- NULL; options(stringsAsFactors=FALSE)
  # if (!is.null(sub.margin)) if (!is.numeric(sub.margin) | length(sub.margin)!=4 | any(sub.margin >= 1) | any(sub.margin < 0)) stop('"sub.margin" must be a 4-length numeric vector between 0 (inclusive) and 1 (exclusive)!')
  # Extract and filter data.
  if (is.vector(data)) {
    vec.na <- make.names(names(data)); if (is.null(vec.na)) stop("Please provide names for the input data!")
    if (!identical(vec.na, names(data))) cat('Syntactically valid column names are made! \n')
    if (any(duplicated(vec.na))) stop('Please make sure data names are unique!')
    form <- grepl("__", vec.na); if (sum(form)==0) { vec.na <- paste0(vec.na, '__', 'con'); con.na <- FALSE } else con.na <- TRUE
    data <- tryCatch({ as.numeric(data) }, warning=function(w) { stop("Please make sure input data are numeric!") }, error=function(e) { stop("Please make sure input data are numeric!") })
    if (is.null(ID)) stop('Please provide a name for the data!')
    gene <- as.data.frame(matrix(data, nrow=1, dimnames=list(ID, vec.na)))

  } else if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'SummarizedExperiment')) {
    id.no <- ID[!ID %in% rownames(data)]
    if (length(id.no)>0) stop(id.no, ': not detected in data! \n')
    # Process data.
    dat.lis <- check_data(data=data, sam.factor=sam.factor, con.factor=con.factor, usage='shm')
    gene <- as.data.frame(dat.lis$dat); con.na <- dat.lis$con.na

  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", or "SummarizedExperiment", except that "spatial_hm" also accepts a "vector".') }

  if (!is.null(trans.scale)) if (trans.scale=='log2') { 
          
      g.min <- min(gene) 
      if (g.min<0) gene <- gene-g.min+1; if (g.min==0) gene <- gene+1; gene <- log2(gene)  

    } else if (trans.scale=='exp2') gene <- 2^gene else if (trans.scale=='row') { 
    gene <- t(scale(t(gene))) } else if (trans.scale=='column') { gene <- scale(gene) }
 
    # Color bar.
    bar.len=1000
    if (col.bar=="all") geneV <- seq(min(gene), max(gene), len=bar.len) else if (col.bar=="selected") geneV <- seq(min(gene[ID, , drop=FALSE]), max(gene[ID, , drop=FALSE]), len=bar.len)
    col <- colorRampPalette(col.com)(length(geneV))
    cs.g <- col_bar(geneV=geneV, cols=col, width=1, bar.title.size=bar.title.size, bar.value.size=bar.value.size)

    # Only take the column names with "__".
    cname <- colnames(gene); form <- grepl('__', cname)
    con <- gsub("(.*)(__)(.*)", "\\3", cname[form]); con.uni <- unique(con)
    sam.uni <- unique(gsub("(.*)(__)(.*)", "\\1", cname))

    # Get SVG names.
    str.lis <- strsplit(svg.path, '/')
    svg.na <- vapply(str.lis, function(x) x[[length(x)]], FUN.VALUE=character(1))
    if (length(svg.na)>1) { 

      shm <- gsub('.*_(shm\\d+).svg$', '\\1', svg.na)
      if (any(duplicated(shm))) stop(paste0('Suffixes of aSVG files are duplicated: ', paste0(shm, collapse='; ')))
      if (!all(grepl('_shm\\d+', svg.na, perl=TRUE))) stop("Suffixes of aSVG files should be indexed as '_shm1', '_shm2', '_shm3', ...") 

    }
    ord <- order(gsub('.*_(shm.*).svg$', '\\1', svg.na))
    svg.path <- svg.path[ord]; svg.na <- svg.na[ord]
    # Coordinates of each SVG are extracted and placed in a list.
    df.attr <- svg.df.lis <- NULL; for (i in seq_along(svg.na)) {
      cat('Coordinates:', svg.na[i], '...', '\n') # '... \n' renders two new lines.
      cores <- deter_core(cores, svg.path[i]); cat('CPU cores:', cores, '\n')
      df_tis <- svg_df(svg.path=svg.path[i], feature=sam.uni, cores=cores)
      if (is.character(df_tis)) stop(paste0(svg.na[i], ': ', df_tis))
      svg.df.lis <- c(svg.df.lis, list(df_tis))
      df.attr0 <- df_tis$df.attr[, c('feature', 'stroke', 'color', 'id', 'element', 'parent', 'index1')]
      df.attr0$SVG <- svg.na[i]; df.attr <- rbind(df.attr, df.attr0)

    }; names(svg.df.lis) <- svg.na; colnames(df.attr)[colnames(df.attr)=='index1'] <- 'order'

    # Extract mapped tissues.
    map.sum <- data.frame(); for (j in svg.na) {
    
      df_tis <- svg.df.lis[[j]] 
      g.df <- df_tis[['df']]; tis.path <- df_tis[['tis.path']]

      not.map <- setdiff(sam.uni, unique(tis.path)); if (verbose==TRUE & length(not.map)>0) cat('Features in data not mapped:', paste0(not.map, collapse=', '), '\n')
      sam.com <- intersect(unique(tis.path), sam.uni) 
      if (length(sam.com)==0) next 

      idx.com <- vapply(cname, function(i) grepl(paste0('^', sam.com, '__', collapse='|'), i), FUN.VALUE=logical(1))
      map.gene <- gene[idx.com]; cna <- colnames(map.gene)
      for (i in ID) {

        df0 <- as.data.frame(t(map.gene[i, ])); colnames(df0) <- 'value'
        featureSVG <- gsub("(.*)(__)(.*)", "\\1", cna)
        rowID <- i; df1 <- data.frame(rowID=rowID, featureSVG=featureSVG)
        if (con.na==TRUE) { condition <- gsub("(.*)(__)(.*)", "\\3", cna); df1 <- cbind(df1, condition=condition) }
        df1 <- cbind(df1, df0); df1$SVG <- j
        map.sum <- rbind(map.sum, df1)

      } 
    
    }; row.names(map.sum) <- NULL

    # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
    w.h.all <- NULL; for (i in seq_along(svg.df.lis)) { w.h.all <- c(w.h.all, svg.df.lis[[i]][['w.h']]); w.h.max <- max(w.h.all) }

    # A set of SHMs (gene_con*) are made for each SVG, and all sets of SHMs under different SVGs are placed in 2 lists in form of ggplots and grobs respectively. Different SHMs of same 'gene_condition' under different SVGs are indexed with suffixed of '_1', '_2', ... E.g. SHMs under shm1.svg: gene_condition1_1, gene_condition2_1; SHMs under shm2.svg: gene_condition1_2, gene_condition2_2; the 4 SHMs are stored in 2 separate lists in form of ggplots and grobs respectively. 
    # The order of ggplot SHMs, grob SHMs, and legends follow the order of SVGs, so all are same.
    gg.lis.all <- gg.all <- grob.all <- lgd.all <- NULL; for (i in seq_along(svg.df.lis)) {
      na0 <- names(svg.df.lis)[i]; cat('ggplots/grobs:', na0, '... \n')
      svg.df <- svg.df.lis[[i]]; g.df <- svg.df[["df"]]; aspect.r <- svg.df[['aspect.r']]
      tis.path <- svg.df[["tis.path"]]; fil.cols <- svg.df[['fil.cols']]
      # if (preserve.scale==TRUE & is.null(sub.margin)) sub.margin <- (1-w.h/w.h.max*0.99)/2
      # SHMs/legend of ggplots under one SVG.
      gg.lis <- gg_shm(gene=gene, con.na=con.na, geneV=geneV, coord=g.df, ID=ID, legend.col=fil.cols, cols=col, tis.path=tis.path, lis.rematch = lis.rematch, ft.trans=ft.trans, sub.title.size=sub.title.size, ft.legend=ft.legend, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.text.size=legend.text.size, legend.plot.title=legend.plot.title, legend.plot.title.size=legend.plot.title.size, line.size=line.size, line.color=line.color, aspect.ratio = aspect.r, ...)
      msg <- paste0(na0, ': no spatial features that have matching sample identifiers in data are detected!'); if (is.null(gg.lis)) stop(msg)

      # Append suffix '_i' for the SHMs of ggplot under SVG[i], and store them in a list.
      ggs <- gg.lis$g.lis.all; names(ggs) <- paste0(names(ggs), '_', i)
      gg.all <- c(gg.all, ggs)

      cores <- deter_core(cores, svg.path[i]); cat('CPU cores:', cores, '\n')
      # Same names with ggplot: append suffix '_i' for the SHMs of grob under each SVG, and store them in a list.
      grobs <- grob_shm(ggs, cores=cores); grob.all <- c(grob.all, grobs)

      # Store SHM ggplots and w.h under each SVG in separate list for use in relative scale.
      lis0 <- list(list(ggs = ggs, w.h = svg.df$w.h)); names(lis0) <- na0
      gg.lis.all <- c(gg.lis.all, lis0)

      # Store legend of ggplot in a list.
      lgd.all <- c(lgd.all, list(gg.lis$g.lgd))
    }; names(lgd.all) <- names(svg.df.lis)

    pat.gen <- paste0(ID, collapse='|'); pat.con <- paste0(con.uni, collapse='|')
    # Use definite patterns and avoid using '.*' as much as possible. Try to as specific as possible.
    pat.all <- paste0('^(', pat.gen, ')_(', pat.con, ')(_\\d+$)')
    # Indexed cons with '_1', '_2', ... at the end.
    con.idxed <- unique(gsub(pat.all, '\\2\\3', names(gg.all)))
    # Layout matrix of SHMs for use in relative scale.
    lay.mat <- lay_shm(lay.shm=lay.shm, con=con.idxed, ncol=ncol, ID.sel=ID, lay.mat = TRUE) 

    # If relative size in multiple SVGs is required, update all SHM ggplots/grobs under each SVG.
    if (is.numeric(relative.scale) & length(svg.path) > 1) if (relative.scale > 0) {
      cat('Applying relative image size ... \n')
      gg.all <- grob.all <- NULL; for (i in seq_along(gg.lis.all)) {
        lis0 <- gg.lis.all[[i]]
        ggs <- rela_size(lis0$w.h['height'], w.h.max, relative.scale, nrow(lay.mat), lis0$ggs)
        gg.all <- c(gg.all, ggs)
        cores <- deter_core(cores, svg.path[i]) # gg.lis.all has the same names with svg.df.lis.
        # Same names with ggplot: append suffix '_i' for the SHMs of grob under each SVG, and store them in a list.
        grobs <- grob_shm(ggs, cores=cores); grob.all <- c(grob.all, grobs)
        
      }
    }; na.all <- names(grob.all)

    # Arrange color scale grob.
    tmp <- normalizePath(tempfile(), winslash='/', mustWork=FALSE)
    png(tmp); cs.grob <- ggplotGrob(cs.g); dev.off()
    if (file.exists(tmp)) do.call(file.remove, list(tmp))
    cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc")) # "mm" is fixed, "npc" is scalable.

    # Arrange SHM grobs.
    na.all <- sort_gen_con(ID.sel=ID, na.all=na.all, con.all=con.idxed, by=lay.shm)
    grob.all <- grob.all[na.all]; gg.all <- gg.all[na.all]
    g.arr <- lay_shm(lay.shm=lay.shm, con=con.idxed, ncol=ncol, ID.sel=ID, grob.list=grob.all, shiny=FALSE)

    if (!is.null(legend.plot)) {
      # Select legend plot to show. 
      if (length(svg.na)>1 & legend.plot!='all') na.lgd <- svg.na[grep(paste0('_(', paste0(legend.plot, collapse='|'), ').svg$'), svg.na)] else na.lgd <- svg.na
      lgd.lis <- lgd.all[na.lgd]
     
      if (legend.2nd==TRUE) {
        gg.all <- gg_2lgd(gg.all=gg.all, sam.dat=sam.uni, ft.trans=ft.trans, position.2nd=position.2nd, legend.nrow.2nd=legend.nrow.2nd, legend.ncol.2nd=legend.ncol.2nd, legend.key.size.2nd=legend.key.size.2nd, add.feature.2nd=add.feature.2nd, legend.text.size.2nd=legend.text.size.2nd, angle.text.key.2nd=angle.text.key.2nd, position.text.key.2nd=position.text.key.2nd)
        grob.all <- lapply(gg.all, ggplotGrob)
      }
      
      # Add labels to target shapes/adjust keys in legend plots.
      lgd.lis <- gg_lgd(gg.all=lgd.lis, angle.text.key=angle.text.key, position.text.key=position.text.key, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, sam.dat=sam.uni, ft.trans=ft.trans)
      grob.lgd.lis <- lapply(lgd.lis, ggplotGrob)
      lgd.tr <- lapply(grob.lgd.lis, grobTree)

      # Number of all rows in SHMs.
      # if (lay.shm=='gene') row.all <- ceiling(length(con.uni)/ncol)*length(ID)
      # if (lay.shm=='con') row.all <- ceiling(length(ID)/ncol)*length(con.uni)
      # In 'arrangeGrob', if numbers in 'layout_matrix' are more than items in 'grobs', there is no difference. The width/height of each subplot is decided by 'widths' and 'heights'.
      #lgd.w <- 0.99; lgd.h <- 0.99/length(na.lgd)/legend.r
      # if (lgd.h*length(na.lgd)>1) { lgd.h <- 0.99/length(na.lgd); lgd.w <- lgd.h*legend.r }

      w.lgd <- (1-bar.width)/(ncol+1); shm.w <- 1-bar.width-w.lgd
      # If legend.r = 0, legend plot size is a square.
      lgd.arr <- arrangeGrob(grobs=lgd.tr, layout_matrix=matrix(seq_along(na.lgd), ncol=1), widths=unit(0.99, "npc"), heights=unit(rep(w.lgd + (0.99 - w.lgd) * legend.r, length(na.lgd)), "npc"))

      # A plot pops up when 'grid.arrange' runs.
      shm <- grid.arrange(cs.arr, g.arr, lgd.arr, ncol=3, widths=unit(c(bar.width-0.005, shm.w, w.lgd), 'npc'))
    } else { shm.w <- 1-bar.width; shm <- grid.arrange(cs.arr, g.arr, ncol=2, widths=unit(c(bar.width-0.005, shm.w), 'npc')) }

    if (!is.null(out.dir)) { 

      for (i in names(gg.all)) {
        g0 <- gg.all[[i]]
        anm.height <- 550 * animation.scale; anm.width <- anm.height / g0$theme$aspect.ratio
        html_ly(gg=gg.all[i], cs.g=cs.g, ft.trans=ft.trans, sam.uni=sam.uni, anm.width=anm.width, anm.height=anm.height, out.dir=out.dir)
      }

      vdo <- video(gg=gg.all, cs.g=cs.g, sam.uni=sam.uni, ft.trans=ft.trans, lgd.key.size=legend.key.size.2nd, lgd.text.size=legend.text.size.2nd, angle.text.key=angle.text.key.2nd, position.text.key=position.text.key.2nd, lgd.row=legend.nrow.2nd, lgd.col=legend.ncol.2nd, legend.value.vdo=legend.value.vdo, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, video.dim=video.dim, res=res, interval=interval, framerate=framerate, out.dir=out.dir)
      if (is.null(vdo)) cat("Video is not generated! \n")

    }
    lis <- list(spatial_heatmap=as.ggplot(shm), mapped_feature=map.sum, feature_attribute=df.attr); invisible(lis)

}

#' Determine the number of CPU cores.
#' @param cores The input number of cores.
#' @param svg.path The path of SVG file.
#' @keywords Internal
#' @noRd
#' @importFrom parallel detectCores

deter_core <- function(cores, svg.path) {
  cores <- as.integer(cores); n.cor <- detectCores(logical=TRUE); fil.size <- file.size(svg.path)
  if (!is.na(n.cor)) { 
    if (fil.size >= 3145728 & n.cor > 2 & is.na(cores)) cores <- 2 else if (fil.size < 3145728 & is.na(cores)) cores <- 1 else if (n.cor > 1 & !is.na(cores)) {
      if (cores >= n.cor) cores <- n.cor - 1
    } 
  } else { if (is.na(cores)) cores <- 1 }; return(cores)
}

#' Adjust relative SHM size for multiple aSVGs
#' @param h The image height in target aSVG file.
#' @param h.max The max image height in all aSVG files.
#' @param scale The coefficient to adjust the relative size.
#' @param lay.row The number of rows in the final SHM layout, returned by \code{lay_shm}.
#' @param gg.lis The list of all SHM ggplots.
#' @keywords Internal
#' @noRd
#' @importFrom ggplot2 theme

rela_size <- function(h, h.max, scale, lay.row, gg.lis) {
  mar.tb <- ((1 - h / h.max) * (0.99 / lay.row)) / 2 * scale + 0.005
  gg.lis <- lapply(gg.lis, function(x) { x + theme(plot.margin = margin(mar.tb, 0.005, mar.tb, 0.005, "npc")) })
  return(gg.lis)
}




