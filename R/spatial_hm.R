#' Plot Spatial Heatmaps
#'
#' The input are a pair of annotated SVG (aSVG) file and formatted data (\code{vector}, \code{data.frame}, \code{SummarizedExperiment}). In the former, spatial features are represented by shapes and assigned unique identifiers, while the latter are numeric values measured from these spatial features and organized in specific formats. In biological cases, aSVGs are anatomical or cell structures, and data are measurements of genes, proteins, metabolites, \emph{etc}. in different samples (\emph{e.g.} cells, tissues). Data are mapped to the aSVG according to identifiers of assay samples and aSVG features. Only the data from samples having matching counterparts in aSVG features are mapped. The mapped features are filled with colors translated from the data, and the resulting images are termed spatial heatmaps. Note, "sample" and "feature" are two equivalent terms referring to cells, tissues, organs \emph{etc.} where numeric values are measured. Matching means a target sample in data and a target spatial feature in aSVG have the same identifier. \cr This function is designed as much flexible as to achieve optimal visualization. For example, subplots of spatial heatmaps can be organized by gene or condition for easy comparison, in multi-layer anotomical structures selected tissues can be set transparent to expose burried features, color scale is customizable to highlight difference among features. This function also works with many other types of spatial data, such as population data plotted to geographic maps.


#' @inheritParams covis
#' @param svg An \code{SVG} object containing one or multiple aSVG instances (see \code{\link{SVG}} and \code{\link{read_svg}}). In the aSVGs, spatial features (tissues, organs, etc) having counterparts with the same identifiers in the `bulk` data will be colored accoording to expression profiles of chosen biomolecules (genes, proteins, etc). 
#' @param lis.rematch The \code{list} for re-matching in SHMs (only bulk data) or matching between single-cell and bulk data in co-visualization.
#' \describe{ 
#'  \item{SHMs}{ A named \code{list} for rematching spatial features between numeric data (ftA, ftB) and aSVGs (ftC, ftD, ftE). In each slot, the slot name is a spatial feature from the data and the corresponding element is one or multiple spatial features from the aSVG. \emph{E.g.} \code{list(ftA = c('ftC', 'ftD'), ftB = c('ftE'))}.
#'  } 
#'  \item{Co-visualization plots}{ Mapping cells to tissues: a named \code{list}, where cell group labels from \code{colData(sce.  dimred)[, 'cell.group']} are the name slots and aSVG features are the corresponding \code{list} elements. Mapping tissues to cells: a named \code{list}, where tissues are the name slots and cells from \code{colData(sce.dimred)[, 'cell.group']} are the corresponding \code{list} elements. Applicable when cell grouping methods are annodation labels, marker genes, clustering, or manual assignments. 
#'  }
#' }
#' @return An image of spatial heatmap(s), a two-component list of the spatial heatmap(s) in \code{ggplot} format and a \code{data.frame} of mapping between assayed samples and aSVG features.

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
#' cnt.chk.simple <- system.file('extdata/shinyApp/data/count_chicken_simple.txt',
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
#' cnt.chk <- system.file('extdata/shinyApp/data/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing samples and conditions is required for toy data2. It should be made
#' # based on the experiment design, which is accessible through the accession number 
#' # "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in this
#' # package and accessed below. 

#' # Access the example targets file. 
#' tar.chk <- system.file('extdata/shinyApp/data/target_chicken.txt', package='spatialHeatmap')
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
#' df.fil.chk <- filter_data(data=df.aggr.chk, pOA=c(0.01, 5), CV=c(0.2, 100))
#' # Filter toy data2.
#' se.fil.chk <- filter_data(data=se.aggr.chk, sam.factor='organism_part', con.factor='age',
#' pOA=c(0.01, 5), CV=c(0.2, 100))
#'
#' ## Spatial heatmaps.
#'
#' # The target chicken aSVG is downloaded from the EBI aSVG repository
#' # (https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg) directly with
#' # function "return_feature". It is included in this package and accessed as below. Details on
#' # how this aSVG is selected are documented in function "return_feature".
#' svg.chk <- system.file("extdata/shinyApp/data", "gallus_gallus.svg",
#' package="spatialHeatmap")
#'
#' # Reading the chicken aSVG file.
#' svg.chk <- read_svg(svg.path=svg.chk)
#' 
#' # Plot spatial heatmaps on gene "ENSGALG00000019846".
#' # Toy data1. 
#' spatial_hm(svg=svg.chk, data=df.fil.chk, ID='ENSGALG00000019846', height=0.4,
#' legend.r=1.9, sub.title.size=7, ncol=3)
#' # Save spaital heatmaps as HTML and video files by assigning "out.dir" "~/test". 
#' \donttest{
#' if (!dir.exists('~/test')) dir.create('~/test')
#' spatial_hm(svg=svg.chk, data=df.fil.chk, ID='ENSGALG00000019846', height=0.4,
#' legend.r=1.9, sub.title.size=7, ncol=3, out.dir='~/test')
#' }
#' # Toy data2.
#' spatial_hm(svg=svg.chk, data=se.fil.chk, ID='ENSGALG00000019846', legend.r=1.9,
#' legend.nrow=2, sub.title.size=7, ncol=3)
#'

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1 \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg 
#' \cr Yu, G., 2020. ggplotify:  Convert Plot to ’grob’ or ’ggplot’ Object. R package version 0.0.5.URLhttps://CRAN.R-project.org/package=ggplotify30
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' \cr Guangchuang Yu (2020). ggplotify: Convert Plot to 'grob' or 'ggplot' Object. R package version 0.0.5. https://CRAN.R-project.org/package=ggplotify
#' \cr Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' Marques A et al. (2016). Oligodendrocyte heterogeneity in the mouse juvenile and adult central nervous system. Science 352(6291), 1326-1329.
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.

#' @name spatial_hm
#' @rdname spatial_hm
#' @aliases spatial_hm,SVG-method
#' @export

setMethod("spatial_hm", c(svg="SVG"), function(svg, data, assay.na=NULL, sam.factor=NULL, con.factor=NULL, ID, charcoal=FALSE, alpha.overlay=1, lay.shm="gene", ncol=2, col.com=c('yellow', 'orange', 'red'), col.bar='selected', thr=c(NA, NA), cores=NA, bar.width=0.08, bar.title=NULL, bar.title.size=0, scale=NULL, ft.trans=NULL, tis.trans=ft.trans, lis.rematch = NULL, legend.r=0.9, sub.title.size=11, sub.title.vjust=2, legend.plot='all', ft.legend='identical', bar.value.size=10, legend.plot.title='Legend', legend.plot.title.size=11, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, angle.text.key=NULL, position.text.key=NULL, legend.2nd=FALSE, position.2nd='bottom', legend.nrow.2nd=NULL, legend.ncol.2nd=NULL, legend.key.size.2nd=0.03, legend.text.size.2nd=10, angle.text.key.2nd=0, position.text.key.2nd='right', add.feature.2nd=FALSE, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, line.width=0.2, line.color='grey70', relative.scale = NULL, verbose=TRUE, out.dir=NULL, animation.scale = 1, selfcontained=FALSE, video.dim='640x480', res=500, interval=1, framerate=1, bar.width.vdo=0.1, legend.value.vdo=NULL, ...) {
  # if ('spatial_hm' %in% as.character(match.call()[[1]])) warning('The function "spatial_hm" is deprecated and replaced by "shm"!', call. = FALSE)
  # if ('spatial_hm' %in% deparse(sys.call(0)[[1]])) warning('The function "spatial_hm" is deprecated and replaced by "shm"!', call. = FALSE)
  warning('The function "spatial_hm" will be deprecated and replaced by "shm"!', call. = FALSE)
  calls <- names(vapply(match.call(), deparse, character(1))[-1])
  if("tis.trans" %in% calls) warning('"tis.trans" is deprecated and replaced by "ft.trans"! \n')
  if("svg.path" %in% calls) warning('"svg.path" is deprecated and replaced by "svg"! \n')

  res <- shm_covis(svg=svg, data=data, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, ID=ID, charcoal=charcoal, alpha.overlay=alpha.overlay, lay.shm=lay.shm, ncol=ncol, col.com=col.com, col.bar=col.bar, thr=thr, cores=cores, bar.width=bar.width, bar.title=bar.title, bar.title.size=bar.title.size, scale=scale, ft.trans=ft.trans, lis.rematch = lis.rematch, legend.r=legend.r, sub.title.size=sub.title.size, sub.title.vjust=sub.title.vjust, legend.plot=legend.plot, ft.legend=ft.legend, bar.value.size=bar.value.size, legend.plot.title=legend.plot.title, legend.plot.title.size=legend.plot.title.size, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.text.size=legend.text.size, angle.text.key=angle.text.key, position.text.key=position.text.key, legend.2nd=legend.2nd, position.2nd=position.2nd, legend.nrow.2nd=legend.nrow.2nd, legend.ncol.2nd=legend.ncol.2nd, legend.key.size.2nd=legend.key.size.2nd, legend.text.size.2nd=legend.text.size.2nd, angle.text.key.2nd=angle.text.key.2nd, position.text.key.2nd=position.text.key.2nd, add.feature.2nd=add.feature.2nd, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, line.width=line.width, line.color=line.color, relative.scale = relative.scale, verbose=verbose, out.dir=out.dir, animation.scale = animation.scale, selfcontained=selfcontained, video.dim=video.dim, res=res, interval=interval, framerate=framerate, bar.width.vdo=bar.width.vdo, legend.value.vdo=legend.value.vdo, ...)
  invisible(res)
})


 

