#' Plot Spatial Heatmaps
#'
#' @description 
#' The input are a pair of annotated SVG (aSVG) file and formatted numeric data (\code{vector}, \code{data.frame}, \code{SummarizedExperiment}). In the aSVG, spatial features (tissues, organs, etc) are represented by shapes and assigned unique identifiers, and the numeric data are measured from these spatial features. In biological applications, aSVGs are anatomical structures, and numeric data are expression values of genes, proteins, metabolites, etc assayed from spatial features. Numeric data are mapped to spatial features in aSVGs according to the same identifiers between the two. The mapped features are filled with colors translated from the numeric data while other features are transparent. The output images are termed spatial heatmaps (SHMs). 
#'
#' This function is designed multi-functional for maximal flexibility. For example, SHMs can be organized horizontally by each gene or each experimental variable to facilitate comparisons. In multi-layer anotomical structures, selected tissues can be set transparent to expose burried features in lower layers. The color scale is customizable to highlight difference across features. This function also works with many other types of spatial data, such as population data plotted to geographic maps.

#' @inheritParams covis

#' @return An `SPHM` object.

#' @section Details:
#' See the package vignette (\code{browseVignettes('spatialHeatmap')}). 

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
#' # The chicken aSVG downloaded from the EBI aSVG repository (https://github.com/ebi-gene-
#' # expression-group/anatomogram/tree/master/src/svg) is included in this package and 
#' # accessed as below.
#' svg.chk <- system.file("extdata/shinyApp/data", "gallus_gallus.svg",
#' package="spatialHeatmap")
#' # Read the chicken aSVG file.
#' svg.chk <- read_svg(svg.path=svg.chk)
#'
#' # Store assay data and aSVG in an "SHM" class.
#' dat.chk <- SPHM(svg=svg.chk, bulk=se.chk.fil)
#' # Plot spatial heatmaps with gene "ENSGALG00000019846".
#' shm(data=dat.chk, ID='ENSGALG00000019846', legend.r=1.9,
#' legend.nrow=5, sub.title.size=7, ncol=3)
#'
#' # Save SHMs as HTML and video files in the "~/test" directory. 
#' \donttest{
#' if (!dir.exists('~/test')) dir.create('~/test')
#' shm(data=dat.chk, ID='ENSGALG00000019846', legend.r=1.9,
#' legend.nrow=5, sub.title.size=7, ncol=3, out.dir='~/test')
#' }
#'
#' @inherit covis author references

#' @name shm
#' @rdname shm
#' @aliases shm,SPHM-method
#' @export

setMethod("shm", c(data="SPHM"), function(data, assay.na=NULL, sam.factor=NULL, con.factor=NULL, ID, charcoal=FALSE, alpha.overlay=1, lay.shm="gene", ncol=2, h=0.99, col.com=c('yellow', 'orange', 'red'), col.bar='selected', thr=c(NA, NA), cores=NA, bar.width=0.08, bar.title=NULL, bar.title.size=0, scale=NULL, ft.trans=NULL, tis.trans=ft.trans, legend.r=0.9, sub.title.size=11, sub.title.vjust=2, legend.plot='all', ft.legend='identical', bar.value.size=10, legend.plot.title='Legend', legend.plot.title.size=11, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, angle.text.key=NULL, position.text.key=NULL, legend.2nd=FALSE, position.2nd='bottom', legend.nrow.2nd=NULL, legend.ncol.2nd=NULL, legend.key.size.2nd=0.03, legend.text.size.2nd=10, angle.text.key.2nd=0, position.text.key.2nd='right', add.feature.2nd=FALSE, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, line.width=0.2, line.color='grey70', relative.scale = NULL, verbose=TRUE, out.dir=NULL, animation.scale = 1, selfcontained=FALSE, video.dim='640x480', res=500, interval=1, framerate=1, bar.width.vdo=0.1, legend.value.vdo=NULL, ...) {
  calls <- names(vapply(match.call(), deparse, character(1))[-1])
  if("tis.trans" %in% calls) warning('"tis.trans" is deprecated and replaced by "ft.trans"! \n')
  if("svg" %in% calls) warning('"svg" is deprecated and replaced by "data"! \n')
  if("svg.path" %in% calls) warning('"svg.path" is deprecated and replaced by "data"! \n')
  # scale=scale: if an argument name is the same with a function name, it may cause errors.
  res <- shm_covis(svg=data@svg, data=data@bulk, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, ID=ID, charcoal=charcoal, alpha.overlay=alpha.overlay, lay.shm=lay.shm, ncol=ncol, h=h, col.com=col.com, col.bar=col.bar, thr=thr, cores=cores, bar.width=bar.width, bar.title=bar.title, bar.title.size=bar.title.size, scale=scale, ft.trans=ft.trans, lis.rematch = data@match, legend.r=legend.r, sub.title.size=sub.title.size, sub.title.vjust=sub.title.vjust, legend.plot=legend.plot, ft.legend=ft.legend, bar.value.size=bar.value.size, legend.plot.title=legend.plot.title, legend.plot.title.size=legend.plot.title.size, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.text.size=legend.text.size, angle.text.key=angle.text.key, position.text.key=position.text.key, legend.2nd=legend.2nd, position.2nd=position.2nd, legend.nrow.2nd=legend.nrow.2nd, legend.ncol.2nd=legend.ncol.2nd, legend.key.size.2nd=legend.key.size.2nd, legend.text.size.2nd=legend.text.size.2nd, angle.text.key.2nd=angle.text.key.2nd, position.text.key.2nd=position.text.key.2nd, add.feature.2nd=add.feature.2nd, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, line.width=line.width, line.color=line.color, relative.scale = relative.scale, verbose=verbose, out.dir=out.dir, animation.scale = animation.scale, selfcontained=selfcontained, video.dim=video.dim, res=res, interval=interval, framerate=framerate, bar.width.vdo=bar.width.vdo, legend.value.vdo=legend.value.vdo, ...)
  data@output <- res; invisible(data)
})


 

