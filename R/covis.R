#' Co-visualizing spatial heatmaps with single cell embedding plots
#'
#' This function is an extension of \code{spatial_hm}. It integrates visualization of single cells and related bulk tissues in a composite plot. The former is in form of an embedding plot (PCA, UMAP, TSNE) while the latter is a spatial heatmap plot returned by \code{spatial_hm}. 

#' @param svg An object of \code{coord} containing one or multiple aSVG instances. See \code{\link{read_svg}} for two to store aSVG files in \code{coord}. 
#' @param charcoal Logical, if \code{TRUE} the raster image will be turned black and white.
#' @param alpha.overlay The opacity of top-layer spatial heatmaps if a raster image is added at the bottom layer. The default is 1.
#' @param sce.dimred A \code{SingleCellExperiment} with reduced dimentionalities such as \code{PCA}, \code{UMAP}, \code{TSNE}.
#' @param dimred One of \code{PCA}, \code{UMAP}, \code{TSNE} in \code{sce.dimred}, specifying which reduced dimensionality to use in co-visualization of bulk tissues and single cells.
#' @param cell.group Applicable in co-visualizing bulk tissues and single cells with annodation-based or manual method. A column name in \code{colData(sce.dimred)}, where one label defines a cell group, and the mapping direction is from cell groups/labels to bulk tissues. 
#' @param tar.cell Applicable in co-visualizing bulk tissues and single cells through annodation-based or manual method. A vector of target cell labels in \code{cell.group}. Cells corresponding to these labels are mapped to bulk tissues through \code{lis.rematch}.  
#' @param tar.bulk A vector of target bulk tissues, which are mapped to single cells through \code{lis.rematch}.
#' @param profile Logical, applicable in co-visualizing bulk tissues and single cells. If \code{TRUE}, one or multiple biological molecule (e.g. gene) identifiers need to be assigned to \code{ID}, and their abundance proifles are included in the co-visualization plot. If \code{FALSE}, abundance proifles are excluded.  
#' @param bar.title.size A numeric of color bar title size. The default is 0. 
#' @param bar.value.size A numeric of value size in the color bar y-axis. The default is 10.

#' @inheritParams filter_data
#' @inheritParams htmlwidgets::saveWidget
#' @inheritParams ggplot2::theme

#' @param position.2nd The position of the secondary legend. One of "top", "right", "bottom", "left", or a two-component numeric vector. The default is "bottom". Applies to the static image and video.
#' @param legend.nrow.2nd An integer of rows of the secondary legend keys. Applies to the static image and video.
#' @param legend.ncol.2nd An integer of columns of the secondary legend keys. Applies to the static image and video.
#' @param legend.key.size.2nd A numeric of legend key size. The default is 0.03. Applies to the static image and video.
#' @param dim.lgd.pos The legend position in dimensionality reduction plot. The default is \code{bottom}.
#' @param dim.lgd.nrow The number of legend rows in dimensionality reduction plot. The default is \code{2}.
#' @param dim.lgd.text.size The size of legend text in dimensionality reduction plot. The default is \code{8}.
#' @param add.feature.2nd Logical TRUE or FALSE. Add feature identifiers to the secondary legend or not. The default is FALSE. Applies to the static image.
#' @param legend.text.size.2nd A numeric of the secondary legend text size. The default is 10. Applies to the static image and video.
#' @param angle.text.key.2nd A value of angle of key text in the secondary legend. Default is 0. Applies to the static image and video.
#' @param position.text.key.2nd The position of key text in the secondary legend, one of "top", "right", "bottom", "left". Default is "right". Applies to the static image and video.
#' @param dim.lgd.pos The legend position in the dimensionality reduction plot. The default is \code{bottom}.
#' @param dim.lgd.nrow The number of legend rows in the dimensionality reduction plot. The default is \code{2}.
#' @param dim.lgd.key.size The size of legend key in the dimensionality reduction plot. The default is \code{4}.
#' @param dim.lgd.text.size The size of legend text in the dimensionality reduction plot. The default is \code{13}.
#' @param dim.capt.size The size of caption text in the dimensionality reduction plot in coclustering. The default is \code{13}.
#' @param angle.text.key A value of key text angle in legend plot. The default is NULL, equivalent to 0.
#' @param position.text.key The position of key text in legend plot, one of "top", "right", "bottom", "left". Default is NULL, equivalent to "right".
#' @param bar.width.vdo The color bar width in video, between 0 and 1.
#' @param legend.value.vdo Logical TRUE or FALSE. If TRUE, the numeric values of matching spatial features are added to video legend. The default is NULL.
#' @param label Logical. If TRUE, spatial features having matching samples are labeled by feature identifiers. The default is FALSE. It is useful when spatial features are labeled by similar colors. 
#' @param label.size The size of spatial feature labels in legend plot. The default is 4.
#' @param label.angle The angle of spatial feature labels in legend plot. Default is 0.
#' @param hjust The value to horizontally adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param vjust The value to vertically adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param opacity The transparency of colored spatial features in legend plot. Default is 1. If 0, features are totally transparent.
#' @param key Logical. The default is TRUE and keys are added in legend plot. If \code{label} is TRUE, the keys could be removed. 

#' @param ft.trans A character vector of tissue/spatial feature identifiers that will be set transparent. \emph{E.g} c("brain", "heart"). This argument is used when target features are covered by  overlapping features and the latter should be transparent.
#' @param lis.rematch \describe{ \item{(1) Spatial heatmap plots of only bulk tissues without single cells.}{ A named \code{list} for rematching between tissues in data (tissue1Data, tissue2Data) and aSVG spatial features (feature1SVG, feature2SVG, feature3SVG). In each slot, the slot name is an tissue identifier in the data and the slot contains one or multiple aSVG features in a vector. \emph{E.g.} \code{list(tissue1Data = c('feature1SVG', 'feature2SVG'), tissue2Data = c('feature3SVG'))}.} \item{(2) Co-visualizing bulk tissues and single cells using annotation-based or manual methods.}{ Mapping cells to bulk tissues: a named \code{list}, where cell labels from \code{colData(sce.dimred)[, 'cell.group']} are the name slots and aSVG features are the corresponding \code{list} elements. Mapping bulk tissues to cells: a named \code{list}, where bulk tissues are the name slots and cells from \code{colData(sce.dimred)[, 'cell.group']} are the corresponding \code{list} elements.}}
#' @param tis.trans This argument is deprecated and replaced by \code{ft.trans}. 
#' @param sub.title.size A numeric of the subtitle font size of each individual spatial heatmap. The default is 11.
#' @param sub.title.vjust A numeric of vertical adjustment for subtitle. The default is \code{2}.
#' @param ft.legend One of "identical", "all", or a character vector of tissue/spatial feature identifiers from the aSVG file. The default is "identical" and all the identical/matching tissues/spatial features between the data and aSVG file are colored in the legend plot. If "all", all tissues/spatial features in the aSVG are shown. If a vector, only the tissues/spatial features in the vector are shown.
#' @param legend.ncol An integer of the total columns of keys in the legend plot. The default is NULL. If both \code{legend.ncol} and \code{legend.nrow} are used, the product of the two arguments should be equal or larger than the total number of shown spatial features.
#' @param legend.nrow An integer of the total rows of keys in the legend plot. The default is NULL. It is only applicable to the legend plot. If both \code{legend.ncol} and \code{legend.nrow} are used, the product of the two arguments should be equal or larger than the total number of matching spatial features.
#' @param legend.key.size A numeric of the legend key size ("npc"), applicable to the legend plot. The default is 0.02. 
#' @param legend.text.size A numeric of the legend label size, applicable to the legend plot. The default is 12.
#' @param line.width The thickness of each shape outline in the aSVG is maintained in spatial heatmaps, \emph{i.e.} the stroke widths in Inkscape. This argument is the extra thickness added to all outlines. Default is 0.2 in case stroke widths in the aSVG are 0.
#' @param line.color A character of the shape outline color. Default is "grey70".
#' @param legend.plot.title The title of the legend plot. The default is 'Legend'. 
#' @param legend.plot.title.size The title size of the legend plot. The default is 11.

#' @param ID A character vector of assyed items (\emph{e.g.} genes, proteins) whose abudance values are used to color the aSVG.
#' @param col.com A character vector of the color components used to build the color scale. The default is c('yellow', 'orange', 'red').
#' @param col.bar One of "selected" or "all", the former uses values of \code{ID} to build the color scale while the latter uses all values from the data. The default is "selected".
#' @param sig.thr A two-numeric vector of the signal thresholds (the range of the color bar). The first and the second element will be the minmun and maximum threshold in the color bar respectively. Signals/values above the max or below min will be assigned the same color as the max or min respectively. The default is \code{c(NA, NA)} and the min and max signals in the data will be used. If one needs to change only max or min, the other should be \code{NA}.    
#' @param cores The number of CPU cores for parallelization, relevant for aSVG files with size larger than 5M. The default is NA, and the number of used cores is 1 or 2 depending on the availability. 
#' @param trans.scale One of "log2", "exp2", "row", "column", or NULL, which means transform the data by "log2" or "2-base expoent", scale by "row" or "column", or no manipuation respectively. This argument should be used if colors across samples cannot be distinguished due to low variance or outliers. 
#' @param bar.width The width of color bar that ranges from 0 to 1. The default is 0.08.

#' @param legend.plot A vector of suffix(es) of aSVG file name(s) such as \code{c('shm1', 'shm2')}. Only aSVG(s) whose suffix(es) are assigned to this arugment will have a legend plot on the right. The default is \code{all} and each aSVG will have a legend plot. If NULL, no legend plot is shown.
#' @param legend.r A numeric (between -1 and 1) to adjust the legend plot size. The default is 0.9.
#' @param legend.2nd Logical, TRUE or FALSE. If TRUE, the secondary legend is added to each spatial heatmap, which are the numeric values of each matching spatial features. The default its FALSE. Only applies to the static image.
#' @param lay.shm One of "gene", "con", or "none". If "gene", spatial heatmaps are organized by genes proteins, or metabolites, \emph{etc.} and conditions are sorted whithin each gene. If "con", spatial heatmaps are organized by the conditions/treatments applied to experiments, and genes are sorted winthin each condition. If "none", spaital heatmaps are organized by the gene order in \code{ID} and conditions follow the order they appear in \code{data}. 
#' @param ncol An integer of the number of columns to display the spatial heatmaps, which does not include the legend plot.

#' @param relative.scale A numeric to adjust the relative sizes between multiple aSVGs. Applicable only if multiple aSVG paths is assigned to \code{svg}. Default is \code{NULL} and all aSVGs have the same size. 
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

#' # Co-visualizing single cells and bulk tissues by mapping data from cells to bulk tissues through manual matching.
#' 
#' # To obtain for examples with randomized data or parameters always the same results, a fixed seed is set.
#' set.seed(10); library(SummarizedExperiment)
#' # Read example single cell data of mouse brain (Marques et al. (2016)).
#' sce.pa <- system.file("extdata/shinyApp/example", "cell_mouse_brain.rds", package="spatialHeatmap")
#' sce <- readRDS(sce.pa)
#' # Pre-process single cell data.
#' sce.dimred.quick <- process_cell_meta(sce, qc.metric=list(subsets=list(Mt=rowData(sce)$featureType=='mito'), threshold=1))
#' colData(sce.dimred.quick)[1:3, 1:2]
#' # Expression values are aggregated by taking means across cell group labels, which are stored in the "colData" slot of sce.dimred.quick.
#' sce.aggr.quick <- aggr_rep(sce.dimred.quick, assay.na='logcounts', sam.factor='label', aggr='mean')
#' # Read the aSVG image into an "SVG" object.
#' svg.mus.brain.pa <- system.file("extdata/shinyApp/example", "mus_musculus.brain.svg", package="spatialHeatmap")
#' svg.mus.brain <- read_svg(svg.mus.brain.pa)
#' tail(attribute(svg.mus.brain)[[1]])[, 1:4]
#' # To map cell group labels to aSVG features, a list with named components is used, where cell labels are in name slots. Note, each cell label can be matched to multiple aSVG features but not vice versa.
#' lis.match.quick <- list(hypothalamus=c('hypothalamus'), cortex.S1=c('cerebral.cortex', 'nose'))
#' # Co-visualization is created on expression values of gene "Eif5b". The target cell groups are "hypothalamus" and "cortex.S1".
#' covis(svg=svg.mus.brain, data=sce.aggr.quick, ID=c('Eif5b'), sce.dimred=sce.dimred.quick, dimred='PCA', cell.group='label', tar.cell=names(lis.match.quick), lis.rematch=lis.match.quick, assay.na='logcounts', bar.width=0.11, dim.lgd.nrow=1, height=0.7, legend.r=1.5, legend.key.size=0.02, legend.text.size=12, legend.nrow=3) 
#' 
#' # More examples of co-visualization are provided in the package vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/spatialHeatmap/inst/doc/covisualize.html


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

#' @name covis
#' @rdname covis
#' @docType methods
#' @export

setMethod("covis", c(svg="SVG"), function(svg, data, assay.na=NULL, sam.factor=NULL, con.factor=NULL, ID, sce.dimred=NULL, dimred='PCA', cell.group=NULL, tar.cell=NULL, tar.bulk=NULL, profile=TRUE, charcoal=FALSE, alpha.overlay=1, lay.shm="gene", ncol=2, col.com=c('yellow', 'orange', 'red'), col.bar='selected', sig.thr=c(NA, NA), cores=NA, bar.width=0.08, bar.title.size=0, trans.scale=NULL, ft.trans=NULL, tis.trans=ft.trans, lis.rematch = NULL, legend.r=0.9, sub.title.size=11, sub.title.vjust=2, legend.plot='all', ft.legend='identical', bar.value.size=10, legend.plot.title='Legend', legend.plot.title.size=11, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, angle.text.key=NULL, position.text.key=NULL, legend.2nd=FALSE, position.2nd='bottom', legend.nrow.2nd=NULL, legend.ncol.2nd=NULL, legend.key.size.2nd=0.03, legend.text.size.2nd=10, angle.text.key.2nd=0, position.text.key.2nd='right', dim.lgd.pos='bottom', dim.lgd.nrow=2, dim.lgd.key.size=4, dim.lgd.text.size=13, dim.capt.size=13, add.feature.2nd=FALSE, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, line.width=0.2, line.color='grey70', relative.scale = NULL, verbose=TRUE, out.dir=NULL, animation.scale = 1, selfcontained=FALSE, video.dim='640x480', res=500, interval=1, framerate=1, bar.width.vdo=0.1, legend.value.vdo=NULL, ...) {
  bulkCell <- NULL
  if ('bulkCell' %in% colnames(sce.dimred)) sce.dimred <- subset(sce.dimred, , bulkCell=='cell')
  res <- shm_covis(svg=svg, data=data, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, ID=ID, sce.dimred=sce.dimred, dimred=dimred, cell.group=cell.group, tar.cell=tar.cell, tar.bulk=tar.bulk, profile=profile, charcoal=charcoal, alpha.overlay=alpha.overlay, lay.shm=lay.shm, ncol=ncol, col.com=col.com, col.bar=col.bar, sig.thr=sig.thr, cores=cores, bar.width=bar.width, bar.title.size=bar.title.size, trans.scale=trans.scale, ft.trans=ft.trans, lis.rematch = lis.rematch, legend.r=legend.r, sub.title.size=sub.title.size, sub.title.vjust=sub.title.vjust, legend.plot=legend.plot, ft.legend=ft.legend, bar.value.size=bar.value.size, legend.plot.title=legend.plot.title, legend.plot.title.size=legend.plot.title.size, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.text.size=legend.text.size, angle.text.key=angle.text.key, position.text.key=position.text.key, legend.2nd=legend.2nd, position.2nd=position.2nd, legend.nrow.2nd=legend.nrow.2nd, legend.ncol.2nd=legend.ncol.2nd, legend.key.size.2nd=legend.key.size.2nd, legend.text.size.2nd=legend.text.size.2nd, angle.text.key.2nd=angle.text.key.2nd, position.text.key.2nd=position.text.key.2nd, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.capt.size=dim.capt.size, add.feature.2nd=add.feature.2nd, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, line.width=line.width, line.color=line.color, relative.scale = relative.scale, verbose=verbose, out.dir=out.dir, animation.scale = animation.scale, selfcontained=selfcontained, video.dim=video.dim, res=res, interval=interval, framerate=framerate, bar.width.vdo=bar.width.vdo, legend.value.vdo=legend.value.vdo, ...)
  invisible(res)
})


#' SHM plots and co-visaulization plots 
#' 
#' @keywords Internal
#' @noRd

#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 ggplot geom_bar aes theme element_blank margin element_rect coord_flip scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs
#' @importFrom rsvg rsvg_ps 
#' @importFrom grImport PostScriptTrace 
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid grobTree unit
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' @importFrom ggplotify as.ggplot
#' @importFrom SingleCellExperiment reducedDimNames

shm_covis <- function(svg, data, assay.na=NULL, sam.factor=NULL, con.factor=NULL, ID, sce.dimred=NULL, dimred='PCA', cell.group=NULL, tar.cell=NULL, tar.bulk=NULL, profile=TRUE, charcoal=FALSE, alpha.overlay=1, lay.shm="gene", ncol=2, col.com=c('yellow', 'orange', 'red'), col.bar='selected', sig.thr=c(NA, NA), cores=NA, bar.width=0.08, bar.title.size=0, trans.scale=NULL, ft.trans=NULL, tis.trans=ft.trans, lis.rematch = NULL, legend.r=0.9, sub.title.size=11, sub.title.vjust=2, legend.plot='all', ft.legend='identical', bar.value.size=10, legend.plot.title='Legend', legend.plot.title.size=11, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, angle.text.key=NULL, position.text.key=NULL, legend.2nd=FALSE, position.2nd='bottom', legend.nrow.2nd=NULL, legend.ncol.2nd=NULL, legend.key.size.2nd=0.03, legend.text.size.2nd=10, angle.text.key.2nd=0, position.text.key.2nd='right', dim.lgd.pos='bottom', dim.lgd.nrow=2, dim.lgd.key.size=4, dim.lgd.text.size=13, dim.capt.size=13, add.feature.2nd=FALSE, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, line.width=0.2, line.color='grey70', relative.scale = NULL, verbose=TRUE, out.dir=NULL, animation.scale = 1, selfcontained=FALSE, video.dim='640x480', res=500, interval=1, framerate=1, bar.width.vdo=0.1, legend.value.vdo=NULL, ...) {
  
  calls <- names(vapply(match.call(), deparse, character(1))[-1])
  if("tis.trans" %in% calls) warning('"tis.trans" is deprecated and replaced by "ft.trans"! \n')
  if("svg.path" %in% calls) warning('"svg.path" is deprecated and replaced by "svg"! \n')
  # save(svg, data, assay.na, sam.factor, con.factor, ID, sce.dimred, dimred, tar.cell, profile, cell.group, tar.bulk, charcoal, alpha.overlay, lay.shm, ncol, col.com, col.bar, sig.thr, cores, bar.width, bar.title.size, trans.scale, ft.trans, tis.trans, lis.rematch, legend.r, sub.title.size, sub.title.vjust, legend.plot, ft.legend, bar.value.size, legend.plot.title, legend.plot.title.size, legend.ncol, legend.nrow, legend.position, legend.direction, legend.key.size, legend.text.size, angle.text.key, position.text.key, legend.2nd, position.2nd, legend.nrow.2nd, legend.ncol.2nd, legend.key.size.2nd, legend.text.size.2nd, angle.text.key.2nd, position.text.key.2nd, dim.lgd.pos, dim.lgd.nrow, dim.lgd.key.size, dim.lgd.text.size, dim.capt.size, add.feature.2nd, label, label.size, label.angle, hjust, vjust, opacity, key, line.width, line.color, relative.scale, verbose, out.dir, animation.scale, selfcontained, video.dim, res, interval, framerate, bar.width.vdo, legend.value.vdo, file='shm.covis.arg')

  x <- y <- color_scale <- tissue <- bulkCell <- NULL; options(stringsAsFactors=FALSE)
  # if (!is.null(sub.margin)) if (!is.numeric(sub.margin) | length(sub.margin)!=4 | any(sub.margin >= 1) | any(sub.margin < 0)) stop('"sub.margin" must be a 4-length numeric vector between 0 (inclusive) and 1 (exclusive)!')
  if (!is(svg, 'SVG')) stop('The "svg" should be a "SVG" object!')
  svgs <- svg; if (missing(ID)) ID <- rownames(data)[1]
  ID <- unique(ID)
  svg.ft.all <- unique(unlist(lapply(seq_along(svgs), function(x) {  unique(attribute(svgs[x])[[1]]$feature)
  })))
  is.sce <- is(sce.dimred, 'SingleCellExperiment')
  profile.no <- profile==FALSE & is.sce 
  # Original lis.rematch used in 'toCell'.
  covis.type <- NULL; lis.rematch.raw <- lis.rematch
  if (is.sce) { # Check coviz related data.
    lis <- check_data_covis(data, dimred, cell.group, sce.dimred, tar.cell, tar.bulk, lis.rematch, svg.ft.all) 
    tar.bulk <- lis$tar.bulk; tar.cell <- lis$tar.cell
    covis.type <- lis$covis.type; lis.rematch <- lis$lis.rematch
    bulk.all <- lis$bulk.all; cell.all <- lis$cell.all
    rm(lis)
  }

  # Extract and filter data.
  if (is.vector(data)) {
    vec.na <- make.names(names(data)); if (is.null(vec.na)) stop("Please provide names for the input data!")
    if (!identical(vec.na, names(data))) cat('Syntactically valid column names are made! \n')
    if (any(duplicated(vec.na))) stop('Please make sure data names are unique!')
    form <- grepl("__", vec.na); if (sum(form)==0) { vec.na <- paste0(vec.na, '__', 'con'); con.na <- FALSE } else con.na <- TRUE
    data <- tryCatch({ as.numeric(data) }, warning=function(w) { stop("Please make sure input data are numeric!") }, error=function(e) { stop("Please make sure input data are numeric!") })
    if (is.null(ID)) stop('Please provide a name for the data!')
    gene <- as.data.frame(matrix(data, nrow=1, dimnames=list(ID, vec.na)))

  } else if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')|is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
    if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) { # Data frame of spatial enrichment.
      cna <- colnames(data)
      if (all(c('gene', 'type', 'total') %in% cna)) {
        data <- subset(data, !duplicated(gene)); rownames(data) <- data$gene
      }
      data <- data[, !colnames(data) %in% c('gene', 'type', 'total', 'metadata', 'edgeR', 'limma', 'DESeq2', 'distinct'), drop=FALSE]
    }
    id.no <- ID[!ID %in% rownames(data)]
    if (length(id.no)>0) stop(paste0(id.no, collapse=' '), ': not detected in data! \n')
    # Process data.
    dat.lis <- check_data(data=data, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, usage='shm')
    gene <- as.data.frame(dat.lis$dat); con.na <- dat.lis$con.na

  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", "dgCMatrix", "SummarizedExperiment", or "SingleCellExperiment" except that "spatial_hm" also accepts a "vector".') }

  if (!is.null(trans.scale)) if (trans.scale=='log2') {      
      g.min <- min(gene) 
      if (g.min<0) gene <- gene-g.min+1; if (g.min==0) gene <- gene+1; gene <- log2(gene)
    } else if (trans.scale=='exp2') gene <- 2^gene else if (trans.scale=='row') { 
    gene <- t(scale(t(gene))) } else if (trans.scale=='column') { gene <- scale(gene) }
 
    # Color bar.
    sig.thr <- dat_fun(sig.thr, as.numeric)
    if (length(sig.thr)!=2) stop('The "sig.thr" must be a two-element vecor and contain as least one numeric!')
    if (any(is.na(sig.thr))) { # Default signal threshold
      if (col.bar=="all") { 
        min.v <- min(gene); max.v <- max(gene) 
      } else if (col.bar=="selected") {
        vals <- gene[ID, , drop=FALSE]
        min.v <- min(vals); max.v <- max(vals)
      }
    }
    # Customized signal threshold.
    thr1 <- sig.thr[1]; thr2 <- sig.thr[2]
    if (is.numeric(thr1) & !is.na(thr1)) min.v <- thr1
    if (is.numeric(thr2) & !is.na(thr2)) max.v <- thr2
     
     if (max.v <= min.v) stop('Make sure the max signal threshold is larger than the min!')
    bar.len=1000; geneV <- seq(min.v, max.v, len=bar.len)

    col <- colorRampPalette(col.com)(length(geneV))
    cs.g <- col_bar(geneV=geneV, cols=col, width=1, bar.title.size=bar.title.size, bar.value.size=bar.value.size)

    # Only take the column names with "__".
    cname <- colnames(gene); form <- grepl('__', cname)
    con <- gsub("(.*)(__)(.*)", "\\3", cname[form]); con.uni <- unique(con)
    sam.uni <- unique(gsub("(.*)(__)(.*)", "\\1", cname))

    svg.pa.na <- img_pa_na(unlist(svgs[, 'svg']))
    # svg.path may not be paths, can be file names, if users provide a SVG class that not includes paths.
    svg.path <- svg.pa.na$path; svg.na <- svg.pa.na$na
    svg.na.cord <- names(svgs)
    # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
    w.h.max <- max(unlist(svgs[, 'dimension']))
   
    # A set of SHMs (gene_con*) are made for each SVG, and all sets of SHMs under different SVGs are placed in 2 lists in form of ggplots and grobs respectively. Different SHMs of same 'gene_condition' under different SVGs are indexed with suffixed of '_1', '_2', ... E.g. SHMs under shm1.svg: gene_condition1_1, gene_condition2_1; SHMs under shm2.svg: gene_condition1_2, gene_condition2_2; the 4 SHMs are stored in 2 separate lists in form of ggplots and grobs respectively. 
    # The order of ggplot SHMs, grob SHMs, and legends follow the order of SVGs, so all are same.
    ft.trans.shm <- NULL; if (is.sce) {
      # The only use of tar.cell, tar.bulk is to get ft.trans.shm, they are not required in dim_color, dim_color2cell, and dim_color_coclus, since the color transfer is data -> SHM -> dim plot and transparent colors are defined in SHM before dim plot.
      lis.covis <- covis_trans(bulk.all=bulk.all, cell.all=cell.all, ft.all=svg.ft.all, tar.bulk=tar.bulk, tar.cell=tar.cell, lis.match=lis.rematch.raw, covis.type=covis.type)
      tar.bulk <- lis.covis$tar.bulk
      tar.cell <- lis.covis$tar.cell
      lis.rematch <- lis.covis$lis.match 
      ft.trans.shm <- lis.covis$ft.trans.shm
    }
    gg.lis.all <- gg.all <- grob.all <- lgd.all <- grob.lgd.all <- gcol.lgd.all <- gcol.all <- NULL
    for (i in seq_along(svgs)) {
      na0 <- svg.na[i]; cat('ggplots/grobs:', na0, '... \n')
      # if (preserve.scale==TRUE & is.null(sub.margin)) sub.margin <- (1-w.h/w.h.max*0.99)/2
      # SHMs/legend of ggplots under one SVG.
      gg.lis <- gg_shm(svg.all=svgs[i], gene=gene, con.na=con.na, geneV=geneV, charcoal=charcoal, alpha.overlay=alpha.overlay, ID=ID, cols=col, covis.type=covis.type, lis.rematch = lis.rematch, ft.trans=ft.trans, ft.trans.shm=ft.trans.shm, sub.title.size=sub.title.size, sub.title.vjust=sub.title.vjust, ft.legend=ft.legend, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.text.size=legend.text.size, legend.plot.title=legend.plot.title, legend.plot.title.size=legend.plot.title.size, line.width=line.width, line.color=line.color, ...)
      msg <- paste0(na0, ': no spatial features that have matching sample identifiers in data are detected!'); if (is.null(gg.lis)) stop(msg)

      # Append suffix '_i' for the SHMs of ggplot under SVG[i], and store them in a list.
      ggs <- gg.lis$g.lis.all; names(ggs) <- paste0(names(ggs), '_', i)
      gg.all <- c(gg.all, ggs)
      cores <- deter_core(cores, svg.path[i]); cat('CPU cores:', cores, '\n')
      # Same names with ggplot: append suffix '_i' for the SHMs of grob under each SVG, and store them in a list. 'i' is equal to SVG.
      grobs <- grob_shm(ggs, cores=cores); grob.all <- c(grob.all, grobs)

      # Store SHM ggplots and w.h under each SVG in separate list for use in relative scale.
      lis0 <- list(list(ggs = ggs, w.h = dimension(svgs)[[1]]))
      names(lis0) <- na0
      gg.lis.all <- c(gg.lis.all, lis0)
    
      # Store legend/colour of ggplot in a list.
      gg.lgd <- gg.lis$g.lgd
      if (profile.no) gg.lgd <- gg.lgd + labs(title=NULL)
      lgd.all <- c(lgd.all, list(gg.lgd))
      grob.lgd.lis <- grob_shm(lgd.all, cores=cores, lgd.pos='bottom')
      grob.lgd.all <- c(grob.lgd.all, grob.lgd.lis)
      gcol.lgd.all <- c(gcol.lgd.all, list(gg.lis$gcol.lgd))

      # Store colors of matching features in each SHM in a list.
      gcols <- gg.lis$gcol.lis.all; names(gcols) <- paste0(names(gcols), '_', i)
      gcol.all <- c(gcol.all, gcols)
    }; names(lgd.all) <- names(grob.lgd.all) <- svg.na
    names(gcol.lgd.all) <- paste0('col_', svg.na)
    # Summary of mapping.
    map.sum <- mapping_sum(svg.na.cord, sam.uni, gcol.all, data=gene, ID, con.na, covis.type=covis.type, lis.match=lis.rematch.raw, verbose)

    pat.gen <- paste0(ID, collapse='|'); pat.con <- paste0(con.uni, collapse='|')
    # Use definite patterns and avoid using '.*' as much as possible. Try to be as specific as possible.
    pat.all <- paste0('^(', pat.gen, ')_(', pat.con, ')(_\\d+$)')
    # Indexed cons with '_1', '_2', ... at the end.
    con.idxed <- unique(gsub(pat.all, '\\2\\3', names(gg.all)))
    # Layout matrix of SHMs for use in relative scale.
    lay.mat <- lay_shm(lay.shm=lay.shm, con=con.idxed, ncol=ncol, ID.sel=ID, lay.mat = TRUE, profile=profile) 

    # If relative size in multiple SVGs is required, update all SHM ggplots/grobs under each SVG.
    if (is.numeric(relative.scale) & length(svg.na) > 1) if (relative.scale > 0) {
      cat('Applying relative image size ... \n')
      gg.all <- grob.all <- NULL; for (i in seq_along(gg.lis.all)) {
        lis0 <- gg.lis.all[[i]]
        ggs <- rela_size(lis0$w.h['height'], w.h.max, relative.scale, nrow(lay.mat), lis0$ggs)
        gg.all <- c(gg.all, ggs)
        cores <- deter_core(cores, svg.path[i]) # gg.lis.all has the same names with svgs.
        # Same names with ggplot: append suffix '_i' for the SHMs of grob under each SVG, and store them in a list.
        grobs <- grob_shm(ggs, cores=cores); grob.all <- c(grob.all, grobs)
        
      }
    }; na.all <- names(grob.all)

    # Arrange color scale grob.
    tmp <- normalizePath(tempfile(), winslash='/', mustWork=FALSE)
    png(tmp); cs.grob <- ggplotGrob(cs.g); dev.off()
    if (file.exists(tmp)) do.call(file.remove, list(tmp))
    cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc")) # "mm" is fixed, "npc" is scalable.
    if (legend.2nd==TRUE) {
      gg.all <- gg_2lgd(gg.all=gg.all, sam.dat=sam.uni, ft.trans=ft.trans, position.2nd=position.2nd, legend.nrow.2nd=legend.nrow.2nd, legend.ncol.2nd=legend.ncol.2nd, legend.key.size.2nd=legend.key.size.2nd, add.feature.2nd=add.feature.2nd, legend.text.size.2nd=legend.text.size.2nd, angle.text.key.2nd=angle.text.key.2nd, position.text.key.2nd=position.text.key.2nd)
      grob.all <- lapply(gg.all, ggplotGrob)
    }

    # Arrange SHM grobs.
    na.all <- sort_gen_con(ID.sel=ID, na.all=na.all, con.all=con.idxed, by=lay.shm)
    grob.all <- grob.all[na.all]; gg.all <- gg.all[na.all]
    gcol.all <- gcol.all[paste0('col_', na.all)]

    if (is.sce) {
      if ('toBulk' %in% covis.type) {
      gg.dim <- plot_dim(sce.dimred, dim=dimred, color.by=cell.group)
      dim.shm.lis <- dim_color(gg.dim=gg.dim, gg.shm.all=gg.all, grob.shm.all=grob.all, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=grob.lgd.all, cell.group=cell.group, tar.cell=tar.cell, con.na=con.na, profile=profile, lis.match=lis.rematch, sub.title.size=sub.title.size, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size)
      } else if (any(c('toBulkAuto', 'toCellAuto') %in% covis.type)) {
      # dim_color_coclus: applicable to toBulkAuto or toCellAuto, since svg features are assigned to cells (one-to-one) as labels, which are the basis for covis.
      # Use the same dim plot through which desired bulk is assigned.
      sce.cell <- subset(sce.dimred, , bulkCell=='cell')
      dimred.pre <- unique(colData(data)$dimred)
      if (!is.null(dimred.pre)) {
        dimred <- dimred.pre
        cat('The reduced dimensionality in "data" is used:', dimred, '.\n')
      }
      gg.dim <- plot_dim(sce.dimred, dim=dimred, color.by='assignedBulk')
      if ('toBulkAuto' %in% covis.type) targ <- tar.cell
      if ('toCellAuto' %in% covis.type) targ <- tar.bulk
      dim.shm.lis <- dim_color_coclus(sce=sce.dimred, targ=targ, profile=profile, gg.dim = gg.dim, gg.shm.all=gg.all, grob.shm.all = grob.all, gg.lgd.all=lgd.all, col.shm.all = gcol.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=grob.lgd.all, con.na=con.na, lis.match=NULL, sub.title.size=sub.title.size, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.capt.size=dim.capt.size)
      } else if ('toCell' %in% covis.type) {
      gg.dim <- plot_dim(sce.dimred, dim=dimred, color.by=cell.group)
      dim.shm.lis <- dim_color2cell(gg.dim=gg.dim, gg.shm.all=gg.all, grob.shm.all=grob.all, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=grob.lgd.all, profile=profile, cell.group=cell.group, con.na=con.na, lis.match=lis.rematch.raw, sub.title.size=sub.title.size, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size)
      }
      grob.all <- dim.shm.lis$dim.shm.grob.lis
      gg.all <- dim.shm.lis$dim.shm.gg.lis
    }
    g.arr <- lay_shm(lay.shm=lay.shm, con=con.idxed, ncol=ncol, ID.sel=ID, grob.list=grob.all, scell=is.sce, profile=profile, shiny=FALSE)
    if (profile.no) legend.plot <- NULL
    if (!is.null(legend.plot)) {
      # Select legend plot to show. 
      if (length(svg.na)>1 & legend.plot!='all') na.lgd <- svg.na[grep(paste0('_(', paste0(legend.plot, collapse='|'), ').svg$'), svg.na)] else na.lgd <- svg.na
      lgd.lis <- lgd.all[na.lgd]
 
      # Add labels to target shapes/adjust keys in legend plots.
      lgd.lis <- gg_lgd(gg.all=lgd.lis, angle.text.key=angle.text.key, position.text.key=position.text.key, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key)
      grob.lgd.lis <- lapply(lgd.lis, ggplotGrob)
      lgd.tr <- lapply(grob.lgd.lis, grobTree)

      w.lgd <- (1-bar.width)/(ncol+1); shm.w <- 1-bar.width-w.lgd
      # If legend.r = 0, legend plot size is a square.
      lgd.arr <- arrangeGrob(grobs=lgd.tr, layout_matrix=matrix(seq_along(na.lgd), ncol=1), widths=unit(0.99, "npc"), heights=unit(rep(w.lgd + (0.99 - w.lgd) * legend.r, length(na.lgd)), "npc"))

      # A plot pops up when 'grid.arrange' runs.
      shm <- grid.arrange(cs.arr, g.arr, lgd.arr, ncol=3, widths=unit(c(bar.width-0.005, shm.w, w.lgd), 'npc'))
    } else { 
      shm.w <- 1-bar.width
      if (!profile.no) shm <- grid.arrange(cs.arr, g.arr, ncol=2, widths=unit(c(bar.width-0.005, shm.w), 'npc')) else shm <- grid.arrange(g.arr, ncol=1, widths=unit(c(1-0.005), 'npc'))
    }

    if (!is.null(out.dir)) { 

      for (i in names(gg.all)) {
        g0 <- gg.all[[i]]
        anm.height <- 550 * animation.scale; anm.width <- anm.height / g0$theme$aspect.ratio
        html_ly(gg=gg.all[i], cs.g=cs.g, ft.trans=ft.trans, sam.uni=sam.uni, anm.width=anm.width, anm.height=anm.height, out.dir=out.dir)
      }

      vdo <- video(gg=gg.all, cs.g=cs.g, bar.width=bar.width.vdo, lgd.key.size=legend.key.size.2nd, lgd.text.size=legend.text.size.2nd, angle.text.key=angle.text.key.2nd, position.text.key=position.text.key.2nd, lgd.row=legend.nrow.2nd, lgd.col=legend.ncol.2nd, legend.value.vdo=legend.value.vdo, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, video.dim=video.dim, res=res, interval=interval, framerate=framerate, out.dir=out.dir)
      if (is.null(vdo)) cat("Video is not generated! \n")

    }
    lis <- list(spatial_heatmap=as.ggplot(shm), mapped_feature=map.sum); return(lis)
}

#' Summary of mapping for data features and colors.
#' @param svg.na.cord The instance names of aSVGs stored in coord.
#' @param sam.uni Data features.
#' @param gcol.all A list of mapped colors.
#' @param data The assay data
#' @param ID The molecule IDs.
#' @param con.na If FALSE, no treatment is included in data.
#' @param covis.type The type of covisualization.
#' @param lis.match The matching list in co-visualization.
#' @param verbose If TRUE, intermediate messages are printed. 
#' @keywords Internal
#' @noRd

mapping_sum <- function(svg.na.cord, sam.uni, gcol.all, data, ID, con.na, covis.type=NULL, lis.match=NULL, verbose) {
  # save(svg.na.cord, sam.uni, gcol.all, data, ID, con.na, verbose, file='mapping.sum.arg')
  map.sum <- data.frame()
  for (j in seq_along(svg.na.cord)) {
    svg.na.cord0 <- svg.na.cord[j] 
    col.na0 <- names(gcol.all) 
    gcol0 <- lis0 <- gcol.all[grepl(paste0('_', j, '$'), col.na0)]
    names(lis0) <- NULL; ft.all0 <- unlist(lis0)
    ft.mapped.col <- grep('NA', ft.all0, value=TRUE, invert=TRUE) 
    ft.mapped <- unique(sub('__\\d+$', '', names(ft.mapped.col)))
    if (!'toBulk' %in% covis.type) {
      not.map <- setdiff(sam.uni, ft.mapped)
      if (verbose==TRUE & length(not.map)>0) cat('Features in data not mapped in', svg.na.cord0, ':', paste0(not.map, collapse=', '), '\n')
    } else {
      ft.mapped <- unique(unlist(lapply(seq_along(lis.match), function(x) { if (any(lis.match[[x]] %in% ft.mapped)) names(lis.match[x]) }))) 
    }
    if (length(ft.mapped)==0) next 
    idx.com <- vapply(colnames(data), function(i) grepl(paste0('^', ft.mapped, '__', collapse='|'), i), FUN.VALUE=logical(1))
    map.gene <- data[idx.com]; cna <- colnames(map.gene)
    for (i in ID) {
      df0 <- as.data.frame(t(map.gene[i, ]))
      colnames(df0) <- 'value'
      feature <- gsub("(.*)(__)(.*)", "\\1", cna)
      ID <- i; df1 <- data.frame(ID=ID, feature=feature)
      if (con.na==TRUE) {
        condition <- gsub("(.*)(__)(.*)", "\\3", cna)
        df1 <- cbind(df1, condition=condition) 
      } else df1 <- cbind(df1, condition='con')
      df1 <- cbind(df1, df0)
      for (k in seq_len(nrow(df1))) {
        col0 <- gcol0[[paste0('col_', df1$ID[k], '_', df1$condition[k], '_', j)]]
        df1$fill[k] <- col0[sub('__\\d+$', '', names(col0)) %in% df1$feature[k]][1]
      } 
      df1$SVG <- svg.na.cord0
      if (con.na==FALSE) df1$condition <- NULL
      map.sum <- rbind(map.sum, df1)
    }  
  }; row.names(map.sum) <- NULL; map.sum
}

#' Get SVG names and order SVG path/name if there are multiple SVGs.
#' @param path.na A vector of image paths or names (SVG or raster).
#' @param raster.ext A vector of raster image extensions.
#' @keywords Internal
#' @noRd

img_pa_na <- function(path.na, raster.ext=c('.jpg', '.JPG', '.png', '.PNG')) {
  str.lis <- strsplit(path.na, '/')
  na <- vapply(str.lis, function(x) x[[length(x)]], FUN.VALUE=character(1))
  ext <- c('.svg', '.SVG', raster.ext); ext <- paste0('(', paste0('\\', ext, collapse='|'), ')')
  pat <- paste0('.*_(shm\\d+)', ext, '$') 
    if (length(na)>1) { 
      shm <- gsub(pat, '\\1', na)
      if (any(duplicated(shm))) stop(paste0('Suffixes of aSVG files are duplicated: ', paste0(shm, collapse='; ')))
      if (!all(grepl('_shm\\d+', na, perl=TRUE))) stop("Suffixes of aSVG files should be indexed as '_shm1', '_shm2', '_shm3', ...") 
    }
    ord <- order(gsub(pat, '\\1', na))
    path <- path.na[ord]; na <- na[ord]; return(list(path=path, na=na))
}

#' Determine the number of CPU cores.
#' @param cores The input number of cores.
#' @param svg.path The path of SVG file.
#' @keywords Internal
#' @noRd
#' @importFrom parallel detectCores

deter_core <- function(cores, svg.path) {
  if (!file.exists(svg.path)) return(cores)
  cores <- as.integer(cores)
  n.cor <- detectCores(logical=TRUE)
  fil.size <- file.size(svg.path)
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

#' Convert vectors.
#' @param from,to,target The \code{target} vector comes from \code{from} vector and is converted to based on mapping between \code{from} and \code{to}. \code{from} and \code{to} is one-to-one mapped. 

#' @keywords Internal
#' @noRd
cvt_vector <- function(from, to, target) { names(to) <- from; to[target] }


#' Check input data in coviz, convert target bulk to SVG bulk.

#' @keywords Internal
#' @noRd

#' @importFrom SummarizedExperiment colData

check_data_covis <- function(data, dimred, cell.group, sce.dimred, tar.cell, tar.bulk, lis.rematch, svg.ft.all) {
  cdat.cell <- colData(sce.dimred)
  coclus.na <- c('assignedBulk', 'similarity')
  bulk.all <- cell.all <- covis.type <- NULL
  if (is.null(tar.cell) & is.null(tar.bulk) & any(!coclus.na %in% colnames(cdat.cell))) stop('In co-visualization, one of "tar.cell" and "tar.bulk" needs to be specified!')
  if (!is.null(tar.cell) & !is.null(tar.bulk)) stop('In co-visualization, only one of "tar.cell" and "tar.bulk" needs to be specified!')
  if (!dimred %in% reducedDimNames(sce.dimred)) stop(paste0(dimred, " is not detected in 'reducedDimNames'!"))
  if (!is.null(cell.group)) if (!cell.group %in% colnames(cdat.cell)) stop(paste0(cell.group, " is not detected in 'colData' slot!"))
  # Check cell2bulk mapping in coclustering. 
  if (all(coclus.na %in% colnames(cdat.cell))) { 
    # Check bulk2cell mapping in coclustering.
    bulk.all <- unique(data$sample)
    cell.all <- setdiff(unique(sce.dimred$assignedBulk), 'none')
    if (is(data, 'SingleCellExperiment') | is(data, 'SummarizedExperiment')) {
      if (unique(data$bulkCell)=='bulk') {
        if (is.null(tar.bulk)) tar.bulk <- bulk.all
        tar.cell <- NULL; covis.type <- 'toCellAuto'
      } else if (unique(data$bulkCell)=='cell'){
        if (is.null(tar.cell)) tar.cell <- cell.all
        tar.bulk <- NULL; covis.type <- 'toBulkAuto'
      }
    }
  } else if (!is.null(lis.rematch) & !is.null(tar.bulk)) {
    if (any(!tar.bulk %in% names(lis.rematch))) stop("Make sure all entries in 'tar.bulk' are in 'names(lis.rematch))'!")
    if (any(!names(lis.rematch) %in% svg.ft.all)) stop("Make sure all entries in 'names(lis.rematch)' are from aSVG features!")
    bulk.all <- unique(colnames(data))
    cell.all <- unique(cdat.cell[, cell.group])
    lis.rematch <- NULL; covis.type <- 'toCell'
  } else if (!is.null(lis.rematch) & !is.null(tar.cell)) {
    bulk.all <- unique(colnames(data))
    cell.all <- unique(cdat.cell[, cell.group])
    covis.type <- 'toBulk'
  }
  return(list(bulk.all=bulk.all, cell.all=cell.all, tar.bulk=unique(tar.bulk), tar.cell=unique(tar.cell), covis.type=covis.type, lis.rematch=lis.rematch))
}


#' Extract transparent aSVG features in co-visualization

#' @keywords Internal
#' @noRd

covis_trans <- function(bulk.all, cell.all, ft.all, tar.bulk='all', tar.cell='all', lis.match, covis.type) {
 #save(bulk.all, ft.all, cell.all, tar.bulk, tar.cell, lis.match, covis.type, file='covis.trans.arg')
  ft.trans.shm <- NULL
  if (covis.type %in% 'toBulk') {
    if (is(lis.match, 'list')) { # Same as rematching in SHM.
      cell.grp.uni <- unique(names(lis.match))
      if ('all' %in% tar.cell) tar.cell <- cell.grp.uni
      # SVG features corresponding to non-target cells are set transparent.
      ft.trans.shm <- unique(unlist(lis.match[setdiff(cell.grp.uni, tar.cell)]))
      # Cell labels in aSVG are set transparent, since they are not expected to overlap.
      cell.ft.inter <- intersect(cell.all, ft.all)
      cell.ft.inter <- cell.ft.inter[!cell.ft.inter %in% unlist(lis.match)]
      ft.trans.shm <- unique(c(ft.trans.shm, cell.ft.inter))
    }
  } else if (covis.type %in% 'toCell') {
    if (!is.null(bulk.all)) {
      if ('all' %in% tar.bulk & is(lis.match, 'list')) tar.bulk <- unique(names(lis.match))
        ft.trans.shm <- unique(setdiff(bulk.all, tar.bulk))
        lis.match <- NULL
      }
  } else if (covis.type %in% 'toBulkAuto') {
    ft.trans.shm <- unique(setdiff(cell.all, tar.cell))
    if ('all' %in% tar.cell) ft.trans.shm <- NULL
  } else if (covis.type %in% 'toCellAuto') {
    ft.trans.shm <- unique(setdiff(bulk.all, tar.bulk))
    if ('all' %in% tar.bulk) ft.trans.shm <- NULL
  }
  return(list(tar.bulk=tar.bulk, tar.cell=tar.cell, lis.match=lis.match, ft.trans.shm=ft.trans.shm))
}
  

