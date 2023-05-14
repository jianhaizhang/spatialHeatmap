#' Co-visualizing spatial heatmaps with single-cell embedding plots
#'
#' This function is an extension of \code{\link{shm}}. It is designed to co-visualize bulk and single-cell data in form of a spatial heatmap (SHM) and embedding plot (PCA, UMAP, TSNE) respectively.

#' @param data An `SHM` class that containing the numeric data and aSVG instances for plotting SHMs or co-visualization plots. See \code{\link{SPHM}}.  
#' @param charcoal Logical, if \code{TRUE} the raster image will be turned black and white.
#' @param alpha.overlay The opacity of the raster image under the SHM when superimposing raster images with SHMs. The default is 1.
#' @param dimred One of \code{PCA}, \code{UMAP}, \code{TSNE} in \code{sce.dimred}, specifying which reduced dimension to use in co-visualization plots.
#' @param cell.group Applicable when co-visualizing bulk and single-cell data with cell grouping methods of annodation labels, marker genes, clustering, or manual assignments. A column name in \code{colData(sce.dimred)}, where one label defines a cell group. 
#' @param tar.cell Applicable when co-visualizing bulk and single-cell data with cell-to-bulk mapping. A vector of cell group labels in \code{cell.group} which are mapped to tissues through \code{lis.rematch}. Cells corresponding to these labels will be colored in the embedding plot while other cells will be grey, and tissues corresponding to these labels will be colored in the SHM while other tissues will be transparent.  
#' @param tar.bulk A vector of tissues of interest, which are mapped to single cells through cell group labels with \code{lis.rematch}. These tissues will be colored in SHMs while other tissues will be transparent, and cells corresponding to these tissues will be colored embedding plots while other cells will be grey. 
#' @param profile Logical, applicable when co-visualizing bulk and single-cell data. If \code{TRUE}, one or multiple biomolecule (e.g. gene) identifiers need to be assigned to \code{ID}, and their abundance proifles will be used for coloring in the co-visualization plots. If \code{FALSE}, constant colors will be used in the co-visualization plot. 
#' @param bar.value.size A numeric of value size in the y-axis of the color bar. The default is 10.

#' @inheritParams filter_data
#' @inheritParams htmlwidgets::saveWidget
#' @inheritParams ggplot2::theme

#' @param position.2nd The position of the secondary legend in SHMs, one of `top`, `right`, `bottom`, `left`, or a two-component numeric vector. The default is `bottom`. Applies to images and videos.
#' @param legend.nrow.2nd An integer of rows of the secondary legend keys in SHMs. Applies to static images and videos.
#' @param legend.ncol.2nd An integer of columns of the secondary legend keys in SHMs. Applies to static images and videos.
#' @param legend.key.size.2nd A numeric of legend key size in SHMs. The default is 0.03. Applies to static images and videos.
#' @param dim.lgd.pos The legend position in dimension reduction plot. The default is \code{bottom}.
#' @param dim.lgd.nrow The number of legend rows in dimension reduction plot. The default is \code{2}.
#' @param add.feature.2nd Logical. If `TRUE`, feature identifiers are added to the secondary legend. The default is FALSE. Applies to static images of SHMs.
#' @param legend.text.size.2nd A numeric of the secondary legend text size in SHMs. The default is 10. Applies to static images and videos.
#' @param angle.text.key.2nd Angle of the key text in the secondary legend in SHMs. Default is 0. Applies to static images and videos.
#' @param position.text.key.2nd The position of key text in the secondary legend in SHMs, one of `top`, `right`, `bottom`, `left`. Default is `right`. Applies to static images and videos.
#' @param dim.lgd.pos The legend position in the dimension reduction plot. The default is \code{bottom}.
#' @param dim.lgd.nrow The number of legend rows in the dimension reduction plot. The default is \code{2}.
#' @param dim.capt.size The size of caption text in the dimension reduction plot in co-clustering (see \code{\link{cocluster}}). The default is \code{13}.
#' @param angle.text.key Key text angle in legend plots. The default is NULL, equivalent to 0.
#' @param position.text.key The position of key text in legend plots, one of `top`, `right`, `bottom`, `left`. Default is NULL, equivalent to `right`.
#' @param bar.width.vdo The color bar width (0-1) in videos.
#' @param legend.value.vdo Logical. If `TRUE`, numeric values of colored spatial features are added to the video legend. The default is NULL.
#' @param label Logical. If `TRUE`, the same spatial features between numeric data and aSVG are labeled by their identifiers. The default is `FALSE`. It is useful when spatial features are labeled by similar colors. 
#' @param label.size The size of spatial feature labels in legend plots. The default is 4.
#' @param label.angle The angle of spatial feature labels in legend plots. Default is 0.
#' @param hjust,vjust The value to horizontally or vertically adjust positions of spatial feature labels in legend plots respectively. Default of both is 0.
#' @param opacity The transparency of colored spatial features in legend plots. Default is 1. If 0, features are totally transparent.
#' @param key Logical. If `TRUE` (default), keys are added in legend plots. If \code{label} is TRUE, the keys could be removed. 
#' @param ft.trans A character vector of spatial features that will be set transparent. When features of interest are covered by overlapping features on the top layers and the latter can be set transparent.
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
#' @param var.cell The column name in the \code{colData} slot of \code{SingleCellExperiment} that indicates the experimental variables. If \code{NULL}, no variables are considered. 
#' @param size.pt,alpha.pt The size and alpha value of points in the embedding plot respectively. 
#' @param shape A named numeric vector of custom shapes for points in the embedding plot, where the names are the cluster labels and values are from \code{c(0, 2:25, 32:127)}.
#' @param col.idp Logical, if \code{TRUE}, each cell in the embedding plot and spatial feature in the SHM are colored independently according the expression values of chosen biomolecules. 
#' @param decon Logical, if \code{TRUE}, the cell data will be considered from bulk deconvolution results.
#' @param h The overall height of embedding plot and SHM, ranging from 0 to 1.
#' @param size.r The scaling ratio (0-1) of tissue section when visualizing the spatially resolved single-cell data.
#' @param bar.title,bar.title.size The title and title size of the color key.
#' @param scale How to scale the assay data. If \code{'no'}, no scaling is applied. If \code{'row'} each row is scaled independently. If \code{'selected'}, the selected rows are scaled together. If \code{'all'}, all rows are scaled together. When \code{col.idp==TRUE}, the cell and bulk data are column-wise combined and scaled together. 
#' @param dim.lgd.pos The position of legends in the embedding plot. One of \code{top}, \code{right}, \code{bottom}, and \code{left}. 
#' @param dim.lgd.nrow The number of rows in the embedding plot legends.
#' @param dim.lgd.key.size,dim.lgd.text.size The key and font size of the embedding plot legends respectively. 
#' @param dim.axis.font.size The size of axis font of the embedding plot.
#' @param dim.lgd.plot.margin The margin of dimension reduction plot in the legend plot, e.g. `margin(t=0.01, r=0.15, b=0.01, l=0.15, unit="npc")`. 
#' @param size.lab.pt The size of point labels, applicable when \code{decon=TRUE}.
#' @param hjust.lab.pt,vjust.lab.pt Numeric values to adjust the horizontal and vertical positions of point labels respectively, applicable when \code{decon=TRUE}. 
#' @param col.com A vector of color components used to build the color scale. The default is `c('yellow', 'orange', 'red')`.
#' @param col.bar One of `selected` or `all`, the former uses expression values of \code{ID} to build the color scale while the latter uses all expression values from the data. The default is `selected`.
#' @param thr A two-numeric vector of expression value thresholds (the range of the color bar). The first and the second element will be the minmun and maximum threshold in the color bar respectively. Values above the max or below min will be assigned the same color as the max or min respectively. The default is \code{c(NA, NA)} and the min and max values in the data will be used. If one needs to change only max or min, the other should be \code{NA}.    
#' @param cores The number of CPU cores for parallelization. The default is `NA`, and the number of used cores is 1 or 2 depending on the availability. 
#' @param scale One of \code{no} (default), \code{selected}, \code{all}, or \code{row}, corresponding to no scaling, scaling selected rows as a whole, scaling all rows as a whole, or scaling each row independently.  
#' @param bar.width The width of color bar that ranges from 0 to 1. The default is 0.08.
#' @param legend.plot A vector of suffix(es) of aSVG file name(s) such as \code{c('shm1', 'shm2')}. Only aSVG(s) whose suffix(es) are assigned to this arugment will have a legend plot on the right. The default is \code{all} and each aSVG will have a legend plot. If NULL, no legend plot is shown.
#' @param legend.r A numeric (-1 to 1) to adjust the legend plot size.
#' @param legend.2nd Logical. If `TRUE`, the secondary legend is added to each SHM, which are the numeric values of each colored spatial features. The default its `FALSE`. Only applies to the static image.
#' @param lay.shm One of `gene`, `con`, or `none`. If `gene`, SHMs are organized horizontally by each biomolecule (gene, protein, or metabolite, \emph{etc.}) and variables are sorted under each biomolecule. If `con`, SHMs are organized horizontally by each experiment vairable and biomolecules are sorted under each variable. If `none`, SHMs are organized by the biomolecule order in \code{ID} and variables follow the order they appear in \code{data}. 
#' @param ncol The number of columns to display SHMs, which does not include the legend plot.
#' @param h The height (0-1) of color key and SHM/co-visualization plots in the middle, not including legend plots.
#' @param relative.scale A numeric to adjust the relative sizes between multiple aSVGs. Applicable only if multiple aSVGs are provided. Default is \code{NULL} and all aSVGs have the same size. 
#' @param verbose Logical. If `TRUE`, spatial features in data not colored in SHMs are printed to R console. Default is `TRUE`.
#' @param out.dir The directory to save SHMs as interactive HTML files and videos. Default is `NULL`, and the HTML files and videos are not saved.
#' @param animation.scale A numeric to scale the SHM size in the HTML files. The default is 1, and the height is 550px and the width is calculated according to the original aspect ratio in the aSVG file.
#' @param aspr The aspect ratio (width to height) in the interative HTML files.
#' @param video.dim A single character of the dimension of video frame in form of 'widthxheight', such as '1920x1080', '1280x800', '320x568', '1280x1024', '1280x720', '320x480', '480x360', '600x600', '800x600', '640x480' (default). The aspect ratio of SHMs are decided by \code{width} and \code{height}. 
#' @param interval The time interval (seconds) between SHM frames in the video. Default is 1.
#' @param framerate An integer of video framerate in frames per seconds. Default is 1. Larger values make the video smoother. 
#' @param res Resolution of the video in dpi.

#' @return A `SPHM` class.

#' @section Details:
#' See the package vignette (\code{browseVignettes('spatialHeatmap')}).  

#' @examples
#' 
#' # To always obtain the same results, a fixed seed is set.
#' set.seed(10); library(SummarizedExperiment); library(ggplot2)
#' # Read example single cell data of mouse brain (Marques et al. (2016)).
#' sce.pa <- system.file("extdata/shinyApp/data", "cell_mouse_brain.rds", 
#' package="spatialHeatmap")
#' # Cell group lables are stored in the "label" column in "colData".
#' sce <- readRDS(sce.pa); colData(sce)[1:3, 1:3]
#' # Read example bulk data of mouse brain (Vacher et al. 2021). 
#' blk.mus.pa <- system.file("extdata/shinyApp/data", "bulk_mouse_cocluster.rds", 
#' package="spatialHeatmap"); blk.mus <- readRDS(blk.mus.pa)
#' assay(blk.mus)[1:3, 1:5]
#'
#' # Jointly normalize bulk and single cell data.
#' cache.pa <- '~/.cache/shm' # Cache path.
#' mus.ann.nor <- read_cache(cache.pa, 'mus.ann.nor') 
#' if (is.null(mus.ann.nor)) {
#'   mus.lis.nor <- norm_cell(sce=sce, bulk=blk.mus, quick.clus=list(min.size = 100, d=15), 
#'   com=FALSE); save_cache(dir=cache.pa, overwrite=TRUE, mus.ann.nor)
#' }
#' # Separate normalized bulk and single cell data.
#' blk.mus.nor <- mus.lis.nor$bulk
#' cell.mus.nor <- mus.lis.nor$cell; colData(cell.mus.nor) <- colData(sce)
#' # Reduce dimensions (PCA, UMAP, TSNE) for single cell data.
#' cell.dim <- reduce_dim(cell.mus.nor, min.dim=5)
#'
#' # Aggregate tissue replicates in bulk data by mean.
#' blk.mus.aggr <- aggr_rep(blk.mus.nor, sam.factor='sample', aggr='mean')
#' assay(blk.mus.aggr)[1:2, ]
#' 
#' # Read the aSVG image into an "SVG" object.
#' svg.mus.brain.pa <- system.file("extdata/shinyApp/data", "mus_musculus.brain.svg", 
#' package="spatialHeatmap")
#' svg.mus.brain <- read_svg(svg.mus.brain.pa)
#' tail(attribute(svg.mus.brain)[[1]])[, 1:4]
#' # Map tissues to cell groups through a "list", where tissue labels are in name slots and
#' # the corresponding elements are cell group labels.
#' lis.match.blk <- list(cerebral.cortex=c('cortex.S1'), hypothalamus=c('corpus.callosum',
#' 'hypothalamus')) 
#' # Store bulk and single-cell data, aSVG, and match list in an "SHM" class.
#' dat.ann.tocell <- SPHM(svg=svg.mus.brain, bulk=blk.mus.aggr, cell=cell.dim, 
#' match=lis.match.blk)
#' # Co-visualization is created on expression values of gene "Eif5b". The target cell
#' # groups are "hypothalamus" and "cortex.S1".
#' covis(data=dat.ann.tocell, ID=c('Cacnb4'), col.idp=TRUE, dimred='TSNE', 
#' cell.group='label', tar.bulk=names(lis.match.blk), bar.width=0.09, dim.lgd.nrow=2, 
#' dim.lgd.text.size=10, h=0.6, legend.r=0.1, legend.key.size=0.01, legend.text.size=10,
#' legend.nrow=2, dim.lgd.plot.margin=margin(t=0.01, r=0.15, b=0.01, l=0.15, unit="npc"))
#'
#' # Save plots to HTML files and videos in the folder "test". 
#' # covis(data=dat.ann.tocell, ID=c('Cacnb4', 'Apod'), col.idp=TRUE, dimred='UMAP', 
#' # cell.group='label', tar.bulk=names(lis.match.blk), bar.width=0.12, bar.value.size=4, 
#' # dim.lgd.nrow=3, dim.lgd.text.size=7, h=0.5, legend.r=1, animation.scale=0.8, aspr=1.1,
#' # legend.key.size=0.01, legend.text.size=7, legend.nrow=3, legend.nrow.2nd=3, 
#' # dim.lgd.plot.margin=margin(t=0.01, r=0.15, b=0.01, l=0.15, unit="npc"), out.dir='test', 
#' # dim.lgd.key.size=2)

#' # More examples of co-visualization are provided in the package vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/spatialHeatmap/inst/doc/covisualize.html
#'
#'
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/
#' https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html
#' http://www.microugly.com/inkscape-quickguide/
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg
#' R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1
#' Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ 
#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra
#' R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/
#' https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg 
#' Yu, G., 2020. ggplotify:  Convert Plot to ’grob’ or ’ggplot’ Object. R package version 0.0.5.URLhttps://CRAN.R-project.org/package=ggplotify30
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' Guangchuang Yu (2020). ggplotify: Convert Plot to 'grob' or 'ggplot' Object. R package version 0.0.5. https://CRAN.R-project.org/package=ggplotify
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' Marques A et al. (2016). Oligodendrocyte heterogeneity in the mouse juvenile and adult central nervous system. Science 352(6291), 1326-1329.
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' Vacher, Claire-Marie, Helene Lacaille, Jiaqi J O’Reilly, Jacquelyn Salzbank, Dana Bakalar, Sonia Sebaoui, Philippe Liere, et al. 2021. “Placental Endocrine Function Shapes Cerebellar Development and Social Behavior.” Nat. Neurosci. 24 (10). Nature Publishing Group: 1392–1401

#' @name covis
#' @rdname covis
#' @aliases covis,SPHM-method

#' @export

setMethod("covis", c(data="SPHM"), function(data, assay.na=NULL, sam.factor=NULL, con.factor=NULL, ID, var.cell=NULL, dimred='PCA', cell.group=NULL, tar.cell=NULL, tar.bulk=NULL, size.pt=1, alpha.pt=0.8, shape=NULL, col.idp=FALSE, decon=FALSE, profile=TRUE, charcoal=FALSE, alpha.overlay=1, lay.shm="gene", ncol=2, h=0.99, size.r=1, col.com=c('yellow', 'orange', 'red'), col.bar='selected', thr=c(NA, NA), cores=NA, bar.width=0.08, bar.title=NULL, bar.title.size=0, scale='no', ft.trans=NULL, tis.trans=ft.trans, legend.r=0, sub.title.size=11, sub.title.vjust=3, legend.plot='all', ft.legend='identical', bar.value.size=10, legend.plot.title='Legend', legend.plot.title.size=11, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, angle.text.key=NULL, position.text.key=NULL, legend.2nd=FALSE, position.2nd='bottom', legend.nrow.2nd=2, legend.ncol.2nd=NULL, legend.key.size.2nd=0.03, legend.text.size.2nd=10, angle.text.key.2nd=0, position.text.key.2nd='right', dim.lgd.pos='bottom', dim.lgd.nrow=2, dim.lgd.key.size=4, dim.lgd.text.size=13, dim.axis.font.size=8, dim.lgd.plot.margin=NULL, dim.capt.size=13, size.lab.pt=5, hjust.lab.pt=0.5, vjust.lab.pt=1.5, add.feature.2nd=FALSE, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, line.width=0.2, line.color='grey70', relative.scale = NULL, verbose=TRUE, out.dir=NULL, animation.scale = 1, selfcontained=FALSE, aspr=1, video.dim='640x480', res=500, interval=1, framerate=1, bar.width.vdo=0.1, legend.value.vdo=NULL, ...) {
  sce.dimred <- data@cell; bulkCell <- NULL
  if (!is(sce.dimred, 'Seurat')) if ('bulkCell' %in% colnames(sce.dimred)) sce.dimred <- subset(sce.dimred, , bulkCell=='cell')
  res <- shm_covis(svg=data@svg, data=data@bulk, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, ID=ID, var.cell=var.cell, sce.dimred=sce.dimred, dimred=dimred, cell.group=cell.group, tar.cell=tar.cell, tar.bulk=tar.bulk, size.pt=size.pt, alpha.pt=alpha.pt, shape=shape, col.idp=col.idp, decon=decon, profile=profile, charcoal=charcoal, alpha.overlay=alpha.overlay, lay.shm=lay.shm, ncol=ncol, h=h, size.r=size.r, col.com=col.com, col.bar=col.bar, thr=thr, cores=cores, bar.width=bar.width, bar.title=bar.title, bar.title.size=bar.title.size, scale=scale, ft.trans=ft.trans, lis.rematch = data@match, legend.r=legend.r, sub.title.size=sub.title.size, sub.title.vjust=sub.title.vjust, legend.plot=legend.plot, ft.legend=ft.legend, bar.value.size=bar.value.size, legend.plot.title=legend.plot.title, legend.plot.title.size=legend.plot.title.size, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.text.size=legend.text.size, angle.text.key=angle.text.key, position.text.key=position.text.key, legend.2nd=legend.2nd, position.2nd=position.2nd, legend.nrow.2nd=legend.nrow.2nd, legend.ncol.2nd=legend.ncol.2nd, legend.key.size.2nd=legend.key.size.2nd, legend.text.size.2nd=legend.text.size.2nd, angle.text.key.2nd=angle.text.key.2nd, position.text.key.2nd=position.text.key.2nd, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.axis.font.size=dim.axis.font.size, dim.lgd.plot.margin=dim.lgd.plot.margin, dim.capt.size=dim.capt.size, size.lab.pt=size.lab.pt, hjust.lab.pt=hjust.lab.pt, vjust.lab.pt=vjust.lab.pt, add.feature.2nd=add.feature.2nd, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, line.width=line.width, line.color=line.color, relative.scale = relative.scale, verbose=verbose, out.dir=out.dir, animation.scale = animation.scale, aspr=aspr, selfcontained=selfcontained, video.dim=video.dim, res=res, interval=interval, framerate=framerate, bar.width.vdo=bar.width.vdo, legend.value.vdo=legend.value.vdo, ...)
  data@output <- res; invisible(data)
})


#' SHM plots and co-visaulization plots 
#' 
#' @param sce.dimred A \code{SingleCellExperiment} with reduced dimentions such as \code{PCA}, \code{UMAP}, \code{TSNE}.
#' @param lis.rematch 
#' \describe{ 
#'  \item{SHMs}{ A named \code{list} for rematching spatial features between numeric data (ftA, ftB) and aSVGs (ftC, ftD, ftE). In each slot, the slot name is a spatial feature from the data and the corresponding element is one or multiple spatial features from the aSVG. \emph{E.g.} \code{list(ftA = c('ftC', 'ftD'), ftB = c('ftE'))}.
#'  } 
#'  \item{Co-visualization plots}{ Mapping cells to tissues: a named \code{list}, where cell group labels from \code{colData(sce.dimred)[, 'cell.group']} are the name slots and aSVG features are the corresponding \code{list} elements. Mapping tissues to cells: a named \code{list}, where tissues are the name slots and cells from \code{colData(sce.dimred)[, 'cell.group']} are the corresponding \code{list} elements. Applicable when cell grouping methods are annodation labels, marker genes, clustering, or manual assignments. 
#'  }
#' }
#' @keywords Internal
#' @noRd

#' @references
#' Hao and Hao et al. Integrated analysis of multimodal single-cell data. Cell (2021) [Seurat V4]
#' Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019) [Seurat V3]
#' Butler et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat Biotechnol (2018) [Seurat V2]
#' Satija and Farrell et al. Spatial reconstruction of single-cell gene expression data. Nat Biotechnol (2015) [Seurat V1]

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

shm_covis <- function(svg, data, assay.na=NULL, sam.factor=NULL, con.factor=NULL, ID, var.cell=NULL, sce.dimred=NULL, dimred='PCA', cell.group=NULL, tar.cell=NULL, tar.bulk=NULL, size.pt=1, alpha.pt=0.8, shape=NULL, col.idp=FALSE, decon=FALSE, profile=TRUE, charcoal=FALSE, alpha.overlay=1, lay.shm="gene", ncol=2, h=0.99, size.r=1, col.com=c('yellow', 'orange', 'red'), col.bar='selected', thr=c(NA, NA), cores=NA, bar.width=0.08, bar.title=NULL, bar.title.size=0, scale='no', ft.trans=NULL, tis.trans=ft.trans, lis.rematch = NULL, legend.r=0, sub.title.size=11, sub.title.vjust=3, legend.plot='all', ft.legend='identical', bar.value.size=10, legend.plot.title='Legend', legend.plot.title.size=11, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.02, legend.text.size=12, angle.text.key=NULL, position.text.key=NULL, legend.2nd=FALSE, position.2nd='bottom', legend.nrow.2nd=2, legend.ncol.2nd=NULL, legend.key.size.2nd=0.03, legend.text.size.2nd=10, angle.text.key.2nd=0, position.text.key.2nd='right', dim.lgd.pos='bottom', dim.lgd.nrow=2, dim.lgd.key.size=4, dim.lgd.text.size=13, dim.axis.font.size=10, dim.lgd.plot.margin=margin(t=0.01, r=0.01, b=0.01, l=0.01, unit="npc"), dim.capt.size=13, size.lab.pt=5, hjust.lab.pt=0.5, vjust.lab.pt=1.5, add.feature.2nd=FALSE, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, line.width=0.2, line.color='grey70', relative.scale = NULL, verbose=TRUE, out.dir=NULL, animation.scale = 1, aspr=1, selfcontained=FALSE, video.dim='640x480', res=500, interval=1, framerate=1, bar.width.vdo=0.1, legend.value.vdo=NULL, ...) {
  
  calls <- names(vapply(match.call(), deparse, character(1))[-1])
  if("tis.trans" %in% calls) warning('"tis.trans" is deprecated and replaced by "ft.trans"! \n')
  if("svg.path" %in% calls) warning('"svg.path" is deprecated and replaced by "svg"! \n')
  # save(svg, data, assay.na, sam.factor, con.factor, ID, sce.dimred, var.cell, dimred, col.idp, decon, tar.cell, profile, cell.group, tar.bulk, size.pt, alpha.pt, shape, charcoal, alpha.overlay, lay.shm, ncol, h, size.r, col.com, col.bar, thr, cores, bar.width, bar.title, bar.title.size, scale, ft.trans, tis.trans, lis.rematch, legend.r, sub.title.size, sub.title.vjust, legend.plot, ft.legend, bar.value.size, legend.plot.title, legend.plot.title.size, legend.ncol, legend.nrow, legend.position, legend.direction, legend.key.size, legend.text.size, angle.text.key, position.text.key, legend.2nd, position.2nd, legend.nrow.2nd, legend.ncol.2nd, legend.key.size.2nd, legend.text.size.2nd, angle.text.key.2nd, position.text.key.2nd, dim.lgd.pos, dim.lgd.nrow, dim.lgd.key.size, dim.lgd.text.size, dim.axis.font.size, dim.lgd.plot.margin, dim.capt.size, add.feature.2nd, label, label.size, label.angle, hjust, vjust, opacity, key, line.width, line.color, relative.scale, verbose, out.dir, animation.scale, aspr, selfcontained, video.dim, res, interval, framerate, bar.width.vdo, legend.value.vdo, file='shm.covis.arg')

  x <- y <- color_scale <- tissue <- bulkCell <- variable <- NULL
  options(stringsAsFactors=FALSE)
  # if (!is.null(sub.margin)) if (!is.numeric(sub.margin) | length(sub.margin)!=4 | any(sub.margin >= 1) | any(sub.margin < 0)) stop('"sub.margin" must be a 4-length numeric vector between 0 (inclusive) and 1 (exclusive)!')
  if (!is(svg, 'SVG')) stop('The "svg" should be a "SVG" object!')
  svgs <- svg; if (missing(ID)) ID <- rownames(data)[1]
  ID <- unique(ID)
  svg.ft.all <- unique(unlist(lapply(seq_along(svgs), function(x) { unique(attribute(svgs[x])[[1]]$feature) })))
  srsc <- is(sce.dimred, 'Seurat')
  sc <- is(sce.dimred, 'SingleCellExperiment') & !srsc & !decon
  profile.no <- profile==FALSE & (sc | srsc | decon)
  if (srsc & decon) stop("If the spatially resolved single-cell data are provided, 'decon' should be FALSE!")
  # Original lis.rematch used in 'toCell'.
  covis.type <- NULL; lis.rematch.raw <- lis.rematch
  if (sc) { # Check coviz related data.
    lis <- check_data_covis(data, dimred, cell.group, sce.dimred, tar.cell, tar.bulk, lis.rematch.raw, col.idp=col.idp, svg.ft.all) 
    tar.bulk <- lis$tar.bulk; tar.cell <- lis$tar.cell
    covis.type <- lis$covis.type; lis.rematch <- lis$lis.rematch
    bulk.all <- lis$bulk.all; cell.all <- lis$cell.all
    cell.group <- lis$cell.group; rm(lis)
  }
  ft.trans.shm <- NULL; if (sc) {
      # The only use of tar.cell, tar.bulk is to get ft.trans.shm, they are not required in dim_color, dim_color2cell, and dim_color_coclus, since the color transfer is data -> SHM -> dim plot and transparent colors are defined in SHM before dim plot.
      lis.covis <- covis_trans(bulk.all=bulk.all, cell.all=cell.all, ft.all=svg.ft.all, tar.bulk=tar.bulk, tar.cell=tar.cell, lis.match=lis.rematch.raw, col.idp=col.idp, covis.type=covis.type)
      tar.bulk <- lis.covis$tar.bulk
      tar.cell <- lis.covis$tar.cell
      lis.rematch <- lis.covis$lis.match 
      ft.trans.shm <- lis.covis$ft.trans.shm
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
    dat.lis <- check_data(data=gene, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, usage='shm')
    gene <- as.data.frame(dat.lis$dat)
  } else if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')|is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) {
    id.no <- ID[!ID %in% rownames(data)]
    if (length(id.no)>0) stop(paste0(id.no, collapse=' '), ': not detected in data! \n')
    # Process data.
    dat.lis <- check_data(data=data, assay.na=assay.na, sam.factor=sam.factor, con.factor=con.factor, usage='shm')
    gene <- as.data.frame(dat.lis$dat); con.na <- dat.lis$con.na
  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", "dgCMatrix", "SummarizedExperiment", or "SingleCellExperiment" except that "spatial_hm" also accepts a "vector".') }
    # Only take the column names with "__".
    # cname <- colnames(gene); form <- grepl('__', cname)
    # con <- gsub("(.*)(__)(.*)", "\\3", cname[form])
    # con <- gsub("(.*)(__)(.*)", "\\3", cname)
    # sam <- gsub("(.*)(__)(.*)", "\\1", cname)
    se <- dat.lis$se; sam.uni <- unique(se$spFeature)
    con.uni <- unique(se$variable)

    con.na.cell <- NULL; if (col.idp==TRUE & sc) {
      sce.lis <- check_sce(sce=sce.dimred, cell.group, var.cell)
      if (is(sce.lis, 'character')) return(sce.lis)
      sce.dimred <- sce.lis$sce; var.cell <- sce.lis$var.cell
      con.na.cell <- sce.lis$con.na.cell
      inter <- intersect(rownames(gene), rownames(sce.dimred))
      gene <- gene[inter, , drop=FALSE]
      sce.dimred <- sce.dimred[inter, ]
      if (length(inter)==0) stop('Bulk and single cell data do not have overlap biomolecules!')
      vars.cell <- unique(colData(sce.dimred)[, var.cell])
      if (!identical(sort(con.uni), sort(vars.cell))) stop('Variables between bulk and single cell data should be the same!')
      gene <- cbind(as.matrix(gene), as.matrix(assay(sce.dimred)))
    }
    if (srsc) {
      lis <- check_srsc(cell=sce.dimred, cell.group=cell.group, var.cell=var.cell, bulk=gene, var.bulk=con.uni)
      if (is(lis, 'character')) return(lis)
      cell <- lis$cell; var.cell <- lis$var.cell
      pkg <- check_pkg('Seurat')
      if (is(pkg, 'character')) stop(pkg)
      sce.sp <- srt2sce(cell=cell, assay=Seurat::DefaultAssay(cell), image=1, x='imagerow', y='imagecol', cell.group=cell.group)
      con.na.cell <- lis$con.na.cell; gene <- lis$bulk
      gene <- cbind(as.matrix(gene), as.matrix(assay(sce.sp)))
    }
    if (decon & !srsc) {
      lis <- check_decon(cell=sce.dimred, cell.group=cell.group, tar.cell=tar.cell, var.cell=var.cell, bulk=gene, var.bulk=con.uni, svg.ft.all=svg.ft.all)
      if (is(lis, 'character')) return(lis)
      sce.dimred <- lis$cell; var.cell <- lis$var.cell
      vars.cell <- lis$vars.cell
      con.na.cell <- lis$con.na.cell; gene <- lis$bulk
      gene <- cbind(as.matrix(gene), as.matrix(assay(sce.dimred)))
    }
    # Color bar.
    thr <- dat_fun(thr, as.numeric)
    if (length(thr)!=2) stop('The "thr" must be a two-element vecor and contain as least one numeric!')
    # Customized signal threshold.
    gene <- thrsd(thr.min=thr[1], thr.max=thr[2], data=gene)
    if (is(gene, 'character')) stop(gene)
    # Scaling. 
    if ('row' %in% scale) { 
      gene <- t(scale(t(gene)))   
    } else if ('all' %in% scale) { 
      gene <- scale_all(gene) 
    } else if ('selected' %in% scale) {
      gene[ID, ] <- scale_all(gene[ID, , drop=FALSE]) 
    }
    if ('selected' %in% col.bar) { 
      min.v <- min(gene[ID, , drop=FALSE])
      max.v <- max(gene[ID, , drop=FALSE]) 
    } else if ('all' %in% col.bar) { 
      min.v <- min(gene); max.v <- max(gene) 
    }; bar.len=1000; geneV <- seq(min.v, max.v, len=bar.len)

    col <- colorRampPalette(col.com)(length(geneV))
    cs.g <- col_bar(geneV=geneV, cols=col, width=1, title=bar.title, bar.title.size=bar.title.size, bar.value.size=bar.value.size)

    svg.pa.na <- img_pa_na(unlist(svgs[, 'svg']))
    # svg.path may not be paths, can be file names, if users provide a SVG class that not includes paths.
    svg.path <- svg.pa.na$path; svg.na <- svg.pa.na$na
    svg.na.cord <- names(svgs)
    # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
    w.h.max <- max(unlist(svgs[, 'dimension']))
   
    # A set of SHMs (gene_con*) are made for each SVG, and all sets of SHMs under different SVGs are placed in 2 lists in form of ggplots and grobs respectively. Different SHMs of same 'gene_condition' under different SVGs are indexed with suffixed of '_1', '_2', ... E.g. SHMs under shm1.svg: gene_condition1_1, gene_condition2_1; SHMs under shm2.svg: gene_condition1_2, gene_condition2_2; the 4 SHMs are stored in 2 separate lists in form of ggplots and grobs respectively. 
    # The order of ggplot SHMs, grob SHMs, and legends follow the order of SVGs, so all are same.
    gg.lis.all <- gg.all <- grob.all <- lgd.all <- grob.lgd.all <- gcol.lgd.all <- gcol.all <- NULL
    for (i in seq_along(svgs)) {
      na0 <- svg.na[i]; cat('ggplots/grobs:', na0, '... \n')
      # if (preserve.scale==TRUE & is.null(sub.margin)) sub.margin <- (1-w.h/w.h.max*0.99)/2
      # SHMs/legend of ggplots under one SVG.
      gg.lis <- gg_shm(svg.all=svgs[i], gene=gene[, seq_len(ncol(se)), drop=FALSE], col.idp=col.idp, con.na=con.na, geneV=geneV, charcoal=charcoal, alpha.overlay=alpha.overlay, ID=ID, cols=col, covis.type=covis.type, lis.rematch = lis.rematch, ft.trans=ft.trans, ft.trans.shm=ft.trans.shm, sub.title.size=sub.title.size, sub.title.vjust=sub.title.vjust, ft.legend=ft.legend, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.text.size=legend.text.size, legend.plot.title=legend.plot.title, legend.plot.title.size=legend.plot.title.size, line.width=line.width, line.color=line.color, ...)
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
    if (!is(h, 'numeric')) stop('"h" should be between 0 and 1!')
    if (h > 1) h <- 1 else if (h < 0) h <- 0
    cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc"), heights=unit(0.8 * h, "npc")) # "mm" is fixed, "npc" is scalable.
    if (legend.2nd==TRUE) {
      gg.all <- gg_2lgd(gg.all=gg.all, sam.dat=sam.uni, ft.trans=ft.trans, position.2nd=position.2nd, legend.nrow.2nd=legend.nrow.2nd, legend.ncol.2nd=legend.ncol.2nd, legend.key.size.2nd=legend.key.size.2nd, add.feature.2nd=add.feature.2nd, legend.text.size.2nd=legend.text.size.2nd, angle.text.key.2nd=angle.text.key.2nd, position.text.key.2nd=position.text.key.2nd)
      grob.all <- lapply(gg.all, ggplotGrob)
    }

    # Arrange SHM grobs.
    na.all <- sort_gen_con(ID.sel=ID, na.all=na.all, con.all=con.idxed, by=lay.shm)
    grob.all <- grob.all[na.all]; gg.all <- gg.all[na.all]
    gcol.all <- gcol.all[paste0('col_', na.all)]
    dim.lgd.lis <- NULL; if (sc) {
      if (FALSE %in% col.idp) if ('toBulk' %in% covis.type) {
      gg.dim <- plot_dim(sce.dimred, dim=dimred, color.by=cell.group)
      dim.shm.lis <- dim_color(gg.dim=gg.dim, gg.shm.all=gg.all, grob.shm.all=grob.all, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=grob.lgd.all, cell.group=cell.group, tar.cell=tar.cell, con.na=con.na, profile=profile, alpha.pt=alpha.pt, lis.match=lis.rematch, sub.title.size=sub.title.size, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.axis.font.size=dim.axis.font.size, shape=shape)
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
      dim.shm.lis <- dim_color_coclus(sce=sce.dimred, targ=targ, profile=profile, alpha.pt=alpha.pt, gg.dim = gg.dim, gg.shm.all=gg.all, grob.shm.all = grob.all, gg.lgd.all=lgd.all, col.shm.all = gcol.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=grob.lgd.all, con.na=con.na, lis.match=NULL, sub.title.size=sub.title.size, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.capt.size=dim.capt.size, dim.axis.font.size=dim.axis.font.size, shape=shape)
      } else if ('toCell' %in% covis.type) {
      gg.dim <- plot_dim(sce.dimred, dim=dimred, color.by=cell.group)
      # Tar bulk is in SHM already.
      dim.shm.lis <- dim_color2cell(gg.dim=gg.dim, gg.shm.all=gg.all, grob.shm.all=grob.all, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=grob.lgd.all, profile=profile, alpha.pt=alpha.pt, cell.group=cell.group, con.na=con.na, lis.match=lis.rematch.raw, sub.title.size=sub.title.size, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.axis.font.size=dim.axis.font.size, shape=shape)
      }
      if (TRUE %in% col.idp & !srsc) {
      if (any(c('toBulkAuto', 'toCellAuto') %in% covis.type)) {
      # sce.cell <- subset(sce.dimred, , bulkCell=='cell')
      dimred.pre <- unique(colData(data)$dimred)
      if (!is.null(dimred.pre)) {
        dimred <- dimred.pre
        cat('The reduced dimensionality in "data" is used:', dimred, '.\n')
      }
      gg.dim <- plot_dim(sce.dimred, dim=dimred, color.by='assignedBulk')
      if ('toBulkAuto' %in% covis.type) targ <- tar.cell
      if ('toCellAuto' %in% covis.type) targ <- tar.bulk
      gg.dim.lis <- list(con=gg.dim)
      } else {
      vars.cell <- unique(sce.lis$sce$variable)
      gg.dim.lis <- NULL; for (i in vars.cell) {
        sce.dimred0 <- subset(sce.dimred, , variable==i) 
        gg.dim <- plot_dim(sce.dimred0, dim=dimred, color.by=cell.group)
        gg.dim.lis <- c(gg.dim.lis, list(gg.dim))
      }; names(gg.dim.lis) <- vars.cell; targ <- NULL
      }
      dim.shm.lis <- dim_color_idp(sce=sce.dimred, covis.type=covis.type, targ=targ, ID=ID, gene=gene[, setdiff(seq_len(ncol(gene)), seq_len(ncol(se))), drop=FALSE], tar.cell=tar.cell, tar.bulk=tar.bulk, con.na.cell=con.na.cell, geneV=geneV, cols=col, gg.dim=gg.dim.lis, gg.shm.all=gg.all, grob.shm.all=grob.all, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=grob.lgd.all, profile=profile, alpha.pt=alpha.pt, cell.group=cell.group, lis.match=lis.rematch.raw, sub.title.size=sub.title.size, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.axis.font.size=dim.axis.font.size, shape=shape, lgd.plot.margin=dim.lgd.plot.margin)
      }
      grob.all <- dim.shm.lis$dim.shm.grob.lis
      gg.all <- dim.shm.lis$dim.shm.gg.lis
      dim.lgd.lis <- dim.shm.lis$dim.lgd.lis
    }
    if (srsc) {
      dim.ovl.shm.lis <- dim_srsc(cell=sce.sp, ID=ID, geneV=geneV, cols=col, tar.cell=tar.cell, con.na.cell=con.na.cell, dimred=dimred, size.pt=size.pt, alpha.pt=alpha.pt, svg.all=svgs[1], profile=profile, gg.shm.all=gg.all, grob.shm.all=grob.all, gg.lgd.all=lgd.all, grob.lgd.all=grob.lgd.all, cell.group=cell.group, sub.title.size=sub.title.size, sub.title.vjust=0, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.axis.font.size=dim.axis.font.size, linewidth=line.width, line.color=line.color, dim.lgd.direc=NULL, size.r=size.r, shape=shape)
      grob.all <- dim.ovl.shm.lis$dim.ovl.shm.grob
      gg.all <- dim.ovl.shm.lis$dim.ovl.shm.gg
    }
    if (decon) {
      dim.shm.lis <- dim_decon(sce=sce.dimred, ID=ID, data=gene[, setdiff(seq_len(ncol(gene)), seq_len(ncol(se))), drop=FALSE], geneV=geneV, cols=col, tar.cell=tar.cell, con.na.cell=con.na.cell, dimred=dimred, profile=profile, size.pt=size.pt, alpha.pt=alpha.pt, gg.shm.all=gg.all, grob.shm.all=grob.all, col.shm.all=gcol.all, gg.lgd.all=lgd.all, col.lgd.all=gcol.lgd.all, grob.lgd.all=grob.lgd.all, cell.group=cell.group, sub.title.size=sub.title.size, dim.lgd.pos=dim.lgd.pos, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.text.size=dim.lgd.text.size, dim.axis.font.size=dim.axis.font.size, size.lab.pt=size.lab.pt, hjust.lab.pt=hjust.lab.pt, vjust.lab.pt=vjust.lab.pt, shape=shape)
      grob.all <- dim.shm.lis$dim.shm.grob
      gg.all <- dim.shm.lis$dim.shm.gg
      dim.lgd.lis <- dim.shm.lis$dim.lgd.lis
    }
    g.arr <- lay_shm(lay.shm=lay.shm, con=con.idxed, ncol=ncol, ID.sel=ID, grob.list=grob.all, scell=sc, decon=decon, srsc=srsc, profile=profile, h=h, shiny=FALSE)
    if (profile.no) legend.plot <- NULL
    if (!is.null(legend.plot)) {
      # Select legend plot to show. 
      if (length(svg.na)>1 & legend.plot!='all') na.lgd <- svg.na[grep(paste0('_(', paste0(legend.plot, collapse='|'), ').svg$'), svg.na)] else na.lgd <- svg.na
      lgd.lis <- lgd.all[na.lgd]
 
      # Add labels to target shapes/adjust keys in legend plots.
      lgd.lis <- gg_lgd(gg.all=lgd.lis, angle.text.key=angle.text.key, position.text.key=position.text.key, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key)
      # Include dim plots in legends: covis color by cell.
      if (!is.null(dim.lgd.lis) & (col.idp==TRUE | decon==TRUE)) {
        for (i in vars.cell) {
          dim.lgd0 <- dim.lgd.lis[[i]]$dim.lgd
          lgd.lis <- c(lgd.lis, setNames(list(dim.lgd0), paste0(i, '_dim.lgd')))
        }
      }; na.lgd <- names(lgd.lis)

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
      message('HTML files ... ') 
      html_ly(gg.all=gg.all, cs.g=cs.g, aspr=aspr, anm.scale=animation.scale, out.dir=out.dir)
      message('Video ...')
      if (!sc) type <- 'shm' else {
        if (FALSE %in% col.idp) type <- 'col.grp' else type <- 'col.idp'
      } 
      vdo <- video(gg=gg.all, cs.g=cs.g, lgd=lgd.lis, lgd.r=legend.r, h=h, type=type, sub.title.size=sub.title.size, bar.width=bar.width, bar.value.size=bar.value.size, lgd.key.size=legend.key.size, lgd.text.size=legend.text.size, lgd.key.size.2nd=legend.key.size.2nd, lgd.text.size.2nd=legend.text.size.2nd, angle.text.key=angle.text.key.2nd, position.text.key=position.text.key.2nd, lgd.row=legend.nrow, lgd.row.2nd=legend.nrow.2nd, lgd.col=legend.ncol.2nd, legend.value.vdo=legend.value.vdo, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, dim.lgd.text.size=dim.lgd.text.size, dim.lgd.key.size=dim.lgd.key.size, dim.lgd.nrow=dim.lgd.nrow, dim.lgd.plot.margin=dim.lgd.plot.margin, video.dim=video.dim, res=res, interval=interval, framerate=framerate, out.dir=out.dir)
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
    map.gene <- data[, idx.com, drop=FALSE]
    cna <- colnames(map.gene)
    for (i in ID) {
      df0 <- as.data.frame(t(map.gene[i, , drop=FALSE]))
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

check_data_covis <- function(data, dimred, cell.group, sce.dimred, tar.cell, tar.bulk, lis.rematch, col.idp=FALSE, svg.ft.all) {
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
    bulk.all <- unique(data$sample); cell.group <- 'assignedBulk'
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
    bulk.all <- unique(sub('__.*', '', colnames(data)))
    cell.all <- unique(cdat.cell[, cell.group])
    if (FALSE %in% col.idp) lis.rematch <- NULL
    covis.type <- 'toCell'
  } else if (!is.null(lis.rematch) & !is.null(tar.cell)) {
    bulk.all <- unique(sub('__.*', '', colnames(data)))
    cell.all <- unique(cdat.cell[, cell.group])
    covis.type <- 'toBulk'
  }
  return(list(bulk.all=bulk.all, cell.all=cell.all, tar.bulk=unique(tar.bulk), tar.cell=unique(tar.cell), covis.type=covis.type, lis.rematch=lis.rematch, cell.group=cell.group))
}


#' Extract transparent aSVG features in co-visualization

#' @keywords Internal
#' @noRd

covis_trans <- function(bulk.all, cell.all, ft.all, tar.bulk='all', tar.cell='all', lis.match, col.idp=FALSE, covis.type) {
 # save(bulk.all, ft.all, cell.all, tar.bulk, tar.cell, lis.match, covis.type, file='covis.trans.arg')
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
        if (FALSE %in% col.idp) lis.match <- NULL
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
  
#' Thresholding assay matrix

#' @keywords Internal
#' @noRd

thrsd <- function(thr.min, thr.max, data) {
  max.v <- max(data); min.v <- min(data) 
  if (is.numeric(thr.min) & !is.na(thr.min)) { 
    lgc.min <- (thr.min < max.v) 
    if (!lgc.min) {
      msg <- paste0('The min threshold should be less than the max value in the assay data ', max.v, '!'); return(msg)
    }
    if (thr.min > min.v) min.v <- thr.min
  }
  if (is.numeric(thr.max) & !is.na(thr.max)) { 
    lgc.max <- (thr.max > min.v)
    if (!lgc.max) {
      msg <- paste0('The max threshold should be larger than the min value in the assay data ', min.v, '!'); return(msg)
    }
    if (thr.max < max.v) max.v <- thr.max 
  }
  if (max.v > min.v) {  
    data[data >= max.v] <- max.v; data[data <= min.v] <- min.v  
  }; return(data)
}


#' When coloring cells and tissues independently, check validity of single cell data.

#' @keywords Internal
#' @noRd

#' @importFrom SummarizedExperiment colData colData<-

check_sce <- function(sce, cell.group, var.cell) { 
  lgc.cna <- length(grep('__', colnames(sce)))==0
  cdat <- colData(sce)
  if (lgc.cna & is.null(cell.group)) { 
    msg <- 'Please specify "cell.group" and optionally "var.cell"!'
    warning(msg); return(msg)
  }
  if (lgc.cna) {
    con.na.cell <- ifelse(!is.null(var.cell), TRUE, FALSE)
    if (!is.null(cell.group) & is.null(var.cell)) { 
      cdat$variable <- 'con'; var.cell <- 'variable'
    } 
    rownames(cdat) <- colnames(sce) <- paste0(cdat[, cell.group], '__', cdat[, var.cell])
    colData(sce) <- cdat
  } else { con.na.cell <- TRUE
    cells <- gsub('(.*)__(.*)', '\\1', colnames(sce))  
    varis <- gsub('(.*)__(.*)', '\\2', colnames(sce)) 
    if (is.null(cell.group)) { 
      cdat$cell <- cells; cell.group <- 'cell'
    }
    if (is.null(var.cell)) { 
      cdat$variable <- varis; var.cell <- 'variable'
    }; colData(sce) <- cdat 
  }
  # if (length(intersect(sp.ft, cdat[, cell.group]))>0) stop('Same identifiers are detected between spatial features in bulk data and single cells!')
  return(list(sce=sce, con.na.cell=con.na.cell, var.cell=var.cell))
}

#' Check spatially resolved single cell data.

#' @keywords Internal
#' @noRd

check_srsc <- function(cell, cell.group, var.cell, bulk, var.bulk) {
  pkg <- check_pkg('Seurat'); if (is(pkg, 'character')) stop(pkg)
  cell <- Seurat::RenameCells(cell, new.names = make.names(colnames(cell)))
  lgc.cna <- length(grep('__', colnames(cell)))==0
  if (lgc.cna & is.null(cell.group)) stop('Please specify "cell.group" and optionally "var.cell"!')
  if (lgc.cna) {
    con.na.cell <- ifelse(!is.null(var.cell), TRUE, FALSE)
    if (!is.null(cell.group) & is.null(var.cell)) { 
      cell[['variable']] <- 'con'; var.cell <- 'variable'
    } 
    cell <- Seurat::RenameCells(cell, new.names = paste0(colnames(cell), '__', cell[[var.cell]][, 1]))
  } else { con.na.cell <- TRUE
    cells <- gsub('(.*)__(.*)', '\\1', colnames(cell))  
    varis <- gsub('(.*)__(.*)', '\\2', colnames(cell)) 
    if (is.null(cell.group)) { 
      cell[['cell']] <- cells; cell.group <- 'cell'
    }
    if (is.null(var.cell)) { 
      cell[['variable']] <- varis; var.cell <- 'variable'
    }
  }
  inter <- intersect(rownames(cell), rownames(bulk))
  if (length(inter)==0) {
    msg <- 'Bulk and single cell data do not have overlap biomolecules!'
    warning(msg); return(msg)
  }
  bulk <- bulk[inter, , drop=FALSE]
  cell <- subset(cell, features=inter) 
  vars.cell <- unique(cell[[var.cell]][, 1])
  if (!identical(sort(var.bulk), sort(vars.cell))) { 
    msg <- 'Variables between bulk and single cell data should be the same!'
    warning(msg); return(msg)
  }
  return(list(cell=cell, bulk=bulk, con.na.cell=con.na.cell, var.cell=var.cell))
}


#' Check single-cell data of deconvolution.

#' @keywords Internal
#' @noRd

check_decon <- function(cell, cell.group, tar.cell, var.cell, bulk, var.bulk, svg.ft.all=svg.ft.all) {
  cdat <- colData(cell)
  if (!'bulk' %in% colnames(cdat)) {
    msg <- "When 'decon' is 'TRUE', a 'bulk' column is required in the 'colData' slot of 'SingleCellExperiment'!"
    warning(msg); return(msg)
  } else {
    if (!any(cdat$bulk %in% svg.ft.all)) {
      msg <- 'At least one bulk tissue in single-cell data (colData) should be the same with an aSVG feature!'
      warning(msg); return()
    }
  }
  if (is.null(tar.cell) | 'all' %in% tar.cell) tar.cell <- unique(cdat[, cell.group]) 
  lgc.cna <- length(grep('__', colnames(cell)))==0
  if (lgc.cna & is.null(cell.group)) stop('Please specify "cell.group" and optionally "var.cell"!')
  if (lgc.cna) {
    con.na.cell <- ifelse(!is.null(var.cell), TRUE, FALSE)
    if (!is.null(cell.group) & is.null(var.cell)) { 
      cell[['variable']] <- 'con'; var.cell <- 'variable'
    }
    colnames(cell) <- paste0(cell[[cell.group]], '__', cell[[var.cell]])
  } else { con.na.cell <- TRUE
    cells <- gsub('(.*)__(.*)', '\\1', colnames(cell))  
    varis <- gsub('(.*)__(.*)', '\\2', colnames(cell)) 
    if (is.null(cell.group)) { 
      cell[['cell']] <- cells; cell.group <- 'cell'
    }
    if (is.null(var.cell)) { 
      cell[['variable']] <- varis; var.cell <- 'variable'
    }
  }
  inter <- intersect(rownames(cell), rownames(bulk))
  if (length(inter)==0) {
    msg <- 'Bulk and single cell data do not have overlap biomolecules!'
    warning(msg); return(msg)
  }
  bulk <- bulk[inter, , drop=FALSE]
  cell <- subset(cell, features=inter) 
  vars.cell <- unique(cell[[var.cell]])
  if (!identical(sort(var.bulk), sort(vars.cell))) { 
    msg <- 'Variables between bulk and single cell data should be the same!'
    warning(msg); return(msg)
  }
  return(list(cell=cell, bulk=bulk, con.na.cell=con.na.cell, var.cell=var.cell, vars.cell=vars.cell))
}


