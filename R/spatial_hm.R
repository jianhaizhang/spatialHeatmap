#' Plot Spatial Heatmap
#'
#'It takes a pair of formatted gene expression data matrix ("SummarizedExperiment") and SVG image as input. In the SVG image, tissue regions are pre-defined and linked with data matrix through tissue names. The expression profiles of the input gene(s) under each condition are mapped to each tissue in the form of different colours. The application is not limited to gene expression data. It is applicable as long as a pair of formatted data matrix and SVG image are provided, such as population data generated in different years across different cities. 

#' @param svg.path The path of the SVG image. Target tissues to be coloured in spaial heatmaps must have identical tissue names with corresponding tissues in the data matrix. \cr E.g.: system.file("extdata/example", "homo_sapiens.brain.svg", package = "spatialHeatmap") (https://www.ebi.ac.uk/gxa/home).

#' @inheritParams filter_data
#' @inheritParams grob_list

#' @param ID A character of gene ID(s) whose expression values are used to colour the spatial heatmaps. It can be a single gene or a vector of multiple genes.

#' @param col.com A character vector of the colour components used to build the colour scale, e.g. the default is c("yellow", "blue", "purple").

#' @param col.bar "selected" or "all", meaning use input genes or whole data matrix to build the colour scale respectively. The default is "selected".
#' @param data.trans "log2", "exp2", or NULL. If colours across tissues cannot distinguish due to low variance or outliers, transform the data matrix by log2 or 2-base expoent (exp2). Default is NULL (data will not be transformed).
#' @param bar.width The width of colour bar. Default if 0.7.
#' @param width A numeric of each subplot width, relative to height. The default is 1.
#' @param height A numeric of each subplot height, relative to width. The default is 1.
#' @param legend.r The ratio of height to width of the legend plot. Default is 1.
#' @param lay.shm "gene" or "con" (condition), the organisation of the subplots.
#' @param ncol Number of columns to display the subplots.
#' @return It generates an image of spatial heatmap(s) along with a colour key.

#' @section Details:
#' Details about how to format an SVG image and a data matrix are provided in the package vignette (\code{browseVignettes('spatialHeatmap')}). The "se" parameter can be the value returned by the function \code{\link{filter_data}}, or built on an expression matrix and a targets file. See examples blow.

#' @examples
#' # The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
#' library(ExpressionAtlas)
#' rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]; assay(rse.hum)[1:3, 1:3]
#'
#' # A targets file describing replicates of samples and conditions is required, which is made based on the "colData" slot in the downloaded "RangedSummarizedExperiment" and available in spatialHeatmap. See the "se" parameter for details. 
#' brain.pa <- system.file('extdata/shinyApp/example/target_brain.txt', package='spatialHeatmap')
#' target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
#' # The "organism_part" and "disease" column describes tissue and condition replicates respectively.  
#' target.hum[c(1:3, 41:42), 4:5]
#' # Place the targets file into "colData" slot as a DataFrame class. 
#' colData(rse.hum) <- DataFrame(target.hum)
#' 
#' # For users with little R expertise, if the gene expression matrix comes as a data frame, it should be placed into "SummarizedExperiment" before proceeding to next step. An example is shown below by borrowing a data matrix from the brain data.
#' # Borrow a data matrix.
#' df <- assay(rse.hum); df[1:2, 1:3]
#' # Place the data matrix and targets file (target.hum) into "SummarizedExperiment".
#' rse.hum <- SummarizedExperiment(assay=df, colData=target.hum, rowData=NULL)
#' 
#' # The count matrix is normalised with estimateSizeFactors (type=‘ratio’).
#' se.nor.hum <- norm_data(data=rse.hum, method.norm='CNF', data.trans='log2')
#'
#' # Average replicates of concatenated sample__condition.
#' se.aggr.hum <- aggr_rep(data=se.nor.hum, sam.factor='organism_part', con.factor='disease', aggr='mean')
#' assay(se.aggr.hum)[49939:49942, ] # The concatenated tissue__conditions are the column names of the output data matrix.
#' 
#' # Genes with low expression level and low variantion are always filtered. 
#' se.fil.hum <- filter_data(data=se.aggr.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.01, 5), CV=c(0.3, 100), dir=NULL)

#' # Formatted SVG image.
#' svg.hum <- system.file("extdata/shinyApp/example", "homo_sapiens.brain.svg", package="spatialHeatmap")
#' # Plot spatial heatmaps of gene ENSG00000268433.
#' spatial_hm(svg=svg.hum, data=se.fil.hum, ID='ENSG00000268433', col.com=c("yellow", "blue", "purple"), width=1, height=0.5, sub.title.size=11, layout="gene", ncol=2, tis.trans=NULL, legend.position=c(0.5, -0.15), legend.nrow=1)


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1 \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Duncan Temple Lang and the CRAN Team (2018). XML: Tools for Parsing and Generating XML Within R and S-Plus. R package version 3.98-1.16. https://CRAN.R-project.org/package=XML \cr Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg 
#' Prudencio, Mercedes, Veronique V Belzil, Ranjan Batra, Christian A Ross, Tania F Gendron, Luc J Pregent, Melissa E Murray, et al. 2015. "Distinct Brain Transcriptome Profiles in C9orf72-Associated and Sporadic ALS." Nat. Neurosci. 18 (8): 1175–82
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom genefilter filterfun genefilter
#' @importFrom ggplot2 ggplot geom_bar aes theme element_blank margin element_rect coord_flip scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs
#' @importFrom rsvg rsvg_ps 
#' @importFrom grImport PostScriptTrace 
#' @importFrom XML addAttributes xmlParse xmlRoot xmlSize xmlSApply xmlAttrs xmlName xmlChildren saveXML xmlApply
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid grobTree unit
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is

spatial_hm <- function(svg.path, data, sam.factor=NULL, con.factor=NULL, ID, col.com=c("yellow", "purple", "blue"), col.bar="selected", bar.width=0.7, data.trans=NULL, tis.trans=NULL, width=1, height=1, legend.r=1, sub.title.size=11, lay.shm="gene", ncol=3, sam.legend='identical', legend.title=NULL, legend.ncol=NULL, legend.nrow=NULL, legend.position='bottom', legend.direction=NULL, legend.key.size=0.5, legend.label.size=8, legend.title.size=8, line.size=0.2, line.color='grey70', ...) {

  x <- y <- color_scale <- tissue <- line.type <- NULL  
  # Extract and filter data.
  options(stringsAsFactors=FALSE) 
  if (is.vector(data)) {

    vec.na <- make.names(names(data)); if (is.null(vec.na)) stop("Please provide names for the input data!")
    if (any(duplicated(vec.na))) stop('Please make sure data names are unique!')
    form <- grepl("__", vec.na); if (sum(form)==0) { vec.na <- paste0(vec.na, '__', vec.na) }
    data <- tryCatch({ as.numeric(data) }, warning=function(w) { stop("Please make sure input data are numeric!") }, error=function(e) { stop("Please make sure input data are numeric!") })
    if (is.null(ID)) stop('Please provide a name for the data!')
    gene <- as.data.frame(matrix(data, nrow=1, dimnames=list(ID, vec.na)))

  } else if (is(data, 'data.frame')|is(data, 'matrix')) {

   gene <- data; cna <- colnames(gene)
   if (any(duplicated(cna))) stop('Please make sure column names are unique!')
   form <- grepl("__", cna); if (sum(form)==0) { colnames(gene) <- paste0(cna, '__', 'con'); con.na=FALSE } else { gene <- gene[, form]; con.na=TRUE }

  } else if (is(data, 'SummarizedExperiment')) {

    gene <- assay(data); r.na <- rownames(gene); gene <- apply(gene, 2, as.numeric) # This step removes rownames of gene2.
    rownames(gene) <- r.na; colnames(gene) <- make.names(colnames(gene))
    col.meta <- as.data.frame(colData(data), stringsAsFactors=FALSE)
    # Factors teated by paste0/make.names are vecters.
    if (!is.null(sam.factor) & !is.null(con.factor)) { colnames(gene) <- paste0(make.names(col.meta[, sam.factor]), '__', make.names(col.meta[, con.factor])); con.na=TRUE } else if (!is.null(sam.factor) & is.null(con.factor)) { sam.na <- make.names(col.meta[, sam.factor]); colnames(gene) <- paste0(sam.na, "__", "con"); con.na=FALSE }

  }
    if (!is.null(data.trans)) if (data.trans=='log2') { 
          
      g.min <- min(gene) 
      if (g.min<0) gene <- gene-g.min+1; if (g.min==0) gene <- gene+1; gene <- log2(gene)  

    } else if (data.trans=='exp2') gene <- 2^gene
 
    # Color bar.
    bar.len=1000
    if (col.bar=="all") geneV <- seq(min(gene), max(gene), len=bar.len) else if (col.bar=="selected") geneV <- seq(min(gene[ID, , drop=FALSE]), max(gene[ID, , drop=FALSE]), len=bar.len)
    col <- colorRampPalette(col.com)(length(geneV))
    cs.g <- col_bar(geneV=geneV, cols=col, width=1, mar=c(3, 0.1, 3, 0.1)); cs.grob <- ggplotGrob(cs.g)    

    df_tis <- svg_df(svg.path=svg.path)
    if (is.character(df_tis)) stop(df_tis)
    g.df <- df_tis[['df']]; tis.path <- df_tis[['tis.path']]
    cname <- colnames(gene); con <- gsub("(.*)(__)(.*)", "\\3", cname); con.uni <- unique(con) 
    grob.lis <- grob_list(gene=gene, con.na=con.na, geneV=geneV, coord=g.df, ID=ID, cols=col, legend.col=df_tis[['fil.cols']], tis.path=tis.path, tis.trans=tis.trans, sub.title.size=sub.title.size, sam.legend=sam.legend, legend.title=legend.title, legend.ncol=legend.ncol, legend.nrow=legend.nrow, legend.position=legend.position, legend.direction=legend.direction, legend.key.size=legend.key.size, legend.label.size=legend.label.size, legend.title.size=legend.title.size, line.size=line.size, line.color=line.color, line.type=line.type, ...)
    g.arr <- lay_shm(lay.shm=lay.shm, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis[['grob.lis']], width=width, height=height, shiny=FALSE)
    cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc")) # "mm" is fixed, "npc" is scalable.
    g.lgd <- grob.lis[['g.lgd']]; grob.lgd <- ggplotGrob(g.lgd)
    # Layout matrix of legend.
    if (lay.shm=='gene') lay.lgd <- matrix(seq_len(ceiling(length(con.uni)/ncol)), byrow=FALSE)
    if (lay.shm=='con') lay.lgd <- matrix(seq_len(ceiling(length(ID)/ncol)), byrow=FALSE)
    lgd.arr <- arrangeGrob(grobs=list(grobTree(grob.lgd)), layout_matrix=lay.lgd, widths=unit(width, "npc"), heights=unit(legend.r, "npc"))
    grid.arrange(cs.arr, g.arr, lgd.arr, ncol=3, widths=c(bar.width, 10/(ncol+1)*ncol, 10/(ncol+1)))

}
