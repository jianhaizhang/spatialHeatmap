#' Spatial Heatmap
#'
#'It takes gene expression profiles and an associated SVG image as input. In the SVG image, tissue regions are pre-defined and associated with expression pfofiles through tissue names. The expression profiles of the input gene(s) under each condition are mapped to each tissue in the form of different colours. Its application is not limited to gene expression data. It can be used as long as a data matrix and an associated SVG image are provided, such as population data generated in different years across different cities. 

#' @param svg.path The path of the SVG image, where different tissues are pre-defined and associated with expression pfofiles through tissue names. \cr E.g.: system.file("extdata/example", "test_final.svg", package = "spatialHeatmap") (Mustroph et al. 2009).

#' @inheritParams filter_data
#' @inheritParams grob_list

#' @param ID The gene IDs used to colour the spatial heatmaps. It can be a single gene or a vector of multiple genes.

#' @param col.com A character vector of the colour components used to build the colour scale, e.g. the default is c("yellow", "purple", "blue").

#' @param col.bar It has two values "selected" and "all", which specifies whether the colour scale is built using the input genes ("selected") or whole data matrix ("all"). The default is "selected".
#' @param data.trans "log2" or "exp.2". Default is NULL. If colours across tissues cannot distinguish due to low variance or outliers, transform the data matrix by log2 or 2-base expoent (exp.2). 
#' @param width The width, relative to height, of each subplot. The default is 1.
#' @param height The height, relative to width, of each subplot. The default is 1.
#' @param lay.shm The organisation of the subplots. The options are "gene" or "con" (condition).
#' @param ncol Number of columns to display the subplots.
#' @return It generates an image of spatial heatmap(s) along with a colour key.

#' @section Details:
#' Details about how to properly format and associate a custom SVG image with a data matrix are provided in a separate SVG_tutorial.html. The "se" parameter can be the value returned by the function "filter_data" or built on the expression matrix. See examples blow.

#' @examples
#' # The example data is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
#' library(ExpressionAtlas)
#' rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]] # Access the RNA-seq raw count matrix.  
#' assay(rse.hum)[1:3, 1:3] # The raw count matrix is in form of "RangedSummarizedExperiment".
#' # A targets file describing replicates of samples and conditions is required, which should be made based on the "colData" slot in "RangedSummarizedExperiment". See the "se" parameter for details. This targets file is available in spatialHeatmap.
#' brain.pa <- system.file('extdata/example_data/target_brain.txt', package='spatialHeatmap')
#' target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
#' The “organism_part” and “disease” column describes tissue and condition replicates respectively. The former includes 2 tissues cerebellum and frontal cortex. They should be identical with respective tissue ids in the SVG image for plotting spatial heatmaps. Note that the replicates of the same tissue or condition should have the identical name.
#' target.hum[c(1:3, 41:42), 4:5]

#' # Place the targets file into "colData" slot. 
#' colData(rse.hum) <- DataFrame(target.hum)

#' # For users with little R expertise, if the gene expression matrix comes as a data frame, it should be placed into "SummarizedExperiment" before proceeding to next step. An example is shown below by borrowing a data frame from the brain data.
#' # Borrow a data matrix.
#' df <- assays(rse.hum)[[1]]; df[1:2, 1:3]
#' # Place the data matrix and targets file (target.brain) into "SummarizedExperiment". The "rowData" slot is optional.
#' rse.hum <- SummarizedExperiment(assays=list(counts=df), colData=target.hum, rowData=NULL)

#' # The count matrix is normalised with estimateSizeFactors (type=‘ratio’) from DESeq2 (Love, Huber, and Anders 2014) and converted to log2 unit.
#' se.nor.hum <- norm_data(se=rse.hum, method.norm='ratio', data.trans='log2')
#' # To aggregate replicates, the tissue and condition column in “colData” slot need to be specified in function aggr_rep. This function concatenates tissue and condition replicates in the targets file with a "double underscore" and treated them as tissue-condition replicates for aggregating. For example, in the above targets file, each of cerebellum and frontal_cortex are concatenated with each of ALS and normal by "__" and cerebellum__ALS, cerebellum__normal, frontal_cortex__ALS, frontal_cortex__normal are used as the replicates for aggregating. The concatenated replicates can be aggregated by mean or median. Here mean is chosen.
#' # In downstream spatial heatmap plotting, the "double underscore" is indispensable as it is the separator when the algorithm recognises tissues and conditions.
#' se.aggr.hum <- aggr_rep(se=se.nor.hum, sam.factor='organism_part', con.factor='disease', aggr='mean')
#' # The concatenated tissue__conditions are the column names of the output data matrix.
#' assay(se.aggr.hum)[49939:49942, ]

#' # Filter genes.
#' se.fil.hum <- filter_data(data=se.aggr.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.01, 5), CV=c(0.6, 100), dir=NULL)

#' # Formatted SVG image.
#' svg.hum <- system.file("extdata/shinyApp/example", "homo_sapiens.brain.svg", package="spatialHeatmap")
#' # Plot spatial heatmaps of gene ENSG00000268433.
#' spatial_hm(svg=svg.hum, se=se.fil.hum, ID='ENSG00000268433', col.com=c("yellow", "blue", "purple"), width=1, height=0.5, sub.title.size=11, layout="gene", ncol=2, tis.trans=NULL)


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1 \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Duncan Temple Lang and the CRAN Team (2018). XML: Tools for Parsing and Generating XML Within R and S-Plus. R package version 3.98-1.16. https://CRAN.R-project.org/package=XML \cr Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr    
#' Prudencio, Mercedes, Veronique V Belzil, Ranjan Batra, Christian A Ross, Tania F Gendron, Luc J Pregent, Melissa E Murray, et al. 2015. “Distinct Brain Transcriptome Profiles in C9orf72-Associated and Sporadic ALS.” Nat. Neurosci. 18 (8): 1175–82
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.” Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8

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

spatial_hm <- function(svg.path, se, sam.factor=NULL, con.factor=NULL, ID, col.com=c("yellow", "purple", "blue"), col.bar="selected", data.trans=NULL, tis.trans=NULL, width=1, height=1, sub.title.size=11, lay.shm="gene", ncol=3, sam.legend='identical', title.legend=NULL, ncol.legend=NULL, nrow.legend=NULL, pos.legend='right', legend.key.size=0.5, legend.label.size=8, legend.title.size=8, line.size=0.2, line.color='grey70', ...) {

    x <- y <- color_scale <- tissue <- NULL
    # Extract and filter data.

    gene <- assay(se); r.na <- rownames(gene); gene <- apply(gene, 2, as.numeric) # This step removes rownames of gene2.
    rownames(gene) <- r.na
    if (!is.null(data.trans)) if (data.trans=='log2') { 
          
      g.min <- min(gene) 
      if (g.min<0) gene <- gene-g.min+1; if (g.min==0) gene <- gene+1; gene <- log2(gene)  

    } else if (data.trans=='exp.2') gene <- 2^gene
 
    if (!is.null(sam.factor) & !is.null(con.factor)) { col.met <- as.data.frame(colData(se)); colnames(gene) <- rownames(col.met) <- paste(col.met[, sam.factor], col.met[, con.factor], sep='__') }

    # Color bar.
    bar.len=1000
    if (col.bar=="all") geneV <- seq(min(gene), max(gene), len=bar.len) else if (col.bar=="selected") geneV <- seq(min(gene[ID, , drop=FALSE]), max(gene[ID, , drop=FALSE]), len=bar.len)
    col <- colorRampPalette(col.com)(length(geneV))
    cs.g <- col_bar(geneV=geneV, cols=col, width=1, mar=c(3, 0.1, 3, 0.1)); cs.grob <- ggplotGrob(cs.g)    

    df_tis <- svg_df(svg.path=svg.path)
    if (is.character(df_tis)) stop(df_tis)
    g.df <- df_tis[['df']]; tis.path <- df_tis[['tis.path']]
    cname <- colnames(gene); con <- gsub("(.*)(__)(.*)", "\\3", cname); con.uni <- unique(con) 
    grob.lis <- grob_list(gene=gene, geneV=geneV, coord=g.df, ID=ID, cols=col, tis.path=tis.path, tis.trans=tis.trans, sub.title.size=sub.title.size, sam.legend=sam.legend, title.legend=title.legend, ncol.legend=ncol.legend, nrow.legend=nrow.legend, pos.legend=pos.legend, legend.key.size=legend.key.size, legend.label.size=legend.label.size, legend.title.size=legend.title.size, line.size=line.size, line.color=line.color, line.type=line.type, ...)
    g.arr <- lay_shm(lay.shm=lay.shm, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis, width=width, height=height, shiny=FALSE)
    cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(25, "mm"))
    grid.arrange(cs.arr, g.arr, ncol=2, widths=c(1.5, 8))

}
