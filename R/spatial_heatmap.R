#' Spatial Heatmap
#'
#'It takes gene expression profiles and an associated SVG image as input. In the SVG image, tissue regions are pre-defined and associated with expression pfofiles through tissue names. The expression profiles of the input gene(s) under each condition are mapped to each tissue in the form of different colours. Its application is not limited to gene expression data. It can be used as long as a data matrix and an associated SVG image are provided, such as population data generated in different years across different cities. 

#' @param svg.path The path of the SVG image, where different tissues are pre-defined and associated with expression pfofiles through tissue names. \cr E.g.: system.file("extdata/example", "test_final.svg", package = "spatialHeatmap") (Mustroph et al. 2009).

#' @param data The gene expression matrix and metadata (optional) in the form of the "SummarizedExperiment" object. In the expression matrix, the row and column names should be gene IDs and sample/conditions, respectively. The sample/condition names MUST be fomatted this way: a sample name is followed by double underscore then the condition, such as "epidermis__140mM_1h", where "epidermis" is the sample and "140mM_1h" is the condition (Geng et al. 2013). In the column names of sample/condition, only letters, digits, single underscore, dots, single space are allowed. Not all samples in the matrix need to be present in the SVG image, and vice versa. Only samples present in the SVG image are recognised and coloured in the spatial heatmap. \cr Example expression matrix: system.file("extdata/example", "root_expr_row_gen.txt", package = "spatialHeatmap").

#' @param pOA It specifies parameters of the filter function "pOverA" from the package "genefilter" (Gentleman et al. 2018). It filters genes according to the proportion "p" of samples where the expression values exceeding a threshold "A". The input is a vector of two numbers, where the first one is the "p" and the second one is "A". The default is c(0, 0), which means no filter is applied. \cr E.g. c(0.1, 2) means genes whose expression values over 2 in at least 10\% of all samples are kept. 

#' @param CV It specifies parameters of the filter function "cv" from the package "genefilter" (Gentleman et al. 2018), which filters genes according to the coefficient of variation (CV). The input is a vector of two numbers, specifying the CV range. The default is c(0, 100), where the range is set from 0 to very large (100) so as not to apply filtering. \cr E.g. c(0.1, 5) means genes with CV between 0.1 and 5 are kept.

#' @param ID The gene IDs used to colour the spatial heatmaps. It can be a single gene or a vector of multiple genes.

#' @param col.com A character vector of the colour components used to build the colour scale, e.g. the default is c("yellow", "purple", "blue").

#' @param col.bar It has two values "selected" and "all", which specifies whether the colour scale is built using the input genes ("selected") or whole data matrix ("all"). The default is "selected".
#' @param tis.trans A vector of tissue names. These tissues cover other tissues and should be set transparent. E.g c("epidermis", "cortex").
#' @param width The width, relative to height, of each subplot. The default is 1.
#' @param height The height, relative to width, of each subplot. The default is 1.
#' @param sub.title.size The size of each subtitle. The default is 11.
#' @param lay.shm The organisation of the subplots. The options are "gene" or "con" (condition).
#' @param ncol Number of columns to display the subplots.
#' @param ...  Other arguments passed to the function "ggplot".
#' @return It generates an image of spatial heatmap(s) along with a colour key.

#' @section Details:
#' Details about how to properly format and associate a custom SVG image with a data matrix are provided in a separate SVG_tutorial.html. The "data" input can be the value returned by the function "filter_data" or built on the expression matrix. See examples blow.

#' @examples
#' # Creat the "SummarizedExperiment" class. Refer to the R package "SummarizedExperiment" for more details.
#' data.path <- system.file("extdata/shinyApp/example", "root_expr_row_gen.txt", package = "spatialHeatmap")   
#' ## The expression matrix, where the row and column names should be gene IDs and sample/conditions, respectively.
#' library(data.table); expr <- fread(data.path, sep='\t', header=TRUE, fill=TRUE)
#' col.na <- colnames(expr)[-ncol(expr)]; row.na <- as.data.frame(expr[, 1])[, 1]
#' expr <- as.matrix(as.data.frame(expr, stringsAsFactors=FALSE)[, -1])
#' rownames(expr) <- row.na; colnames(expr) <- col.na
#' library(SummarizedExperiment); expr <- SummarizedExperiment(assays=list(expr=expr)) # Metadata is not necessary.  

#' # The svg image path.
#' svg.path <- system.file("extdata/shinyApp/example", "root_cross_final.svg", package="spatialHeatmap")
#' # The expression profiles of gene "PSAC" and "NDHG" under different conditions are mapped to tissues defined in the SVG image in the form of different colours. 
#' spatial_hm(svg.path=svg.path, data=expr, pOA=c(0, 0), CV=c(0, 100), ID=c("PSAC", "NDHG"), col.com=c("yellow", "blue", "purple"), width=1, height=1, sub.title.size=11, layout="gene", ncol=3)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1 \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Duncan Temple Lang and the CRAN Team (2018). XML: Tools for Parsing and Generating XML Within R and S-Plus. R package version 3.98-1.16. https://CRAN.R-project.org/package=XML \cr Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr    
#' Geng Y, Wu R, Wee CW, Xie F et al. A spatio-temporal understanding of growth regulation during the salt stress response in Arabidopsis. Plant Cell 2013 Jun;25(6):2132-54. PMID: 23898029 \cr Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8
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

spatial_hm <- function(svg.path, data, pOA=c(0, 0), CV=c(0, 10000), ID, col.com=c("yellow", "purple", "blue"), col.bar="selected", tis.trans=NULL, width=1, height=1, sub.title.size=11, lay.shm="gene", ncol=3, ...) {

    x <- y <- color_scale <- tissue <- NULL
    # Extract and filter data.

    gene <- assay(data); r.na <- rownames(gene); gene <- apply(gene, 2, as.numeric) # This step removes rownames of gene2.
    rownames(gene) <- r.na; ffun <- filterfun(pOverA(pOA[1], pOA[2]), cv(CV[1], CV[2]))
    filtered <- genefilter(gene, ffun); gene <- gene[filtered, ]
    id.in <- ID %in% rownames(gene); if (!all(id.in)) stop(paste0(ID[!id.in], " is filtered.")) 

    # Color bar.
    bar.len=1000
    if (col.bar=="all") geneV <- seq(min(gene), max(gene), len=bar.len) else if (col.bar=="selected") geneV <- seq(min(gene[ID, , drop=FALSE]), max(gene[ID, , drop=FALSE]), len=bar.len)
    col <- colorRampPalette(col.com)(length(geneV))
    cs.g <- col_bar(geneV=geneV, cols=col, width=1, mar=c(3, 0.1, 3, 0.1)); cs.grob <- ggplotGrob(cs.g)    

    df_tis <- svg_df(svg.path=svg.path); g.df <- df_tis[['df']]; tis.path <- df_tis[['tis.path']]
    cname <- colnames(gene); con <- gsub("(.*)(__)(.*)", "\\3", cname); con.uni <- unique(con) 
    grob.lis <- grob_list(gene=gene, geneV=geneV, coord=g.df, ID=ID, cols=col, tis.path=tis.path, tis.trans=tis.trans, sub.title.size=11)
    g.arr <- lay_shm(lay.shm=lay.shm, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis, width=width, height=height, shiny=FALSE)
    cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(25, "mm"))
    grid.arrange(cs.arr, g.arr, ncol=2, widths=c(1.5, 8))

}
