#' Spatial Heatmap
#'
#'It takes gene expression profiles and an associated SVG image as input. In the SVG image, tissue regions are pre-defined and associated with expression pfofiles through tissue names. The expression profiles of the input gene(s) under each condition are mapped to each tissue in the form of different colours. Its application is not limited to gene expression data. It can be used as long as a data matrix and an associated SVG image are provided, such as population data generated in different years across different cities. 

#' @param svg The path of the SVG image, where different tissues are pre-defined and associated with expression pfofiles through tissue names. \cr E.g.: system.file("extdata/example", "test_final.svg", package = "spatialHeatmap") (Mustroph et al. 2009).

#' @param data The gene expression matrix and metadata (optional) in the form of the "SummarizedExperiment" object. In the expression matrix, the row and column names should be gene IDs and sample/conditions, respectively. The sample/condition names MUST be fomatted this way: a sample name is followed by double underscore then the condition, such as "epidermis__140mM_1h", where "epidermis" is the sample and "140mM_1h" is the condition (Geng et al. 2013). In the column names of sample/condition, only letters, digits, single underscore, dots, single space are allowed. Not all samples in the matrix need to be present in the SVG image, and vice versa. Only samples present in the SVG image are recognised and coloured in the spatial heatmap. \cr Example expression matrix: system.file("extdata/example", "root_expr_row_gen.txt", package = "spatialHeatmap").

#' @param pOA It specifies parameters of the filter function "pOverA" from the package "genefilter" (Gentleman et al. 2018). It filters genes according to the proportion "p" of samples where the expression values exceeding a threshold "A". The input is a vector of two numbers, where the first one is the "p" and the second one is "A". The default is c(0, 0), which means no filter is applied. \cr E.g. c(0.1, 2) means genes whose expression values over 2 in at least 10\% of all samples are kept. 

#' @param CV It specifies parameters of the filter function "cv" from the package "genefilter" (Gentleman et al. 2018), which filters genes according to the coefficient of variation (CV). The input is a vector of two numbers, specifying the CV range. The default is c(0, 10000), where the range is set from 0 to very large (10000) so as not to apply filtering. \cr E.g. c(0.1, 5) means genes with CV between 0.1 and 5 are kept.

#' @param ID The gene IDs used to colour the spatial heatmaps. It can be a single gene or a vector of multiple genes.

#' @param col.com A character vector of the colour components used to build the colour scale, e.g. the default is c("yellow", "blue", "purple").

#' @param col.bar It has two values "selected" and "all", which specifies whether the colour scale is built using the input genes ("selected") or whole data matrix ("all"). The default is "selected".
#' @param tis.trans A vector of tissue names. These tissues cover other tissues and should be set transparent. E.g c("epidermis", "cortex").
#' @param width The width, relative to height, of each subplot. The default is 1.

#' @param height The height, relative to width, of each subplot. The default is 1.

#' @param sub.title.size The size of each subtitle. The default is 11.

#' @param layout The organisation of the subplots. The options are "gene" or "con" (condition).

#' @param ncol Number of columns to display the subplots.
#' @param ...  Other arguments passed to the function "ggplot".

#' @return It generates an image of spatial heatmap(s) along with a colour key.

#' @section Details:
#' Details about how to properly format and associate a custom SVG image with a data matrix are provided in a separate SVG_tutorial.html. The "data" input can be the value returned by the function "filter.data" or built on the expression matrix. See examples blow.

#' @examples
#' # Creat the "SummarizedExperiment" class. Refer to the R package "SummarizedExperiment" for more details.
#' data.path <- system.file("extdata/example", "root_expr_row_gen.txt", package = "spatialHeatmap")   
#' ## The expression matrix, where the row and column names should be gene IDs and sample/conditions, respectively.
#' library(data.table); expr <- fread(data.path, sep='\t', header=TRUE, fill=TRUE)
#' col.na <- colnames(expr)[-ncol(expr)]; row.na <- as.data.frame(expr[, 1])[, 1]
#' expr <- as.matrix(as.data.frame(expr, stringsAsFactors=FALSE)[, -1])
#' rownames(expr) <- row.na; colnames(expr) <- col.na
#' library(SummarizedExperiment); expr <- SummarizedExperiment(assays=list(expr=expr)) # Metadata is not necessary.  

#' # The svg image path.
#' svg.path <- system.file("extdata/example", "root_cross_final.svg", package="spatialHeatmap")
#' # The expression profiles of gene "PSAC" and "NDHG" under different conditions are mapped to tissues defined in the SVG image in the form of different colours. 
#' spatial.hm(svg=svg.path, data=expr, pOA=c(0, 0), CV=c(0, 10000), ID=c("PSAC", "NDHG"), col.com=c("yellow", "blue", "purple"), width=1, height=1, sub.title.size=11, layout="gene", ncol=3)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1 \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Duncan Temple Lang and the CRAN Team (2018). XML: Tools for Parsing and Generating XML Within R and S-Plus. R package version 3.98-1.16. https://CRAN.R-project.org/package=XML \cr Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr Simon Urbanek and Jeffrey Horner (2015). Cairo: R graphics device using cairo graphics library for creating high-quality bitmap (PNG, JPEG, TIFF), vector (PDF, SVG, PostScript) and display (X11 and Win32) output. R package version 1.5-9. https://CRAN.R-project.org/package=Cairo \cr    
#' Geng Y, Wu R, Wee CW, Xie F et al. A spatio-temporal understanding of growth regulation during the salt stress response in Arabidopsis. Plant Cell 2013 Jun;25(6):2132-54. PMID: 23898029 \cr Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggplot geom_bar aes theme element_blank margin element_rect coord_flip scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual theme element_blank ggtitle  element_rect margin element_text labs 
#' @importFrom rsvg rsvg_ps 
#' @importFrom grImport PostScriptTrace 
#' @importFrom XML addAttributes xmlParse xmlRoot xmlSize xmlSApply xmlAttrs xmlName xmlChildren saveXML xmlApply
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid grobTree unit
#' @importFrom grDevices colorRampPalette
#' @importFrom Cairo Cairo
#' @importFrom methods is
#' @importFrom genefilter filterfun genefilter

spatial.hm <- function(svg, data, pOA=c(0, 0), CV=c(0, 10000), ID, col.com=c("yellow", "blue", "purple"), col.bar="selected", tis.trans=NULL, width=1, height=1, sub.title.size=11, layout="gene", ncol=3, ...) {

    x <- y <- color_scale <- tissue <- NULL
    # Extract and filter data.
    gene <- assay(data); r.na <- rownames(gene)
    gene <- apply(gene, 2, as.numeric) # This step removes rownames of gene2.
    rownames(gene) <- r.na
    ffun <- filterfun(pOverA(pOA[1], pOA[2]), cv(CV[1], CV[2]))
    filtered <- genefilter(gene, ffun); gene <- gene[filtered, ]
    id.in <- ID %in% rownames(gene)
    if (!all(id.in)) stop(paste0(ID[!id.in], " is filtered.")) 

    # Color bar.
    bar.len=1000
    if (col.bar=="all") geneV <- seq(min(gene), max(gene), len=bar.len) else if (col.bar=="selected") geneV <- seq(min(gene[ID, , drop=FALSE]), max(gene[ID, , drop=FALSE]), len=bar.len)
    col <- colorRampPalette(col.com)(length(geneV))
    cs.df <- data.frame(color_scale=geneV, y=1)
    cs.g <- ggplot()+geom_bar(data=cs.df, aes(x=color_scale, y=y), fill=col, stat="identity", width=((max(geneV)-min(geneV))/bar.len)*1)+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=margin(3, 0.1, 3, 0.1, "cm"), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"))+coord_flip()
    if (max(geneV)<10000) cs.g <- cs.g+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand = c(0,0))
    if (max(geneV)>=10000) cs.g <- cs.g+scale_y_continuous(expand=c(0,0))+scale_x_continuous(labels=function(x) format(x, scientific=TRUE)); cs.grob <- ggplotGrob(cs.g)

    # Make sure the style is correct. If the stroke width is not the same across polygons such as '0.0002px', '0.216px', some stroke outlines cannot be recognised by 'PostScriptTrace'. Then some polygons are missing. Since the ggplot is based on 'stroke' not 'fill'.
    tmp <- system.file("extdata/tmp", package="spatialHeatmap")
    xmlfile <- xmlParse(svg); xmltop <- xmlRoot(xmlfile); ply <- xmltop[[xmlSize(xmltop)]]
    style <- 'stroke:#000000;stroke-width:5.216;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1' # 'fill' is not necessary.
    # Change 'style' of all polygons.
    for (i in seq_len(xmlSize(ply))) {                      
         
      addAttributes(ply[[i]], style=style) 
      if (xmlSize(ply[[i]])>=1) for (j in seq_len(xmlSize(ply[[i]]))) { addAttributes(ply[[i]][[j]], style=style) }
        
    }; svg.inter <- paste0(tmp, '/internal.svg'); saveXML(doc=xmlfile, file=svg.inter)
    # SVG file conversion. 
    rsvg_ps(svg.inter, file=sub("svg$", "ps", svg.inter))
    p1 <- sub("svg$", "ps", svg.inter); p2 <- paste0(sub("svg$", "ps", svg.inter), ".xml")  
    if (length(grep("~", svg.inter))) {

        wd1 <- getwd(); setwd("~"); hm <- getwd(); setwd(wd1)
        p1 <- sub("~", hm, p1); p2 <- sub("~", hm, p2)

    }; PostScriptTrace(p1, p2) 
    grml <- xmlParse(p2); top <- xmlRoot(grml) # Use internal svg to get coordinates.
    xml <- xmlParse(svg); xmltop <- xmlRoot(xml); size <- xmlSize(xmltop) # Use original not internal svg to get path ids. Otherwise, errors can come up.
    do.call(file.remove, list(svg.inter, p1, p2))

        lis.ma <- xmlApply(xmltop[[size]], xmlAttrs)
        if (is(lis.ma, "matrix")) { id.xml <- lis.ma["id", ] } else if (is(lis.ma, "list")) {

            id.xml <- NULL
            for (i in seq_len(length(lis.ma))) { id.xml <- c(id.xml, lis.ma[[i]][["id"]]) }

        }

        xml.na <- NULL
        for (i in seq_len(xmlSize(xmltop[[size]]))) { xml.na <- c(xml.na, xmlName(xmltop[[size]][[i]])) } 

        for (j in seq_len(length(xml.na))) {

            if (xml.na[j]=="g") {

                len.dif <- length(id.xml)-length(xml.na); g.size <- xmlSize(xmltop[[size]][[j]])
                if ((j+1+len.dif) <= length(id.xml)) {

                if (j==1) { 

                    id.xml <- c(paste0(id.xml[j+len.dif], "_", seq_len(g.size)), 
                    id.xml[(j+1+len.dif):length(id.xml)]) 

                } else if (j>1) {

                    id.xml <- c(id.xml[seq_len(j-1+len.dif)], paste0(id.xml[j+len.dif], "_", 
                    seq_len(g.size)), id.xml[(j+1+len.dif):length(id.xml)])

                }

                } else if ((j+1+len.dif) >= length(id.xml)) { 

                    id.xml <- c(id.xml[seq_len(j-1+len.dif)], paste0(id.xml[j+len.dif], "_", seq_len(g.size))) 

                }

            }

        }

        tis.path <- gsub("_\\d+$", "", id.xml)

        k <-0; df <- NULL; 
        for (i in seq_len(xmlSize(top)-1)) {

            if (xmlAttrs(top[[i]])['type']=='stroke') {

                k <- k+1; chil <- xmlChildren(top[[i]])
                xy <- chil[grep("move|line", names(chil))]
                coor <- matrix(NA, nrow=length(xy), ncol=2, dimnames=list(NULL, c("x", "y")))

            for (j in seq_len(length(xy))) {

                coor[j, "x"] <- as.numeric(xmlAttrs(xy[[j]])["x"])
                coor[j, "y"] <- as.numeric(xmlAttrs(xy[[j]])["y"])

            }

            df0 <- cbind(tissue=id.xml[k], data.frame(coor, stringsAsFactors=FALSE), stringsAsFactors=TRUE)
            df <- rbind(df, df0, stringsAsFactors=TRUE)

            }

        }; g.df <- df

    # Map colours to samples according to expression level.
    cname <- colnames(gene); con <- gsub("(.*)(__)(.*)", "\\3", cname)

        g.list <- function(j) {

            g.col <- NULL; con.idx <- grep(paste0("^", j, "$"), con)
            tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]

            for (i in tis.path) {

                tis.idx <- which(tis.col1 %in% i)
                if (length(tis.idx)==1) { g.col <- c(g.col, scol1[tis.idx])
                } else if (length(tis.idx)==0) { g.col <- c(g.col, "white") }

            }
            names(g.col) <- tis.df <- unique(g.df[, 'tissue']) # The colors might be internally re-ordered alphabetically during mapping, so give them names to fix the match with tissues. E.g. c('yellow', 'blue') can be re-ordered to c('blue', 'yellow'), which makes tissue mapping wrong. Correct: colours are not re-ordered. The 'tissue' in 'data=g.df' are internally re-ordered according to a factor. Therfore, 'tissue' should be a factor with the right order. Otherwise, disordered mapping can happen.
            # Make selected tissues transparent by setting their colours as NA.
            if (!is.null(tis.trans)) for (i in tis.df) { if (sub('_\\d+$', '', i) %in% tis.trans) g.col[i] <- NA }
            g <- ggplot(...)+geom_polygon(data=g.df, aes(x=x, y=y, fill=tissue), color="black")+scale_fill_manual(values=g.col, guide=FALSE)+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=sub.title.size))+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))+ggtitle(paste0(k, "_", j)); return(g)

        }

        grob.na <- grob.lis <- NULL; for (k in ID) {

            scol <- NULL; for (i in gene[k, ]) { 

                ab <- abs(i-geneV); col.ind <- which(ab==min(ab))[1]; scol <- c(scol, col[col.ind])

            }

            idx <- grep("__", cname); c.na <- cname[idx]
            tis.col <- gsub("(.*)(__)(.*)", "\\1", c.na); g.lis <- NULL
            con.uni <- unique(con); grob.na0 <- paste0(k, "_", con.uni)
            g.lis <- lapply(con.uni, g.list); grob <- lapply(g.lis, ggplotGrob)
            names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob) 

        }; grob.all <- c(list(cs=cs.grob), grob.lis)

    cs.arr <- arrangeGrob(grobs=list(grobTree(grob.all[[1]])), layout_matrix=cbind(1), widths=unit(25, "mm"))

    # Organise layout.
    if (layout=="gene") {

        all.cell <- ceiling(length(unique(con))/as.numeric(ncol))*as.numeric(ncol)
        cell.idx <- c(seq_len(length(unique(con))), rep(NA, all.cell-length(unique(con))))
        m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
        lay <- NULL; for (i in seq_len(length(ID))) { lay <- rbind(lay, m+(i-1)*length(unique(con))) }
        g.tr <- lapply(grob.all[2:length(grob.all)], grobTree)
        n.col <- ncol(lay); n.row <- nrow(lay)
        g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc"))
        grid.arrange(cs.arr, g.arr, ncol=2, widths=c(1.5, 8)) 

    } else if (layout=="con") {

        grob.all.na <- names(grob.all); grob.all.idx <- NULL
        for (i in unique(con)) {

            grob.all.idx <- c(grob.all.idx, grep(paste0(i, "$"), grob.all.na))

        }

        grob.all.con <- grob.all[grob.all.idx]
        all.cell <- ceiling(length(ID)/as.numeric(ncol))*as.numeric(ncol)
        cell.idx <- c(seq_len(length(ID)), rep(NA, all.cell-length(ID)))
        m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
        lay <- NULL; for (i in seq_len(length(unique(con)))) { lay <- rbind(lay, m+(i-1)*length(ID)) }
        g.tr <- lapply(grob.all.con, grobTree); g.tr <- g.tr[names(grob.all.con)]
        n.col <- ncol(lay); n.row <- nrow(lay)
        g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc"))
        grid.arrange(cs.arr, g.arr, ncol=2, widths=c(1.5, 8)) 

    }

}
