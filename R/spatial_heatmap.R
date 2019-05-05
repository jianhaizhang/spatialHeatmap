#' Spatial Heatmap
#'
#'It takes a gene expression matrix and an associated svg image as input, where each tissue region is defined. The expression profiles of the input gene(s) under each condition are mapped to each tissue in the form of different colours. Its application is not limited to gene expression data, since it can be used as long as a data matrix and an associated svg image are provided, such as population data generated in different years across different cities. 

#' @param svg The path of the svg image, where different tissues are defined and can be labeled with different colors. \cr E.g.: system.file("extdata/example", "test_final.svg", package = "spatialHeatmap")

#' @param data The path of the data matrix. In the gene expression matrix, the row and column names are usually gene IDs and sample/conditions, respectively. The sample/condition names MUST be fomatted this way: a sample name is followed by double underscore then the condition, such as “stele__140mM_48h" (Mustroph et al. 2009), where stele is the sample and 140mM_48h is the condition. One column of meta data (e.g. gene annotation) can also be included in parallel with sample/condition at the end. In the column names of sample/condition and meta data, only letters, digits, single underscore, dots are allowed. Not all samples in the matrix necessarily need to be present in the svg image, and vice versa. Only samples present in the svg image are recognised and coloured. The rows and columns can be swaped and the function identifies gene IDs by "isRowGen". \cr E.g.: system.file("extdata/example", "root_expr_ann_row_gen.txt", package = "spatialHeatmap").

#' @param sep The seprator of the data matrix, e.g. ",", "\\t", ";".

#' @param isRowGene If row names are gene IDs it is "TRUE", otherwise it is "FALSE". If the data matrix is not gene expression, this argument is set according to if the row is equivalent to gene ID.

#' @param pOA It specifies parameters of the filter function "pOverA" from the package "genefilter" (Gentleman et al. 2018). It filters genes according to the proportion "p" of samples where the expression values exceeding a threshold "A". The input is a two-component vector, where the first one is the "p" and the second one is "A". The default is c(0, 0), which means no filter is applied. \cr E.g. c(0.1, 2) means genes whose expression values over 2 in at least 10\% of all samples are kept.                           

#' @param CV It specifies parameters of the filter function "cv" from the package "genefilter" (Gentleman et al. 2018), which filters genes according to the coefficient of variation (CV). The input is a two-component vector, specifying the CV range. The default is c(0, 10000), where the range is set from very samll (0) to very large (10000) so as to not apply filtering. \cr E.g. c(0.1, 5) means genes with CV between 0.1 and 5 are kept.

#' @param ID The gene IDs used to display the spatial heatmaps.

#' @param col.com The colour components used to build the colour scale, which must be separated with comma and no space, e.g. the default is "green,blue,purple,yellow,red".

#' @param col.bar It has two values "selected" and "all", which specifies whether the colour scale is built using the input genes ("selected") or whole data matrix ("all"). The default is "selected".

#' @param width The width, relative to height, of each subplot. The default is 1.

#' @param height The height, relative to width, of each subplot. The default is 1.

#' @param sub.title.size The size of each subtitle. The default is 11.

#' @param layout The layout of the subplots. The options are "gene" or "con" (condition). For example, in gene expression matrix the spatial tissue heatmaps can be displayed by gene or condition (con).

#' @param ncol Number of columns to display the subplots.

#' @return It generates an image of spatial heatmap(s) along with a colour key.

#' @section Details:
#' Details about how to properly format and associate a custom svg image with a data matrix are provided in the vignette "SVG_tutorial.html". \cr The data matrix can be first filtered with the function "filter.data" then provide the resulting matrix path to this function.

#' @examples
#' # The expression matrix path.
#' data.path <- system.file("extdata/example", "root_expr_ann_row_gen.txt", package = "spatialHeatmap")
#' # The svg image path.
#' svg.path <- system.file("extdata/example", "root_cross_final.svg", package = "spatialHeatmap")
#' # The expression profiles of gene "PSAC" and "NDHG" under different conditions are mapped to tissues defined in the svg image in the form of different colours. 
#' spatial.hm(svg=svg.path, data=data.path, sep="\t", isRowGene=TRUE, pOA=c(0.1, 3), CV=c(0.05, 1000), ID=c("PSAC", "NDHG"), col.com=c("green", "blue", "purple", "yellow", "red"), width=1, height=1, sub.title.size=11, layout="gene", ncol=3)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8
#' R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1
#' @export
#' @importFrom ggplot2 ggplot geom_bar aes theme element_blank margin element_rect coord_flip scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual theme element_blank ggtitle  element_rect margin element_text labs 
#' @importFrom rsvg rsvg_ps 
#' @importFrom grImport PostScriptTrace 
#' @importFrom XML xmlParse xmlRoot xmlSize xmlSApply xmlAttrs xmlName xmlChildren 
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid grobTree unit
#' @importFrom shiny runApp
#' @importFrom grDevices colorRampPalette
#' @importFrom Cairo Cairo
#' @importFrom methods is

spatial.hm <- function(svg, data, sep, isRowGene, pOA=c(0, 0), CV=c(0, 10000), ID, col.com=c("green", "blue", "purple", "yellow", "red"), col.bar="selected", width=1, height=1, sub.title.size=11, layout, ncol) {

    # require(grImport); require(rsvg); require(ggplot2); require(gridExtra); require(Cairo); require(grid); require(XML); require(data.table); require(genefilter)
    x <- y <- color_scale <- tissue <- NULL
    # Import and filter data.
    gene.f <- fread(data, header=TRUE, sep=sep, fill=TRUE)
    c.na <- colnames(gene.f)[-ncol(gene.f)]; r.na <- as.data.frame(gene.f[, 1])[, 1]
    df <- as.data.frame(gene.f, stringsAsFactors=FALSE)[, -1]
    rownames(df) <- r.na; colnames(df) <- c.na

    if(isRowGene==FALSE) df <- t(df); r.na <- rownames(df)
    idx <- grep("__", colnames(df)); gene2 <- df[, idx, drop=FALSE]
    gene2 <- apply(gene2, 2, as.numeric) # This step removes rownames of gene2.
    rownames(gene2) <- r.na

    ffun <- filterfun(pOverA(pOA[1], pOA[2]), cv(CV[1], CV[2]))
    filtered <- genefilter(gene2, ffun); gene2 <- gene2[filtered, ]
    id.in <- ID %in% rownames(gene2)
    if (!all(id.in)) stop(paste0(ID[!id.in], " is filtered.")) 

    # Color bar.
    if (col.bar=="all") geneV <- seq(min(gene2), max(gene2), len=1000) else if (col.bar=="selected") geneV <- seq(min(gene2[ID, , drop=FALSE]), max(gene2[ID, , drop=FALSE]), len=1000)
    col <- colorRampPalette(col.com)(length(geneV))
    cs.df <- data.frame(color_scale=geneV, y=1)
    cs.g <- ggplot()+geom_bar(data=cs.df, aes(x=color_scale, y=y), fill=col, stat="identity", width=0.2)+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=margin(3, 0.1, 3, 0.1, "cm"), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"))+coord_flip()+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand = c(0,0)); cs.grob <- ggplotGrob(cs.g)

    # SVG file conversion.
    rsvg_ps(svg, file=sub("svg$", "ps", svg))
    p1 <- sub("svg$", "ps", svg); p2 <- paste0(sub("svg$", "ps", svg), ".xml")  
    if (length(grep("~", svg))) {

        wd1 <- getwd(); setwd("~"); hm <- getwd(); setwd(wd1)
        p1 <- sub("~", hm, p1); p2 <- sub("~", hm, p2)

    }; PostScriptTrace(p1, p2) 

    grml <- xmlParse(paste0(sub("svg$", "ps", svg), ".xml")); top <- xmlRoot(grml)

    xml <- xmlParse(svg); xmltop <- xmlRoot(xml); size <- xmlSize(xmltop)

        lis.ma <- xmlSApply(xmltop[[size]], xmlAttrs)
        if (is(lis.ma, "matrix")) { id.xml <- lis.ma["id", ] } 
        else if (is(lis.ma, "list")) {

            id.xml <- NULL
            for (i in seq_len(length(lis.ma))) { id.xml <- c(id.xml, lis.ma[[i]][["id"]]) }

        }

        xml.na <- NULL
        for (i in seq_len(xmlSize(xmltop[[size]]))) { xml.na <- c(xml.na, 
        xmlName(xmltop[[size]][[i]])) } 

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

            if (xmlAttrs(top[[i]])[1]=="fill") {

                k <- k+1; chil <- xmlChildren(top[[i]])
                xy <- chil[grep("move|line", names(chil))]
                coor <- matrix(NA, nrow=length(xy), ncol=2, dimnames=list(NULL, c("x", "y")))

            for (j in seq_len(length(xy))) {

                coor[j, "x"] <- as.numeric(xmlAttrs(xy[[j]])["x"])
                coor[j, "y"] <- as.numeric(xmlAttrs(xy[[j]])["y"])

            }

            df0 <- cbind(tissue=id.xml[k], data.frame(coor), stringsAsFactors=TRUE)

            }

            df <- rbind(df, df0)

        }; g.df <- df

    # Map colours to samples according to expression level.
    cname <- colnames(gene2); con <- gsub("(.*)(__)(\\w+$)", "\\3", cname)

        g.list <- function(j) {

            g.col <- NULL; con.idx <- grep(paste0("^", j, "$"), con)
            tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]

            for (i in tis.path) {

                tis.idx <- which(tis.col1 %in% i)
                if (length(tis.idx)==1) { g.col <- c(g.col, scol1[tis.idx])
                } else if (length(tis.idx)==0) { g.col <- c(g.col, "white") }

            }

        g <- ggplot()+geom_polygon(data=g.df, aes(x=x, y=y, fill=tissue), color="black")+
        scale_fill_manual(values=g.col, guide=FALSE)+theme(axis.text=element_blank(), 
        axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=
        element_rect(fill="white", colour="grey80"), plot.margin=
        margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), 
        plot.title=element_text(hjust=0.5, size=sub.title.size))+labs(x="", y="")+
        scale_y_continuous(expand=c(0.01,0.01))+scale_x_continuous(expand=c(0.01,0.01))+
        ggtitle(paste0(k, "_", j)); return(g)

        }

        grob.na <- grob.lis <- NULL
        for (k in ID) {

            scol <- NULL
            for (i in gene2[k, ]) { 

                ab <- abs(i-geneV); col.ind <- which(ab==min(ab))[1]
                scol <- c(scol, col[col.ind])

            }

            idx <- grep("__", cname); c.na <- cname[idx]
            tis.col <- gsub("(.*)(__)(\\w+$)", "\\1", c.na); g.lis <- NULL
            con.uni <- unique(con); grob.na0 <- paste0(k, "_", con.uni)
            g.lis <- lapply(con.uni, g.list); grob <- lapply(g.lis, ggplotGrob)
            names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob) 

        }; grob.all <- c(list(cs=cs.grob), grob.lis)

    cs.arr <- arrangeGrob(grobs=list(grobTree(grob.all[[1]])), layout_matrix=cbind(1), widths=unit(15, "mm"))

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

    do.call(file.remove, list(list.files(".", "capture.*.ps")))
    do.call(file.remove, list(list.files("./tmp", ".ps$", full.names=TRUE)))
    do.call(file.remove, list(list.files("./tmp", ".ps.xml$", full.names=TRUE)))

    }
