\name{spatial.hm}
\alias{spatial.hm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Spatial heatmap.
}
\description{
 It takes a data matrix and an associated svg image to display corresponding spatial heatmaps. Take the gene expression matrix and an associated tissue svg image as an example. It can display spatial tissue heatmaps for mutltiple genes under multiple conditions. In each such tissue heatmap, the gene expression levels are represented by colours for each sample under each condition.
}
\usage{
 spatial.hm(svg, data, sep, isRowGene, pOA, CV, ID, colour, width, height, sub.title.size, layout, ncol)
}
\arguments{
  \item{svg}{
   The path of the svg image, where different regions (e.g. tissues) are labeled with 
   different colors.
   E.g.: system.file("extdata/example", "test_final.svg", package = "spatialHeatmap")
}
  \item{data}{
   The path of the data matrix. In the example of gene expression matrix, the dimension 
   names are gene IDs and sample/conditions. The sample/condition names MUST be fomatted 
   this way: a sample name is followed by double underscore then the condition, such as 
   "sample name__condition name". The meta data (e.g. gene annotation) can also be included 
   in parallel with sample/condition. In the names of sample/condition and meta data, only 
   letters, digits, single underscore, dots are allowed.   
   E.g.: system.file("extdata/example", "gene_expr_test.txt", package = "spatialHeatmap")
}
  \item{sep}{
   The seprator of the data matrix, e.g. ",", "\t", ";".
}
  \item{isRowGene}{
   It specifies if the row names are used to display spatial heatmaps. The options are 
   "TRUE" or "FALSE". For example, in a gene expression matrix genes are used to display 
   heatmaps and the gene IDs are rows names, then the option is "TRUE".
}
  \item{pOA}{
   It is a filter function that filters according to the proportion of elements exceeding A.
   The input is a two-component vector, where the first one is the proportion and the 
   second one is A, e.g.: c(0.1, 2). The default is c(0, 0), which means no filter is 
   applied. Refer to "pOverA" from the package "genefilter".
}
  \item{CV}{
   A filter function that filters according to the coefficient of variation (CV). The input
   is a two-component vector, where the first and second mean the lower and upper bound of 
   CVs used to filter, e.g.: cv(0.1, 5). The default is cv(0, 10000), which tries to aviod 
   filtering.  Refer to "cv" from the package "genefilter".
}
  \item{ID}{
   The IDs used to display the spatial heatmaps, e.g. gene IDs.
}
  \item{colour}{
   The colours components used to make the colour scale, which must be separated with comma
   and no space. The default is "green,blue,purple,yellow,red".
}
  \item{width}{
   The width of each subplot. The default is 1.
}
  \item{height}{
   The height of each subplot. The default is 1.
}
  \item{sub.title.size}{
   The size of each subtitle. The default is 11. 
}
  \item{layout}{
   The layout of the subplots. The options are "gene" or "con". For example, in gene 
   expression matrix the spatial tissue heatmaps can be organised by gene or condition (con).
}
  \item{ncol}{
   Number of columns to organise the subplots.
}
}
\details{
 Details about how to properly format and associate custom SVG images with data matrices 
 are provided here: http://biocluster.ucr.edu/~jzhan067/shiny_HM_tutorial/shiny_heatmap_tutorial.html.
}
\value{
 It generates an image of spatial heatmap(s) along with a color key.
}
\references{
 1. https://www.gimp.org/tutorials/.
 2. https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html.
 3. http://www.microugly.com/inkscape-quickguide/.
}
\author{
 Jianhai Zhang <jzhan067@ucr.edu; zhang.jianhai@hotmail.com>.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 data.path <- system.file("extdata/example", "gene_expr_test.txt", package = "spatialHeatmap")
 svg.path <- system.file("extdata/example", "test_final.svg", package = "spatialHeatmap")
 spatial.hm(svg=svg.path, data=data.path, sep="\t", isRowGene=T, pOA=c(0.1, 3), 
 CV=c(0.1, 1000), ID=c("ATMG00580", "ATCG01060"), colour=c("green", "blue", "purple", 
 "yellow", "red"), width=1, height=1, sub.title.size=11, layout="gene", ncol=3)

## The function is currently defined as
function (svg, data, sep, isRowGene, pOA, CV, ID, colour, width, 
    height, sub.title.size, layout, ncol) 
{
    library(grImport)
    library(rsvg)
    library(ggplot2)
    library(gridExtra)
    library(Cairo)
    library(grid)
    library(XML)
    library(data.table)
    library(genefilter)
    gene.f <- fread(data, header = T, sep = sep, fill = T)
    c.na <- colnames(gene.f)[-ncol(gene.f)]
    r.na <- as.data.frame(gene.f[, 1])[, 1]
    df <- as.data.frame(gene.f, stringsAsFactors = F)[, -1]
    rownames(df) <- r.na
    colnames(df) <- c.na
    if (isRowGene == F) 
        df <- t(df)
    r.na <- rownames(df)
    idx <- grep("__", colnames(df))
    gene2 <- df[, idx, drop = F]
    gene2 <- apply(gene2, 2, as.numeric)
    rownames(gene2) <- r.na
    if (missing(pOA)) 
        pOA <- pOverA(0, 0)
    if (missing(CV)) 
        CV <- cv(0, 10000)
    ffun <- filterfun(pOverA(pOA[1], pOA[2]), cv(CV[1], CV[2]))
    filtered <- genefilter(gene2, ffun)
    gene2 <- gene2[filtered, ]
    geneV <- seq(min(gene2), max(gene2), len = 1000)
    if (missing(colour)) 
        colour <- c("green", "blue", "purple", "yellow", "red")
    col <- colorRampPalette(colour)(length(geneV))
    cs.df <- data.frame(color_scale = geneV, y = 1)
    cs.g <- ggplot() + geom_bar(data = cs.df, aes(x = color_scale, 
        y = y), fill = col, stat = "identity", width = 0.2) + 
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), plot.margin = margin(3, 
                0.1, 3, 0.1, "cm"), panel.grid = element_blank(), 
            panel.background = element_rect(fill = "white", colour = "grey80")) + 
        coord_flip() + scale_y_continuous(expand = c(0, 0)) + 
        scale_x_continuous(expand = c(0, 0))
    cs.grob <- ggplotGrob(cs.g)
    rsvg_ps(svg, file = sub("svg$", "ps", svg))
    p1 <- sub("svg$", "ps", svg)
    p2 <- paste0(sub("svg$", "ps", svg), ".xml")
    if (length(grep("~", svg))) {
        wd1 <- getwd()
        setwd("~")
        hm <- getwd()
        setwd(wd1)
        p1 <- sub("~", hm, p1)
        p2 <- sub("~", hm, p2)
    }
    PostScriptTrace(p1, p2)
    grml <- xmlParse(paste0(sub("svg$", "ps", svg), ".xml"))
    top <- xmlRoot(grml)
    xml <- xmlParse(svg)
    xmltop <- xmlRoot(xml)
    size <- xmlSize(xmltop)
    lis.ma <- xmlSApply(xmltop[[size]], xmlAttrs)
    if (class(lis.ma) == "matrix") {
        id.xml <- lis.ma["id", ]
    }
    else if (class(lis.ma) == "list") {
        id.xml <- NULL
        for (i in 1:length(lis.ma)) {
            id.xml <- c(id.xml, lis.ma[[i]][["id"]])
        }
    }
    xml.na <- NULL
    for (i in 1:xmlSize(xmltop[[size]])) {
        xml.na <- c(xml.na, xmlName(xmltop[[size]][[i]]))
    }
    for (j in 1:length(xml.na)) {
        if (xml.na[j] == "g") {
            len.dif <- length(id.xml) - length(xml.na)
            g.size <- xmlSize(xmltop[[size]][[j]])
            if ((j + 1 + len.dif) <= length(id.xml)) {
                if (j == 1) {
                  id.xml <- c(paste0(id.xml[j + len.dif], "_", 
                    1:g.size), id.xml[(j + 1 + len.dif):length(id.xml)])
                }
                else if (j > 1) {
                  id.xml <- c(id.xml[1:(j - 1 + len.dif)], paste0(id.xml[j + 
                    len.dif], "_", 1:g.size), id.xml[(j + 1 + 
                    len.dif):length(id.xml)])
                }
            }
            else if ((j + 1 + len.dif) >= length(id.xml)) {
                id.xml <- c(id.xml[1:(j - 1 + len.dif)], paste0(id.xml[j + 
                  len.dif], "_", 1:g.size))
            }
        }
    }
    tis.path <- gsub("_\\d+$", "", id.xml)
    k <- 0
    df <- NULL
    for (i in 1:((xmlSize(top) - 1))) {
        if (xmlAttrs(top[[i]])[1] == "fill") {
            k <- k + 1
            chil <- xmlChildren(top[[i]])
            xy <- chil[grep("move|line", names(chil))]
            coor <- matrix(NA, nrow = length(xy), ncol = 2, dimnames = list(NULL, 
                c("x", "y")))
            for (j in 1:length(xy)) {
                coor[j, "x"] <- as.numeric(xmlAttrs(xy[[j]])["x"])
                coor[j, "y"] <- as.numeric(xmlAttrs(xy[[j]])["y"])
            }
            df0 <- cbind(tissue = id.xml[k], data.frame(coor), 
                stringsAsFactors = T)
        }
        df <- rbind(df, df0)
    }
    g.df <- df
    cname <- colnames(gene2)
    con <- gsub("(.*)(__)(\\w+$)", "\\3", cname)
    if (missing(sub.title.size)) 
        sub.title.size <- 11
    g.list <- function(j) {
        g.col <- NULL
        con.idx <- grep(paste0("^", j, "$"), con)
        tis.col1 <- tis.col[con.idx]
        scol1 <- scol[con.idx]
        for (i in tis.path) {
            tis.idx <- which(tis.col1 \%in\% i)
            if (length(tis.idx) == 1) {
                g.col <- c(g.col, scol1[tis.idx])
            }
            else if (length(tis.idx) == 0) {
                g.col <- c(g.col, "white")
            }
        }
        g <- ggplot() + geom_polygon(data = g.df, aes(x = x, 
            y = y, fill = tissue), color = "black") + scale_fill_manual(values = g.col, 
            guide = F) + theme(axis.text = element_blank(), axis.ticks = element_blank(), 
            panel.grid = element_blank(), panel.background = element_rect(fill = "white", 
                colour = "grey80"), plot.margin = margin(0.1, 
                0.1, 0.1, 0.3, "cm"), axis.title.x = element_text(size = 16, 
                face = "bold"), plot.title = element_text(hjust = 0.5, 
                size = sub.title.size)) + labs(x = "", y = "") + 
            scale_y_continuous(expand = c(0.01, 0.01)) + scale_x_continuous(expand = c(0.01, 
            0.01)) + ggtitle(paste0(k, "_", j))
        return(g)
    }
    grob.na <- grob.lis <- NULL
    for (k in ID) {
        scol <- NULL
        for (i in gene2[k, ]) {
            ab <- abs(i - geneV)
            col.ind <- which(ab == min(ab))[1]
            scol <- c(scol, col[col.ind])
        }
        idx <- grep("__", cname)
        c.na <- cname[idx]
        tis.col <- gsub("(.*)(__)(\\w+$)", "\\1", c.na)
        g.lis <- NULL
        con.uni <- unique(con)
        grob.na0 <- paste0(k, "_", con.uni)
        g.lis <- lapply(con.uni, g.list)
        grob <- lapply(g.lis, ggplotGrob)
        names(grob) <- grob.na0
        grob.lis <- c(grob.lis, grob)
    }
    grob.all <- c(list(cs = cs.grob), grob.lis)
    cs.arr <- arrangeGrob(grobs = list(grobTree(grob.all[[1]])), 
        layout_matrix = cbind(1), widths = unit(15, "mm"))
    if (missing(width)) 
        width <- 1
    if (missing(height)) 
        height <- 1
    if (layout == "gene") {
        all.cell <- ceiling(length(unique(con))/as.numeric(ncol)) * 
            as.numeric(ncol)
        cell.idx <- c(1:length(unique(con)), rep(NA, all.cell - 
            length(unique(con))))
        m <- matrix(cell.idx, ncol = as.numeric(ncol), byrow = T)
        lay <- NULL
        for (i in 1:length(ID)) {
            lay <- rbind(lay, m + (i - 1) * length(unique(con)))
        }
        g.tr <- lapply(grob.all[2:length(grob.all)], grobTree)
        n.col <- ncol(lay)
        n.row <- nrow(lay)
        g.arr <- arrangeGrob(grobs = g.tr, layout_matrix = lay, 
            widths = unit(rep(width/n.col, n.col), "npc"), heights = unit(rep(height/n.row, 
                n.row), "npc"))
        grid.arrange(cs.arr, g.arr, ncol = 2, widths = c(1.5, 
            8))
    }
    else if (layout == "con") {
        grob.all.na <- names(grob.all)
        grob.all.idx <- NULL
        for (i in unique(con)) {
            grob.all.idx <- c(grob.all.idx, grep(paste0(i, "$"), 
                grob.all.na))
        }
        grob.all.con <- grob.all[grob.all.idx]
        all.cell <- ceiling(length(ID)/as.numeric(ncol)) * as.numeric(ncol)
        cell.idx <- c(1:length(ID), rep(NA, all.cell - length(ID)))
        m <- matrix(cell.idx, ncol = as.numeric(ncol), byrow = T)
        lay <- NULL
        for (i in 1:length(unique(con))) {
            lay <- rbind(lay, m + (i - 1) * length(ID))
        }
        g.tr <- lapply(grob.all.con, grobTree)
        g.tr <- g.tr[names(grob.all.con)]
        n.col <- ncol(lay)
        n.row <- nrow(lay)
        g.arr <- arrangeGrob(grobs = g.tr, layout_matrix = lay, 
            widths = unit(rep(width/n.col, n.col), "npc"), heights = unit(rep(height/n.row, 
                n.row), "npc"))
        grid.arrange(cs.arr, g.arr, ncol = 2, widths = c(1.5, 
            8))
    }
    do.call(file.remove, list(list.files(".", "capture.*.ps")))
    do.call(file.remove, list(list.files("./tmp", ".ps$", full.name = T)))
    do.call(file.remove, list(list.files("./tmp", ".ps.xml$", 
        full.name = T)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
