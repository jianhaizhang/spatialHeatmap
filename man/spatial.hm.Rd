\name{spatial.hm}
\alias{spatial.hm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Spatial heatmap
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
   The path of the data matrix. In the example of gene expression matrix, the dimension names are gene IDs and sample/conditions. The sample/condition names MUST be fomatted this way: a sample name is followed by double underscore then the condition, such as "sample name__condition name". The meta data (e.g. gene annotation) can also be included in parallel with sample/condition. In the names of sample/condition and meta data, only letters, digits, single underscore, dots are allowed.   
   E.g.: system.file("extdata/example", "gene_expr_test.txt", package = "spatialHeatmap")
}
  \item{sep}{
   The seprator of the data matrix, e.g. ",", "\\t", ";".
}
  \item{isRowGene}{
   It specifies if the row names are used to display spatial heatmaps. The options are "TRUE" or "FALSE". For example, in a gene expression matrix genes are used to display heatmaps and the gene IDs are rows names, then the option is "TRUE".
}
  \item{pOA}{
   It is a filter function that filters according to the proportion of elements exceeding A. The input is a two-component vector, where the first one is the proportion and the second one is A, e.g.: c(0.1, 2). The default is c(0, 0), which means no filter is applied. Refer to "pOverA" from the package "genefilter".
}
  \item{CV}{
   A filter function that filters according to the coefficient of variation (CV). The input is a two-component vector, where the first and second mean the lower and upper bound of CVs used to filter, e.g.: cv(0.1, 5). The default is cv(0, 10000), which tries to aviod filtering.  Refer to "cv" from the package "genefilter".
}
  \item{ID}{
   The IDs used to display the spatial heatmaps, e.g. gene IDs.
}
  \item{colour}{
   The colours components used to make the colour scale, which must be separated with comma and no space. The default is "green,blue,purple,yellow,red".
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
   The layout of the subplots. The options are "gene" or "con". For example, in gene expression matrix the spatial tissue heatmaps can be organised by gene or condition (con).
}
  \item{ncol}{
   Number of columns to organise the subplots.
}
}
\details{
 Details about how to properly format and associate custom SVG images with data matrices are provided here: http://biocluster.ucr.edu/~jzhan067/shiny_HM_tutorial/shiny_heatmap_tutorial.html.
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
 data.path <- system.file("extdata/example", "gene_expr_ann_row_gen.txt", package = "spatialHeatmap")
 svg.path <- system.file("extdata/example", "test_final.svg", package = "spatialHeatmap")
 spatial.hm(svg=svg.path, data=data.path, sep="\t", isRowGene=TRUE, pOA=c(0.1, 3), 
 CV=c(0.05, 1000), ID=c("244902_at", "244903_at"), colour=c("green", "blue", "purple", "yellow", "red"), width=1, height=1, sub.title.size=11, layout="gene", ncol=3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
