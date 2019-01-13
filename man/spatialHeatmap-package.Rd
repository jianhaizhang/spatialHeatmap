\name{spatialHeatmap-package}
\alias{spatialHeatmap-package}
\alias{spatialHeatmap}
\docType{package}
\title{
\packageTitle{spatialHeatmap}
 Interactive data matrix, spatial heatmap, matrix heatmap, and network
}
\description{
\packageDescription{spatialHeatmap}
 The Spatial Heatmap is an integrated app, which includes interactive data matrix, interactive spatial heatmap, interactive matrix heatmap, and interactive network.
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{spatialHeatmap}
\packageIndices{spatialHeatmap}
 This Shiny App based spatialHeatmap can be used for interactive visualisation as long as a data matrix and an associated SVG image are available, This tool has three main functionality, which makes it informative, interactive, and versatile in visualizing large-scale data. First, it generates spatial heatmaps where user defined tissue regions are coloured by the expression profile of a gene of interest. The gene expression information is uploaded as a table matrix and the associated tissue image is uploaded as an SVG format. Second, the app computes gene network modules based on gene expression profiles across tissue samples. It uses a matrix heatmap to visualize the expression of a chosen gene in the context of the corresponding gene network module where the chosen gene is contained. Third, the app displays the same network module in the matrix heatmap in the form of an interactive network graph. The network module identification is computationally intense. To maintain performance, the large gene expression matrix (e.g.: > 10,000 genes) should be computed by the "Compute locally" mode, while small ones (e.g.: < 10,000 genes) can be computed by the "Compute online" mode. If the app times out at some time, then users want to refresh the page in their web browser.
}
\author{
\packageAuthor{spatialHeatmap}

Maintainer: \packageMaintainer{spatialHeatmap}
  Jianhai Zhang <jzhan067@ucr.edu; zhang.jianhai@hotmail.com>.
}
\references{
  https://www.w3schools.com/graphics/svg_intro.asp  
  https://shiny.rstudio.com/tutorial/  
  https://shiny.rstudio.com/articles/datatables.html  
  https://rstudio.github.io/DT/010-style.html  
  https://plot.ly/r/heatmaps/  
  https://www.gimp.org/tutorials/  
  https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html  
  http://www.microugly.com/inkscape-quickguide/  
  https://cran.r-project.org/web/packages/visNetwork/vignettes/Introduction-to-visNetwork.html  

  Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017).
  shiny: Web Application Framework for R. R package version 1.0.3.
  https://CRAN.R-project.org/package=shiny

  Winston Chang and Barbara Borges Ribeiro (2017). shinydashboard: Create
  Dashboards with 'Shiny'. R package version 0.6.1.
  https://CRAN.R-project.org/package=shinydashboard

  Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R.
  Journal of Statistical Software, 30(4), 1-37. URL
  http://www.jstatsoft.org/v30/i04/.

  Jeroen Ooms (2017). rsvg: Render SVG Images into PDF, PNG, PostScript, or
  Bitmap Arrays. R package version 1.1. https://CRAN.R-project.org/package=rsvg

  H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New
  York, 2016.

  Yihui Xie (2016). DT: A Wrapper of the JavaScript Library 'DataTables'. R
  package version 0.2. https://CRAN.R-project.org/package=DT

  Baptiste Auguie (2016). gridExtra: Miscellaneous Functions for "Grid" Graphics.
  R package version 2.2.1. https://CRAN.R-project.org/package=gridExtra

  Andrie de Vries and Brian D. Ripley (2016). ggdendro: Create Dendrograms and
  Tree Diagrams Using 'ggplot2'. R package version 0.1-20.
  https://CRAN.R-project.org/package=ggdendro

  Langfelder P and Horvath S, WGCNA: an R package for weighted correlation
  network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559

  Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
  Correlations and Hierarchical Clustering. Journal of Statistical Software,
  46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.

  Simon Urbanek and Jeffrey Horner (2015). Cairo: R graphics device using cairo
  graphics library for creating high-quality bitmap (PNG, JPEG, TIFF), vector
  (PDF, SVG, PostScript) and display (X11 and Win32) output. R package version
  1.5-9. https://CRAN.R-project.org/package=Cairo

  R Core Team (2017). R: A language and environment for statistical computing. R
  Foundation for Statistical Computing, Vienna, Austria. URL
  https://www.R-project.org/.

  Duncan Temple Lang and the CRAN Team (2017). XML: Tools for Parsing and
  Generating XML Within R and S-Plus. R package version 3.98-1.9.
  https://CRAN.R-project.org/package=XML

  Carson Sievert, Chris Parmer, Toby Hocking, Scott Chamberlain, Karthik Ram,
  Marianne Corvellec and Pedro Despouy (NA). plotly: Create Interactive Web
  Graphics via 'plotly.js'. https://plot.ly/r,
  https://cpsievert.github.io/plotly_book/, https://github.com/ropensci/plotly.

  Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R
  package version 1.10.4. https://CRAN.R-project.org/package=data.table

  R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter:
  methods for filtering genes from high-throughput experiments. R package version
  1.58.1.

  Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
  Correlations and Hierarchical Clustering. Journal of Statistical Software,
  46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.

  Almende B.V., Benoit Thieurmel and Titouan Robert (2017). visNetwork: Network
  Visualization using 'vis.js' Library. R package version 2.0.1.
  https://CRAN.R-project.org/package=visNetwork
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%%~~ documentation directory ~~
\keyword{ package }
\seealso{
 \code{\link{filter.data}}, \code{\link{adj.mod}}, \code{\link{spatial.hm.all}}, \code{\link{spatial.hm}}
}
\examples{
#Get the path of SVG file and data matrix.
svg.path <- system.file("extdata/example", "root_cross_final.svg", package = "spatialHeatmap")
data.path <- system.file("extdata/example", "root_expr_ann_row_gen.txt", package = "spatialHeatmap")

#Local work mode for large data matrix.

## Filter the data matrix. The resulting data matrix "processed_data.txt" is automatically saved in the directory "local_mode_result", which should be uploaded to the app on web browser.
exp <- filter.data(data=data.path, sep="\t", isRowGen=TRUE, c(0, 0), c(0.1, 10000),
 "processed_data") 

## Compute the adjacency matrix and identify modules. The resulting adjacency file "adj.txt" and "mod.txt" are also automatically saved in "local_mode_result", which should be uploaded to the app on web browser.
adj_mod <- adj.mod(data=exp, type="signed", minSize=20)

# Lauch of the app.
\donttest{ spatial.hm.all() }

# Plot the spatial heatmaps on R console.
spatial.hm(svg=svg.path, data=data.path, sep="\t", isRowGene=TRUE, pOA=c(0.1, 3), 
CV=c(0.05, 1000), ID=c("PSAC", "NDHE"), colour=c("green", "blue", "purple", "yellow", "red"), width=1, height=1, sub.title.size=11, layout="gene", ncol=3)

}
