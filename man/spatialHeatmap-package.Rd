\name{spatialHeatmap-package}
\alias{spatialHeatmap-package}
\alias{spatialHeatmap}
\docType{package}
\title{
\packageTitle{spatialHeatmap}
 Spatial heatmap, matrix heatmap, network
}
\description{
\packageDescription{spatialHeatmap}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{spatialHeatmap}
\packageIndices{spatialHeatmap}

 This spatialHeatmap is designed for visualising gene expression data such as sequencing, microarray, etc.. The core feature is to map expression values of target genes under different conditions to different tissues, which are called spatial heatmaps. To make it more versatile, the matrix heatmap and network features were also developed, which display the target gene in the context of corresponding gene module. All the utilities implemented in this R package are also combined as a web-browser based Shiny app (Chang et al., n.d.; Chang and Borges Ribeiro 2018).

The application of this package is not limited to gene expression data. It can be used as long as a data matrix and an associated SVG image are provided, such as population data generated in different years across different cities. The tutorial of how to properly format and associate a custom SVG image with a data matrix is provided in an independent vignette “SVG_tutorial.html”.

For a quick test, users can launch the combined Shiny app then either select the default examples from the left menu or download the example data matrix and svg image in the insturction page.
}
\author{
\packageAuthor{spatialHeatmap}
Author: \packageAuthor{spatialHeatmap}
  Jianhai Zhang (PhD student at Genetics, Genomics and Bioinformatics, University of California, Riverside), Dr. Thomas Girke (professor at Department of Botany and Plant Sciences, University of California, Riverside)
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
\keyword{ spatial heatmap}
\seealso{
 \code{\link{filter.data}}, \code{\link{adj.mod}}, \code{\link{spatial.hm.all}}, \code{\link{spatial.hm}}, \code{\link{matrix.heatmap}}, \code{\link{network}}
}
\examples{
# Get the path of svg image and associated data matrix.
svg.path <- system.file("extdata/example", "root_cross_final.svg", package = "spatialHeatmap")
data.path <- system.file("extdata/example", "root_expr_ann_row_gen.txt", package = "spatialHeatmap")

# Filter the data matrix. The directory "local_mode_result/" is automatically created in the location specified by "dir" and is used to save filtered data matrix "processed_data.txt".
exp <- filter.data(data=data.path, sep="\t", isRowGen=TRUE, pOA=c(0, 0), CV=c(0.1, 10000), dir="./") 

# Plot the spatial heatmaps for gene "PSAC" and "NDHG".
spatial.hm(svg=svg.path, data=data.path, sep="\t", isRowGene=TRUE, pOA=c(0.1, 3), CV=c(0.05, 1000), ID=c("PSAC", "NDHG"), col.com=c("green", "blue", "purple", "yellow", "red"), width=1, height=1, sub.title.size=11, layout="gene", ncol=3)

# Compute the adjacency matrix and identify modules. The directory "local_mode_result/" is automatically created in the location specified by "dir" and is used to save the resulting adjacency file "adj.txt" and module definition file "mod.txt".
adj_mod <- adj.mod(data=exp, type="signed", minSize=20, dir="./")

# The gene "PSAC" is represented in the context of its network module in the form of an interactive matrix heatmap.
\donttest{ matrix.heatmap(geneID="PSAC", data=exp, adj.mod=adj_mod, ds="2", scale="row") }

# The gene "PSAC" is represented in the context of its network module in the form of an interactive network.
\donttest{ network(geneID="PSAC", data=exp, adj.mod=adj_mod) }

# Lauch of the integrated Shiny app.
\donttest{ spatial.hm.all() }

}
