\name{spatialHeatmap-package}
\alias{spatialHeatmap-package}
\alias{spatialHeatmap}
\docType{package}
\title{
\packageTitle{spatialHeatmap}
 Spatial Heatmap, Matrix Heatmap, Network
}
\description{
\packageDescription{spatialHeatmap}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{spatialHeatmap}
\packageIndices{spatialHeatmap}

This R package spatialHeatmap is designed to allow users to visualise custom data matrix on SVG images in an intuitive way, and therefore promote hypothesis generation. In biological research, examples include but not limit to gene expression profiling by RNA-seq or microarray, qPCR, protein subcellular localisation, etc., regardless of species.

The core feature “Spatial Heatmap” is to map expression profile of a target gene under different conditions to different cells/tissues/organs (samples) on an associated SVG image, where different samples are predefined. After mapping, the expression profile is represented as different colours across samples in the image, which are called spatial heatmaps. In addition, the accessory features of “Matrix Heatmap” and “Network” display the target gene in the context of corresponding gene module, which make this package more versatile. All the utilities implemented in this R package are also combined as a web-browser based Shiny app, which can be used locally or online on the Shiny server.

The application of this package is not limited to gene expression data. It can be used as long as a data matrix and a configured SVG image are provided. Other exmaples include population data collected in different years across different cities, health data of an individual under different conditions, ect.

}
\author{
\packageAuthor{spatialHeatmap}
Author: \packageAuthor{spatialHeatmap}
  Jianhai Zhang (PhD candidate Genetics, Genomics and Bioinformatics, University of California, Riverside), Dr. Thomas Girke (professor at Department of Botany and Plant Sciences, University of California, Riverside)
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
 \code{\link{filter_data}}, \code{\link{adj_mod}}, \code{\link{spatial_hm_all}}, \code{\link{spatial_hm}}, \code{\link{matrix_heatmap}}, \code{\link{network}}
}
\examples{
# Creat the "SummarizedExperiment" class. Refer to the R package "SummarizedExperiment" for more details.
data.path <- system.file("extdata/shinyApp/example", "root_expr_row_gen.txt", package = "spatialHeatmap")   
## The expression matrix, where the row and column names should be gene IDs and sample/conditions respectively. This data matrix is truncated from a GEO dataset (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46205), which is already normalised.
library(data.table); expr <- fread(data.path, sep='\t', header=TRUE, fill=TRUE)
col.na <- colnames(expr)[-ncol(expr)]; row.na <- as.data.frame(expr[, 1])[, 1]
expr <- as.matrix(as.data.frame(expr, stringsAsFactors=FALSE)[, -1])
rownames(expr) <- row.na; colnames(expr) <- col.na
col.met.path <- system.file("extdata/shinyApp/example", "col_metadata.txt", package = "spatialHeatmap") 
## Condition metadata is data frame. It has a column of tissues and a column of contidions, which correspond to columns of the data matrix "expr".
col.metadata <- read.table(col.met.path, header=TRUE, row.names=NULL, sep='\t', stringsAsFactors=FALSE)
row.met.path <- system.file("extdata/shinyApp/example", "row_metadata.txt", package = "spatialHeatmap")
## Row metadata is a data frame. It has a column of row (gene) annotations, which correspond to rows of the data matrix "expr".
row.metadata <- read.table(row.met.path, header=TRUE, row.names=1, sep='\t', stringsAsFactors=FALSE)
## The expression matrix, row metadata, and column metadata are stored in a "SummarizedExperiment" object. The row metadata is optional while column metadata is mandatory. The column names in the expression matrix are not important, since they are ultimately renewed by column metadata.
library(SummarizedExperiment); expr <- as.matrix(expr); se <- SummarizedExperiment(assays=list(expr=expr), rowData=row.metadata, colData=col.metadata)  
exp <- filter_data(data=se, pOA=c(0, 0), CV=c(0.1, Inf), ann='ann', samples='sample', conditions='condition', dir=NULL) 

# The svg image path.
svg.path <- system.file("extdata/shinyApp/example", "root_cross_final.svg", package="spatialHeatmap")
# The expression profiles of gene "PSAC" and "NDHG" under different conditions are mapped to tissues defined in the SVG image in the form of different colours. 
spatial_hm(svg=svg.path, data=exp, samples='sample', conditions='condition', pOA=c(0, 0), CV=c(0, 10000), ID=c("PSAC", "NDHG"), col.com=c("yellow", "blue", "purple"), width=1, height=1, sub.title.size=11, layout="gene", ncol=3)

# Compute adjacency matrix and identify modules
adj.mod <- adj_mod(data=exp, type="signed", minSize=15, dir=NULL)

# The gene "PSAC" is represented in the context of its gene module in the form of a static matrix heatmap.
matrix_heatmap(geneID="PSAC", data=exp, adj.mod=adj.mod, ds="2", scale="row", static=TRUE)

# The gene "PSAC" is represented in the context of its network module in the form of an interactive network.
network(geneID="PSAC", data=exp, ann='ann', adj.mod=adj.mod, static=TRUE)

# Launch the integrated Shiny app.
\donttest{ spatial_hm_all() }

}
