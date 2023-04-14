\name{spatialHeatmap-package}
\alias{spatialHeatmap-package}
\alias{spatialHeatmap}
\docType{package}
\title{
\packageTitle{spatialHeatmap}
 Spatial Heatmap, Spatial Enrichment, Data Mining, Co-visualization
}
\description{
\packageDescription{spatialHeatmap}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{spatialHeatmap}
\packageIndices{spatialHeatmap}

The spatialHeatmap package provides functionalities for visualizing cell-, tissue- and organ-specific data of biological assays by coloring the corresponding spatial features defined in anatomical images according to a numeric color key. The color scheme used to represent the assay values can be customized by the user. This core functionality is called a spatial heatmap plot. It is enhanced with nearest neighbor visualization tools for groups of measured items (e.g. gene modules) sharing related abundance profiles, including matrix heatmaps combined with hierarchical clustering dendrograms and network representations. The functionalities of spatialHeatmap can be used either in a command-driven mode from within R or a graphical user interface (GUI) provided by a Shiny App that is also part of this package. While the R-based mode provides flexibility to customize and automate analysis routines, the Shiny App includes a variety of convenience features that will appeal to many biologists. Moreover, the Shiny App has been designed to work on both local computers as well as server-based deployments (e.g. cloud-based or custom servers) that can be accessed remotely as a centralized web service for using spatialHeatmap's functionalities with community and/or private data.

As anatomical images the package supports both tissue maps from public repositories and custom images provided by the user. In general any type of image can be used as long as it can be provided in SVG (Scalable Vector Graphics) format, where the corresponding spatial features have been defined (see aSVG below). The numeric values plotted onto a spatial heatmap are usually quantitative measurements from a wide range of profiling technologies, such as microarrays, next generation sequencing (e.g. RNA-Seq and scRNA-Seq), proteomics, metabolomics, or many other small- or large-scale experiments. For convenience, several preprocessing and normalization methods for the most common use cases are included that support raw and/or preprocessed data. Currently, the main application domains of the spatialHeatmap package are numeric data sets and spatially mapped images from biological and biomedical areas. Moreover, the package has been designed to also work with many other spatial data types, such a population data plotted onto geographic maps. This high level of flexibility is one of the unique features of spatialHeatmap. Related software tools for biological applications in this field are largely based on pure web applications (Winter et al. 2007; Waese et al. 2017) or local tools (Maag 2018; Muschelli, Sweeney, and Crainiceanu 2014) that typically lack customization functionalities. These restrictions limit users to utilizing pre-existing expression data and/or fixed sets of anatomical image collections. To close this gap for biological use cases, we have developed spatialHeatmap as a generic R/Bioconductor package for plotting quantitative values onto any type of spatially mapped images in a programmable environment and/or in an intuitive to use GUI application.

}
\author{
\packageAuthor{spatialHeatmap}
Author: \packageAuthor{spatialHeatmap}
  Jianhai Zhang (PhD candidate at Genetics, Genomics and Bioinformatics, University of California, Riverside), Dr. Thomas Girke (Professor at Department of Botany and Plant Sciences, University of California, Riverside)
Maintainer: \packageMaintainer{spatialHeatmap}
  Jianhai Zhang <jzhan067@@ucr.edu>.
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

  https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg 

  Winter, Debbie, Ben Vinegar, Hardeep Nahal, Ron Ammar, Greg V Wilson, and Nicholas J Provart. 2007. 
  "An 'Electronic Fluorescent Pictograph' Browser for Exploring and Analyzing Large-Scale Biological Data Sets." PLoS One 2 (8): e718

  Waese, Jamie, Jim Fan, Asher Pasha, Hans Yu, Geoffrey Fucile, Ruian Shi, Matthew Cumming, et al. 2017. 
  "EPlant: Visualizing and Exploring Multiple Levels of Data for Hypothesis Generation in Plant Biology." Plant Cell 29 (8): 1806–21

  Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angelica Liechti, et al. 2019.   
  "Gene Expression Across Mammalian Organ Development." Nature 571 (7766): 505–9

  Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas

  Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8

  McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97

  Maag, Jesper L V. 2018. "Gganatogram: An R Package for Modular Visualisation of Anatograms and Tissues Based on Ggplot2." F1000Res. 7 (September): 1576

  Muschelli, John, Elizabeth Sweeney, and Ciprian Crainiceanu. 2014. "BrainR: Interactive 3 and 4D Images of High Resolution Neuroimage Data." R J. 6 (1): 41–48

  Morgan, Martin, Valerie Obenchain, Jim Hester, and Hervé Pagès. 2018. SummarizedExperiment: SummarizedExperiment Container

  Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017).
  shiny: Web Application Framework for R. R package version 1.0.3.
  https://CRAN.R-project.org/package=shiny

  Winston Chang and Barbara Borges Ribeiro (2017). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.6.1.
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

  Lori Shepherd and Martin Morgan (2020). BiocFileCache: Manage Files Across 
  Sessions. R package version 1.12.1.
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%%~~ documentation directory ~~
\keyword{ spatial heatmap}
\seealso{
 \code{\link{norm_data}}, \code{\link{aggr_rep}}, \code{\link{filter_data}}, \code{\link{spatial_hm}}, \code{\link{submatrix}}, \code{\link{adj_mod}}, \code{\link{matrix_hm}}, \code{\link{network}}, \code{\link{return_feature}}, \code{\link{update_feature}}, \code{\link{shiny_shm}}, \code{\link{custom_shiny}}
}
\examples{

#' ## The example data included in this package come from an RNA-seq analysis on 
#' ## development of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019). 
#' ## The complete raw count data are downloaded using the R package ExpressionAtlas
#' ## (Keays 2019) with the accession number "E-MTAB-6769". 
#'
#' # Access example count data. 
#' count.chk <- read.table(system.file('extdata/shinyApp/data/count_chicken.txt', 
#' package='spatialHeatmap'), header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing spatial features and variables is made based on the 
#' # experiment design.
#' target.chk <- read.table(system.file('extdata/shinyApp/data/target_chicken.txt', 
#' package='spatialHeatmap'), header=TRUE, row.names=1, sep='\t')
#' # Every column in example data 2 corresponds with a row in the targets file. 
#' target.chk[1:5, ]
#' # Store example data in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#'
#' # Normalize data.
#' se.chk.nor <- norm_data(data=se.chk, norm.fun='CNF', log2.trans=TRUE)
#'
#' # Aggregate replicates of "spatialFeature_variable", where spatial features are organs
#' # and variables are ages.
#' se.chk.aggr <- aggr_rep(data=se.chk.nor, sam.factor='organism_part', con.factor='age',
#' aggr='mean')
#' assay(se.chk.aggr)[1:3, 1:3]
#'
#' # Genes with experssion values >= 5 in at least 1% of all samples (pOA), and coefficient
#' # of variance (CV) between 0.2 and 100 are retained.
#' se.chk.fil <- filter_data(data=se.chk.aggr, sam.factor='organism_part', con.factor='age', 
#' pOA=c(0.01, 5), CV=c(0.2, 100), file=NULL)
#' 
#' # The chicken aSVG downloaded from the EBI aSVG repository (https://github.com/ebi-gene-
#' # expression-group/anatomogram/tree/master/src/svg) is included in this package and 
#' # accessed as below.
#' svg.chk <- system.file("extdata/shinyApp/data", "gallus_gallus.svg",
#' package="spatialHeatmap")
#' # Read the chicken aSVG file.
#' svg.chk <- read_svg(svg.path=svg.chk)
#'
#' # Store assay data and aSVG in an "SHM" class.
#' dat.chk <- SPHM(svg=svg.chk, bulk=se.chk.fil)
#' # Plot spatial heatmaps with gene "ENSGALG00000019846".
#' shm(data=dat.chk, ID='ENSGALG00000019846', legend.r=1.9,
#' legend.nrow=5, sub.title.size=7, ncol=3)
#'
#' # Save SHMs as HTML and video files in the "~/test" directory. 
#' \donttest{
#' if (!dir.exists('~/test')) dir.create('~/test')
#' shm(data=dat.chk, ID='ENSGALG00000019846', legend.r=1.9,
#' legend.nrow=5, sub.title.size=7, ncol=3, out.dir='~/test')
#' }
#'
}
