#' Integrated Shiny App
#' 

#'  
 
#' @return A web browser based Shiny app.

#' @section Details:
#' No argument is required, this function launches the Shiny app directly. 

#' @examples
#' \donttest{ desired_bulk_shiny() }

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.w3schools.com/graphics/svg_intro.asp  
#'
#' https://shiny.rstudio.com/tutorial/  
#'
#' https://shiny.rstudio.com/articles/datatables.html  
#'
#' https://rstudio.github.io/DT/010-style.html  
#'
#' https://plot.ly/r/heatmaps/  
#'
#' https://www.gimp.org/tutorials/  
#'
#' https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html  
#'
#' http://www.microugly.com/inkscape-quickguide/  
#'
#' https://cran.r-project.org/web/packages/visNetwork/vignettes/Introduction-to-visNetwork.html  
#'
#' Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017). shiny: Web Application Framework for R. R package version 1.0.3. https://CRAN.R-project.org/package=shiny  
#'
#' Winston Chang and Barbara Borges Ribeiro (2017). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.6.1. https://CRAN.R-project.org/package=shinydashboard
#'
#' Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/
#'
#' Jeroen Ooms (2017). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.1. https://CRAN.R-project.org/package=rsvg
#'
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#'
#' Yihui Xie (2016). DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.2. https://CRAN.R-project.org/package=DT
#'
#' Baptiste Auguie (2016). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.2.1. https://CRAN.R-project.org/package=gridExtra
#'
#' Andrie de Vries and Brian D. Ripley (2016). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1-20. https://CRAN.R-project.org/package=ggdendro
#'
#' Langfelder P and Horvath S, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559
#'
#' Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. Journal of Statistical Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/
#'
#' Simon Urbanek and Jeffrey Horner (2015). Cairo: R graphics device using cairo graphics library for creating high-quality bitmap (PNG, JPEG, TIFF), vector (PDF, SVG, PostScript) and display (X11 and Win32) output. R package version 1.5-9. https://CRAN.R-project.org/package=Cairo
#'
#' R Core Team (2017). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#'
#' Duncan Temple Lang and the CRAN Team (2017). XML: Tools for Parsing and Generating XML Within R and S-Plus. R package version 3.98-1.9. https://CRAN.R-project.org/package=XML
#'
#' Carson Sievert, Chris Parmer, Toby Hocking, Scott Chamberlain, Karthik Ram, Marianne Corvellec and Pedro Despouy (NA). plotly: Create Interactive Web Graphics via 'plotly.js'. https://plot.ly/r, https://cpsievert.github.io/plotly_book/, https://github.com/ropensci/plotly
#'
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table
#'
#' R. Gentleman, V. Carey, W. Huber and F. Hahne (2017). genefilter: genefilter: methods for filtering genes from high-throughput experiments. R package version 1.58.1.
#'
#' Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. Journal of Statistical Software, 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/
#'
#' Almende B.V., Benoit Thieurmel and Titouan Robert (2017). visNetwork: Network Visualization using 'vis.js' Library. R package version 2.0.1. https://CRAN.R-project.org/package=visNetwork

#' @export desired_bulk_shiny
#' @importFrom shiny runApp

desired_bulk_shiny <- function() {
  # if (as.character(match.call()[[1]])=="shiny_all") warning("'shiny_all()' is deprecated and replaced by 'shiny_shm()'!", call. = FALSE)
  path <- system.file("extdata/shinyApp/desired_bulk_app.R", package="spatialHeatmap"); runApp(path)
}

