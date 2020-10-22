#' Integrated Shiny App
#' 
#' In additon to generating spatial heatmaps and corresponding item (genes, proteins, metabolites, \emph{etc.}) context plots from R, spatialHeatmap includes a Shiny App (\url{https://shiny.rstudio.com/}) that provides access to the same functionalities from an intuitive-to-use web browser interface. Apart from being very user-friendly, this App conveniently organizes the results of the entire visualization workflow in a single browser window with options to adjust the parameters of the individual components interactively. Upon launched, the app automatically displays a pre-formatted example. 

#' To use this app, the data matrix (\emph{e.g.} gene expression matrix) and aSVG image are uploaded as tabular text (\emph{e.g.} in CSV or TSV format) and SVG file, respectively. To also allow users to upload data matrix stored in \code{SummarizedExperiment} objects, one can export them from R to a tabular file with the \code{\link{filter_data}} function. In this function call, the user sets a desired directory path under \code{dir}. Within this directory the tabular file will be written to "customComputedData/sub_matrix.txt" in TSV format. The column names in the exported tabular file preserve the experimental design information from the \code{colData} slot by concatenating the corresponding sample and condition information separated by double underscores. To interactively view functional descriptions by moving the cursor over network nodes, the corresponding annotation column needs to be present in the \code{rowData} slot and its column name assigned to the \code{ann} argument. In the exported tabular file the extra annotation column is appended to the expression matrix. See function \code{\link{filter_data}} for details.

#' If the subsetted data matrix in the Matrix Heatmap is too large, \emph{e.g.} >10,000 rows, the "customComputedData" under "Step 1: data sets" is recommended. Since this subsetted matrix is fed to the Network, and the internal computation of adjacency matrix and module identification would be intensive. In order to protect the app from crash, the intensive computation should be performed outside the app, then upload the results under "customComputedData". When using "customComputedData", the data matrix to upload is the subsetted matrix "sub_matrix.txt" generated with \code{\link{submatrix}}, which is a TSV-tabular text file. The adjacency matrix and module assignment to upload are "adj.txt" and "mod.txt" generated in function \code{\link{adj_mod}} respectively. Note, "sub_matrix.txt", "adj.txt", and "mod.txt" are downstream to the same call on \code{\link{filter_data}}, so the three files should not be mixed between different filtering when uploading. See the instruction page in the app for details.

#' The large matrix issue could be resolved by increasing the subsetting strigency to get smaller matrix in \code{\link{submatrix}} in most cases. Only in rare cases users cannot avoid very large subsetted matrix, the "customComputedData" is recommended. 
 
#' @return A web browser based Shiny app.

#' @section Details:
#' No argument is required, this function launches the Shiny app directly. 

#' @examples
#' \donttest{ shiny_all() }

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

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

#' @export
#' @importFrom shiny runApp

shiny_all <- function() {

    path <- system.file("extdata/shinyApp", package="spatialHeatmap")
    runApp(path)

}

