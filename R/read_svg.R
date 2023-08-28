#' Parsing annotated SVG (aSVG) files
#'
#' Parse one or multiple aSVG files and store their coordinates and related attributes in an \code{SVG} container, which will be used for creating spatial heatmap (SHM) plots.

#' @param svg.path A vector of one or multiple paths of aSVG files. If multiple aSVGs, such as aSVGs depicting organs development across mutiple times, the aSVGs should be indexed with suffixes "_shm1", "_shm2", ..., such as "arabidopsis.thaliana_organ_shm1.svg", "arabidopsis.thaliana_organ_shm2.svg". 
#' @param raster.path \itemize{ \item A vector of one or multiple paths of raster images in form of jpg or png, which are usually used as templates for creating aSVG images in \code{svg.path}. Optional (default is \code{NULL}), only applicable when superimposing raster images with SHM plots that are created from aSVG images. \item Matching raster and aSVG images is indicated by identical base names such as imageA.png and imageA.svg. The layout order in SHMs composed of multiple independent images is controlled by numbering the corresponding file pairs accordingly such as imageA_1.png and imageA_1.svg, imageB_2.png and imageB_2.svg, etc.}
#' @param cores Number of CPUs to parse the aSVG files (default is \code{1}). 
#' @param srsc Logical. If `TRUE`, the aSVG is considered for co-visualizing spatially resolved single-cell and bulk data, and the rotation angle of the tissue section for spatial assays will be recorded.   
#' @return An object of \code{SVG} class, containing one or multiple aSVG instances. 

#' @seealso
#' \code{\link{SVG}}: the \code{SVG} class. 

#' @examples
#' 
#' # The first raste image used as a template to create an aSVG. 
#' raster.pa1 <- system.file('extdata/shinyApp/data/maize_leaf_shm1.png',
#' package='spatialHeatmap')
#' # The first aSVG created with the first template. 
#' svg.pa1 <- system.file('extdata/shinyApp/data/maize_leaf_shm1.svg',
#' package='spatialHeatmap')
#' # The second raster image used as a template to create an aSVG. 
#' raster.pa2 <- system.file('extdata/shinyApp/data/maize_leaf_shm2.png',
#' package='spatialHeatmap')
#' # The second aSVG created with the second template. 
#' svg.pa2 <- system.file('extdata/shinyApp/data/maize_leaf_shm2.svg',
#' package='spatialHeatmap')
#'
#' # Parse these two aSVGs without association with raster images.
#' svgs <- read_svg(svg.path=c(svg.pa1, svg.pa2), raster.path=NULL)

#' # Parse these two aSVGs. The raster image paths are provide so as to
#' # be associated with respective aSVGs, which will be used when 
#' # superimposing raster images with SHM plots.
#' svgs <- read_svg(svg.path=c(svg.pa1, svg.pa2), raster.path=c(raster.pa1, raster.pa2))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2

#' @export
#' @importFrom xml2 read_html

read_svg <- function(svg.path, raster.path=NULL, cores=1, srsc=FALSE) {
    # Get SVG/raster names and order their path/name if there are multiple images.
    svg.pa.na <- img_pa_na(path.na=svg.path); svg.path <- svg.pa.na$path; svg.na <- svg.pa.na$na
    raster.na <- NULL; if (length(raster.path)==0) raster.path <- NULL
    if (is.character(raster.path)) { 
      raster.pa.na <- img_pa_na(raster.path); raster.pa <- raster.pa.na$path; raster.na <- raster.pa.na$na
      # Ensure aSVGs and rasters are paired. After this step, the aSVG and raster are paired, and if multiple aSVGs exist, aSVGs and templates are ordered by shm1, shm2, etc.
      msg <- svg_raster(c(svg.na, raster.na), raster.ext=c('.jpg', '.JPG', '.png', '.PNG'), shiny=FALSE)
      if (!is.null(msg)) stop(msg)
    }
    # Coordinates of each SVG are extracted and placed in a list.
    cordn <- attrb <- dimen <- svg.l <- raster <- angle <- NULL
    for (i in seq_along(svg.na)) {
      cat('Parsing:', svg.na[i], '...', '\n') # '... \n' renders two new lines.
      cores <- deter_core(cores, svg.pa=svg.path[i]); cat('CPU cores:', cores, '\n')
      svg.lis <- svg_df(svg.path=svg.path[i], feature=NULL, cores=cores, srsc=srsc)
      if (is.character(svg.lis)) {
        msg <- paste0(svg.na[i], ': ', svg.lis)
        warning(msg); return(msg) 
      }
      cordn <- c(cordn, list(svg.lis$coordinate))
      attrb <- c(attrb, list(svg.lis$attribute))
      dimen <- c(dimen, list(svg.lis$dimension))
      raster <- c(raster, list(raster.path[i]))
      angle <- c(angle, list(svg.lis$angle))
      svg.l <- c(svg.l, list(svg.lis$svg.obj))
    }
    names(cordn) <- names(attrb) <- names(dimen) <- names(raster) <- names(svg.l) <- names(angle) <- svg.na
    svg.all <- SVG(coordinate=cordn, attribute=attrb, dimension=dimen, svg=svg.l, raster=raster, angle=angle); svg.all
}



  

