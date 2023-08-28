#' Exporting each spatial heatmap to a separate SVG file
#' 
#' This function exports each spatial heatmap (associated with a specific gene and condition) as separate SVG files. In contrast to the original aSVG file, spatial features in the output SVG files are assigned heat colors.   

#' @param input The output returned by \link{shm} or \link{covis}. 
#' @param out.dir The directory path where the colored aSVG file will be saved.
#' @return Nothing is returned. 

#' @examples
#' 
#' # Read the aSVG file.
#' svg.hum.pa <- system.file("extdata/shinyApp/data", 'homo_sapiens.brain.svg', 
#' package="spatialHeatmap")
#' svg.hum <- read_svg(svg.hum.pa)
#' # Attributes of spatial features.
#' feature.hum <- attribute(svg.hum)[[1]]
#' set.seed(20) # To obtain reproducible results, a fixed seed is set.
#' # Spatial features.
#' unique(feature.hum$feature)[1:10]
#' # Create a random numeric vector.
#' my_vec <- setNames(sample(1:100, 4), c('substantia.nigra', 'putamen', 
#' 'prefrontal.cortex', 'notMapped'))
#' my_vec
#' # Plot spatial heatmaps with the imported aSVG and random vector.
#' dat.quick <- SPHM(svg=svg.hum, bulk=my_vec)
#' shm.res <- shm(data=dat.quick, ID='testing', ncol=1, sub.title.size=20, legend.nrow=3,
#' bar.width=0.1)
#' # Export each spatial heatmap (under a certain gene and condition) to a separate SVG 
#' # file in a temporary directory.
#' write_svg(input=shm.res, out.dir=tempdir(check = TRUE))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references 
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2 

#' @export
#' @importFrom xml2 write_xml xml_unserialize 

write_svg <- function(input, out.dir) {
  ID <- condition <- feature <- index.sub <- NULL
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  svg.obj.lis <- svg_obj(svg(input))
  svg.na <- names(svg.obj.lis)[1]
  svg.obj <- svg.obj.lis[[1]]
  df.att <- attribute(svg(input))[[1]]
  map <- output(input)$mapped_feature
  # Genes and variables.
  ids <- unique(map$ID); vars <- unique(map$condition)
  len.var <- length(vars)>0
  if (!len.var) vars <- 'con' # In order to run the for loop.
  # Coloring SVG for each gene and each variable.
  for (i in ids) {
  for (j in vars) {
    # One gene + one variable.
    if (len.var) df0 <- subset(map, ID %in% i & condition %in% j)
    if (!len.var) df0 <- subset(map, ID %in% i)
    # Attributes of mapped features.
    df.att0 <- subset(df.att, feature %in% df0$feature)
    df.att0 <- df.att0[, c('feature', 'stroke', 'fill', 'id', 'element', 'parent', 'index.sub')]
    df.att0 <- subset(df.att0, !duplicated(index.sub))
    colnames(df.att0)[colnames(df.att0)=='index.sub'] <- 'order'
    # Output SVG file name.
    if (len.var) na <- paste0(i, '_', j, '_', svg.na)
    if (!len.var) na <- paste0(i, '_', svg.na)
    write_xml(xml_unserialize(svg.obj), file=file.path(out.dir, na)) 
    df.new0 <- cbind(colorNew=df0$fill, df.att0, SVG=na)
    # Color output SVG file.
    update_feature(df.new=df.new0, dir=out.dir) 
  }
  }
}



