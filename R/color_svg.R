#' Coloring spatial features in aSVG files.
#' 
#' This function copies the original aSVG file into a specified directory, and colors the spatial features in the new aSVG file with heat colors in spatial heatmaps. 

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
#' # Color the spatial features in the original aSVG file with heat colors in the spatial
#' # heatmap. The output aSVG file is saved in a temporary directory.
#' color_svg(input=shm.res, out.dir=tempdir(check = TRUE))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @export

color_svg <- function(input, out.dir) {
  ID <- condition <- feature <- index.sub <- NULL
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  svg.path <- svg_pa(svg(input))[[1]] 
  svg.na <- basename(svg.path)
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
    file.copy(svg.path, file.path(out.dir, na))
    df.new0 <- cbind(colorNew=df0$fill, df.att0, SVG=na)
    # Color output SVG file.
    update_feature(df.new=df.new0, dir=out.dir) 
  }
  }
}



