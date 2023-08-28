#' Create Customized spatialHeatmap Shiny Apps
#'
#' This function creates customized spatialHeatmap Shiny Apps with user-provided data, aSVG files, and default parameters by using the spatialHeatmap Shiny App as the template. 

#' @param data A nested `list` of custom data and aSVG files that will be the pre-included examples in the custom App. File paths of each data-aSVG pair should be included in a `list` that have four slots: `name`, `display`, `data`, and `svg`, e.g. `lis1 <- list(name='mus.brain', display='Mouse brain (SHM)', data='./mus_brain.txt', svg='./mus_brain.svg')`. The `name` will be syntactically valide for R code, while the `display` will be shown on the user interface, so the latter can include special characters. The `data` and `svg` is the file path of numeric data and corresponding aSVG respectively. If multiple aSVGs (e.g. growth stages) correspond to a single numeric data set, the respective paths will be stored in a vector in the `svg` slot (see example below). After store the data-aSVG pairs in separate `list`s, store these `list`s in another `list` (nested `list`), e.g. `list(lis1, lis2)`. The data and aSVGs in the nested `list` will be copied to the `data` folder in the custom App.
#' @param db Paired assay data and aSVGs in a tar file that will be used as a backend database. See \code{\link{write_hdf5}}.  
#' @param lis.par A `list` of custom parameters of the Shiny App that will be the default when the custom App is launched. See \code{par.tmp}. Default is `NULL` and the default parameters in the spatialHeatmap Shiny App will be inherited. 
#' @param par.tmp Logical. If `TRUE` the default paramters in the spatialHeatmap Shiny App are returned in a `list`, and users can edit these settings then assign the `list` back to \code{lis.par}. Note, only the existing values in the `list` can be edited and the hierarchy of the `list` should be preserved. Otherwise, it cannot be recognized by the internal program. 
#' @param dld.sgl A `list` of paired data matrix and single aSVG file, which would be downloadable example dataset on the App for testing. The `list` consists of paths of the data matrix and aSVG file with name slots of `data` and `svg` respectively, e.g. \code{list(data='./data_download.txt', svg='./root_download_shm.svg')}. The specified data and aSVG will copied to the `data` folder in the App. 
#' @param dld.mul A `list` of paired data matrix and multiple aSVG files, which would be downloadable example dataset on the App for testing. It is the same as `dld.sgl` except that in the `svg` slot, multiple aSVG file paths are stored, e.g. `list(data='./data_download.txt', svg=c('./root_young_download_shm1.svg', './root_old_download_shm2.svg'))`. 
#' @param dld.vars A `list` of paired data matrix and single aSVG file, which would be downloadable data.tmp dataset on the App for testing. It is the same as `dld.sgl` except that multiple experimental variables are combined in the data matrix (see package viengettes for details.)
#' @param data.tmp Logical. If `TRUE` (default), both example data sets in the spatialHeatmap Shiny App and custom data sets in `data` will be included in the custom App.
#' @param ldg,gallery,about The paths of ".html" files that will be used for the landing, gallery, and about pages in the custom Shiny App respectively. The default is `NULL`, indicating default ".html" files will be used.  
#' @param app.dir The directory to create the Shiny App. Default is current work directory \code{.}.

#' @return If \code{par.tmp==TRUE}, the default paramters in spatialHeatmap Shiny App are returned in a `list`. Otherwise, a customized Shiny App is generated in the path of \code{app.dir}. 

#' @examples
#'
#' # The data sets in spatialHeatmap are used for demonstrations.
#' 
#' ## Below are demonstrations of simple usage.
#' # File paths of one data matrix and one aSVG.
#' data.path1 <- system.file('extdata/shinyApp/data/expr_mouse.txt', package='spatialHeatmap')
#' svg.path1 <- system.file('extdata/shinyApp/data/mus_musculus.male.svg', 
#' package='spatialHeatmap')
#' # Save the file paths in a list with name slots of "name", "display", "data", and "svg".
#' lis.dat1 <- list(name='Mouse', display='Mouse (SHM)', data=data.path1, svg=svg.path1)
#' # Create custom Shiny Apps.
#' \donttest{
#' if (!dir.exists('~/test_shiny')) dir.create('~/test_shiny')
#' # Create custom Shiny app by feeding this function these datasets and parameters.
#' custom_shiny(data=list(lis.dat1), app.dir='~/test_shiny')
#' # Lauch the app.
#' shiny::runApp('~/test_shiny/shinyApp') 
#' }
#'
#' ## Below are demonstrations of advanced usage. 
#'
#' # Paths of one data matrix and two aSVGs (two growth stages).
#' data.path2 <- system.file('extdata/shinyApp/data/random_data_multiple_aSVGs.txt', 
#' package='spatialHeatmap')
#' svg.path2.1 <- system.file('extdata/shinyApp/data/arabidopsis.thaliana_organ_shm1.svg', 
#' package='spatialHeatmap')
#' svg.path2.2 <- system.file('extdata/shinyApp/data/arabidopsis.thaliana_organ_shm2.svg', 
#' package='spatialHeatmap')
#' # Save the file paths in a list with name slots of "name", "display", "data", and "svg".
#' lis.dat2 <- list(name='growthStage', display='Multiple aSVGs (SHM)',data=data.path2, svg=c(svg.path2.1, svg.path2.2))
#'
#' # Paths of one data matrix with combined variables and one aSVG.
#' data.path.vars <- system.file('extdata/shinyApp/data/mus_brain_vars_se_shiny.rds', 
#' package='spatialHeatmap')
#' svg.path.vars <- system.file('extdata/shinyApp/data/mus_musculus.brain.svg', 
#' package='spatialHeatmap')
#' # Save the file paths in a list with name slots of "name", "display", "data", and "svg".
#' lis.dat.vars <- list(name='multiVariables', display='Multiple variables (SHM)', data=data.path.vars, svg=svg.path.vars)
#'
#' # Paths of one data matrix and one aSVGs for creating downloadable example data sets.
#' data.path.dld1 <- system.file('extdata/shinyApp/data/expr_mouse.txt', 
#' package='spatialHeatmap')
#' svg.path.dld1 <- system.file('extdata/shinyApp/data/mus_musculus.male.svg', 
#' package='spatialHeatmap')
#' # Save the file paths in a list with name slots of "name", "display", "data", and "svg".
#' dld.sgl <- list(data=data.path.dld1, svg=svg.path.dld1)
#'
#' # For demonstration purpose, the same data and aSVGs are used to make the list for creating 
#' # downloadable example dataset of two growth stages. 
#' dld.mul <- list(data=data.path2, svg=c(svg.path2.1, svg.path2.2))
#'
#' # For demonstration purpose, the same multi-variable data and aSVG are used to create the
#' # downloadable example multi-variable dataset.
#' dld.vars <- list(data=data.path.vars, svg=svg.path.vars)
#'
#' # Retrieve the default parameters in the App template.
#' lis.par <- custom_shiny(par.tmp=TRUE)
#' # Change default setting of color scheme in the color key.
#' lis.par$shm.img['color', ] <- 'yellow,orange,blue'
#' # The default dataset to show when the app is launched.
#' lis.par$default.dataset <- 'shoot'
#' 
#' # Store all default data sets in a nested list.
#' dat.all <- list(lis.dat1, lis.dat2, lis.dat.vars)
#'
#' # Create custom Shiny Apps.
#' \donttest{
#' if (!dir.exists('~/test_shiny')) dir.create('~/test_shiny')
#' custom_shiny(data=dat.all, lis.par=lis.par, dld.sgl=dld.sgl, dld.mul=dld.mul, dld.vars=dld.vars, 
#' app.dir='~/test_shiny')
#' # Lauch the App.
#' shiny::runApp('~/test_shiny/shinyApp') 
#' }
#'
#' # The customized Shiny App is able to take database backend as well. Examples are 
#' # demonstrated in the function "write_hdf5".
#'
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Jeremy Stephens, Kirill Simonov, Yihui Xie, Zhuoer Dong, Hadley Wickham, Jeffrey Horner, reikoch, Will Beasley, Brendan O'Connor and Gregory R. Warnes (2020). yaml: Methods to Convert R Data to YAML and Back. R package version 2.2.1. https://CRAN.R-project.org/package=yaml
#' Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017). shiny: Web Application Framework for R. R package version 1.0.3. https://CRAN.R-project.org/package=shiny 

#' @export
#' @importFrom grDevices colors

custom_shiny <- function(data=NULL, db=NULL, lis.par=NULL, par.tmp=FALSE, dld.sgl=NULL, dld.mul=NULL, dld.vars=NULL, data.tmp=TRUE, ldg=NULL, gallery=NULL, about=NULL, app.dir='.') {
  options(stringsAsFactors=FALSE)
  pkg <- check_pkg('yaml'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
  if (is.null(data) & is.null(db) & par.tmp==FALSE) stop("Both 'data' and 'db' are 'NULL'!")
  if (!dir.exists(app.dir)) dir.create(app.dir) 
  # Default config file.
  cfg.def <- yaml::yaml.load_file(system.file('extdata/shinyApp/config/config.yaml', package='spatialHeatmap'))
  # Default parameters.
  lis.par.def <- cfg.def[!grepl('^dataset\\d+|download_single|download_multiple|download_multiple_variables|download_covisualization|download_batched_data_aSVGs', names(cfg.def))]
  # Return parameter template.
  if (par.tmp==TRUE) {
    for (i in seq_along(lis.par.def)) {
      lis0 <- lis.par.def[[i]]; if (length(lis0)>1) { 
        name <- default <- NULL; for (j in seq_along(lis0)) {
          pair <- strsplit(lis0[j], ':')[[1]]
          name <- c(name, pair[1]); default <- c(default, pair[2])
        }; df0 <- data.frame(default=default)
        rownames(df0) <- name; lis.par.def[[i]] <- df0
      } 
    }; return(lis.par.def)
  }
  app.dir0 <- normalizePath(app.dir, winslash="/", mustWork=FALSE)
  app.dir <- file.path(app.dir0, 'shinyApp')
  cp_app(app.dir=app.dir0) # Copy Shiny App from spatialHeatmap to the target directory.
  if (!is.null(ldg)) {
    if (!grepl('\\.html$', ldg)) stop('The landing page should be an ".html" file!') 
    file.copy(ldg, file.path(app.dir, 'www/html/landing.html'), overwrite=TRUE)
  }
  if (!is.null(gallery)) { 
    if (!grepl('\\.html$', gallery)) stop('The gallery page should be an ".html" file!')
    file.copy(gallery, file.path(app.dir, 'www/html/gallery.html'), overwrite=TRUE)
  }
  if (!is.null(about)) { 
    if (!grepl('\\.html$', about)) stop('The About page should be an ".html" file!')
    file.copy(about, file.path(app.dir, 'www/html/about.html'), overwrite=TRUE)
  } 
  # Load default parameter list.
  if (is.null(lis.par)) { lis.par <- lis.par.def } else {  
    for (i in seq_along(lis.par)) {
      # Names and values are concatenated by ':', thus ':' cannot be used in values.
      df0 <- lis.par[[i]]; if (is.data.frame(df0)) {
      pair <- paste0(row.names(df0), ':', df0$default); lis.par[[i]] <- pair
      }
    }
  }

  # Include default examples or not.
  if (data.tmp==FALSE) { 
    # If exclude default examples and no download files are provided, default download files are retained.
    # pat.dld <- 'dummyfile|expr_mouse.txt|arabidopsis.thaliana_root.cross_shm.svg|random_data_multiple_aSVGs.txt|arabidopsis.thaliana_organ_shm1.svg|arabidopsis.thaliana_organ_shm2.svg'
    # file.remove(grep(pat.dld, list.files(paste0(app.dir, '/data'), '*', full.names=TRUE), invert=TRUE, value=TRUE))
    # if (!is.null(dld.sgl)) file.remove(list.files(paste0(app.dir, '/data'), 'expr_mouse.txt|arabidopsis.thaliana_root.cross_shm.svg', full.names=TRUE))
    # if (!is.null(dld.mul)) file.remove(list.files(paste0(app.dir, '/data'), 'random_data_multiple_aSVGs.txt|arabidopsis.thaliana_organ_shm1.svg|arabidopsis.thaliana_organ_shm2.svg', full.names=TRUE))
    exp <- NULL
  } else {
    # Include all default examples.
    cfg.def <- yaml::yaml.load_file(system.file('extdata/shinyApp/config/config.yaml', package='spatialHeatmap'))
    lis.dat.def <- cfg.def[grepl('^dataset\\d+', names(cfg.def))]
    idx.rm <- NULL; for (i in seq_along(lis.dat.def)) {
    if (any(lis.dat.def[[i]]$svg=='none')) idx.rm <- c(idx.rm, i)
    }; exp <- lis.dat.def[-idx.rm]
  }
  # Validate colours.
  col_check('shm.img', lis.par$shm.img)
  col_check('network', lis.par$network)
  dat.db <- NULL; app.dat.pa <- file.path(app.dir, 'data')
  if (!is.null(db)) if (grepl('\\.tar$', db)) {
    # Copy data/svg in tar file to the App.
    file.copy(db, app.dat.pa, overwrite=TRUE, recursive=TRUE)
    # Check overlap between data lists from "data" and "db".
    check.dat <- ovl_dat_db(data=data, db=db)
    data <- check.dat$data; dat.db <- check.dat$dat.db 
  }
  # Use default download files.
  if (!is.null(dld.sgl)) lis.dld1 <- cp_file(dld.sgl, app.dir, 'data') else lis.dld1 <- cfg.def$download_single
  if (!is.null(dld.mul)) lis.dld2 <- cp_file(dld.mul, app.dir, 'data') else lis.dld2 <- cfg.def$download_multiple 
  if (!is.null(dld.vars)) lis.dld3 <- cp_file(dld.vars, app.dir, 'data') else lis.dld3 <- cfg.def$download_multiple_variables 
  lis.dld <- list(download_single=lis.dld1, download_multiple=lis.dld2, download_multiple_variables=lis.dld3)
  # Custom options are always included.
  # lis.cus1 <- lis.cus2 <- NULL
  # if (custom==TRUE) 
  lis.cus1 <- list(name='customBulkData', display='customBulkData', data='none', svg='none')
  # if (custom.computed==TRUE) 
  lis.cus2 <- list(name='customCovisData', display='customCovisData', data='none', svg='none')
  # All data sets.
  if (is(data, 'list')) for (i in seq_along(data)) {
    dat.lis0 <- data[[i]]
    dat.pa0 <- dat.lis0$data; svg.pa0 <- dat.lis0$svg
    # Change data path based on the App.
    data[[i]]$data <- file.path('data', basename(dat.pa0))
    data[[i]]$svg <- file.path('data', basename(svg.pa0))
    # Copy data/aSVG to the App. 
    file.copy(dat.pa0, app.dat.pa); file.copy(svg.pa0, app.dat.pa) 
  }
  data <- c(list(lis.cus1, lis.cus2), data, dat.db, exp)
  data <- data[!vapply(data, is.null, logical(1))]
  # Name the complete list.
  names(data) <- paste0('dataset', seq_along(data))
  lis.all <- c(data, lis.par, lis.dld)
  yaml::write_yaml(lis.all, paste0(app.dir, '/config/config.yaml')); cat('Done! \n')
}


#' Check overlap between separate data sets and database
#'
#' @keywords Internal
#' @noRd

ovl_dat_db <- function(data, db) { 
  dat.db <- read_hdf5(db, name='match')
  for (i in seq_along(dat.db)) dat.db[[i]]$data <- dat.db[[i]]$svg <- 'data/data_shm.tar'
  if (is.null(data)) return(list(data=data, dat.db=dat.db))
  nas <- lapply(data, function(x) x$name) 
  nas.db <- lapply(dat.db, function(x) x$name)
  # Data in tar file take precedence over in list.
  data <- data[!nas %in% nas.db]
  return(list(data=data, dat.db=dat.db))
}  

#' Copy Shiny App from spatialHeatmap to the target directory
#'
#' @keywords Internal
#' @noRd

cp_app <- function(app.dir) { 
  app.dir <- normalizePath(app.dir, winslash="/", mustWork=FALSE)
  app.dir <- file.path(app.dir, 'shinyApp')
  if (!dir.exists(app.dir)) dir.create(app.dir)
  app.path <- system.file('extdata/shinyApp', package='spatialHeatmap')
  # Remove residues from spatialHeatmap.
  for (i in c('app.R', 'config', 'data', 'R', 'www')) {
    file.copy(file.path(app.path, i), app.dir, recursive=TRUE, overwrite=TRUE) 
  }
  file.remove(list.files(file.path(app.dir, 'www/video'), '*.mp4$', full.names=TRUE))
  file.remove(list.files(file.path(app.dir, 'www/html_shm'), '*.html$', full.names=TRUE))
  file.remove(list.files(file.path(app.dir, 'www/html_shm/lib'), '*', full.names=TRUE))
}
  

#' Check validity of color indredients in the yaml file
#'
#' @keywords Internal
#' @noRd

col_check <- function(element, vec.all) {

  col0 <- vec.all[grepl('^color:', vec.all)]
  color <- gsub('.*:(.*)', '\\1', col0)
  color <- gsub(' |\\.|-|;|,|/', '_', color)
  color <- strsplit(color, '_')[[1]]
  color <- color[color!='']; color1 <- color[!color %in% colors()]
  if (length(color1)>0) stop(paste0('Colors in ', element, ' not valid: ', paste0(color1, collapse=', '), '!'))

}


#' Convert the data-aSVG pair from data frame to list
#'
#' @keywords Internal
#' @noRd

pair2lis <- function(df.pair, db=FALSE) { 

  na.all <- as.vector(df.pair$name); dat.all <- as.vector(df.pair$data)
  svg.all <- df.pair[, 'aSVG']
  lis.dat0 <- NULL; for (i in seq_len(nrow(df.pair))) {
    # Separate multiple aSVGs.
    svg <- svg.all[i]; if (grepl(';| |,', svg)) {

      strs <- strsplit(svg, ';| |,')[[1]]; svg <- strs[strs!='']

    }
  # Data of txt file and db are distinguished by 'data/' in the dataset list.
  if (db==FALSE) lis.dat0 <- c(lis.dat0, list(list(name=na.all[i], data=paste0('data/', dat.all[i]), svg=paste0('data/', svg)))) else lis.dat0 <- c(lis.dat0, list(list(name=na.all[i], data=dat.all[i], svg=svg)))

  }; return(lis.dat0)

}


#' Copy user-provided files, and change data/svg path.
#'
#' @keywords Internal
#' @noRd

cp_file <- function(lis, app.dir, folder) {
  if (is.null(lis)) return()
  for (i in seq_along(lis)) { 
    lis0 <- lis[[i]]; for (k in seq_along(lis0)) {
      # Copy files.
      vec <- lis0[[k]]; if (!all(file.exists(vec))) next
      files <- NULL; for (v in vec) files <- c(files, v)
      # cat('Copying files: \n'); print(files)
      file.copy(files, paste0(app.dir, '/', folder), overwrite=TRUE)
      # Shorten paths.
      if (length(vec)==1) { 
        str <- strsplit(vec, '/')[[1]]
        lis[[i]][[k]] <- paste0(folder, '/', str[length(str)])
      } else if (length(vec)>1) {
        svgs <- NULL; for (j in seq_along(vec)) {
          str <- strsplit(vec[j], '/')[[1]]
          svgs <- c(svgs, paste0(folder, '/', str[length(str)]))
        }; lis[[i]][[k]] <- svgs
      }
    }
  }; return(lis)
}
