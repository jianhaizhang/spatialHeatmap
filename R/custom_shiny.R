#' Create Customized Shiny App of Spaital Heatmap
#'
#' This function creates customized Shiny App with user-provided data, aSVG files, and default parameters. Default settings are defined in the "config.yaml" file in the "config" folder of the app, and can be edited directly in a yaml file editor.  

#' @param ... Separate lists of paired data matrix and aSVG files, which are included as default datasets in the Shiny app. Each list must have three elements with name slots of "name", "data", and "svg" respectively. For example, list(name='dataset1', data='./data1.txt', svg='./root_shm.svg'). The "name" element (\emph{e.g.} 'dataset1') is listed under "Step 1: data sets" in the app, while "data" and "svg" are the paths of data matrix and aSVG files. If multiple aSVGs (\emph{e.g.} growth stages) are included in one list, the respective paths are stored in a vector in the "svg" slot (see example below). After calling this function, the data and aSVGs are copied to the "example" folder in the app. See detailed examples below.

#' @param lis.par A list of default parameters of the Shiny app. See \code{ lis.par.tmp }. Default is NULL, which means default parameters are adopted.

#' @param lis.par.tmp Logical, TRUE or FALSE. Default is FALSE. If TRUE the template of default paramter list is returned, and users can set customized default values then assign this list to \code{ lis.par }. Note, only the default values in the list can be changed while the hierarchy of the list should be preserved. Otherwise, it cannot be recognized by the internal program. 
#' @param lis.dld.single A list of paired data matrix and single aSVG file, which would be downloadable on the app for testing. The list should have two elements with name slots of "data" and "svg" respectively, which are the paths of the data matrix and aSVG file repectively. After the function call, the specified data and aSVG are copied to the "example" folder in the app. Note the two name slots should not be changed. \emph{E.g.} \code{list(data='./data_download.txt', svg='./root_download_shm.svg')}.
#' @param lis.dld.mul A list of paired data matrix and multiple aSVG files, which would be downloadable on the app for testing. The multiple aSVG files could be multiple growth stages of a plant. The list should have two elements with name slots of "data" and "svg" respectively, which are the paths of the data matrix and aSVG files repectively. After the function call, the specified data and aSVGs are copied to the "example" folder in the app. Note the two name slots should not be changed. \emph{E.g.} list(data='./data_download.txt', svg=c('./root_young_download_shm.svg', './root_old_download_shm.svg')).
#' @param custom Logical, TRUE or FALSE. If TRUE (default), the "customData" option under "Step 1: data sets" is included, which allows to upload datasets from local computer.
#' @param custom.computed Logical, TRUE or FALSE. If TRUE (default), the "customComputdData" option under "Step 1: data sets" is included, which allows to upload computed datasets from local computer. See \code{\link{adj_mod}}. 
#' @param example Logical, TRUE or FALSE. If TRUE (default), the default examples in "spatialHeatmap" package are included in the app as well as those provided to \code{...} by users.
#' @param app.dir The directory to create the Shiny app. Default is current work directory \code{'.'}.

#' @return If \code{lis.par.tmp==TRUE}, the template of default paramter list is returned. Otherwise, a customized Shiny app is generated in the path of \code{app.dir}. 

#' @examples

#' # The pre-packaged examples are used for illustration purpose.
#' # Get one data path and one aSVG path and assembly them into a list for creating default dataset.

#' data.path1 <- system.file('extdata/shinyApp/example/expr_arab.txt', package='spatialHeatmap')
#' svg.path1 <- system.file('extdata/shinyApp/example/arabidopsis_thaliana.shoot_shm.svg', 
#' package='spatialHeatmap')
#' # The list with name slots of "name", "data", and "svg".
#' lis.dat1 <- list(name='shoot', data=data.path1, svg=svg.path1)

#' # Get one data path and two aSVG paths and assembly them into a list for creating default 
#' # dataset, which include two growth stages.
#' data.path2 <- system.file('extdata/shinyApp/example/random_data_multiple_aSVGs.txt', 
#' package='spatialHeatmap')
#' svg.path2.1 <- system.file('extdata/shinyApp/example/arabidopsis_thaliana.organ_shm1.svg', 
#' package='spatialHeatmap')
#' svg.path2.2 <- system.file('extdata/shinyApp/example/arabidopsis_thaliana.organ_shm2.svg', 
#' package='spatialHeatmap')
#' # The list with name slots of "name", "data", and "svg", where the two aSVG paths are stored
#' # in a vector in "svg".
#' lis.dat2 <- list(name='growthStage', data=data.path2, svg=c(svg.path2.1, svg.path2.2))

#' # Get one data path and one aSVG path and assembly them into a list for creating downloadable
#' # dataset.
#' data.path.dld1 <- system.file('extdata/shinyApp/example/expr_arab.txt', 
#' package='spatialHeatmap')
#' svg.path.dld1 <- system.file('extdata/shinyApp/example/arabidopsis_thaliana.organ_shm.svg', 
#' package='spatialHeatmap')
#' # The list with name slots of "data", and "svg".
#' lis.dld.single <- list(name='organ', data=data.path.dld1, svg=svg.path.dld1)

#' # For demonstration purpose, the same data and aSVGs are used to make the list for creating 
#' # downloadable dataset, which include two growth stages. 
#' lis.dld.mul <- list(data=data.path2, svg=c(svg.path2.1, svg.path2.2))

#' # Retrieve the default parameters.
#' lis.par <- custom_shiny(lis.par.tmp=TRUE)
#' # Change default values.
#' lis.par$shm.img['color', ] <- 'yellow,orange,blue'
#' # The default dataset to show upon the app is launched.
#' lis.par$default.dataset <- 'shoot'

#' \donttest{
#' if (!dir.exists('~/test_shiny')) dir.create('~/test_shiny')
#' # Create custom Shiny app by feeding this function these datasets and parameters.
#' custom_shiny(lis.dat1, lis.dat2, lis.par=lis.par, lis.dld.single=lis.dld.single, 
#' lis.dld.mul=lis.dld.mul, app.dir='~/test_shiny')
#' # Lauch the app.
#' shiny::runApp('~/test_shiny/shinyApp') 
#' }
#'
#' # The customized Shiny app is able to take database backend as well. Examples are demonstrated in the function "write_hdf5".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Jeremy Stephens, Kirill Simonov, Yihui Xie, Zhuoer Dong, Hadley Wickham, Jeffrey Horner, reikoch, Will Beasley, Brendan O'Connor and Gregory R. Warnes (2020). yaml: Methods to Convert R Data to YAML and Back. R package version 2.2.1. https://CRAN.R-project.org/package=yaml
#' \cr Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017). shiny: Web Application Framework for R. R package version 1.0.3. https://CRAN.R-project.org/package=shiny 

#' @export custom_shiny
#' @importFrom yaml yaml.load_file write_yaml
#' @importFrom grDevices colors

custom_shiny <- function(..., lis.par=NULL, lis.par.tmp=FALSE, lis.dld.single=NULL, lis.dld.mul=NULL, custom=TRUE, custom.computed=TRUE, example=TRUE, app.dir='.') {

  options(stringsAsFactors=FALSE)
  # Default config file.
  cfg.def <- yaml.load_file(system.file('extdata/shinyApp/config/config.yaml', package='spatialHeatmap'))
  # Default parameters.
  lis.par.def <- cfg.def[!grepl('^dataset\\d+|download_single|download_multiple', names(cfg.def))]
  # Return parameter template.
  if (lis.par.tmp==TRUE) {

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
  app.dir <- paste0(app.dir0, '/shinyApp')
  app.path <- system.file('extdata/shinyApp', package='spatialHeatmap')
  # Remove residues from spatialHeatmap.
  file.copy(app.path, app.dir0, recursive=TRUE, overwrite=TRUE) 
  file.remove(list.files(paste0(app.dir, '/www/video'), '*.mp4$', full.names=TRUE))
  file.remove(list.files(paste0(app.dir, '/www/ggly'), '*.html$', full.names=TRUE))
  file.remove(list.files(paste0(app.dir, '/www/ggly/lib'), '*', full.names=TRUE))
  unlink(paste0(app.dir, '/rsconnect'), recursive=TRUE) 
  file.remove(list.files(paste0(app.dir, '/html_shm'), '*.html$', full.names=TRUE))
  file.remove(list.files(paste0(app.dir, '/html_shm/lib/'), '*', full.names=TRUE))
  lis.dat <- list(...)
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
  if (example==FALSE) { 

    # If exclude default examples and no download files are provided, default download files are retained.
    pat.dld <- 'dummyfile|expr_arab.txt|arabidopsis_thaliana.root.cross_shm.svg|random_data_multiple_aSVGs.txt|arabidopsis_thaliana.organ_shm1.svg|arabidopsis_thaliana.organ_shm2.svg'
    file.remove(grep(pat.dld, list.files(paste0(app.dir, '/example'), '*', full.names=TRUE), invert=TRUE, value=TRUE))
    if (!is.null(lis.dld.single)) file.remove(list.files(paste0(app.dir, '/example'), 'expr_arab.txt|arabidopsis_thaliana.root.cross_shm.svg', full.names=TRUE))
    if (!is.null(lis.dld.mul)) file.remove(list.files(paste0(app.dir, '/example'), 'random_data_multiple_aSVGs.txt|arabidopsis_thaliana.organ_shm1.svg|arabidopsis_thaliana.organ_shm2.svg', full.names=TRUE))
    exp <- NULL

  } else {

    # Include all default examples.
    cfg.def <- yaml.load_file(system.file('extdata/shinyApp/config/config.yaml', package='spatialHeatmap'))
    lis.dat.def <- cfg.def[grepl('^dataset\\d+', names(cfg.def))]
    idx.rm <- NULL; for (i in seq_along(lis.dat.def)) {
    if (any(lis.dat.def[[i]]$svg=='none')) idx.rm <- c(idx.rm, i)
    }; exp <- lis.dat.def[-idx.rm]

  }
  # Validate colours.
  col_check('shm.img', lis.par$shm.img)
  col_check('network', lis.par$network)
  # Copy data/svg in list of two tar files and isolate data/svg list of tar file. 
  idx <- lis.dat1 <- NULL; for (i in seq_along(lis.dat)) {

    lis0 <- lis.dat[[i]]; if (length(lis0)==2) { 
  
      idx <- c(idx, i); for (j in lis0) {

        j <- normalizePath(j)
        if (grepl('data_shm.tar$', j)) lis.dat1 <- pair2lis(read_hdf5(j, 'df_pair')[[1]], db=TRUE)
        if (grepl('\\.tar$', j)) system(paste0('cp ', j, ' ', app.dir, '/example')) else stop('Compressed data and aSVGs should be two independent ".tar" files respectively!')

      }

    }

  }
  if (is.null(idx)) lis.dat <- cp_file(lis.dat, app.dir, 'example') else { 
    lis.dat <- cp_file(lis.dat[-idx], app.dir, 'example')
    # Data in tar file take precedence over in list.
    na.all <- vapply(lis.dat, function(x) { na <- NULL; na <- c(na, x$name) }, character(1)) 
    na.all1 <- vapply(lis.dat1, function(x) { na <- NULL; na <- c(na, x$name) }, character(1))
    lis.dat <- lis.dat[!(na.all %in% na.all1)]
  }
  if (!is.null(lis.dld.single)) lis.dld1 <- cp_file(lis.dld.single, app.dir, 'example') else {
    # Use default download files.
    lis.dld1 <- list(data="example/expr_arab.txt", svg="example/arabidopsis_thaliana.root.cross_shm.svg")
  }
  if (!is.null(lis.dld.mul)) lis.dld2 <- cp_file(lis.dld.mul, app.dir, 'example') else {
    # Use default download files. 
    lis.dld2 <- list(data="example/random_data_multiple_aSVGs.txt", svg=c('example/arabidopsis_thaliana.organ_shm1.svg', 'example/arabidopsis_thaliana.organ_shm2.svg'))

  }; lis.dld <- list(download_single=lis.dld1, download_multiple=lis.dld2)
  lis.cus1 <- lis.cus2 <- NULL
  if (custom==TRUE) lis.cus1 <- list(name='customData', data='none', svg='none')
  if (custom.computed==TRUE) lis.cus2 <- list(name='customComputedData', data='none', svg='none')
  # All data sets.
  lis.dat <- c(list(list(name='none', data='none', svg='none'), lis.cus1, lis.cus2), lis.dat, lis.dat1, exp)
  # Name the complete list.
  names(lis.dat) <- paste0('dataset', seq_along(lis.dat))
  lis.all <- c(lis.dat, lis.par, lis.dld)
  write_yaml(lis.all, paste0(app.dir, '/config/config.yaml')); cat('Done! \n')

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
  # Data of txt file and db are distinguished by 'example/' in the dataset list.
  if (db==FALSE) lis.dat0 <- c(lis.dat0, list(list(name=na.all[i], data=paste0('example/', dat.all[i]), svg=paste0('example/', svg)))) else lis.dat0 <- c(lis.dat0, list(list(name=na.all[i], data=dat.all[i], svg=svg)))

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
