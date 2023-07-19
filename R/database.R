#' Creat databases for the Shiny App
#'
#' The function \link{write_hdf5} is designed to construct the database backend for the Shiny App (\link{shiny_shm}), while \link{read_hdf5} is designed to query the database. 

#' @name Database
#' @rdname Database
#' @aliases write_hdf5 read_hdf5

#' @return 
#' `write_hdf5`: a file 'data_shm.tar' saved on disk. \cr
#' `read_hdf5`: a nested `list` of file paths corresponding to data-aSVG pairs. 

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
#' R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Hervé Pagès (2020). HDF5Array: HDF5 backend for DelayedArray objects. R package version 1.16.1.
#' Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8 
#' Davis, Sean, and Paul Meltzer. 2007. “GEOquery: A Bridge Between the Gene Expression Omnibus (GEO) and BioConductor.” Bioinformatics 14: 1846–7
#' Gautier, Laurent, Leslie Cope, Benjamin M. Bolstad, and Rafael A. Irizarry. 2004. “Affy—analysis of Affymetrix GeneChip Data at the Probe Level.” Bioinformatics 20 (3). Oxford, UK: Oxford University Press: 307–15. doi:10.1093/bioinformatics/btg405
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Huber, W., V. J. Carey, R. Gentleman, S. An ders, M. Carlson, B. S. Carvalho, H. C. Bravo, et al. 2015. “Orchestrating High-Throughput Genomic Analysis Wit H Bioconductor.” Nature Methods 12 (2): 115–21. http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9 

#' @examples
#'
#' ## Example of a single aSVG instance.
#'
#' # The example data included in this package come from an RNA-seq analysis on 
#' # development of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019). 
#' # The complete raw count data are downloaded using the R package ExpressionAtlas
#' # (Keays 2019) with the accession number "E-MTAB-6769". 
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
#' # Indicate spatial features and experiment variables with "spFeature" and "variable"
#' # in the "colData" slot respectively.
#' colnames(colData(se.chk))[8] <- 'spFeature'
#' colnames(colData(se.chk))[6] <- 'variable'
#' colData(se.chk)[1:2, ]
#'
#' # Create a temporary directory "db_shm".
#' dir.db <- file.path(tempdir(check=TRUE), 'db_shm')
#' if (!dir.exists(dir.db)) dir.create(dir.db)
#' 
#' # Save the assay data in the temporary directory.
#' chk.dat.pa <- file.path(dir.db, 'test_chicken.rds')
#' saveRDS(se.chk, file=chk.dat.pa)
#'
#' # The chicken aSVG downloaded from the EBI aSVG repository (https://github.com/ebi-gene-
#' # expression-group/anatomogram/tree/master/src/svg) is included in this package and 
#' # accessed as below.
#' svg.chk <- system.file("extdata/shinyApp/data", "gallus_gallus.svg",
#' package="spatialHeatmap")
#' # Store file paths of data and aSVG pair in a list. 
#' dat1 <- list(name='test_chicken', display='Test (chicken SHM)', data=chk.dat.pa, svg=svg.chk)
#'
#' ## Example of multiple aSVG and raster images.
#'
#' # Two aSVG instances of maize leaves represent two growth stages and each aSVG has a raster 
#' # image counterpart.
#' # Random numeric data.
#' pa.leaf <- system.file("extdata/shinyApp/data", 'dat_overlay.txt',
#' package="spatialHeatmap")
#' dat.leaf <- read_fr(pa.leaf); dat.leaf[1:2, ]
#' # Save the assay data in the temporary directory.
#' leaf.dat.pa <- file.path(dir.db, 'test_leaf.rds')
#' saveRDS(dat.leaf, file=leaf.dat.pa)
#' 
#' # Paths of the two aSVG files.
#' svg.leaf1 <- system.file("extdata/shinyApp/data", 'maize_leaf_shm1.svg',
#' package="spatialHeatmap")
#' svg.leaf2 <- system.file("extdata/shinyApp/data", 'maize_leaf_shm2.svg',
#' package="spatialHeatmap") 
#' # Paths of the two corresponsing raster images.
#' rst.leaf1 <- system.file("extdata/shinyApp/data", 'maize_leaf_shm1.png',
#' package="spatialHeatmap")
#' rst.leaf2 <- system.file("extdata/shinyApp/data", 'maize_leaf_shm2.png',
#' package="spatialHeatmap") 
#'
#' # Store file paths of data, aSVG, raster images in a list. 
#' dat2 <- list(name='test_leaf', display='Test (maize leaf SHM)', data=leaf.dat.pa, svg=c(svg.leaf1, svg.leaf2, rst.leaf1, rst.leaf2))
#' 
#' # Store these two data sets in a nested list. 
#' dat.lis <- list(dat1, dat2) 
#'
#' \donttest{
#' # Save the database in the temporary directory. 
#' write_hdf5(data=dat.lis, dir=dir.db, replace=TRUE)
#' # Read data-aSVG pair of "test_leaf" from "data_shm.tar".
#' dat.q <- read_hdf5(file=file.path(dir.db, 'data_shm.tar'), name='test_leaf')
#' # Numeric data.
#' dat.leaf <- readRDS(dat.q[[1]]$data)
#' # aSVG.
#' svg1 <- dat.q[[1]]$svg[1]; svg2 <- dat.q[[1]]$svg[2]
#' rst1 <- dat.q[[1]]$svg[3]; rst2 <- dat.q[[1]]$svg[4]
#' svg.leaf <- read_svg(svg=c(svg1, svg2), raster=c(rst1, rst2))
#' }
#'
#' \donttest{  
#' # Create customized Shiny Apps with the database.
#' custom_shiny(db=file.path(dir.db, 'data_shm.tar'), app.dir='~/test_shiny')
#' # Create customized Shiny Apps with the data sets in a nested list.
#' custom_shiny(data=dat.lis, app.dir='~/test_shiny')
#' # Run the app.
#' shiny::runApp('~/test_shiny/shinyApp') 
#' }
NULL

#' @rdname Database
#' @param data A nested `list` of each numeric data and aSVG pair. Each pair consists of four slots: `name`, `display`, `data`, and `svg`, e.g. `lis1 <- list(name='mus.brain', display='Mouse brain (SHM)', data='./mus_brain.txt', svg='./mus_brain.svg')`. The `name` contains a syntactically valid entry for R while the `display` contains an entry for display on the user interface. The `data` and `svg` contain paths of numeric data and aSVGs that will be included in the data base respectively. The supported data containers include `data.frame` and `SummarizedExperiment`. See the `data` argument in \link{filter_data} for data formats. After formatted, the data should be saved as ".rds" files with the function \link{saveRDS}, then assign the corresponding paths the `data` slot. By contrast, the `svg` slot contains one or multiple aSVG/raster image file paths. Each numeric data set is first saved in an independent DHF5-based `SummarizedExperiment` (SE) object, then all these saved SE objects and aSVG/raster images are compressed together into a file "data_shm.tar", i.e. the backend database. In addition, the nested `list` assigned to `data` is also included in the database for querying purpose (\link{read_hdf5}). 
#' @param dir The directory to save the database "data_shm.tar".
#' @param replace If data set with the same `name` in \code{data} are already saved in \code{dir}, should the \code{dir} be emptied? Default is `FALSE`. If `TRUE`, the existing content in \code{dir} will be lost.

#' @inheritParams HDF5Array::saveHDF5SummarizedExperiment

#' @export
#' @importFrom SummarizedExperiment assay 
#' @importFrom utils untar tar

write_hdf5 <- function(data, dir='./data_shm', replace=FALSE, chunkdim=NULL, level=NULL, verbose=FALSE) {
  pkg <- check_pkg('HDF5Array'); if (is(pkg, 'character')) stop(pkg)
  options(stringsAsFactors=FALSE)
  dir <- normalizePath(dir, winslash="/", mustWork=FALSE)
  # Check data list format.
  dat.vld <- check_dat_db(data)
  if (!is.null(dat.vld)) return(dat.vld)
  cat('Data in progress... \n')
  dir.tmp <- file.path(tempdir(check=TRUE), 'shm_db') 
  hd.nas <- NULL
  msg.var <- 'Experimental variables are not detected: '
  for (i in seq_along(data)) {
    # Check "sample_condition" format, and stored data in SE.
    dat.pa <- data[[i]]$data
    if (!file.exists(dat.pa)) {
      msg <- 'Data do not exist!'; message(msg); return(msg) 
    }; dat <- readRDS(dat.pa)
    if (is(dat, 'data.frame')|is(dat, 'matrix')|is(dat, 'DFrame')|is(dat, 'dgCMatrix')) {# Convert data to SE.
      form <- grepl('__', colnames(dat))
      if (sum(form)==0) message(msg.var, 'data set ', i)
      dat <- check_data(data=dat)$se
    }
    # Check "spFeature" and "variable" in colData.
    if (is(dat, 'SummarizedExperiment')) {
      dat <- check_se(dat); if (is(dat, 'character')) { message(data[[i]]$name); return(dat) }
      dat <- dat$se
      hd.na <- paste0(data[[i]]$name, '_') 
      hd.nas <- c(hd.na, hd.nas)
      assay(dat) <- as(assay(dat), 'matrix')
      HDF5Array::saveHDF5SummarizedExperiment(dat, dir=dir.tmp, prefix=hd.na, replace=replace, chunkdim=chunkdim, level=level, verbose=verbose)
    } else cat(paste0('Accepted data containers are: data.frame, matrix, DFrame, dgCMatrix, SummarizedExperiment! This data is not saved: ', i), '\n')
  }
  file.tar <- file.path(dir.tmp, 'data_shm.tar')
  # dat.pat <- paste(paste0('^', hd.nas, 'assays.h5$', collapse='|'), paste0('^', hd.nas, 'se.rds$', collapse='|'), sep='|')
  # hd.all <- list.files(dir, dat.pat, full.names=FALSE)
  hd.all <- c(paste0(hd.nas, 'assays.h5'), paste0(hd.nas, 'se.rds'))
  svg.all <- unlist(lapply(data, function(x) x$svg))
  # Copy all SVGs to the same directory with data.
  file.copy(svg.all, dir.tmp, overwrite=TRUE)
  svg.na <- basename(svg.all)
  mat.na <- 'match.rds'
  saveRDS(data, file=file.path(dir.tmp, mat.na))  
  file.all <- c(hd.all, svg.na, mat.na)
  # This step is necessary to avoid recursive paths in tar files.
  wd <- getwd(); setwd(dir.tmp)
  # Data, svg, and matching list in the same tar file.
  tar(file.tar, files=file.all, tar="tar", extra_flags='--absolute-names')
  file.rename(file.tar, file.path(dir, 'data_shm.tar'))
  setwd(wd); unlink(dir.tmp, recursive=TRUE); message('Done!')
}

#' Check the format of SummarizedExperiment for use in Shiny App
#'
#' @return A scaled data frame.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan07@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

check_se <- function(se) {
  metadata <- NULL
  cdat <- colData(se); cna <- colnames(cdat)
  if ('spFeature' %in% cna & 'variable' %in% cna) { 
    colnames(se) <- paste0(se$spFeature, '__', se$variable); con.na <- TRUE
  } else if ('spFeature' %in% cna & !'variable' %in% cna) {
    se$variable <- 'con'
    colnames(se) <- paste0(se$spFeature, '__', se$variable); con.na <- FALSE
  } else if (!'spFeature' %in% cna) {
    msg <- 'Missing columns in the "colData" slot: 1. "spFeature"-spatial features such as tissues, organs, etc; 2. "variable"-experimental variables such as treatments, time points, etc.'
    message(msg); return(msg)
  }; df.meta <- metadata(se)$df.meta 
  if (!is.null(df.meta)) {
    if (ncol(data.frame(df.meta))<2) { msg <- 'The "df.meta" in the "metadata" slot should be a "data.frame" with at least two columns!'; message(msg); return(msg) }
  }; return(list(se=se, con.na=con.na))
}


#' Check the format of data list
#'
#' @return A scaled data frame.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan07@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

check_dat_db <- function(data) {
  for (i in data) { # List formats.
    if (!all(c('name', 'display', 'data', 'svg') %in% names(i))) {
      msg <- "Each data set need to have four slots 'name', 'display', 'data', and 'svg'!"
      warning(msg); print(i); return(msg)
    }
  }
  nas <- unlist(lapply(data, function(x) x$name))
  dup.na <- duplicated(nas)
  if (any(dup.na)) { 
    msg <- paste0("The 'names' are duplicated: ", paste0(nas[dup.na], collapse=', '), '!')
    warning(msg); return(msg)
  } 
}



#' @rdname Database
#' @param file The path of "data_shm.tar" generated by \link{write_hdf5}.
#' @param name One or multiple data set names (see the `data` argument in \link{write_hdf5}) in a vector, such as \code{c('test_leaf', 'test_chicken')}. If `match`, the matching `list` between assay data and aSVGs will be returned.   

#' @export
#' @importFrom utils untar

read_hdf5 <- function(file, name=NULL) {
  pkg <- check_pkg('HDF5Array'); if (is(pkg, 'character')) stop(pkg)
  options(stringsAsFactors=FALSE)
  dir <- paste0(tempdir(check=TRUE), '/data_shm')
  if (!dir.exists(dir)) dir.create(dir)
  # Extract matching list.
  untar(file, 'match.rds', exdir=dir, tar='tar')
  mat.lis <- readRDS(file.path(dir, 'match.rds'))
  for (i in seq_along(mat.lis)) {
    mat.lis[[i]]$data <- basename(mat.lis[[i]]$data)
    mat.lis[[i]]$svg <- basename(mat.lis[[i]]$svg)
  }; if ('match' %in% name) return(mat.lis)
  lis <- list(); for (i in name) {
    # A single data set (data and svg).
    lis0 <- mat.lis[lapply(mat.lis, function(x) x$name) %in% i]
    if (length(lis0)==0) {
      msg <- paste0('Data set not found: ', i); warning(msg)
      return(msg)
    }; lis0 <- lis0[[1]]
    hds <- paste0(i, c('_assays.h5', '_se.rds'))
    svgs <- basename(lis0$svg)
    # Extract a pair of data and svg.
    untar(file, c(hds, svgs), exdir=dir, tar='tar')
    # Load data in DHF5.
    dat0 <- HDF5Array::loadHDF5SummarizedExperiment(dir=dir, prefix=paste0(i, '_'))
    dat.pa <- file.path(dir, paste0(i, '.rds'))
    # Save data as .rds file.
    saveRDS(dat0, file=dat.pa)
    lis0$data <- dat.pa; lis0$svg <- file.path(dir, svgs)
    lis <- c(lis, list(lis0))
  }; return(lis)
}
