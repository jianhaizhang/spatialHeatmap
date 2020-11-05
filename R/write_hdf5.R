#' Construct Database for the Shiny App
#'
#' This is a convenience function for constructing the database backend in the Shiny app (\link{shiny_all}). The data to store in the database should be in the class of "data.frame" or "SummarizedExperiment" and should be formatted according to the conventions in the "data" argument of \link{spatial_hm}. After formatted, all these data should be arranged in a list and each data slot should have a unique name such as "expr_arab", "expr_chicken", \emph{etc.}. \cr In addition, a pairing data frame describing the matching relationship between the data and aSVG files must also be included in the list with the exclusive slot name "df_pair". This data frame should contain at least three columns: name, data, aSVG. The name column includes concise description of each data-aSVG pair, and entries in this column will be listed under "Step 1: data sets" on the Shiny app. The data column contains slot names of all data in the list ("expr_arab", "expr_chicken", \emph{etc.}), and the aSVG column includes the aSVG file names corresponding to each data respectively such as "gallus_gallus.svg", \emph{etc.} If one data is related to multiple aSVG files (\emph{e.g.} multiple development  stages), these aSVGs should be concatenated by comma, space, or semicolon, \emph{e.g.} "arabidopsis_thaliana.organ_shm1.svg;arabidopsis_thaliana.organ_shm2.svg". Inclusion of other columns providing metadata of the data and aSVGs are optional, which is up to the users. \cr After calling this function, all the data including "df_pair" in the list are saved into independent DHF5 databases, and all the DHF5 databases are finally compressed in the file "data_shm.tar". Accordingly, all the corresponding aSVG files listed in the "df_pair" should be compressed in another "tar" file such as "aSVG.tar". If the directory path containing the aSVG files are assigned to \code{svg.dir}, all the SVG files in the diretory are compressed in "aSVGs.tar" automatically. The two tar files compose the database in the Shiny app and should be placed in the "example" folder in the app or uploaded on the user interface. 

#' @param dat.lis A list of data of class "data.frame" or "SummarizedExperiment", where every data should have a unique slot name such as "expr_arab", "expr_chicken", \emph{etc.}. In addition to the data, a pairing data frame describing pairing between the data and aSVG files must be included under the exclusive slot name "df_pair". This data frame has three required columns: the "name" column includes concise names of the data-aSVG pair, the "data" column contains all slot names of the data ("expr_arab", "expr_chicken", \emph{etc.}) and the "aSVG" column contains the aSVG file names corresponding to each data. If one data is related to multiple aSVG files (\emph{e.g.} multiple development  stages), these aSVGs should be concatenated by comma, space, or semicolon, \emph{e.g.} "arabidopsis_thaliana.organ_shm1.svg;arabidopsis_thaliana.organ_shm2.svg". The metadata of data and aSVGs could be optionally included in extra columns.

#' @param dir The directory path to save the "data_shm.tar" file. Default is \code{./data_shm}.
#' @param replace If data with the same slot names in \code{dat.lis} are already saved in \code{dir}, should the \code{dir} be emptied? Default is FALSE. If TRUE, the existing content in \code{dir} will be lost.
#' @param svg.dir The directory path of aSVG files listed in "df_pair". If provded, all SVG files in the directory are compressed in "aSVGs.tar" and saved in \code{dir}. Default is NULL, which requires users to compress the aSVGs in a tar file.

#' @inheritParams HDF5Array::saveHDF5SummarizedExperiment

#' @return A file of "data_shm.tar" is save in \code{dir}. If \code{svg.dir} is assigned a valid value, all relevant SVG files are compressed in "aSVGs.tar" in \code{dir}.

#' @examples

#' ## The examples below demonstrate 1) how to dump Expression Atlas data set into the Shiny database; 2) how to dump GEO data set into the Shiny database; 3) how to include aSVGs of multiple development stages; 4) how to read the database; 5) how to create customized Shiny app with the database. 
#'
#' # 1. Dump data from Expression Atlas into "data_shm.tar" using ExpressionAtlas package (Keays 2019).
#'
#' # The chicken data derived from an RNA-seq analysis on developments of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019) is chosen as example.
#' # The following searches the Expression Atlas for expression data from ‘heart’ and ‘gallus’.
#' library(ExpressionAtlas)
#' cache.pa <- '~/.cache/shm' # The path of cache.
#' all.chk <- read_cache(cache.pa, 'all.chk') # Retrieve data from cache.
#' if (is.null(all.chk)) { # Save downloaded data to cache if it is not cached.
#'  all.chk <- searchAtlasExperiments(properties="heart", species="gallus")
#'  save_cache(dir=cache.pa, overwrite=TRUE, all.chk)
#' }
#'
# Among the matching entries, accession ‘E-MTAB-6769’ is downloaded.
#' all.chk[3, ]
#' rse.chk <- read_cache(cache.pa, 'rse.chk') # Read data from cache.
#' if (is.null(rse.chk)) { # Save downloaded data to cache if it is not cached.
#'   rse.chk <- getAtlasData('E-MTAB-6769')[[1]][[1]]
#'   save_cache(dir=cache.pa, overwrite=TRUE, rse.chk)
#' }

#' # The downloaded data is stored in "SummarizedExperiment" by default (SE, M. Morgan et al. 2018). The experiment design is described in the "colData" slot. The following returns first three rows.
#' colData(rse.chk)[1:3, ]

#' # In the "colData" slot, it is required to define the "sample" and "condition" columns respectively. Both "sample" and "condition" are general terms. The former refers to entities where the numeric data are measured such as cell organelles, tissues, organs, ect. while the latter denotes experimental treatments such as drug dosages, gender, trains, time series, PH values, ect. In the downloaded data, the two columns are not explicitly defined, so "organism_part" and "age" are selected and renamed as "sample" and "condition" respectively.
#' colnames(colData(rse.chk))[c(6, 8)] <- c('condition', 'sample'); colnames(colData(rse.chk))

#' # The raw RNA-Seq count are preprocessed with the following steps: (1) normalization, (2) aggregation of replicates, and (3) filtering of reliable expression data. The details of these steps are explained in the pacakge vignette.
#' \donttest{ browseVignettes('spatialHeatmap') }
#' se.nor.chk <- norm_data(data=rse.chk, norm.fun='ESF', data.trans='log2') # Normalization
#' se.aggr.chk <- aggr_rep(data=se.nor.chk, sam.factor='sample', con.factor='condition', aggr='mean') # Replicate agggregation using mean 
#' se.fil.chk <- filter_data(data=se.aggr.chk, sam.factor='sample', con.factor='condition', pOA=c(0.01, 5), CV=c(0.6, 100), dir=NULL) # Genes are filtered out if not meet these criteria: expression values are at least 5 in at least 1% of all samples, coeffient of variance is between 0.6 and 100.

#' # The aSVG file corresponding with the data is pre-packaged and copied to a temporary directory.
#' dir.svg <- paste0(tempdir(check=TRUE), '/svg_shm') # Temporary directory. 
#' if (!dir.exists(dir.svg)) dir.create(dir.svg)
#' svg.chk <- system.file("extdata/shinyApp/example", 'gallus_gallus.svg', package="spatialHeatmap") # Path of the aSVG file.
#' file.copy(svg.chk, dir.svg, overwrite=TRUE) # Copy the aSVG file.
#'
#' # 2. Dump data from GEO into "data_shm.tar" using GEOquery package (S. Davis and Meltzer 2007).
#'
#' # The Arabidopsis thaliana (Arabidopsis) data from an microarray assay of hypoxia treatment on Arabidopsis root and shoot cell types (Mustroph et al. 2009) is selected as example.
#' # The data set is downloaded with the accession number "GSE14502". It is stored in ExpressionSet container (W. Huber et al. 2015) by default, and then converted to a SummarizedExperiment object.
#' library(GEOquery)
#' gset <- read_cache(cache.pa, 'gset') # Retrieve data from cache.
#' if (is.null(gset)) { # Save downloaded data to cache if it is not cached.
#'   gset <- getGEO("GSE14502", GSEMatrix=TRUE, getGPL=TRUE)[[1]]
#'   save_cache(dir=cache.pa, overwrite=TRUE, gset)
#' }
#' se.sh <- as(gset, "SummarizedExperiment") # Converted to SummarizedExperiment

#' # The gene symbol identifiers are extracted from the rowData component to be used as row names.
#' rownames(se.sh) <- make.names(rowData(se.sh)[, 'Gene.Symbol'])

#' # A slice of the experimental design in colData slot is shown. Both the samples and conditions are contained in the "title" column. The samples are indicated by promoters: pGL2 (root atrichoblast epidermis), pCO2 (root cortex meristematic zone), pSCR (root endodermis), pWOL (root vasculature), etc., and conditions are control and hypoxia.
#' colData(se.sh)[60:63, 1:4]

#' # Since the samples and conditions need to be listed in two independent columns, like the the chicken data above, a targets file is recommended to separate samples and conditions. The main reason to choose this Arabidopdis data is to illusrate the usage of targets file when necessary. A pre-packaged targets file is accessed and partially shown below.
#' sh.tar <- system.file('extdata/shinyApp/example/target_arab.txt', package='spatialHeatmap')
#' target.sh <- read_fr(sh.tar); target.sh[60:63, ]
#' # Load custom the targets file into colData slot.
#' colData(se.sh) <- DataFrame(target.sh)

#' # This data set was already normalized with the RMA algorithm (Gautier et al. 2004). Thus, the pre-processing steps are restricted to aggregation of replicates and filtering of reliably expressed genes. 

#' se.aggr.sh <- aggr_rep(data=se.sh, sam.factor='sample', con.factor='condition', aggr='mean') # Replicate agggregation using mean
#' se.fil.arab <- filter_data(data=se.aggr.sh, sam.factor='sample', con.factor='condition', pOA=c(0.03, 6), CV=c(0.30, 100), dir=NULL) # Filtering of genes with low intensities and variance
#'
#' # Similarly, the aSVG file corresponding to this data is pre-packaged and copied to the same temporary directory.
#' svg.arab <- system.file("extdata/shinyApp/example", 'arabidopsis_thaliana.organ_shm.svg', package="spatialHeatmap") 
#' file.copy(svg.arab, dir.svg, overwrite=TRUE)
#'
#' # 3. The random data and aSVG files of two development stages of Arabidopsis organs.
#' 
#' # The gene expression data is randomly generated and pre-packaged. 
#' pa.growth <- system.file("extdata/shinyApp/example", 'random_data_multiple_aSVGs.txt', package="spatialHeatmap")
#' dat.growth <- read_fr(pa.growth); dat.growth[1:3, ]
#' # Paths of the two corresponsing aSVG files.
#' svg.arab1 <- system.file("extdata/shinyApp/example", 'arabidopsis_thaliana.organ_shm1.svg', package="spatialHeatmap") 
#' svg.arab2 <- system.file("extdata/shinyApp/example", 'arabidopsis_thaliana.organ_shm2.svg', package="spatialHeatmap") 
#' # Copy the two aSVG files to the same temporary directory.
#' file.copy(c(svg.arab1, svg.arab2), dir.svg, overwrite=TRUE)
#'
#' # Make the pairing table, which describes matchings between the data and aSVG files.
#' df.pair <- data.frame(name=c('chicken', 'arab', 'growth'), data=c('expr_chicken', 'expr_arab',  'random_data_multiple_aSVGs'), aSVG=c('gallus_gallus.svg', 'arabidopsis_thaliana.organ_shm.svg', 'arabidopsis_thaliana.organ_shm1.svg;arabidopsis_thaliana.organ_shm2.svg'))
#' # Note that multiple aSVGs should be concatenated by comma, semicolon, or single space.
#' df.pair
#'
#' # Organize the data and pairing table in a list, and create the database.
#' dat.lis <- list(df_pair=df.pair, expr_chicken=se.fil.chk, expr_arab=se.fil.arab, random_data_multiple_aSVGs=dat.growth)
#' # Create the database in a temporary directory "db_shm".
#' dir.db <- paste0(tempdir(check=TRUE), '/db_shm') # Temporary directory. 
#' \donttest{ 
#' if (!dir.exists(dir.db)) dir.create(dir.db)
#' write_hdf5(dat.lis=dat.lis, dir=dir.db, svg.dir=dir.svg, replace=TRUE)
#'
#' # 4. Read data and/or pairing table from "data_shm.tar".
#' dat.lis1 <- read_hdf5(paste0(dir.db, '/data_shm.tar'), names(dat.lis))
#' }
#'
#' # 5. Create customized Shiny app with the database.
#' \donttest{ 
#' if (!dir.exists('~/test_shiny')) dir.create('~/test_shiny')
#' lis.tar <- list(data=paste0(dir.db, '/data_shm.tar'), svg=paste0(dir.db, '/aSVGs.tar'))
#' custom_shiny(lis.tar, app.dir='~/test_shiny')
#' # Run the app.
#' shiny::runApp('~/test_shiny/shinApp') 
#' }
#'
#' # Except "SummarizedExperiment", the database also accepts data in form of "data.frame". In that case, the columns should follow the naming scheme "sample__condition", i.e. a sample and a condition are concatenated by double underscore. The details are seen in the "data" argument of the function "spatial_hm".
#' # The following takes the Arabidopsis data as example.
#' df.arab <- assay(se.fil.arab); df.arab[1:3, 1:3]
#' # The new data list.
#' dat.lis2 <- list(df_pair=df.pair, expr_chicken=se.fil.chk, expr_arab=df.arab, random_data_multiple_aSVGs=dat.growth)
#' 
#' # If the data does not have an corresponding aSVG or vice versa, in the pairing table the slot of missing data or aSVG should be filled with "none". In that case, on the Shiny user interface, users will be prompted to select an aSVG for the unpaired data or select a data for the unpaired aSVG. 
#' # For example, if the aSVG "arabidopsis_thaliana.organ_shm.svg" has no matching data, the pairing table should be made like below.
#' df.pair1 <- data.frame(name=c('chicken', 'arab', 'growth'), data=c('expr_chicken', 'none',  'random_data_multiple_aSVGs'), aSVG=c('gallus_gallus.svg', 'arabidopsis_thaliana.organ_shm.svg', 'arabidopsis_thaliana.organ_shm1.svg;arabidopsis_thaliana.organ_shm2.svg'))
#' df.pair1
#' # The new data list.
#' dat.lis3 <- list(df_pair=df.pair, expr_chicken=se.fil.chk, none='none', random_data_multiple_aSVGs=dat.growth)
#'

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Hervé Pagès (2020). HDF5Array: HDF5 backend for DelayedArray objects. R package version 1.16.1.
#' Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8
#' Davis, Sean, and Paul Meltzer. 2007. “GEOquery: A Bridge Between the Gene Expression Omnibus (GEO) and BioConductor.” Bioinformatics 14: 1846–7
#' Gautier, Laurent, Leslie Cope, Benjamin M. Bolstad, and Rafael A. Irizarry. 2004. “Affy—analysis of Affymetrix GeneChip Data at the Probe Level.” Bioinformatics 20 (3). Oxford, UK: Oxford University Press: 307–15. doi:10.1093/bioinformatics/btg405
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Huber, W., V. J. Carey, R. Gentleman, S. An ders, M. Carlson, B. S. Carvalho, H. C. Bravo, et al. 2015. “Orchestrating High-Throughput Genomic Analysis Wit H Bioconductor.” Nature Methods 12 (2): 115–21. http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export write_hdf5
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment
#' @importFrom HDF5Array saveHDF5SummarizedExperiment

write_hdf5 <- function(dat.lis, dir='./data_shm', replace=FALSE, chunkdim=NULL, level=NULL, verbose=FALSE, svg.dir=NULL) {

  options(stringsAsFactors=FALSE); dir <- normalizePath(dir)
  na.all <- names(dat.lis)
  if (!'df_pair' %in% na.all) stop('The "df_pair" data frame is required!')
  if (is.null(na.all)) stop('Every name slot of the data list should be assigned a value!')
  if ('' %in% na.all) stop('Every name slot of the data list should be assigned a value!')
  df.pair <- dat.lis$df_pair; if (!all(c('name', 'data', 'aSVG') %in% colnames(df.pair))) stop('The "df_pair" table should include at least three columns: name, data, aSVG!')
  if (any(duplicated(df.pair$name))) stop('The entries in "name" column of "df_pair" should be unique!')
  na.all1 <- na.all[na.all!='df_pair']
  if (!all(sort(as.vector(df.pair$data))==sort(na.all1))) stop('The name slots of data in "dat.lis" must be the same with the "data" column in "df_pair"!')
  cat('Data in progress... \n')
  for (i in seq_along(dat.lis)) {

    # Check "sample_condition" format, and stored data in SE.
    data <- dat.lis[[i]]
    if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')) { 
      form <- grepl('__', colnames(data))
      if (sum(form)==0 & na.all[i]!='df_pair') stop(paste0(na.all[i], ': the colnames of assay slot should follow the naming scheme "sample__condition"!'))
      if (na.all[i]=='df_pair') { data <- SummarizedExperiment(assays=list(expr=data)) } else {
        
        dat.lis.all <- check_data(data=data); dat <- dat.lis.all$dat
        row.meta <- dat.lis.all$row.meta
        data <- SummarizedExperiment(assays=list(merge=dat), rowData=row.meta) 

      }

    }
    # Check "sample" and "condition" in colData.
    if (is(data, 'SummarizedExperiment')) {
      form <- grepl('__', colnames(assay(data)))
      cold <- colData(data)
      if (sum(form)==0 & na.all[i]!='df_pair') {
        if (!any(c('sample', 'condition') %in% colnames(colData(data)))) stop(paste0(na.all[i], ': one of the following requirements should be fulfilled: 1) the colnames of assay slot follow the naming scheme "sample__condition"; 2) the "sample", "condition" columns are defined in the colData slot!'))
        if (all(c('sample', 'condition') %in% colnames(cold))) { 
          cna <- paste0(cold$sample, '__', cold$condition)
          if (any(duplicated(cna))) stop(paste0(na.all[i], ': duplicated "sample__condition" replicates are detected! Please use "aggr_rep" to aggregate replicates.'))
        } else if ('sample' %in% colnames(cold)) {
          if (any(duplicated(cold$sample))) stop(paste0(na.all[i], ': the "sample" should not be duplicated in the absence of "condition"! Please use "aggr_rep" to aggregate replicates.'))
        } else if (!('sample' %in% colnames(cold))) stop(paste0(na.all[i], ': "sample" column must be present in the colData slot!'))
     } else if (any(form) & na.all[i]!='df_pair') {
       if (any(duplicated(colnames(assay(data))))) stop(paste0(na.all[i], ': duplicated "sample__condition" replicates are detected! Please use "aggr_rep" to aggregate replicates.'))
     }
     saveHDF5SummarizedExperiment(data, dir=dir, prefix=paste0(na.all[i], '_'), replace=replace, chunkdim=chunkdim, level=level, verbose=verbose)
    } else cat(paste0('Accepted data classes are: data.frame, matrix, DFrame, SummarizedExperiment! This data is not saved: ', na.all[i]), '\n')

  }; wd <- getwd(); setwd(dir)
  file.tar <- paste0(tempdir(check=TRUE), '/data_shm.tar')
  system(paste0('tar -cf ', file.tar, ' *assays.h5 *se.rds'))
  file.remove(list.files(dir, pattern='\\.h5$|\\.rds$|\\.tar$', full.names=TRUE))
  system(paste0('mv ', file.tar, ' .'))
  if (!is.null(svg.dir)) if (dir.exists(svg.dir)) {
    setwd(svg.dir); cat('aSVGs in progress... \n')
    system(paste0('tar -cf ', paste0(dir, '/aSVGs.tar'), ' *.svg'))
  }; setwd(wd); cat('Done! \n')

}






