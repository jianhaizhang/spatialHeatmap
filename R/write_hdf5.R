#' Construct DHF5 Database Backend for the Shiny App
#'
#' This is a convenience function for constructing the database backend in the Shiny app (\link{shiny_all}). The data to store in the database should be in the class of "data.frame" or "SummarizedExperiment" and should be formatted according to the conventions in the "data" argument of \link{spatial_hm}. After formatted, all these data should be arranged in a list and each data slot should have a unique name such as "data1", "data2", \emph{etc.}. In addition, a data frame describing the matching relationship between the data and aSVG files must also be included in the list with the exclusive slot name "df_pair". This data frame should contain at least two columns: data, aSVG. The data column includes the slot names of all data in the list ("data1", "data2", \emph{etc.}), and the aSVG column includes the aSVG file names corresponding to each data respectively such as "gallus_gallus.svg", \emph{etc.} If one data is related to multiple aSVG files (\emph{e.g.} multiple development  stages), these aSVGs should be concatenated by comma, space, or semicolon, \emph{e.g.} "arabidopsis_thaliana.organ_shm1.svg;arabidopsis_thaliana.organ_shm2.svg". Inclusion of other columns providing metadata of the data and aSVGs are optional, which is up to the users. After calling this function, all the data including "df_pair" in the list are saved into independent DHF5 databases, and all the DHF5 databases are finally compressed in the file "data_shm.tar". Accordingly, all the corresponding aSVG files listed in the "df_pair" should be compressed in another "tar" file such as "aSVG.tar". The two tar files compose the database backend in the Shiny app.

#' @param dat.lis A list of data of class "data.frame" or "SummarizedExperiment", where every data should have a unique slot name such as "data1", "data2", \emph{etc.}. In addition to the data, a data frame describing pairing between the data and aSVG files must be included under the exclusive slot name "df_pair", where the "data" column contains all slot names of the data ("data1", "data2", \emph{etc.}) and the "aSVG" column contains the aSVG file names corresponding to each data. If one data is related to multiple aSVG files (\emph{e.g.} multiple development  stages), these aSVGs should be concatenated by comma, space, or semicolon, \emph{e.g.} "arabidopsis_thaliana.organ_shm1.svg;arabidopsis_thaliana.organ_shm2.svg". The metadata of data and aSVGs could be optionally included in extra columns.

#' @param dir The directory path to save the "data_shm.tar" file. Default is \code{./data_shm}.
#' @param replace If data with the same slot names in \code{dat.lis} are already saved in \code{dir}, should the \code{dir} be emptied? Default is FALSE. If TRUE, the existing content in \code{dir} will be lost.

#' @inheritParams HDF5Array::saveHDF5SummarizedExperiment

#' @return A file of "data_shm.tar" is save in \code{dir}.

#' @examples

#' ## In the following examples, the 2 toy data come from an RNA-seq analysis on developments of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are included in this package. The complete raw count data are downloaded using the R package ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as a "data frame" input to exemplify data with simple samples/conditions, while toy data2 as "SummarizedExperiment" to illustrate data involving complex samples/conditions.   
#'
#' ## Set up toy data.
#'
#' ## In the following examples, the 2 toy data come from an RNA-seq analysis on development of 7 chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, they are included in this package. The complete raw count data are downloaded using the R package ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769". Toy data1 is used as a "data frame" input to exemplify data of simple samples/conditions, while toy data2 as "SummarizedExperiment" to illustrate data involving complex samples/conditions.   
#' 
#' ## Set up toy data.
#' 
#' # Access toy data1.
#' cnt.chk.simple <- system.file('extdata/shinyApp/example/count_chicken_simple.txt', package='spatialHeatmap')
#' df.chk <- read.table(cnt.chk.simple, header=TRUE, row.names=1, sep='\t', check.names=FALSE)
#' # Columns follow the namig scheme "sample__condition", where "sample" and "condition" stands for organs and time points respectively.
#' df.chk[1:3, ]
#'
#' # A column of gene annotation can be appended to the data frame, but is not required.  
#' ann <- paste0('ann', seq_len(nrow(df.chk))); ann[1:3]
#' df.chk <- cbind(df.chk, ann=ann)
#' df.chk[1:3, ]
#'
#' # Access toy data2. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing samples and conditions is required for toy data2. It should be made based on the experiment design, which is accessible through the accession number "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in this package and accessed below. 

#' # Access the example targets file. 
#' tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
#' target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
#' # Every column in toy data2 corresponds with a row in targets file. 
#' target.chk[1:5, ]
#' # Store toy data2 in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene annotation, but not required.
#' rowData(se.chk) <- DataFrame(ann=ann)
#'
#' # Aggregate "sample_condition" replicates in toy data1.
#' df.aggr.chk <- aggr_rep(data=df.chk, aggr='mean')
#' df.aggr.chk[1:3, ]
#'
#' # Aggregate "sample_condition" replicates in toy data2, where "sample" is "organism_part" and "condition" is "age". 
#' se.aggr.chk <- aggr_rep(data=se.chk, sam.factor='organism_part', con.factor='age', aggr='mean')
#' assay(se.aggr.chk)[1:3, 1:3]

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Hervé Pagès (2020). HDF5Array: HDF5 backend for DelayedArray objects. R package version 1.16.1.

#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. "Differential Expression Analysis of Multifactor RNA-Seq Experiments with Respect to Biological Variation." Nucleic Acids Research 40 (10): 4288–97
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export write_hdf5
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
library(HDF5Array); library(SummarizedExperiment)

write_hdf5 <- function(dat.lis, dir='./data_shm', replace=FALSE, chunkdim=NULL, level=NULL, verbose=FALSE) {

  options(stringsAsFactors=FALSE); na.all <- names(dat.lis)
  if (! 'df_pair' %in% na.all) stop('The "df_pair" data frame is required!')
  if (is.null(na.all)) stop('Every name slot of the data list should be assigned a value!')
  if ('' %in% na.all) stop('Every name slot of the data list should be assigned a value!')
  cat('In progress... \n')
  for (i in seq_along(dat.lis)) {

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
    } else stop(paste0('Accepted data classes are: data.frame, matrix, DFrame, SummarizedExperiment! Please check data: ', na.all[i]))

  }; wd <- getwd(); setwd(dir)
  file.tar <- paste0(tempdir(check=TRUE), '/data_shm.tar')
  system(paste0('tar -cf ', file.tar, ' *assays.h5 *se.rds'))
  file.remove(list.files(dir, pattern='\\.h5$|\\.rds$|\\.tar$', full.names=TRUE))
  system(paste0('mv ', file.tar, ' .')); setwd(wd)
  cat('Done! \n')

}


dat.lis <- list(df_pair=df.pair, arab=dat.arab, chicken=se.aggr.chk, growth=dat.growth)
write_hdf5(dat.lis, dir='~/test5', replace=T)

df.pair <- data.frame(row.names=c('arab', 'chicken', 'growth'), data=c('expr_arab', 'expr_chicken', 'random_data_multiple_aSVGs'), aSVG=c('arabidopsis_thaliana.organ_shm.svg', 'gallus_gallus.svg', 'arabidopsis_thaliana.organ_shm1.svg;arabidopsis_thaliana.organ_shm2.svg'))
write.table(df.pair, '~/test3/df_pair.txt', sep='\t', row.names=T, col.names=T)
read_fr('~/test3/df_pair.txt', header=T)

pa.arab <- system.file("extdata/shinyApp/example", 'expr_arab.txt', package="spatialHeatmap")
dat.arab <- read_fr(pa.arab)
pa.growth <- system.file("extdata/shinyApp/example", 'random_data_multiple_aSVGs.txt', package="spatialHeatmap")
dat.growth <- read_fr(pa.growth)
 
     # Access toy data2. 
     cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
     count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
     count.chk[1:3, 1:5]
     
     # A targets file describing samples and conditions is required for toy data2. It should be made based on the experiment design, which is accessible through the accession number "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in this package and accessed below. 
     # Access the example targets file. 
     tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
     target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
     # Every column in toy data2 corresponds with a row in targets file. 
     target.chk[1:5, ]
     colnames(target.chk)[c(6, 8)] <- c('condition', 'sample')
     # Store toy data2 in "SummarizedExperiment".
     library(SummarizedExperiment)
     se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
     # The "rowData" slot can store a data frame of gene annotation, but not required.
     ann <- paste0('ann', seq_len(nrow(se.chk))); ann[1:3]
     rowData(se.chk) <- DataFrame(ann=ann)
     se.aggr.chk <- aggr_rep(data=se.chk, sam.factor='sample', con.factor='condition', aggr='mean')

read_hdf5 <- function(file, prefix) {

  options(stringsAsFactors=FALSE)
  dir <- paste0(tempdir(check=TRUE), '/data_shm')
  if (!dir.exists(dir)) dir.create(dir)
  lis <- list(); for (i in prefix) {

    system(paste0('tar -xf ', file, ' -C ', dir, ' ', paste0(i, c('_assays.h5', '_se.rds'), collapse=' ')))
    dat <- loadHDF5SummarizedExperiment(dir=dir, prefix=paste0(i, '_'))
    assay(dat) <- as.data.frame(assay(dat))
    assay.na <- names(assays(dat))
    if (!is.null(assay.na)) if (assay.na[1]=='merge') dat <- cbind(assay(dat), rowData(dat))
    if (i=='df_pair') dat <- assay(dat)
    dat <- list(dat); names(dat) <- i; lis <- c(lis, dat)

  }; return(lis)

}

dat.lis1 <- read_hdf5('~/test5/data_shm.tar', names(dat.lis))



