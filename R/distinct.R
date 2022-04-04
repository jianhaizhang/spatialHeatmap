#' Perform All Pairwise Comparisons for Selected Spatial Features
#'
#' Perform all pairwise comparisons for all spatial features in the selected column in \code{colData} slot.

#' @param se A \code{SummarizedExperiment} object, which is returned by \code{sub_data}. The \code{colData} slot is required to contain at least two columns of "features" and "factors" respectively. The \code{rowData} slot can optionally contain a column of discriptions of each gene and the column name should be \code{metadata}. 
#' @param com.factor The column name of spatial features to compare in \code{colData} slot.
#' @param return.all Logical. If \code{TRUE}, all the comparison results are returned in a data frame. The default is \code{FALSE}.
#' @param log2.fc The log2-fold change cutoff. The default is 1.
#' @param fdr The FDR cutoff. The default is 0.05.

#' @inheritParams norm_data

#' @return If \code{return.all=TRUE}, all comparison results are returned in a data frame. If \code{return.all=FALSE}, the up and down genes are returned in data frames for each feature, and the data frames are organized in a nested list.

#' @details
#' Requires three columns in colData: sample, replicates, condition. All condition pairs are compared for each sample.
#' Compare samples across all conditions: con.factor=NULL, which means actual samples are treated as conditions internally, and the sample column is assumed to include only one sample. Conversely, compare conditions across all samples: sam.factor=NULL, which means actual samples are treated as one sample internally.
#' log2FC of con1 VS con2 for gene1 in the output approximates to log2(mean(con1)/mean(con2)) in the un-normalized count matrix for gene1.

#' @keywords Internal
#' @noRd

#' @examples

#' ## In the following examples, the toy data come from an RNA-seq analysis on development of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, it is
#' ## included in this package. The complete raw count data are downloaded using the R package
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769".   
#'
#' ## Set up toy data.
#' 
#' # Access toy data. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing samples and conditions is required for toy data. It should be made
#' # based on the experiment design, which is accessible through the accession number 
#' # "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in this
#' # package and accessed below. 

#' # Access the count table. 
#' cnt.chk <- system.file('extdata/shinyApp/example/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]

#' # Access the example targets file. 
#' tar.chk <- system.file('extdata/shinyApp/example/target_chicken.txt', package='spatialHeatmap')
#' target.chk <- read.table(tar.chk, header=TRUE, row.names=1, sep='\t')
#' # Every column in toy data corresponds with a row in targets file. 
#' target.chk[1:5, ]
#' # Store toy data in "SummarizedExperiment".
#' library(SummarizedExperiment)
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene metadata, but not required. Only the 
#' # column named "metadata" will be recognized. 
#' # Pseudo row metadata.
#' metadata <- paste0('meta', seq_len(nrow(count.chk))); metadata[1:3]
#' rowData(se.chk) <- DataFrame(metadata=metadata)
#'
#' ## As conventions, raw sequencing count data should be normalized and filtered to
#' ## reduce noise. Since normalization will be performed in spatial enrichment, only filtering
#' ## is required before subsetting the data.  
#'
#' # Filter out genes with low counts and low variance. Genes with counts over 5 in
#' # at least 10% samples (pOA), and coefficient of variance (CV) between 3.5 and 100 are 
#' # retained.
#' se.fil.chk <- filter_data(data=se.chk, sam.factor='organism_part', con.factor='age',
#' pOA=c(0.1, 5), CV=c(3.5, 100), dir=NULL)
#' # Subset the data.
#' data.sub <- sub_data(data=se.fil.chk, feature='organism_part', features=c('brain', 'heart', 'kidney'),
#' factor='age', factors=c('day10', 'day12'), com.by='feature', target='brain')
#' # Perform all pairwise comparisons in "organism_part" in "colData" slot.
#' dis <- distinct(se=data.sub, com.factor='organism_part')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
#' \cr Simone Tiberi and Mark D. Robinson. (2020). distinct: distinct: a method for differential analyses via hierarchical permutation tests. R package version 1.2.0. https://github.com/SimoneTiberi/distinct

#' @importFrom SummarizedExperiment assay colData colData<-
#' @importFrom distinct distinct_test log2_FC
#' @importFrom utils combn
#' @importFrom stats model.matrix

distt <- function(se, norm.fun='CNF', parameter.list=NULL, log2.trans=TRUE, com.factor, return.all=FALSE, log2.fc=1, fdr=0.05) {
  filtered <- NULL
  if (!is.null(norm.fun)) { cat('Normalizing counts ... \n')
  se <- norm_data(data=se, norm.fun=norm.fun, parameter.list=parameter.list, log2.trans=log2.trans) }
   
  # Replicate names are the rownames of design matrix, and should not be duplicated, so avoid duplicated replicate names. 
  # cdat <- as.data.frame(colData(se)) # Duplicated row names are numbered.
  cdat <- colData(se); cdat$OneSample <- 'sample'; sam.factor <- 'OneSample'
  con.factor <- com.factor

  if (any(duplicated(rownames(cdat)))) rownames(cdat) <- paste0(rownames(cdat), '.', seq_along(rownames(cdat)))
  cdat$rep.distt <- rownames(cdat); colData(se) <- DataFrame(cdat)
  sam.all <- unique(cdat[, sam.factor])
  # Compare each pair of conditions for each sample.
  lis.all <- list(); for (i in sam.all) {
    idx <- cdat[, sam.factor] %in% i
    cdat0 <- cdat[idx, , drop=FALSE]
    con.all0 <- cdat0[, con.factor]; con.uni0 <- unique(con.all0)
    # Two or more conditions are required.
    if (length(unique(con.all0))==1) { cat(i, 'is skipped, since only one condition is detected! \n'); next }
    com <- combn(con.uni0, 2); rna <- rownames(se) 
    df.all1 <- data.frame(rm=rep(NA, length(rna)), row.names=rna)
    # Compare each condition pair of a sample.
    for (j in seq_len(ncol(com))) {
      con.all1 <- com[, j]
      idx1 <- (cdat[, con.factor] %in% con.all1) & idx
      se0 <- se[, idx1]
      fct <- factor(colData(se0)[, con.factor]); dsg <- model.matrix(~fct)
      vs <- paste0(levels(fct), collapse='_VS_')
      cat(i, ':', vs, '... \n')
      rownames(dsg) <- colData(se0)[, 'rep.distt']
      # set.seed(61217)
      # Genes with values <=0 in less than min_non_zero_cells are disgarded.
      res0 <- distinct_test(x=se0, name_assays_expression=assayNames(se0)[1], name_cluster=sam.factor, name_sample='rep.distt', design=dsg, column_to_test=2, min_non_zero_cells=ceiling((min(table(fct))+1)/2), n_cores=2)
      # Log2 FC requires count matrix.
      SummarizedExperiment::assay(se0) <- 2^SummarizedExperiment::assay(se0)
      res1 <- log2_FC(res=res0, x=se0, name_assays_expression=assayNames(se0)[1], name_cluster=sam.factor, name_group=con.factor)
      cna.res <- colnames(res1)
      colnames(res1)[cna.res=='p_adj.loc'] <- paste0(vs, '_FDR') 
      colnames(res1)[grepl('^log2FC', cna.res)] <- paste0(vs, '_log2FC')
      rownames(res1) <- res1$gene; res1 <- subset(res1, filtered==FALSE)
      res1 <- res1[, grepl('_FDR$|_log2FC$|^filtered$', colnames(res1))]
      # Avoid NA in rownames.
      int0 <- intersect(rownames(df.all1), rownames(res1))
      df.all1 <- cbind(df.all1[int0, ], res1[int0, ])
    }; df.all1 <- df.all1[, -1]
    if (return.all==TRUE) return(df.all1)
    UD <- up_dn(sam.all=con.uni0, df.all=df.all1, log.fc=abs(log2.fc), fdr=fdr, log.na='log2FC', fdr.na='FDR'); cat('\n')
    lis0 <- list(UD); names(lis0) <- i; lis.all <- c(lis.all, lis0)
  }; return(lis.all)
}



