#' Perform All Pairwise Comparisons for Selected Spatial Features
#'
#' Perform all pairwise comparisons for all spatial features in the selected column in \code{colData} slot.

#' @param se A \code{SummarizedExperiment} object, which is returned by \code{tar_ref}. The \code{colData} slot is required to contain at least two columns of "features" and "factors" respectively. The \code{rowData} slot can optionally contain a column of discriptions of each gene and the column name should be \code{metadata}. 
#' @param com.factor The column name of spatial features to compare in \code{colData} slot.
#' @param method.adjust The method to adjust p values in multiple testing. The default is \code{BH}.
#' @param return.all Logical. If \code{TRUE}, all the comparison results are returned in a data frame. The default is \code{FALSE}.
#' @param log2.fc The log2-fold change cutoff. The default is 1.
#' @param fdr The FDR cutoff. The default is 0.05.

#' @return If \code{return.all=TRUE}, all comparison results are returned in a data frame. If \code{return.all=FALSE}, the up and down genes are returned in data frames for each feature, and the data frames are organized in a nested list.

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
#' cnt.chk <- system.file('extdata/shinyApp/data/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]
#'
#' # A targets file describing samples and conditions is required for toy data. It should be made
#' # based on the experiment design, which is accessible through the accession number 
#' # "E-MTAB-6769" in the R package ExpressionAtlas. An example targets file is included in this
#' # package and accessed below. 

#' # Access the count table. 
#' cnt.chk <- system.file('extdata/shinyApp/data/count_chicken.txt', package='spatialHeatmap')
#' count.chk <- read.table(cnt.chk, header=TRUE, row.names=1, sep='\t')
#' count.chk[1:3, 1:5]

#' # Access the example targets file. 
#' tar.chk <- system.file('extdata/shinyApp/data/target_chicken.txt', package='spatialHeatmap')
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
#' pOA=c(0.1, 5), CV=c(3.5, 100))
#' # Subset the data.
#' data.sub <- tar_ref(data=se.fil.chk, feature='organism_part', ft.sel=c('brain', 'heart', 'kidney'),
#' variable='age', var.sel=c('day10', 'day12'), com.by='feature', target='brain')
#' # Perform all pairwise comparisons in "organism_part" in "colData" slot.
#' dsq <- deseq2(se=data.sub, com.factor='organism_part')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
#' \cr Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

#' @importFrom SummarizedExperiment assay colData
#' @importFrom utils combn 

deseq2 <- function(se, com.factor, method.adjust='BH', return.all=FALSE, log2.fc=1, fdr=0.05, outliers=outliers, pairwise = FALSE, verbose=TRUE) {
  # save(se, com.factor, method.adjust, return.all, log2.fc, fdr, outliers,  pairwise, verbose, file='deseq2.arg')
  pkg <- check_pkg('DESeq2'); if (is(pkg, 'character')) stop(pkg)
  expr <- assay(se); cdat <- colData(se)
  fct <- cdat[, com.factor] <- as.factor(cdat[, com.factor])
  # Change the column name of target comparison elements to a uniform name, since user-provided names are different.   
  colnames(cdat)[colnames(cdat)==com.factor] <- 'Factor'
  
  com <- combn(x=levels(fct), m=2)
  df.all <- data.frame(rm=rep(NA, nrow(expr)))
  if (pairwise==FALSE) { 
    # Alternatively, if use a design matrix for design, then use a numeric contrast vector to contrasts in "results", such as c(-1,  0, 0, 1). DESeq2 operates on raw counts.
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=expr, colData=cdat, design=~Factor) # "design" does not impact "rlog" and "varianceStabilizingTransformation". Rownames of colData do not need to be identical with countData.
    dds <- DESeq2::DESeq(dds)
    for (i in seq_len(ncol(com))) {
      contr <- c('Factor', com[2, i], com[1, i])
      if (verbose==TRUE) message(contr[2], '-', contr[3])
      df0 <- as.data.frame(DESeq2::results(object=dds, contrast=contr, pAdjustMethod=method.adjust)[, c('log2FoldChange', 'padj')])
      colnames(df0) <- paste0(contr[2], '_VS_', contr[3], '_', colnames(df0)); df.all <- cbind(df.all, df0)
    }; sam.all <- levels(fct) 
  } else if (pairwise==TRUE) { 
    if (grepl('__', cdat[, 'Factor'][1])) wng("If compare by 'feature_variable', please use 'pairwise=FALSE'.") 
    if (all(unique(cdat[, 'Factor'])==unique(cdat[, 'feature'])))  { 
      vari <- cdat$feature; ft <- cdat$variable; com.by <- 'feature'  
    } # Compare by feature.
    if (all(unique(cdat[, 'Factor'])==unique(cdat[, 'variable']))) {
      ft <- cdat$feature; vari <- cdat$variable; com.by <- 'variable'
    } # Compare by variable.
    com <- combn(x=sort(unique(vari)), m=2)
    
    for (i in seq_len(ncol(com))) { # Pairwise comparison.
      w0 <- vari %in% com[, i]; se0 <- se[, w0]
      se0$feature <- factor(se0$feature, levels=unique(se0$feature))
      se0$variable <- factor(se0$variable)
      if (com.by=='feature') {
        # "factor" implies reference level (control), if "contrast" in "results" is not specified
        # the default comparison is treatment VS reference. Regardless of reference, the comparison 
        # can be specified by "contrast" in "results" such as "contrast=c('condition', 'treatment', 'control')" 
        # is "treatment VS control" even though reference is "treatment".
        dds0 <- DESeqDataSet(se0, design = ~ variable + feature)
        contr <- c('feature', com[2, i], com[1, i])
      }
      if (com.by=='variable') { 
        dds0 <- DESeqDataSet(se0, design = ~ feature + variable) 
        contr <- c('variable', com[2, i], com[1, i])
      }; if (verbose==TRUE) message(contr[2], '-', contr[3])
      dds0 <- DESeq(dds0)
      df0 <- as.data.frame(DESeq2::results(object=dds0, contrast=contr, pAdjustMethod=method.adjust)[, c('log2FoldChange', 'padj')])
      colnames(df0) <- paste0(contr[2], '_VS_', contr[3], '_', colnames(df0)); df.all <- cbind(df.all, df0)
    }; sam.all <- unique(vari)
  }; df.all <- df.all[, -1]; if (return.all==TRUE) return(df.all)
  UD <- up_dn(sam.all=sam.all, df.all=df.all, log.fc=abs(log2.fc), fdr=fdr, log.na='log2FoldChange', fdr.na='padj', method='DESeq2', outliers=outliers); return(UD)
  
}


