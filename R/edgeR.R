#' Perform All Pairwise Comparisons for Selected Spatial Features
#'
#' Perform all pairwise comparisons for all spatial features in the selected column in \code{colData} slot.

#' @param se A \code{SummarizedExperiment} object, which is returned by \code{sub_data}. The \code{colData} slot is required to contain at least two columns of "features" and "factors" respectively. The \code{rowData} slot can optionally contain a column of discriptions of each gene and the column name should be \code{metadata}. 
#' @param method.norm The normalization method in edgeR (Robinson et al, 2010). The default is \code{TMM}.
#' @param com.factor The column name of spatial features to compare in \code{colData} slot.
#' @param method.adjust The method to adjust p values in multiple testing. The default is \code{BH}.
#' @param return.all Logical. If \code{TRUE}, all the comparison results are returned in a data frame. The default is \code{FALSE}.
#' @param log2.fc The log2-fold change cutoff. The default is 1.
#' @param fdr The FDR cutoff. The default is 0.05.

#' @return If \code{return.all=TRUE}, all comparison results are returned in a data frame. If \code{return.all=FALSE}, the up and down genes are returned in data frames for each feature, and the data frames are organized in a nested list.

#' @details
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
#' edg <- edgeR(se=data.sub, com.factor='organism_part')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
# \cr Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140

#' @importFrom SummarizedExperiment assay colData
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT topTags
#' @importFrom limma makeContrasts
#' @importFrom stats model.matrix 
#' @importFrom utils combn

edgeR <- function(se, method.norm='TMM', com.factor, method.adjust='BH', return.all=FALSE, log2.fc=1, fdr=0.05) {

  df.cnt <- assay(se); y <- DGEList(counts=df.cnt)
  cat('Normalizing ...', '\n') # norm.factors are used in glmFit.
  # To store normalized counts in 'se' and set method.norm='none' is not right, since the 'norm.factors' are essentially used but they are 1 if method.norm='none'.
  y <- calcNormFactors(y, method=method.norm)
  fct <- factor(colData(se)[, com.factor]); design <- model.matrix(~0+fct)
  colnames(design) <- levels(fct); rownames(design) <- colnames(df.cnt) # Duplicated row names in design are allowed.
  y <- estimateDisp(y, design); fit <- glmFit(y, design)

  com <- combn(x=colnames(design), m=2); con <- paste(com[1,], com[2,], sep="-")
  con.ma <- makeContrasts(contrasts=con, levels=design)
  cna.con <- colnames(con.ma); cna.con1 <- sub('-', '_VS_', cna.con)
  cat('Computing DEGs ...', '\n')
  df.all <- data.frame(rm=rep(NA, nrow(df.cnt))); for (i in seq_len(length(cna.con))) {

    lrt <- glmLRT(fit, contrast=con.ma[, i])
    tag <- as.data.frame(topTags(lrt, n=nrow(lrt), adjust.method=method.adjust, sort.by="none", p.value=1))
    cat(cna.con1[i], '\n')
    colnames(tag) <- paste0(cna.con1[i], '_', colnames(tag)); df.all <- cbind(df.all, tag)

  }; df.all <- df.all[, -1]
  if (return.all==TRUE) return(df.all)
  UD <- up_dn(sam.all=levels(fct), df.all=df.all, log.fc=abs(log2.fc), fdr=fdr, log.na='logFC', fdr.na='FDR'); return(UD)

}



