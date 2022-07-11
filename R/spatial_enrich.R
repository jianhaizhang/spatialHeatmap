#' Identify Spatial Feature-Specifcally Expressed Genes 
#'
#' This functionality is an extension of the spatial heatmap. It identifies spatial feature-specifically expressed genes and thus enables the spatial heatmap to visualize feature-specific profiles. The spatial features include cellular compartments, tissues, organs, \emph{etc}. The function compares the target feature with all other selected features in a pairwise manner. The genes significantly up- or down-regulated in the target feature across all pairwise comparisons are denoted final target feature-specifcally expressed genes. The underlying methods include edgeR (Robinson et al, 2010), limma (Ritchie et al, 2015), DESeq2 (Love et al, 2014), distinct (Tiberi et al, 2020). The feature-specific genes are first detected with each method and can be summarized across methods. \cr In addition to feature-specific genes, this function is also able to identify genes specifically expressed in certain condition or in composite factor. The latter is a combination of multiple expermental factors. \code{E.g.} the spatiotemporal factor is a combination of feature and time points.

#' @param data A \code{SummarizedExperiment} object, which is returned by \code{sub_data}. The \code{colData} slot is required to contain at least two columns of "features" and "factors" respectively. The \code{rowData} slot can optionally contain a column of discriptions of each gene and the column name should be \code{metadata}. 
#' @param methods One or more of \code{edgeR}, \code{limma}, \code{DESeq2}, \code{distinct}. The default is \code{c('edgeR', 'limma')}.
#' @param norm The normalization method (\code{TMM}, \code{RLE}, \code{upperquartile}, \code{none}) in edgeR. The default is \code{TMM}. Details: https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors

#' @param log2.trans.dis Logical, only applicable when \code{distinct} is in \code{methods}. The default is \code{TRUE}, and the count data is transformed to log-2 scale.

#' @param aggr One of \code{mean} (default), \code{median}. The method to aggregated replicates in the data frame of feature-specific genes.

#' @param log2.trans.aggr Logical. If \code{TRUE} (default), the aggregated data (see \code{aggr}) is transformed to log2-scale, included in the returned data frame of feature-specific genes, and would be further used in the spatial heatmaps.
 
#' @param p.adjust The method (\code{holm}, \code{hochberg}, \code{hommel}, \code{bonferroni}, \code{BH}, \code{BY}, \code{fdr}, \code{none}) to adjust p values in multiple hypothesis testing. The default is \code{BH}.
#' @param log2.fc The log2-fold change cutoff. The default is 1.
#' @param fdr The FDR cutoff. The default is 0.05.

#' @return A nested list containing the feature-specific genes summarized across methods within methods.  

#' @examples

#' ## In the following examples, the toy data come from an RNA-seq analysis on development of 7
#' ## chicken organs under 9 time points (Cardoso-Moreira et al. 2019). For conveninece, it is
#' ## included in this package. The complete raw count data are downloaded using the R package
#' ## ExpressionAtlas (Keays 2019) with the accession number "E-MTAB-6769".   
#'
#' library(SummarizedExperiment)
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
#' se.chk <- SummarizedExperiment(assay=count.chk, colData=target.chk)
#' # The "rowData" slot can store a data frame of gene metadata, but not required. Only the 
#' # column named "metadata" will be recognized. 
#' # Pseudo row metadata.
#' metadata <- paste0('meta', seq_len(nrow(count.chk))); metadata[1:3]
#' rowData(se.chk) <- DataFrame(metadata=metadata)
#'
#' # Subset the data by selected features (brain, heart, kidney) and factors (day10, day12).
#' data.sub <- sub_data(data=se.chk, feature='organism_part', features=c('brain', 'heart', 
#' 'kidney'), factor='age', factors=c('day10', 'day12'), com.by='feature', target='brain')
#' 
#' ## As conventions, raw sequencing count data should be normalized and filtered to
#' ## reduce noise. Since normalization will be performed in spatial enrichment, only filtering
#' ## is required.  
#'
#' # Filter out genes with low counts and low variance. Genes with counts over 5 in
#' # at least 10% samples (pOA), and coefficient of variance (CV) between 3.5 and 100 are 
#' # retained.
#' data.sub.fil <- filter_data(data=data.sub, sam.factor='organism_part', con.factor='age',
#' pOA=c(0.1, 5), CV=c(0.7, 100), dir=NULL)
#' # Identify brain-specifically expressed genes relative to heart and kidney, where day10 and
#' # day12 are treated as replicates. 
#' deg.lis <- spatial_enrich(data.sub.fil)
#' # All up- and down-regulated genes in brain across methods. On the right is the data after
#' # replicates aggregated, and will be used in the spatial heatmaps.
#' deg.lis$deg.table[1:3, ]
#' # Up-regulated genes detected by edgeR.
#' deg.lis$lis.up.down$up.lis$edgeR.up[1:5]

#' # The aSVG path.
#' svg.chk <- system.file("extdata/shinyApp/example", "gallus_gallus.svg", 
#' package="spatialHeatmap")
#' # Plot one brain-specific gene in spatial heatmap.
#' spatial_hm(svg.path=svg.chk, data=deg.lis$deg.table, ID=deg.lis$deg.table$gene[1], legend.r=1.9, legend.nrow=2, sub.title.size=7, ncol=2, bar.width=0.11)

#' # Overlap of up-regulated brain-specific genes across methods.
#' deg_ovl(deg.lis$lis.up.down, type='up', plot='upset')
#' deg_ovl(deg.lis$lis.up.down, type='up', plot='matrix')
#' # Overlap of down-regulated brain-specific genes across methods.
#' deg_ovl(deg.lis$lis.up.down, type='down', plot='upset')
#' deg_ovl(deg.lis$lis.up.down, type='down', plot='matrix')
#' # Line graph of gene expression profile.
#' profile_gene(deg.lis$deg.table[1, ])

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9
#' \cr Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' \cr Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
#' \cr Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#' \cr Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
#' \cr Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' \cr Simone Tiberi and Mark D. Robinson. (2020). distinct: distinct: a method for differential analyses via hierarchical permutation tests. R package version 1.2.0. https://github.com/SimoneTiberi/distinct

#' @export spatial_enrich
#' @importFrom SummarizedExperiment colData rowData

spatial_enrich <- function(data, methods=c('edgeR', 'limma'), norm='TMM', log2.trans.dis=TRUE, log2.fc=1, p.adjust='BH', fdr=0.05, aggr='mean', log2.trans.aggr=TRUE) {
  cdat <- colData(data); com.factor <- colnames(cdat)[1]
  edg <- dsq <- lim <- dis <- NULL
  if ('edgeR' %in% methods) { cat('edgeR ... \n')
    edg <- edgeR(data, method.norm=norm, com.factor=com.factor, method.adjust=p.adjust, return.all=FALSE, log2.fc=log2.fc, fdr=fdr)
    cat('Done! \n')
  }
  if ('limma' %in% methods) { cat('limma ... \n')
    lim <- limma(data, m.array=FALSE, method.norm=norm, com.factor=com.factor, method.adjust=p.adjust, return.all=FALSE, log2.fc=log2.fc, fdr=fdr)
    cat('Done! \n') 
  }
  if ('DESeq2' %in% methods) { cat('DESeq2 ... \n')
    dsq <- deseq2(data, com.factor=com.factor, method.adjust=p.adjust, return.all=FALSE, log2.fc=log2.fc, fdr=fdr)
    cat('Done! \n')
  }
  if ('distinct' %in% methods) { cat('distinct ... \n')
    dis <- distt(data, norm.fun='CNF', parameter.list=list(method=norm), log2.trans=log2.trans.dis, com.factor=com.factor, return.all=FALSE, log2.fc=log2.fc, fdr=fdr)
    cat('Done! \n')
  }

  lis <- list(edgeR=edg, limma=lim, DESeq2=dsq, distinct=dis)[c('edgeR', 'limma', 'DESeq2', 'distinct') %in% methods]
  tar <- unique(cdat[, 1][cdat$target=='yes'])[1]
  deg.lis1 <- deg_lis(lis, sam=tar, 'up')
  deg.lis2 <- deg_lis(lis, sam=tar, 'down')
  if (length(deg.lis1) == 1) { # Only one method is selected.
    deg.lis1 <- c(deg.lis1[1], deg.lis1[1])
    names(deg.lis1) <- paste0(names(deg.lis1), c('', '.1')) 
  }
  if (length(deg.lis2) == 1) { # Only one method is selected.
    deg.lis2 <- c(deg.lis2[1], deg.lis2[1])
    names(deg.lis2) <- paste0(names(deg.lis2), c('', '.1'))
  }
  lis.up.dn <- list(up.lis = deg.lis1, down.lis = deg.lis2)
  cat('DEG table ... \n')
  # venn_inter returns matrix, since matrix accepts duplicate row names while data frame not. Some genes might be up in one method while down in other methods, so there could be duplicated rows in df.deg.
  # DataFrame allows duplicated row names.
  df.deg <- DataFrame(rbind(venn_inter(deg.lis1), venn_inter(deg.lis2)))
  cna.deg <- colnames(df.deg)
  for (i in c('total', 'edgeR', 'limma', 'DESeq2', 'distinct')) {
    if (i %in% cna.deg) df.deg[, i] <- as.numeric(df.deg[, i])
  }
  na.deg <- rownames(df.deg)
  if (nrow(df.deg)==0) warning('No up or down regulated genes are detected!')
  df.met <- rowData(data); if (ncol(df.met)>0) {
    idx.met <- grep('^metadata$', colnames(df.met), ignore.case=TRUE)
    if (length(idx.met)>0) df.deg <- cbind(df.deg, df.met[na.deg, idx.met[1], drop = FALSE]) 
  }; 
  # df.deg <- cbind(gene=rownames(df.deg), df.deg)
  dat.nor <- norm_data(assay(data), norm.fun='CNF', parameter.list=list(method=norm), log2.trans=log2.trans.aggr)
  df.aggr <- DataFrame(aggr_rep(dat.nor, sam.factor=NULL, con.factor=NULL, aggr=aggr))
  # Keep DFrame all the time, and include a gene column in the end. 
  df.deg <- cbind(gene=na.deg, df.deg, df.aggr[na.deg, , drop=FALSE])
  cat('Done! \n')
  # The "feature__factor" columns in df.deg will be extacted in 'spatial_hm'.
  return(list(deg.table=df.deg, lis.up.down=lis.up.dn))
}


#' Extract up or down DEGs for provided sample, which are fed to UpsetR.
#'
#' @param lis The list containing up and down DEGs of each method. \emph{E.g.} list(edg=deg.edg, dsq=deg.dsq, lim=deg.lim.r).
#' @param sam The target sample/feature.
#' @param deg One of \code{up} or \code{down}, which means to extract up or down DEGs.

#' @return A list of extract up or down DEGs.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

deg_lis <- function(lis, sam, deg='up') {

  sam.all <- NULL; for (i in lis) { sam.all <- c(sam.all, names(i)) }
  sam.all <- unique(sam.all)
  meths <- names(lis); if ('' %in% meths) stop('Each element in "lis" should have a name!')
  lis.all <- NULL; for (i in meths) {
    lis0 <- lis[[i]]
    rna <- rownames(lis0[[sam]][[paste0(sam, '_', deg)]])
    if (length(rna)>0) { lis1 <- list(rna); names(lis1) <- paste0(i, '.', deg); lis.all <- c(lis.all, lis1) } else cat('No genes detected in', i, '! \n')
  }; return(lis.all)

}


#' Translate overlap up/down genes from different methods in a list to a data frame
#'
#' @param lis.all The list of either extracted up or down DEGs returned by \code{deg_lis}

#' @return A data frame of up/down DEGs across methods.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber, Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2020). gplots: Various R Programming Tools for Plotting Data. R package version 3.1.1. https://CRAN.R-project.org/package=gplots
 
#' @importFrom gplots venn 

venn_inter <- function(lis.all) {
  if (is.null(lis.all)) return (data.frame())
  gen.all <- unique(unlist(lis.all))
  # Create an empty data frame, where rows are all genes and columns are methods.
  zero <- rep(0, length(gen.all))
  df.all <- as.data.frame(matrix(rep(0, length(gen.all)*length(lis.all)), ncol=length(lis.all)))
  cna <- names(lis.all)
  suf <- unique(gsub('(.*)\\.(.*)', '\\2', cna))
  meth <- unique(gsub('(.*)\\.(.*)', '\\1', cna))
  names(lis.all) <- colnames(df.all) <- meth
  rownames(df.all) <- gen.all
  # Retrieve all overlaps.
  lis <- venn(lis.all, show.plot=FALSE) 
  inter <- attributes(lis)$intersections
  # Translate all overlaps into a data frame.
  for (i in seq_along(inter)) {
    lis0 <- inter[[i]]; na0 <- strsplit(names(inter[i]), ':')[[1]]
    w1 <- which(gen.all %in% lis0); w2 <- which(meth %in% na0)  
    df.all[w1, w2] <- 1
  }
  df.all <- cbind(type=suf, total=rowSums(df.all), df.all)
  # Matrix accepts dupliated rows. Some genes might be up in one method while down in other methods, so when combine up and down table in a single table, there could be duplicated row names. In data frame, duplicated row names are appended 1.
  df.all <- as.matrix(df.all[order(df.all$total, decreasing=TRUE), ])
  return(df.all)
}



