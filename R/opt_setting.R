#' Bar plots of co-clustering optimization results.
#'
#' @param df.res A \code{data.frame} of validation results in co-clustering optimization.
#' @param nas The target metrics, which are column names of \code{df.res}.
#' @param summary The method to summarize ranks of \code{nas} across datasets.

#' @return An object of ggplot.

#' @examples

#' set.seed(10)
#' dimred <- c('PCA', 'UMAP'); dims <- seq(5, 80, 15) 
#' graph <- c('knn', 'snn'); cluster <- c('wt', 'fg', 'le') 
#' df.para <- expand.grid(dataset=c('dataset1', 'dataset2'), norm='FCT', fil='fil1', dimred=dimred, dims=dims, graph=graph, cluster=cluster, stringsAsFactors = FALSE) 
#' df.para$auc <- round(runif(nrow(df.para), 0, 1), 2)
#' df.para$accuracy <- round(runif(nrow(df.para), 0, 1), 2)
#' df.para[1:5, ] 
#' opt_setting(df.para, nas=c('auc', 'accuracy'))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @export 

opt_setting <- function(df.res, nas, summary='mean') { 
  dataset <- NULL 
  dat.na <- unique(df.res$dataset)   
  dat1 <- subset(df.res, dataset=='dataset1') 
  dat1$na.pas <- paste0(dat1$norm, dat1$fil, dat1$dimred, dat1$dims, dat1$graph, dat1$cluster)  
  # Results from each dataset will be sorted in the same way.  
  ord <- order(dat1$na.pas); dat1 <- dat1[ord, ]  
 
  df.rank <- NULL 
  for (i in dat.na) { 
    dat0 <- subset(df.res, dataset==i)  
    dat0$na.pas <- paste0(dat0$norm, dat0$fil, dat0$dimred, dat0$dims, dat0$graph, dat0$cluster) 
    if (!all(sort(dat1$na.pas)==sort(dat0$na.pas))) {   
      message('Different column names are detected!'); return()    
    }
    dat0 <- dat0[ord, ]; df.sel <- dat0[, nas]
    df.rank0 <- sapply(df.sel, rank)  
    # Combine ranks from each dataset in a table. 
    df.rank <- cbind(df.rank, df.rank0)   
  } # Summary ranks across datasets.  
    if (summary=='mean') summ <- rowMeans(df.rank)   
    if (summary=='median') summ <- Biobase::rowMedians(df.rank) 
    data.rank <- cbind(dat1, rank=summ)[, c('rank', 'norm', 'fil', 'dimred', 'dims', 'graph', 'cluster')] 
    data.rank <- data.rank[order(data.rank[, 'rank'], decreasing=FALSE), ]  
    data.rank[, 'rank'] <- order(data.rank[, 'rank'], decreasing=TRUE) 
    data.rank <- data.rank[data.rank[, 'rank'], ]; data.rank 
}  


 
