#' Plotting the clusters returned by K-means clustering
#'
#' This function is designed to plot the clusters returned by K-means clustering.

#' @param data A numeric \code{data.frame} scaled with \code{scale}.
#' @param res The clustering results returned by \code{kmeans}.
#' @param query A rowname in \code{data}. The cluster containing this rowname will be labeled.
#' @param dimred The dimension reduction method, one of \code{'PCA'}, \code{'UMAP'}, or \code{'TSNE'}.
#' @param title The plot title. 

#' @return A ggplot.

#' @examples 

#' set.seed(10)
#' data <- iris[, 1:4]
#' rownames(data) <- paste0('gene', seq_len(nrow(data)))
#' dat.scl <- scale(data)
#' clus <- kmeans(dat.scl, 6)
#' plot_kmeans(data=dat.scl, res=clus, query='gene1', dimred='TSNE') 

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' Melville J (2022). _uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction_. R package version 0.1.14, <https://CRAN.R-project.org/package=uwot>
#' Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation, URL: https://github.com/jkrijthe/Rtsne
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @export
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse labs xlab ylab ggtitle theme element_text 
#' @importFrom stats princomp

plot_kmeans <- function(data, res, query, dimred='TSNE', title='Kmeans Clustering Results') {
  dim1 <- dim2 <- cluster <- NULL
  # Dimension reduction.  
  if (dimred=='PCA') dim.df <- princomp(data)$scores[, 1:2]
  if (dimred=='UMAP') {
    pkg <- check_pkg('uwot'); if (is(pkg, 'character')) stop(pkg) 
    dim.df <- uwot::umap(data)
  }
  if (dimred=='TSNE') {
    pkg <- check_pkg('Rtsne'); if (is(pkg, 'character')) stop(pkg) 
    dim.df <- Rtsne::Rtsne(data, check_duplicates = FALSE)$Y; rownames(dim.df) <- row.names(data)
  }
  # Combine dimensions and clusters.
  colnames(dim.df) <- c('dim1', 'dim2')
  dim.cl <- data.frame(dim.df, cluster=res$cluster)
  clus.all <- dim.cl$cluster <- paste0('clus', dim.cl$cluster)
  # Indicating the cluster containing the query gene.
  w <- which(rownames(dim.cl)==query)
  cl.query <- clus.all[w]
  cl.qu.w <- which(clus.all==cl.query)
  dim.cl$cluster[cl.qu.w] <- paste0(cl.query, '.query')
  xylab <- paste0(dimred, seq_len(2))

  ggplot(dim.cl, aes(x = dim1, y = dim2, color = as.factor(cluster), fill = as.factor(cluster))) + geom_point() + 
  stat_ellipse(type = "t", geom = "polygon", alpha = 0.4) + labs(color='Clusters', fill='Clusters') + xlab(xylab[1]) + 
  ylab(xylab[2]) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
}


