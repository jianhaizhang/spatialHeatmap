#' Using the elbow method to find the optimal number of clusters in K-means clustering
#'
#' This function is designed to find the optimal number of clusters in K-means clustering using the elbow method. 

#' @param data A numeric \code{data.frame} scaled with \code{scale}.
#' @param max.k The max number of clusters to consider.
#' @param title A character string of the elbow plot title.

#' @return A ggplot.

#' @examples 

#' data <- iris[, 1:4]
#' dat.scl <- scale(data)
#' optimal_k(dat.scl, 5)

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}
#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_x_continuous xlab ylab ggtitle theme element_text
#' @importFrom stats kmeans

optimal_k <- function(data, max.k, title='Elbow Method for Finding the Optimal k') {
  wss <- k <- total.withinss <- NULL
  wss <- function(data, k) {
    cluster <- kmeans(data, k); return (cluster$tot.withinss)
  }
  ks <- seq(2, max.k, by=1)
  # Run algorithm over a range of k 
  wss.op <- unlist(lapply(ks, wss, data=data))
  elbow <- data.frame(k=ks, total.withinss= wss.op)
  ggplot(elbow, aes(x = k, y = total.withinss)) + geom_point() + geom_line() + scale_x_continuous(breaks = ks) + xlab('Max         clusters') + ylab('Total withinss') + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
}


