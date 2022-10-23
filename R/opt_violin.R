#' Violin plots of co-clustering validation results
#'
#' @param data A \code{data.frame} of co-clustering validation results.
#' @param para.na Targe parameters, which are one or multiple column names in \code{data}.
#' @param bar.width Width of the bar.
#' @param thr A y-axis threshold, which will be used to draw a horizontal line.
#' @param title,title.size The plot title and its size.
#' @param xlab,ylab The x and y axis label respectively.
#' @param axis.title.size The size of x and y axis labels.
#' @param x.text.size,y.text.size The size of x and y axis text.
#' @param x.agl,x.vjust The angle and vertical position of x-axis text.


#' @return An object of ggplot.

#' @examples

#' set.seed(10)
#' data <- data.frame(auc=round(runif(30, 0, 1), 2), accuracy=round(runif(30, 0, 1), 2))
#' opt_violin(data=data, para.na=c('auc', 'accuracy'))

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @export 
#' @importFrom ggplot2 ggplot geom_violin geom_boxplot labs theme element_text geom_hline

opt_violin <- function(data, para.na, bar.width=0.1, thr=NULL, title=NULL, title.size=25, xlab=NULL, ylab=NULL, axis.title.size=25, x.text.size=25, y.text.size=25, x.agl=0, x.vjust=0.6) {
  x <- y <- NULL
  df <- NULL; for (i in seq_len(length(para.na))) {
    na0 <- para.na[i] 
    df0 <- data.frame(x=na0, y=data[, na0, drop=TRUE])   
    df <- rbind(df, df0)
  } 
  gg <- ggplot(df, aes(x = x, y = y)) + geom_violin(trim = FALSE) +
  geom_boxplot(width=bar.width) + labs(title=title, x=xlab, y=ylab) +
  theme(plot.title = element_text(size=title.size, hjust = 0.5, vjust=1.5)) +
  theme(axis.text.x=element_text(angle=x.agl, vjust=x.vjust, size = x.text.size), axis.title = element_text(size=axis.title.size), axis.text.y = element_text(size = y.text.size))
  if (is.numeric(thr)) gg <- gg + geom_hline(yintercept=thr, linetype="dashed", color = "black", size=1)  
  gg  
} 
 
