#' Bar plots of co-clustering optimization results.
#'
#' @param df.res A \code{data.frame} of co-clustering optimization results.
#' @param para.na The targe parameter, which is a column name in \code{df.res}.
#' @param bar.width Width of a single bar.
#' @param thr A y-axis threshold, which will be used to draw a horizontal line.
#' @param title,title.size The plot title and its size.
#' @param xlab,ylab The x and y axis label respectively.
#' @param axis.title.size The size of x and y axis labels.
#' @param x.text.size,y.text.size The size of x and y axis text.
#' @param x.agl,x.vjust The angle and vertical position of x-axis text.
#' @param fill The color of bars.


#' @return An object of ggplot.

#' @examples

#' set.seed(10)
#' df.res <- data.frame(dimred=sample(c('PCA', 'UMAP'), 50, replace=TRUE), cluster=sample(c('wt', 'fg', 'le'), 50, replace=TRUE))
#' opt_bar(df.res=df.res, para.na='cluster', ylab='Remaining outcomes')

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @export 
#' @importFrom ggplot2 ggplot geom_bar position_dodge2 labs theme element_text geom_hline

opt_bar <- function(df.res, para.na, bar.width=0.8, thr=NULL, title=NULL, title.size=25, xlab=NULL, ylab=NULL, axis.title.size=25, x.text.size=25, y.text.size=25, x.agl=80, x.vjust=0.6, fill='#FF6666') {
  runs <- NULL 
  tab <- table(df.res[, para.na]); df <- data.frame(tab)  
  colnames(df) <- c('para.na', 'runs')  
  if (is.null(xlab)) xlab <- para.na  
  gg <- ggplot(df, aes(x=para.na, y=runs)) + 
  geom_bar(color=NA, width=bar.width, fill=fill, position=position_dodge2(width=bar.width, preserve="single"), stat="identity") + labs(title=title, x=xlab, y=ylab)+theme(plot.title = element_text(size=title.size, hjust = 0.5, vjust=1.5)) + 
  theme(axis.text.x=element_text(angle=x.agl, vjust=x.vjust, size = x.text.size), axis.title = element_text(size=axis.title.size), axis.text.y = element_text(size = y.text.size)) 
  if (is.numeric(thr)) gg <- gg + geom_hline(yintercept=thr, linetype="dashed", color = "black", size=1)   
  gg    
}    

 
