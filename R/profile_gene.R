#' Plot Gene Expression Profiles in a Data Frame
#'
#' @param data A data frame, where rows are genes and columns are features/conditions.
#' @param scale The way to to scale the data. If \code{none} (default), no scaling. If \code{row}, the data is scalaed independently. If \code{all}, all the data is scaled as a whole.
#' @param x.title,y.title X-axis title and Y-axis title respectively.
#' @param text.size The size of axis title and text.
#' @param text.angle The angle of axis text.

#' @return An image of ggplot.

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @examples

#' # See examples in the function "spatial_enrich".

#' @seealso \code{spatial_enrich}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' \cr Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

#' @export profile_gene
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_line theme labs element_text element_rect element_line

profile_gene <- function(data, scale='none', x.title='Sample/conditions', y.title='Value', text.size=15, text.angle=45) {                                                                                   
  if (all(c('gene', 'type', 'total') %in% colnames(data))) { # Data frame of spatial enrichment.
    data <- subset(data, !duplicated(gene)); rownames(data) <- data$gene                        
    data <- data[, !colnames(data) %in% c('gene', 'type', 'total', 'metadata', 'edgeR', 'limma', 'DESeq2', 'distinct'), drop=FALSE]            
  }                                                                                                                                
  if (scale=='row') data <- t(scale(t(data))) else if (scale=='all') data <- scale_all(data)                                       
  # convert to long format                                                                                                         
  df.long <- reshape2::melt(as.matrix(data))                                                                                       
  colnames(df.long) <- c('Genes', 'Samples', 'Value')                                                                              
  # The levels in melted data frame has the same order (left to right) with row order (top to bottom) in original data frame       before melted.                                                                                                                     
  # Colours map to variables in original data frame before melted.                                                                 
  # Possible: the colour order (left to right) matches with the row order (top to bottom) in original data frame before melted,    but the coloured lined is plotted in the order of levels (left to right) in melted data frame.                                     
  # if (length(cols)<nrow(data)) cols <- diff_cols(nrow(data))                                                                     
  # Custom colours: scale_color_manual(values=cols)                                                                                
  g <- ggplot(data=df.long, aes(x=Samples, y=Value, colour=Genes, group=Genes))+geom_line()+labs(title="", x=x.title,  y=y.title)+theme(legend.position="right", axis.text=element_text(size=text.size), axis.title=element_text(size=text.size, face="bold"), axis.text.x=element_text(angle=text.angle, hjust=1), panel.background = element_rect(fill = "gray95", colour = "gray95", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour ="white"))
  return(g)
}


#' Scale the data as a whole
#'
#' @param dat A data frame or matrix of of numeric data. 

#' @return A scaled data frame.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jianhai.zhang@@email.ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

scale_all <- function(dat) {
  vals <- scale(data.frame(value=unlist(dat)))[, 1]
  dat.s <- as.data.frame(matrix(vals, nrow=nrow(dat), byrow=FALSE))
  dimnames(dat.s) <- dimnames(dat); return(dat.s)
}



