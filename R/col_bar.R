#' Build Color Bar for the Spatial Heatmaps
#'
#' @param geneV The gene expression values used for building the color bar.
#' @param cols The color codes corresponding to geneV.
#' @param width A numeric between 0 and 1, used to reduce the widths of bins in the color bar so as to avoid overlaps.
#' @param bar.title.size A numeric of color bar title size. The default is 0. 
#' @param bar.value.size A numeric of value size in the color bar y-axis. The default is 10.
#' @param mar A 4-length numeric vector specifying margins of the color bar, in the "top", "right", "bottom", "left" order. The unit is "cm". Default is "c(3, 0.1, 3, 0.1)".

#' @return An image of ggplot.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @importFrom ggplot2 ggplot geom_bar theme element_blank margin coord_flip scale_y_continuous scale_x_continuous 

col_bar <- function(geneV, cols, width, bar.title.size=0, bar.value.size=10, mar=c(0.07, 0.01, 0.05, 0.1)) {        

  color_scale <- y <- NULL
  cs.df <- data.frame(color_scale=geneV, y=1)
  cs.g <- ggplot(data=cs.df)+geom_bar(aes(x=color_scale, y=y), fill=cols, orientation='x', stat="identity", width=((max(geneV)-min(geneV))/length(geneV))*width)+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=bar.value.size), plot.margin=margin(t=mar[1], r=mar[2], b=mar[3], l=mar[4], "npc"), panel.grid=element_blank(), panel.background=element_blank(), plot.title= element_text(hjust=0.1, size=bar.title.size))+coord_flip()+labs(title="value", x=NULL, y=NULL)
  if (max(abs(geneV))<=10000 & max(abs(geneV))>=0.0001) cs.g <- cs.g+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0)) else cs.g <- cs.g+scale_x_continuous(expand=c(0,0), labels=function(x) format(x, scientific=TRUE))+scale_y_continuous(expand=c(0,0))+theme(axis.text.y=element_text(angle=45))
  return(cs.g)

}





