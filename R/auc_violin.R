#' Plot extracted AUCs by parameter settings 
#'
#' In coclustering optimization, visualize extracted AUCs by each parameter settings in violin plots.

#' @param df.lis The nested list of extracted aucs returned by \code{auc_stat}.
#' @param xlab,ylab The x and y axis labels in the violin plots.
#' @param nrow The numbers of rows of all the violin plots.
#' @param title The title of composite violin plots.
#' @param key.title The title of legend.
#' @param lgd.key.size The size of legend keys.

#' @return An object of ggplot.

#' @examples

#' See function "coclus_opt" by calling "?coclus_opt".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra

#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme element_text unit ggplotGrob 

#' @export auc_violin

auc_violin <- function(df.lis, xlab, ylab='AUC', nrow=3, title=NULL, key.title=NULL, lgd.key.size=0.05) {
  gg.lis <- NULL; df.na.all <- names(df.lis)
  df.nas <- grep('^auc', df.na.all, value=TRUE)
  df.par.na <- names(df.lis[[df.nas[1]]])
  par.nas <- sub('^df\\.', '', df.par.na)
  names(par.nas) <- df.par.na
  for (i in df.nas) {
    df.all <- NULL; for (j in seq_along(par.nas)) {  
    df0 <- df.lis[[i]][[df.par.na[j]]]
    if (is.null(df0)) next
    df0$parameter <- par.nas[j]; df.all <- rbind(df.all, df0)
    }
    gg <- ggplot(df.all, aes(x=parameter, y=auc, fill=parameter)) +
    geom_violin() + geom_boxplot(width=0.1, show.legend=FALSE) +
    labs(title=title, x=xlab, y=ylab, fill=key.title)+theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.key.size = unit(lgd.key.size, 'npc'))
    gg.lis <- c(gg.lis, list(gg))
  }
  gg.grob <- lapply(gg.lis, ggplotGrob)
  grid.arrange(grobs=gg.grob, nrow=nrow)
}
