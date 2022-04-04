#' Violin plot of extracted AUCs by top spd.sets
#'
#' In coclustering optimization, visualize extracted AUCs by top spd.sets ranked by frequency in violin plots.

#' @param df.lis The nested list of extracted aucs returned by \code{auc_stat}.
#' @param n Number of top \code{spd.set} ranked by frequencies to plot.
#' @param xlab,ylab The x and y axis labels in the violin plots.
#' @param x.agl,x.vjust Angle and vertical position to adjust x-axis text.
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

#' @export spd_auc_violin

spd_auc_violin <- function(df.lis, n=5, ylab='AUC', xlab, x.agl=45, x.vjust=0.6, nrow=3, title=NULL, key.title=NULL, lgd.key.size=0.03) {
  df.nas <- names(df.lis); df.frq.na='df.frq'; df.all.na='df.all'
  gg.lis <- NULL; for (i in df.nas) {
    df.frq0 <- df.lis[[i]][[df.frq.na]]
    df.frq0$spd.set <- paste0('s', df.frq0$sim, 'p', df.frq0$sim.p, 'd', df.frq0$dim)
    df.all0 <- df.lis[[i]][[df.all.na]]
    df.all0$spd.set <- paste0('s', df.all0$sim, 'p', df.all0$sim.p, 'd', df.all0$dim)
    df0 <- subset(df.all0, spd.set %in% df.frq0[seq_len(n), 'spd.set'])
    gg <- ggplot(df0, aes(x=spd.set, y=auc, fill=spd.set)) +
    geom_violin() + geom_boxplot(width=0.1, show.legend=FALSE) +
    labs(title=title, x=xlab, y=ylab, fill=key.title)+theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.key.size = unit(lgd.key.size, 'npc'), axis.text.x=element_text(angle=x.agl, vjust=x.vjust))
    gg.lis <- c(gg.lis, list(gg))
  }
  gg.grob <- lapply(gg.lis, ggplotGrob)
  grid.arrange(grobs=gg.grob, nrow=nrow)
}
