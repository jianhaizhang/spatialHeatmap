#' Plot mean of extracted AUCs by parameter settings
#'
#' In coclustering optimization, visualize means of extracted AUCs of each parameter settings by each AUC cutoff in bar plots.

#' @param df.auc The \code{data.frame} of mean AUCs in the nested list of extracted aucs returned by \code{auc_stat}.
#' @param parameter The coloumn name of parameters in \code{df.auc}.
#' @param auc.thr The coloumn name of AUC cutoffs in \code{df.auc}.
#' @param mean The coloumn name of mean AUCs in \code{df.auc}.
#' @param bar.width Width of a single bar.
#' @param title The title of composite violin plots.
#' @param key.title The title of legend.
#' @param lgd.key.size The size of legend keys.

#' @return An object of ggplot.

#' @examples

#' See function "coclus_opt" by calling "?coclus_opt".

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @importFrom ggplot2 ggplot aes geom_bar position_dodge2 labs theme element_text unit

#' @export mean_auc_bar

mean_auc_bar <- function(df.auc, parameter='parameter', auc.thr='auc.thr', mean='mean', bar.width=0.8, title=NULL, key.title=NULL, lgd.key.size=0.05) { 
  cna <- colnames(df.auc)
  cna[cna==parameter] <- 'parameter'
  cna[cna==auc.thr] <- 'auc.thr'; cna[cna==mean] <- 'mean'
  colnames(df.auc) <- cna
  gg <- ggplot(df.auc, aes(x=auc.thr, y=mean, fill=parameter)) +
  geom_bar(color='black', width=bar.width, position=position_dodge2(width=bar.width, preserve="single"), stat="identity") +
  labs(title=title, x="AUC threshold", y='Mean AUC', fill=key.title)+theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.key.size = unit(lgd.key.size, 'npc'))
  gg
}
