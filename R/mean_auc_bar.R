#' Plot mean of extracted AUCs by parameter settings
#'
#' In coclustering optimization, visualize means of extracted AUCs of each parameter settings by each AUC cutoff in bar plots.

#' @param df.auc The \code{data.frame} of mean AUCs in the nested list of extracted aucs returned by \code{auc_stat}.
#' @param parameter The coloumn name of parameters in \code{df.auc}.
#' @param auc.thr The coloumn name of AUC cutoffs in \code{df.auc}.
#' @param mean The coloumn name of mean AUCs in \code{df.auc}.
#' @param bar.width Width of a single bar.
#' @param title The title of composite violin plots.
#' @param title.size The title size. Default is 20.
#' @param x.lab,ylab The x and y axis label respectively.
#' @param axis.title.size The size of x and y axis labels.
#' @param x.text.size,y.text.size The size of x and y axis text.
#' @param key.title The title of legend.
#' @param lgd.key.size The size of legend keys.
#' @param lgd.text.size The size of legend text.

#' @return An object of ggplot.

#' @examples

#' # To obtain reproducible results, always start a new R session and set a fixed seed for Random Number Generator at the beginning, which is required only once in each R session.  
#' set.seed(10)
#' 
#' # Example bulk data of Arabidopsis thaliana (Arabidopsis) root for coclustering optimization (Li et al 2016).
#' blk <- readRDS(system.file("extdata/cocluster/data", "bulk_cocluster.rds", package="spatialHeatmap"))
#'
#' # Example single cell data of Arabidopsis thaliana (Arabidopsis) root for coclustering optimization (Shahan et al 2020).
#' sc10 <- readRDS(system.file("extdata/cocluster/data", "sc10_cocluster.rds", package="spatialHeatmap"))
#'sc11 <- readRDS(system.file("extdata/cocluster/data", "sc11_cocluster.rds", package="spatialHeatmap"))
#'
#' \donttest{

#' # These example data are already pre-processed. To demonstrate the optimization process the pre-processing steps are perfomed again with few genes or cells removed.
#'
#' # Inital filtering before normalization. 
#' blk <- filter_data(data=blk, pOA=c(0.2, 15), CV=c(1.5, 100)); dim(blk)
#' 
#' fil.init <- filter_cell(lis=list(sc10=sc10, sc11=sc11), bulk=blk, gen.rm='^ATCG|^ATCG', min.cnt=1, p.in.cell=0.3, p.in.gen=0.1)
#' 
#' # Normalization.
#' # sum.factor.
#' norm.fct <- norm_multi(dat.lis=fil.init, cpm=FALSE)
#' # sum.factor + CPM.
#' norm.cpm <- norm_multi(dat.lis=fil.init, cpm=TRUE)
#'
#' # Secondary filtering.
#' # Filtering parameter sets.
#' df.par.fil <- data.frame(p=c(0.1, 0.2, 0.3, 0.4), A=rep(1, 4), cv1=c(0.1, 0.2, 0.3, 0.4), cv2=rep(100, 4), min.cnt=rep(1, 4), p.in.cell=c(0.1, 0.25, 0.3, 0.35), p.in.gen=c(0.01, 0.05, 0.1, 0.15))
#' df.par.fil
#'
#' # Filtered results are saved in "opt_res".
#' if (!dir.exists('opt_res')) dir.create('opt_res')
#' fct.fil.all <- filter_iter(bulk=norm.fct$bulk, cell.lis=list(sc10=norm.fct$sc10, sc11=norm.fct$sc11), df.par.fil=df.par.fil, gen.rm='^ATCG|^ATCG', wk.dir='opt_res', norm.meth='fct')
#'
#' cpm.fil.all <- filter_iter(bulk=norm.cpm$bulk, cell.lis=list(sc10=norm.cpm$sc10, sc11=norm.cpm$sc11), df.par.fil=df.par.fil, gen.rm='^ATCG|^ATCG', wk.dir='opt_res', norm.meth='cpm')
#'
#' # Matching table between bulk tissues and single cells.
#' match.pa <- system.file("extdata/cocluster/data", "match_arab_root_cocluster.txt", package="spatialHeatmap")
#' df.match.arab <- read.table(match.pa, header=TRUE, row.names=1, sep='\t')
#' df.match.arab[1:3, ]
#'
#' # Optimization. 
#' # Check parallelization guide.
#' coclus_opt(wk.dir='opt_res', parallel.info=TRUE, dimred=c('PCA', 'UMAP'), graph.meth=c('knn', 'snn'), sim=seq(0.2, 0.4, by=0.1), sim.p=seq(0.2, 0.4, by=0.1), dim=seq(5, 7, by=1))
#' 
#' # The first-level parallel computing relies on the slurm scheduler (https://slurm.schedmd.com/documentation.html), so if it is available the whole optimization process could be parallelized at two levels. 
# Copy slurm template to current directory. Edit the template according to the parallelization guide and available computing resources.
#' file.copy(system.file("extdata/cocluster", "slurm.tmpl", package="spatialHeatmap"), './slurm.tmpl')
#' 
#' # The first- and second-level parallelizations are set 3 and 2 respectively.
#' library(BiocParallel)
#' opt <- coclus_opt(wk.dir='opt_res', dimred=c('PCA', 'UMAP'), graph.meth=c('knn', 'snn'), sim=seq(0.2, 0.4, by=0.1), sim.p=seq(0.2, 0.4, by=0.1), dim=seq(5, 7, by=1), df.match=df.match.arab, batch.par=BatchtoolsParam(workers=3, cluster="slurm", template='slurm.tmpl'), multi.core.par=MulticoreParam(workers=2))
#'
#' # If slurm is not available, parallelize the optimization only at the second-level through 2 workers. 
#' opt <- coclus_opt(wk.dir='opt_res', dimred=c('PCA', 'UMAP'), graph.meth=c('knn', 'snn'), sim=seq(0.2, 0.4, by=0.1), sim.p=seq(0.2, 0.4, by=0.1), dim=seq(5, 7, by=1), df.match=df.match.arab, batch.par=NULL, multi.core.par=MulticoreParam(workers=2))
#'
#' # The performaces of parameter settings are measured by AUC values in ROC curve. The following demonstrates how to visualize the AUCs and select optimal parameter settings.
# If one AUC is larger than 0.5, 0.6, 0.7, 0.8, or 0.9, and has total bulk assignments >= 500, true assignments >= 500, the corresponding parameter settings are extracted. 
#'
#' # Extract AUCs and other parameter settings for filtering parameter sets.
#' df.lis.fil <- auc_stat(wk.dir='opt_res', tar.par='filter', total.min=500, true.min=300, aucs=round(seq(0.5, 0.9, 0.1), 1))
#' df.lis.fil$df.auc.mean[1:3, ]
#' 
#' # Mean AUCs by each filtering settings and AUC cutoff.
#' mean_auc_bar(df.lis.fil[[1]], bar.width=0.07, title='Mean AUCs by filtering settings')
#' }



#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#' @importFrom ggplot2 ggplot aes geom_bar position_dodge2 labs theme element_text unit

#' @export mean_auc_bar

mean_auc_bar <- function(df.auc, parameter='parameter', auc.thr='auc.thr', mean='mean', bar.width=0.8, title=NULL, title.size=20, xlab="AUC threshold", ylab='Mean AUC', axis.title.size=11, x.text.size=11, y.text.size=11, key.title=NULL, lgd.key.size=0.05, lgd.text.size=10) { 
  cna <- colnames(df.auc)
  cna[cna==parameter] <- 'parameter'
  cna[cna==auc.thr] <- 'auc.thr'; cna[cna==mean] <- 'mean'
  colnames(df.auc) <- cna
  gg <- ggplot(df.auc, aes(x=auc.thr, y=mean, fill=parameter)) +
  geom_bar(color='black', width=bar.width, position=position_dodge2(width=bar.width, preserve="single"), stat="identity") +
  labs(title=title, x=xlab, y=ylab, fill=key.title)+theme(plot.title = element_text(size=title.size, hjust = 0.5, vjust=1.5)) +
  theme(legend.key.size = unit(lgd.key.size, 'npc'), legend.text=element_text(size=lgd.text.size), axis.title = element_text(size=axis.title.size), axis.text.x = element_text(size = x.text.size), axis.text.y = element_text(size = y.text.size))
  gg
}
