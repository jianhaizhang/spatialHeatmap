#' Embedding plots of single cells before and after co-clustering optimization
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

#' @keywords Internal                                                                                                              
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Xavier Robin, Natacha Turck, Alexandre Hainard, Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller     (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics, 12, p. 77.  DOI: 10.1186/ 1471-2105-12-77 <http://www.biomedcentral.com/1471-2105/12/77/>
#' Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. https://www.nature.com/articles/s41592-019-0654-x.
#' SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.  org/package=gridExtra

#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot ggplot aes geom_point theme_classic theme element_text margin element_rect element_blank unit labs scale_colour_manual guide_legend geom_text
#' @importFrom gridExtra grid.arrange
 
dim_opt <- function(input, dimred='TSNE', tit1='Before co-clustering', tit2='After co-clustering', title.size=20, pt.size=1, lgd.pos='bottom', lgd.text.size=13, lgd.nrow=4, lgd.key.size=3, lgd.spa.y=-0.02, lgd.spa.x=-0.005, alp=0.8, axis.font.size=15, ann.size=8, adj.col=0) {
  bulkCell <- x <- y <- NULL
  sce <- input$sce.all; auc <- round(pROC::auc(input$roc.obj), 2)
  sce <- subset(sce, , bulkCell %in% 'cell')
  labxy <- paste0(dimred, seq_len(2))

  # Data frame for plotting.
  dat <- reducedDim(sce, dimred); colnames(dat) <- c('x', 'y')
  dat <- data.frame(cbind(dat, colData(sce)))
  sams <- unique(dat$sample)
  cols <- setNames(diff_cols(length(sams)+adj.col)[seq_along(sams)], sams)
  # All cells.
  br <- sams; lgd.show=sams 
  # Add AUC value.
  min.x <- min(dat$x); min.y <- min(dat$y)
  ann.pos <- c(min.x + abs(min.x)*0.05, min.y + abs(min.y)*0.05)
  # Original dim plot.
  g1 <- ggplot(dat, aes(x=x, y=y)) + geom_point(size=pt.size, shape=19, alpha=alp, aes(colour=sample)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=title.size), legend.position=lgd.pos, legend.text=element_text(size=lgd.text.size), legend.margin = margin(t=0.001, l=0.001, r=0.001, unit='npc'), legend.title=element_text(size=lgd.text.size+1), legend.background =element_rect(fill='transparent'), axis.text = element_blank(), axis.ticks = element_blank(), axis.title=element_text(size=axis.font.size), legend.spacing.y = unit(lgd.spa.y, 'npc'), legend.spacing.x = unit(lgd.spa.x, 'npc')) + labs(title=tit1, x=labxy[1], y=labxy[2]) + scale_colour_manual(name='', values=cols, breaks=br, labels=lgd.show, guide=guide_legend(title='', nrow=lgd.nrow, byrow = TRUE, override.aes = list(size=lgd.key.size))) + scale_y_continuous(expand=expansion(mult=c(0.01, 0.01)))+scale_x_continuous(expand=expansion(mult=c(0.01, 0.01))) 
  # Incorrect/no assignments have color gray80.
  fal <- dat$response %in% c(FALSE, 'none')
  dat$sample[fal] <- paste0(dat$sample[fal], '.false')
  sams.new <- unique(dat$sample)
  idx <- grepl('\\.false$', sams.new)
  col.fal <- setNames(rep('gray80', sum(idx)), sams.new[idx])
  cols <- c(cols, col.fal)
  # Legend to show.
  br <- c(sams.new[!idx], sams.new[idx][1]); lgd.show=c(sams.new[!idx], 'other') 
  # Plot after co-clustering.
  g2 <- ggplot(dat, aes(x=x, y=y)) + geom_point(size=pt.size, shape=19, alpha=alp, aes(colour=sample)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=title.size), legend.position=lgd.pos, legend.text=element_text(size=lgd.text.size), legend.margin = margin(t=0.001, l=0.001, r=0.001, unit='npc'), legend.title=element_text(size=lgd.text.size+1), legend.background =element_rect(fill='transparent'), axis.text = element_blank(), axis.ticks = element_blank(), axis.title=element_text(size=axis.font.size), legend.spacing.y = unit(lgd.spa.y, 'npc'), legend.spacing.x = unit(lgd.spa.x, 'npc')) + labs(title=tit2, x=labxy[1], y=labxy[2]) + scale_colour_manual(name='', values=cols, breaks=br, labels=lgd.show, guide=guide_legend(title='', nrow=lgd.nrow, byrow = TRUE, override.aes = list(size=lgd.key.size)))+ geom_text(x=ann.pos[1], y=ann.pos[2], label=paste0('AUC: ', auc), size=ann.size, hjust=-0) + scale_y_continuous(expand=expansion(mult=c(0.01, 0.01)))+scale_x_continuous(expand=expansion(mult=c(0.01, 0.01)))
  grid.arrange(g1, g2, nrow = 1)
}
