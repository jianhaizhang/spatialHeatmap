#' Plot single-cell quality control statistics
#'
#' @param sce.unfil Unfiltered \code{SingleCellExperiment} object with quality control statistics stored in \code{colData}.
#' @param ercc Logical, indicating whether ERCC spike-in controls are available. The default is \code{FALSE}.

#' @return A composite ggplot.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). "Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R." _Bioinformatics_, *33*, 1179-1186. doi: 10.1093/bioinformatics/btw777 (URL: https://doi.org/10.1093/bioinformatics/btw777).
#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra.

#' @importFrom ggplot2 ggplot theme_void
#' @importFrom scater plotColData
#' @importFrom gridExtra grid.arrange

sc_qc_plot <- function(sce.unfil, ercc = FALSE) {
  g.sum <- plotColData(sce.unfil, y="sum", colour_by="discard") + scale_y_log10() + ggtitle("Total count") 
  g.det <- plotColData(sce.unfil, y="detected", colour_by="discard") + scale_y_log10() + ggtitle("Detected features")              
  g.mt.per <- plotColData(sce.unfil, y="subsets_Mt_percent", colour_by="discard") + ggtitle("Mito percent") 
  if (ercc) g.ercc.per <- plotColData(sce.unfil, y="altexps_ERCC_percent", colour_by="discard") + ggtitle("ERCC percent") else g.ercc.per <- ggplot() + theme_void()
  cat('Single cell: mitochondrial counts vs total, spike-in ... \n')
  g.mt.log <- plotColData(sce.unfil, x="sum", y="subsets_Mt_percent", colour_by="discard") + scale_x_log10()
  if (ercc) g.mt.ercc <- plotColData(sce.unfil, x="altexps_ERCC_percent", y="subsets_Mt_percent", colour_by="discard") else g.mt.ercc <- ggplot() + theme_void()
  if (ercc) grid.arrange(
    g.sum, g.det, g.mt.per, g.ercc.per, g.mt.log, g.mt.ercc,                                                                       
    ncol = 2
  ) else grid.arrange(g.sum, g.det, g.mt.per, g.mt.log, ncol = 2)                                                                  
}                                                      
     
