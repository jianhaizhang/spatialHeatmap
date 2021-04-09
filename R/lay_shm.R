#' Organise Spatial Heatmaps by Gene or Condition
#'
#' @param con A charater vector of all conditions.
#' @param ncol An integer specifying number of columns in the layout.
#' @param ID.sel A vector of target genes for spatial heatmaps.
#' @param grob.list A list of spatial heatmaps in the form of ggplot2 plot grob, returned by \code{\link{grob_list}}.
#' @param lay.mat Logical. If true, only the layout matrix is returned.
#' @param shiny Logical, TRUE or FALSE, apply this function to shiny app or not respectively.
#' @inheritParams spatial_hm

#' @return Organised spatial heatmaps according to the provided layout.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr    
#' Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @importFrom ggplot2 ggplot aes theme element_blank margin element_rect scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs 
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom grid grobTree unit

lay_shm <- function(lay.shm, con, ncol, ID.sel, grob.list = NULL, lay.mat = FALSE, shiny = FALSE) {
  
  # save(lay.shm, con, ncol, ID.sel, grob.list, shiny, file = 'all.lay')

  ncol <- as.numeric(ncol); con <- unique(con); grob.all.na <- names(grob.list)
  if (lay.shm=="gene"|lay.shm=="none") {

    all.cell <- ceiling(length(con)/ncol)*ncol
    cell.idx <- c(seq_len(length(con)), rep(NA, all.cell-length(con)))
    m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
    lay <- NULL; for (i in seq_len(length(ID.sel))) { lay <- rbind(lay, m+(i-1)*length(con)) }
    if (lay.mat == TRUE) return(lay)
    # Sort conditions under each gene.
    #na.sort <- sort_gen_con(ID.sel=ID.sel, na.all=grob.all.na, con.all=con, by=lay.shm)
    # grob.list <- grob.list[na.sort]
    if (shiny==TRUE & length(grob.list)>=1) return(list(shm=grid.arrange(grobs=grob.list, layout_matrix=lay, newpage=TRUE), lay=lay))
       
    g.tr <- lapply(grob.list[seq_len(length(grob.list))], grobTree)
    n.col <- ncol(lay); n.row <- nrow(lay)
    g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(0.99/n.col, n.col), "npc"), heights=unit(rep(0.99/n.row, n.row), "npc"))

  } else if (lay.shm=="con") {

    all.cell <- ceiling(length(ID.sel)/ncol)*ncol
    cell.idx <- c(seq_len(length(ID.sel)), rep(NA, all.cell-length(ID.sel)))
    m <- matrix(cell.idx, ncol=ncol, byrow=TRUE)
    lay <- NULL; for (i in seq_len(length(con))) { lay <- rbind(lay, m+(i-1)*length(ID.sel)) }
    if (lay.mat == TRUE) return(lay)
    # na.sort <- sort_gen_con(ID.sel=ID.sel, na.all=grob.all.na, con.all=con, by='con')
    # grob.list <- grob.list[na.sort]
    if (shiny==TRUE & length(grob.list)>=1) return(list(shm=grid.arrange(grobs=grob.list, layout_matrix=lay, newpage=TRUE), lay=lay))
    g.tr <- lapply(grob.list, grobTree); g.tr <- g.tr[names(grob.list)]
    n.col <- ncol(lay); n.row <- nrow(lay)
    g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(0.99/n.col, n.col), "npc"), heights=unit(rep(0.99/n.row, n.row), "npc")) 

  }; return(g.arr)

}
