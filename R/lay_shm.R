#' Organise spatial heatmaps according to user-specified layout
#'

#' @param lay.shm The layout of spatial heatmaps, "gene" or "con", which means organise spatial heatmaps by genes or by conditions.
#' @param con A charater vector of all conditions.
#' @param ncol An integer specifying columns of layout.
#' @param ID.sel A vector of target genes for spatial heatmaps.
#' @param grob.list A list of spatial heatmaps in the form of ggplot2 plot grob, returned by "grob_list()".
#' @param shiny "TRUE" or "FALSE", apply this function to shiny app or not respectively.
#' @inheritParams spatial_hm

#' @return Organised spatial heatmaps according to the layout.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016. \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @importFrom ggplot2 ggplot aes theme element_blank margin element_rect scale_y_continuous scale_x_continuous ggplotGrob geom_polygon scale_fill_manual ggtitle element_text labs 

#' @references
#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr    
#' Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. “Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.” Proc Natl Acad Sci U S A 106 (44): 18843–8

#' @importFrom gridExtra arrangeGrob grid.arrange


lay_shm <- function(lay.shm, con, ncol, ID.sel, grob.list, width, height, shiny) {

    if (lay.shm=="gene") {

        all.cell <- ceiling(length(unique(con))/as.numeric(ncol))*as.numeric(ncol)
        cell.idx <- c(seq_len(length(unique(con))), rep(NA, all.cell-length(unique(con))))
        m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
        lay <- NULL; for (i in seq_len(length(ID.sel))) { lay <- rbind(lay, m+(i-1)*length(unique(con))) }
      if (shiny==TRUE & length(grob.list)>=1) return(grid.arrange(grobs=grob.list, layout_matrix=lay, newpage=TRUE))
       
      g.tr <- lapply(grob.list[seq_len(length(grob.list))], grobTree)
      n.col <- ncol(lay); n.row <- nrow(lay)
      g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc"))

    } else if (lay.shm=="con") {

        grob.all.na <- names(grob.list); grob.all.idx <- NULL
        for (i in unique(con)) {

            grob.all.idx <- c(grob.all.idx, grep(paste0(i, "$"), grob.all.na))

        }

        grob.all.con <- grob.list[grob.all.idx]
        all.cell <- ceiling(length(ID.sel)/as.numeric(ncol))*as.numeric(ncol)
        cell.idx <- c(seq_len(length(ID.sel)), rep(NA, all.cell-length(ID.sel)))
        m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
        lay <- NULL; for (i in seq_len(length(unique(con)))) { lay <- rbind(lay, m+(i-1)*length(ID.sel)) }
 
        if (shiny==TRUE & length(grob.all.con)>=1) return(grid.arrange(grobs=grob.all.con, layout_matrix=lay, newpage=TRUE))
        g.tr <- lapply(grob.all.con, grobTree); g.tr <- g.tr[names(grob.all.con)]
        n.col <- ncol(lay); n.row <- nrow(lay)
        g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc")) 

    }; return(g.arr)

}

