#' Organise Spatial Heatmaps by Gene or Condition
#'
#' @param lay.shm The layout of spatial heatmaps, "gene" or "con", which means organise spatial heatmaps by genes or by conditions.
#' @param con A charater vector of all conditions.
#' @param ncol An integer specifying number of columns in the layout.
#' @param ID.sel A vector of target genes for spatial heatmaps.
#' @param grob.list A list of spatial heatmaps in the form of ggplot2 plot grob, returned by \code{\link{grob_list}}.
#' @param shiny Logical, TRUE or FALSE, apply this function to shiny app or not respectively.
#' @inheritParams spatial_hm

#' @return Organised spatial heatmaps according to the provided layout.
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

  # Sort vector of letter and number mixture.
  sort_mix <- function(vec) {

    w <- is.na(as.numeric(gsub('\\D', '', vec)))
    let <- vec[w]; num <- vec[!w]
    let.num <- as.numeric(gsub('\\D', '', vec))
    vec[!w] <- num[order(let.num[!w])]
    vec[w] <- let[order(let.num[w])]; return(vec)
        
  }        

  width <- as.numeric(width); height <- as.numeric(height); ncol <- as.numeric(ncol); con <- unique(con)
  grob.all.na <- names(grob.list)
  if (lay.shm=="gene") {

    all.cell <- ceiling(length(con)/ncol)*ncol
    cell.idx <- c(seq_len(length(con)), rep(NA, all.cell-length(con)))
    m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
    lay <- NULL; for (i in seq_len(length(ID.sel))) { lay <- rbind(lay, m+(i-1)*length(con)) }
    # Sort conditions under each gene.
    con.pat1 <- paste0('.*_(', paste0(con, collapse='|'), ')$')
    na.sort <- NULL; for (i in ID.sel) {
        
      na0 <- grob.all.na[grepl(paste0('^', i, '_'), grob.all.na)]
      con1 <- gsub(con.pat1, '\\1', na0)
      na.sort <- c(na.sort, paste0(i, '_', sort_mix(con1)))

    }; grob.list <- grob.list[na.sort]
    if (shiny==TRUE & length(grob.list)>=1) return(grid.arrange(grobs=grob.list, layout_matrix=lay, newpage=TRUE))
       
    g.tr <- lapply(grob.list[seq_len(length(grob.list))], grobTree)
    n.col <- ncol(lay); n.row <- nrow(lay)
    g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc"))

  } else if (lay.shm=="con") {

    all.cell <- ceiling(length(ID.sel)/ncol)*ncol
    cell.idx <- c(seq_len(length(ID.sel)), rep(NA, all.cell-length(ID.sel)))
    m <- matrix(cell.idx, ncol=ncol, byrow=TRUE)
    lay <- NULL; for (i in seq_len(length(con))) { lay <- rbind(lay, m+(i-1)*length(ID.sel)) }
    # Sort conditions and genes.
    gen.pat1 <- paste0('^(', paste0(ID.sel, collapse='|'), ')_.*')
    na.sort <- NULL; for (i in sort_mix(con)) {
      
      na0 <- grob.all.na[grepl(paste0('_', i, '$'), grob.all.na)]
      gen1 <- gsub(gen.pat1, '\\1', na0)
      na.sort <- c(na.sort, paste0(sort_mix(gen1), '_', i))

    }; grob.list <- grob.list[na.sort]
    if (shiny==TRUE & length(grob.list)>=1) return(grid.arrange(grobs=grob.list, layout_matrix=lay, newpage=TRUE))
    g.tr <- lapply(grob.list, grobTree); g.tr <- g.tr[names(grob.list)]
    n.col <- ncol(lay); n.row <- nrow(lay)
    g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc")) 

  }; return(g.arr)

}

