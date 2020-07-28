#' Make Spatial Heatmap Video
#'
#' @param gg A list of spatial heatmaps of ggplot.
#' @param cs.g The color key of ggplot.
#' @param sam.uni A vector of unique samples extracted from data matrix.
#' @param tis.trans A vector of tissues to be transparent.
#' @param lgd.key.size The size of legend key (including text). Default is 0.02.
#' @param lgd.text.size The size of legend text. Default is 8.
#' @param lgd.row An integer of legend rows.
#' @inheritParams spatial_hm
#' @return A video is saved in \code{out.dir}.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Yihui Xie (2013). animation: An R Package for Creating Animations and Demonstrating Statistical Methods. Journal of Statistical Software, 53(1), 1-27. URL http://www.jstatsoft.org/v53/i01/.

#' @importFrom animation saveVideo ani.options

video <- function(gg, cs.g, sam.uni, tis.trans, lgd.key.size=0.02, lgd.text.size=8, lgd.row=2, width=1, height=1, video.dim='640x480', interval=1, out.dir) {

  ffm <- tryCatch({ system('which ffmpeg', intern=TRUE) }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
  if (!grepl('ffmpeg', ffm)) return()

    na <- names(gg)
    cat('Video: adjust legend size/rows... \n')
    gg1 <- gg_lgd(gg.all=gg, size.key=lgd.key.size, size.text=lgd.text.size, row=lgd.row, sam.dat=sam.uni, tis.trans=tis.trans)
    lay <- rbind(c(NA, NA), c(1, 2), c(NA, NA))
    cat('Saving video... \n')
    saveVideo(
    for (i in na) { print(grid.arrange(cs.g, gg1[[i]], widths=unit(c(0.08, width*0.92), 'npc'), heights=unit(c(0.05, height*0.99, 0.05), 'npc'), layout_matrix=lay)) 
    ani.options(interval=interval)
    }, video.name=paste0(normalizePath(out.dir), "/shm.mp4"), other.opts=paste0('-pix_fmt yuv420p -b 300k -s:v ', video.dim), img.name='shm', ani.res=1000, verbose=FALSE)
    
}


