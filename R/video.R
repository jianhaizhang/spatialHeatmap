#' Make Spatial Heatmap Video
#'
#' @param gg A list of spatial heatmaps of ggplot.
#' @param cs.g The color key of ggplot.
#' @param sam.uni A vector of unique samples extracted from data matrix.
#' @param tis.trans A vector of tissues to be transparent.
#' @param lgd.key.size The size of legend key (including text). Default is 0.02.
#' @param lgd.text.size The size of legend text. Default is 8.
#' @param lgd.row An integer of legend rows.
#' @param width The image width in video in "npc", ranging from 0 to 0.92. Default is 0.92.
#' @param height The image height in video in "npc", ranging from 0 to 0.99. Default is 0.99.
#' @inheritParams spatial_hm
#' @inheritParams gg_lgd
#' @inheritParams col_bar
#' @return A video is saved in \code{out.dir}.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Jeroen Ooms (2020). av: Working with Audio and Video in R. R package version 0.5.0. https://CRAN.R-project.org/package=av
#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra

#' @importFrom av av_capture_graphics
#' @importFrom gridExtra grid.arrange

video <- function(gg, cs.g, sam.uni, tis.trans, sub.title.size=NULL, bar.value.size=NULL, lgd.key.size=0.02, lgd.text.size=8, angle.text.key=NULL, position.text.key=NULL, lgd.row=2, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, width=0.92, height=0.99, video.dim='640x480', res=500, interval=1, framerate=1, out.dir) {

  # Test if "av" works.
  test <- function() {
    av_capture_graphics(expr=for (i in seq_along(1:2)) plot(i), output=paste0(tempdir(check=TRUE), '/tmp.mp4'))
  }; try(test())
  ffm <- tryCatch({ test() }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
  if (grepl('error|warning', ffm)) return()

  if (!is.null(bar.value.size)) cs.g <- cs.g+theme(axis.text.y=element_text(size=bar.value.size))
  na <- names(gg)
  cat('Video: adjust legend size/rows... \n')
  gg1 <- gg_lgd(gg.all=gg, size.key=lgd.key.size, size.text.key=lgd.text.size, angle.text.key=angle.text.key, position.text.key=position.text.key, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, sub.title.size=sub.title.size, row=lgd.row, sam.dat=sam.uni, tis.trans=tis.trans)
  lay <- rbind(c(NA, NA), c(1, 2), c(NA, NA))
  cat('Saving video... \n')
  res.r=res/144; w.h <- round(as.numeric(strsplit(video.dim, 'x')[[1]])*res.r)
  if (w.h[1] %% 2!=0) w.h[1] <- w.h[1]+1
  if (w.h[2] %% 2!=0) w.h[2] <- w.h[2]+1
  av_capture_graphics(expr=for (i in na) { print(grid.arrange(cs.g, gg1[[i]],widths=unit(c(0.08, width), 'npc'), 
  heights=unit(c(0.05, height, 0.05), 'npc'), layout_matrix=lay)) }, 
  output=paste0(normalizePath(out.dir), "/shm.mp4"), width=w.h[1], height=w.h[2], res=res, vfilter=paste0('framerate=fps=', framerate))

}



