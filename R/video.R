#' Make Spatial Heatmap Video
#'
#' @param gg A list of spatial heatmaps of ggplot.
#' @param cs.g The color key of ggplot.
#' @param sam.uni A vector of unique samples extracted from data matrix.
#' @param ft.trans A vector of tissues to be transparent.
#' @param lgd.key.size The size of legend key (including text). Default is 0.02.
#' @param lgd.text.size The size of legend text. Default is 8.
#' @param lgd.row An integer of legend rows.
#' @param lgd.col An integer of legend columns.
#' @param width The image width in video in "npc", ranging from 0 to 0.92. Default is 0.92.
#' @param height The image height in video in "npc", ranging from 0 to 0.99. Default is 0.99.
#' @inheritParams spatial_hm

#' @param angle.text.key A value of key text angle in legend plot. The default is NULL, equivalent to 0.
#' @param position.text.key The position of key text in legend plot, one of "top", "right", "bottom", "left". Default is NULL, equivalent to "right".
#' @param legend.value.vdo Logical TRUE or FALSE. If TRUE, the numeric values of matching spatial features are added to video legend. The default is NULL.
#' @param sub.title.size The title size of ggplot.
#' @param label Logical. If TRUE, spatial features having matching samples are labeled by feature identifiers. The default is FALSE. It is useful when spatial features are labeled by similar colors. 
#' @param label.size The size of spatial feature labels in legend plot. The default is 4.
#' @param label.angle The angle of spatial feature labels in legend plot. Default is 0.
#' @param hjust The value to horizontally adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param vjust The value to vertically adjust positions of spatial feature labels in legend plot. Default is 0.
#' @param opacity The transparency of colored spatial features in legend plot. Default is 1. If 0, features are totally transparent.
#' @param key Logical. The default is TRUE and keys are added in legend plot. If \code{label} is TRUE, the keys could be removed. 


#' @return A video is saved in \code{out.dir}.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Jeroen Ooms (2020). av: Working with Audio and Video in R. R package version 0.5.0. https://CRAN.R-project.org/package=av
#' Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra

#' @importFrom av av_capture_graphics
#' @importFrom gridExtra grid.arrange

video <- function(gg, cs.g, sam.uni, ft.trans, sub.title.size=NULL, bar.value.size=NULL, lgd.key.size=0.02, lgd.text.size=8, angle.text.key=NULL, position.text.key=NULL, lgd.row=2, lgd.col=NULL, legend.value.vdo=NULL, label=FALSE, label.size=4, label.angle=0, hjust=0, vjust=0, opacity=1, key=TRUE, width=0.92, height=0.99, video.dim='640x480', res=500, interval=1, framerate=1, out.dir) {

  try(test_ffm()); ffm <- tryCatch({ test_ffm() }, error=function(e){ return('error') }, warning=function(w) { return('warning') } )
  if (grepl('error|warning', ffm)) return()

  if (!is.null(bar.value.size)) cs.g <- cs.g+theme(axis.text.y=element_text(size=bar.value.size))
  na <- names(gg)
  cat('Video: adjust legend size/rows... \n')
  gg1 <- gg_lgd(gg.all=gg, size.key=lgd.key.size, size.text.key=lgd.text.size, angle.text.key=angle.text.key, position.text.key=position.text.key, legend.value.vdo=legend.value.vdo, label=label, label.size=label.size, label.angle=label.angle, hjust=hjust, vjust=vjust, opacity=opacity, key=key, sub.title.size=sub.title.size, row=lgd.row, col=lgd.col, sam.dat=sam.uni, ft.trans=ft.trans)
  lay <- rbind(c(NA, NA), c(1, 2), c(NA, NA))
  cat('Saving video... \n')
  res.r=res/144; w.h <- round(as.numeric(strsplit(video.dim, 'x')[[1]])*res.r)
  if (w.h[1] %% 2!=0) w.h[1] <- w.h[1]+1
  if (w.h[2] %% 2!=0) w.h[2] <- w.h[2]+1
  av_capture_graphics(expr=for (i in na) { print(grid.arrange(cs.g, gg1[[i]],widths=unit(c(0.08, width), 'npc'), 
  heights=unit(c(0.05, height, 0.05), 'npc'), layout_matrix=lay)) }, 
  output=paste0(normalizePath(out.dir, winslash="/", mustWork=FALSE), "/shm.mp4"), width=w.h[1], height=w.h[2], res=res, vfilter=paste0('framerate=fps=', framerate))

}

#' # Test if "av" works
#'
#' @keywords Internal
#' @noRd

test_ffm <- function() {
  av_capture_graphics(expr=for (i in seq_len(2)) plot(i), output=paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/tmp.mp4'))
}



