shms_rela_size <- function(input, svg.df, shm, lay.shm, svg.path) {

  observeEvent(list(input$relaSize), {
    relaSize <- input$relaSize; svg.df <- svg.df(); grob.gg.all <- shm$grob.gg.all
    cat('Adjust relative plot size in multiple aSVGs ... \n')
    if (!is.numeric(relaSize) | !is.list(svg.df) | length(svg.df) == 1 | !is.list(grob.gg.all) | is.null(lay.shm()) | is.null(svg.path())) return()
    if (relaSize < 0) return()
    # Get max width/height of multiple SVGs, and dimensions of other SVGs can be set relative to this max width/height.
    w.h.all <- NULL; for (i in seq_along(svg.df)) { w.h.all <- c(w.h.all, svg.df[[i]][['w.h']]); w.h.max <- max(w.h.all) }
    gg.all <- grob.all <- NULL
    for (i in seq_along(grob.gg.all)) {
      gg.lis <- grob.gg.all[[i]]$gg.lis 
      # Also update the central shm$grob.gg.all
      grob.gg.all[[i]]$gg.lis <- gg.lis <- rela_size(svg.df[[i]]$w.h['height'], w.h.max, relaSize, nrow(lay.shm()), gg.lis)
      gg.all <- c(gg.all, gg.lis)
      # Also update the central shm$grob.gg.all
      grob.gg.all[[i]]$grob.lis <- grob.lis <- grob_shm(gg.lis, cores = deter_core(2, svg.path()$svg.path[i]))
      grob.all <- c(grob.all, grob.lis) 
    }; shm$grob.all <- grob.all; shm$gg.all <- gg.all; shm$grob.gg.all <- grob.gg.all
  })
return(shm)
}

