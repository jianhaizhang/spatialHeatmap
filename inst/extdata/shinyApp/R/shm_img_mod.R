
shm_img_ui <- function(id) {
  plotOutput(NS(id, "shm"), height='100%', width = '100%')
}


shm_img_server <- function(id, ipt, shm, col.reorder, gID, color, arg.lis) {
  moduleServer(id, function(input, output, session) {
  observe({

  lay <- arg.lis$lay.shm; svg.df <- arg.lis$svg.df; geneIn <- arg.lis$geneIn; pat.all <- arg.lis$pat.all; shm.bar <- arg.lis$shm.bar; sam <- arg.lis$sam
  scale.shm <- ipt$scale.shm

    if (is.null(lay)|!is.numeric(scale.shm)) return()
    if (scale.shm <= 0) return()
    # subplot: height 300, width 250 
    # Avoid: if one column has all NAs in the layout matrix, the aspect ratio is distroyed. So only take the columns not containing all NAs.
    col.vld <- sum(unlist(lapply(seq_len(ncol(lay)), function(x) !all(is.na(lay[, x])))))
    # width/height relate to scrolling in box. 
    output$shm <- renderPlot(width = col.vld * 250 * scale.shm, height = nrow(lay) * 300 * scale.shm, { 
    if (col.reorder$col.re=='N') return()
    if.con <-  is.null(ipt$dt_rows_selected)|is.null(svg.df)|gID$geneSel[1]=="none"|is.null(shm$grob.all1)
    if (length(if.con==FALSE)==0) if (length(if.con)==0) return(); if (is.na(if.con)|if.con==TRUE) return(NULL)
    if (is.na(color$col[1])|length(color$col=="none")==0|ipt$color=="") return(NULL)
    r.na <- rownames(geneIn[["df.aggr.tran"]])
    grob.na <- names(shm$grob.all1)
    # Select target grobs.
    # Use definite patterns and avoid using '.*' as much as possible. Try to as specific as possible.
    pat.all <- paste0('^', pat.all, '(_\\d+$)')
    grob.lis.p <- shm$grob.all1[grepl(pat.all, grob.na)] # grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    # Indexed cons with '_1', '_2', ... at the end.
    con <- unique(gsub(pat.all, '\\2\\3', names(grob.lis.p))); if (length(con)==0) return()
    cat('Plotting spatial heatmaps... \n')
    lay <- ipt$gen.con; ID <- gID$geneSel; ncol <- ipt$col.n
    # This step is plotting.
    shm.lay <- lay_shm(lay.shm=lay, con=con, ncol=ncol, ID.sel=ID, grob.list=grob.lis.p, shiny=TRUE); shm.arr <- shm.lay$shm
    # Adjust the dimension in chicken example. 
    #svg.df <- svg.df(); w.all <- h.all <- 0
    # for (i in svg.df) { w.all <- w.all+i$w.h['width']; h.all <- h.all+i$w.h['height'] }
    #shm.row <- nrow(shm.lay$lay)
    #if (length(svg.df)==1) updateNumericInput(session, inputId="height", label="Overall height:", value=as.numeric(h.all/w.all*shm.row*input$width), min=0.1, max=Inf, step=NA)
    if (ipt$ext!='NA') {
      validate(need(try(ipt$res>0), 'Resolution should be a positive numeric!'))
      validate(need(try(ipt$lgd.w>=0 & ipt$lgd.w <1), 'Legend width should be between 0 to 1!'))
      validate(need(try(ipt$lgd.ratio>0), 'Legend aspect ratio should be a positive numeric!'))
      cs.grob <- ggplotGrob(shm.bar)
      cs.arr <- arrangeGrob(grobs=list(grobTree(cs.grob)), layout_matrix=cbind(1), widths=unit(1, "npc"))
      # Legend size in downloaded SHM is reduced.
      lgd.lis <- shm$lgd.all; lgd.lis <- gg_lgd(gg.all=lgd.lis, sam.dat=sam, ft.trans=ipt$tis, label=FALSE)
      lgd.lis <- gg_lgd(gg.all=lgd.lis, size.key=ipt$lgd.key.size*0.5, size.text.key=NULL, label.size=ipt$lgd.lab.size, row=ipt$lgd.row, sam.dat=sam, ft.trans=ipt$tis, position.text.key='right', label=(ipt$lgd.label=='Yes'))
      if (ipt$lgd.w>0) {
        grob.lgd.lis <- lapply(lgd.lis, ggplotGrob)
        lgd.tr <- lapply(grob.lgd.lis, grobTree)
    # In 'arrangeGrob', if numbers in 'layout_matrix' are more than items in 'grobs', there is no difference. The width/height of each subplot is decided by 'widths' and 'heights'.
    lgd.arr <- arrangeGrob(grobs=lgd.tr, layout_matrix=matrix(seq_along(lgd.lis), ncol=1), widths=unit(1, "npc"), heights=unit(rep(1/length(lgd.lis)/ipt$lgd.ratio, length(lgd.lis)), "npc"))
    w.lgd <- (1-0.08)/(ncol+1)*ipt$lgd.w # Legend is reduced.
    png(paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/tmp.png')); shm1 <- grid.arrange(cs.arr, shm.arr, lgd.arr, ncol=3, widths=unit(c(0.08-0.005, 1-0.08-w.lgd, w.lgd), 'npc')); dev.off() } else { png(paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/tmp.png')); shm1 <- grid.arrange(cs.arr, shm.arr, ncol=2, widths=unit(c(0.08-0.005, 1-0.08), 'npc')); dev.off() }
    ggsave(paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/shm.', ipt$ext), plot=shm1, device=ipt$ext, width=ipt$width/72, height=ipt$height/72, dpi=ipt$res, unit='in') 
    }
    
    })

  })

  })
}
