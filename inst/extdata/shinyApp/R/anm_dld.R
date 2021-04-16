# Generate html files for animation.
anm_dld <- function(ipt, gly.url, color, shm, svg.df, geneIn, pat.all, shm.bar, sam) {
  anm.dld <- reactive({
    input <- ipt
    if (is.null(input$ggly.but)) return()
    if (input$ggly.but=='No'|is.null(gly.url())) return()
    if (is.null(svg.df())|is.null(geneIn())|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL) 
    withProgress(message="Downloading animation: ", value=0, {
    incProgress(0.1, detail="in progress ...")
    gg.all <- shm$gg.all1; na <- names(gg.all)
    gg.na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]
    gg <- gg.all[gg.na]
    pro <- 0.1; for (i in seq_along(gg.na)) {
    incProgress(pro+0.2, detail=paste0('preparing ', gg.na[i], '.html...'))
    html_ly(gg=gg[i], cs.g=shm.bar(), ft.trans=input$tis, sam.uni=sam(), anm.width=input$width.ly, anm.height=input$height.ly, out.dir='.') }
   })
  }); return(anm.dld)
}

anm_dld_handler <- function(anm.dld, output, out.na) {
 # This step leaves 'fil.na' in 'output$dld.anm' being a global variable.
  output[[out.na]] <- downloadHandler( 
    # The rest code will run only after 'anm.dld()' is done.
    filename=function(){ anm.dld(); "html_shm.zip" },
    fil.na <- paste0(normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE), '/html_shm.zip'),
    content=function(fil.na){ cat('Downloading animation... \n'); zip(fil.na, 'html_shm/') }
  ); return(output)
}
