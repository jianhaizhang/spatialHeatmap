# Change animation dimension according to user input.
anm_dim_ui <- function(id, tran.t = FALSE, anm.w = FALSE, anm.h = FALSE, dld.anm.but = FALSE) {
  tagList(
    if (tran.t == TRUE) htmlOutput(NS(id, "tran.t")),
    if (anm.w == TRUE) htmlOutput(NS(id, "anm.w")),
    if (anm.h == TRUE) htmlOutput(NS(id, "anm.h")),
    if (dld.anm.but == TRUE) htmlOutput(NS(id, "dld.anm.but"))
  )
}

anm_dim_server <- function(id, ipt, color, shm, cfg, arg.lis) {
  moduleServer(id, function(input, output, session) {
  input <- ipt; gly.url <- arg.lis$gly.url; svg.df <- arg.lis$svg.df; geneIn <- arg.lis$geneIn
  pat.all <- arg.lis$pat.all; # anm.dld <- arg.lis$anm.dld
  # Variables in 'observe' are accessible anywhere in the same 'observe'.
  observe({
  
    if (is.null(input$ggly.but)|is.null(input$fm)) return() 
    if (input$ggly.but=='No'|is.null(gly.url)) return()
    if (is.null(svg.df)|is.null(geneIn)|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
    
    gg.all <- shm$gg.all1; na <- names(gg.all)
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all, '_\\d+$'), na)]; na1 <- na[as.integer(input$fm)]
    na2 <- list.files('www/ggly', pattern=paste0(na1, '\\.html$')); if (length(na2)==0) return(NULL)
    gg <- gg.all[[na1]]
    dat <- layer_data(gg); x.max <- max(dat$x); y.max <- max(dat$y)
    w <- cfg$lis.par$shm.anm['width', 'default']
    if (w!='NA') w <- as.numeric(w) else w <- NA
    h <- cfg$lis.par$shm.anm['height', 'default']
    if (h!='NA') h <- as.numeric(h) else h <- NA
    if (!is.na(w)) {
      h <- y.max/x.max*w; if (h>550) { h <- 550; w <- x.max/y.max*h }
    } else if (is.na(w) & !is.na(h)) { w <- x.max/y.max*h }
   output$tran.t <- renderUI({
      numericInput(inputId='t', label='Transition time (s)', value=as.numeric(cfg$lis.par$shm.anm['transition', 'default']), min=0.1, max=Inf, step=NA, width=270)
    }) 
   output$anm.w <- renderUI({
      numericInput(inputId='width.ly', label='Width', value=w, min=1, max=Inf, step=NA, width=170)
    })
   output$anm.h <- renderUI({
      numericInput(inputId='height.ly', label='Height', value=h, min=1, max=Inf, step=NA, width=170)
    })
  output$dld.anm.but <- renderUI({ downloadButton("dld.anm", "Download") })
  })

})

}
