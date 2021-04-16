# Access html files of targt genes in the animation.
gly_url <- function(ipt, shm, gID, pat.all) {
  # As long as the variable of 'reactive' is used in the 'ui.R', changes of elements in 'reactive' would cause chain change all the way to 'ui.R'. E.g. the change in "input$ggly.but=='No'" leads to changes in 'output$ggly' and 'ui.R', not necessarily changes in 'output$ggly' call changes in 'gly.url'.
  gly.url <- reactive({
    input <- ipt
    if (is.null(input$ggly.but)) return() 
    if (is.null(shm$gg.all1)|input$ggly.but=='No'|gID$geneSel[1]=='none'|is.null(pat.all())) return(NULL)
    gg.all <- shm$gg.all1; na <- names(gg.all)
    # Only take the selected genes.
    na <- na[grepl(paste0('^', pat.all(), '_\\d+$'), na)]; na1 <- na[as.integer(input$fm)]
    na2 <- list.files('www/ggly', pattern=na1)
    if (length(na2)==0|is.na(na2)) return(NULL)
    cat('Animation: access', na2, 'path \n')
    paste0('ggly/', na2) 
  }); return(gly.url)
}
