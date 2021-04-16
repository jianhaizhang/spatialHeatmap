# Display SHMs by frames (animation).
anm_ui <- function(id) { htmlOutput(NS(id, "ggly")) }

anm_server <- function(id, ipt, color, arg.lis) {
  moduleServer(id, function(input, output, session) {

  gly.url <- arg.lis$gly.url; svg.df <- arg.lis$svg.df
  geneIn <- arg.lis$geneIn; input <- ipt
  observeEvent(list(log=input$log, tis=input$tis, col.but=input$col.but, cs.v=input$cs.v, preScale=input$preScale, ggly.but=input$ggly.but, fm=input$fm), {
  
  output$ggly <- renderUI({ 

    if (input$ggly.but=='No'|is.null(gly.url)) return()
    if (is.null(svg.df)|is.null(geneIn)|is.null(input$dt_rows_selected)|color$col[1]=='none') return(NULL)
    withProgress(message="Animation: ", value=0, {
    incProgress(0.75, detail="plotting ...")
    gly.url <- gly.url; cat('Animation: plotting', gly.url, '\n')
    tags$iframe(src=gly.url, height=input$height.ly, width=input$width.ly, scrolling='auto')
  
    })

  })

  })

  })
}

