# The Shiny modules (e.g. search_ui) are temporarily placed in this file only for debugging purpose, and will be moved to independent files in the R folder after the App development is completed.

library(shiny); library(shinydashboard); library(shinydashboardPlus); library(yaml); library(plotly); library(visNetwork); library(DT); library(shinyWidgets); library(shinyBS); library(shinyjs)
requireNamespace('DESeq2'); requireNamespace('av'); requireNamespace('BiocGenerics'); requireNamespace('distinct')
requireNamespace('dendextend'); requireNamespace('HDF5Array'); requireNamespace('magick'); requireNamespace('DT');
requireNamespace('pROC'); requireNamespace('shinyWidgets'); requireNamespace('shinyjs'); requireNamespace('htmltools');
requireNamespace('shinyBS'); requireNamespace('sortable')


js <- "function openFullscreen(elem) {
  if (elem.requestFullscreen) {
    elem.requestFullscreen();
  } else if (elem.mozRequestFullScreen) { /* Firefox */
    elem.mozRequestFullScreen();
  } else if (elem.webkitRequestFullscreen) { /* Chrome, Safari and Opera */
    elem.webkitRequestFullscreen();
  } else if (elem.msRequestFullscreen) { /* IE/Edge */
    elem.msRequestFullscreen();
  }
}"

# A module can only have a "ui" element without the "server".

ui <- function(request) {
  dashboardPage( 
    # includeCSS("style.css"),
    dashboardHeader(title = NULL, titleWidth = 0),
    dashboardSidebar(collapsed = TRUE, disable = TRUE, width = 0, sidebarMenu() ),
    controlbar = dashboardControlbar(id = "right.bar", collapsed = FALSE, overlay = FALSE, width = 51),
    dashboardBody(
   # tags$head(HTML('<title>spatialHeatmap</title>')),
     useShinyjs(), 
     includeCSS("www/style.css"),
     tags$head(tags$link(rel="stylesheet", type="text/css", href="style.css")),
   
     tags$head(tags$script(src = "javascript.js"), tags$script(HTML(js))),
     includeScript(path = "www/javascript.js"),
     tags$script(src="javascript.js"),
     tags$head(HTML("<script type='text/javascript' src='javascript.js'></script>")),
     HTML('<head>
         <link rel="stylesheet" type="text/css" href="style.css">
         <script type="text/javascript" src="www/javascript.js"></script>
         </head>
     '),

   # Place title on the right of dashboard header.
   tags$script(HTML('
     $(document).ready(function() {
       $("header").find("nav").append(\'<img src="logo.png" width="175">\');
     })
   ')),
    fluidRow( 
      do.call(tabsetPanel, append(list(type = "pills", id = 'shm.sup', selected="landing", upload_ui('upl'), shm_ui('shmAll', data_ui('dat', dim_ui('datDim')), search_ui('sear')), deg_ui('deg'), scell_ui('scell')),
        list(tabPanel("About", value='about',
          if (0) box(width = 12, title = "", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
          ),
        includeHTML("instruction/about.html")
        ))
    )) # do.call append

    )
 ) # dashboardBody
)

}
