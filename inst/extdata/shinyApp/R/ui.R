# The Shiny modules (e.g. search_ui) are temporarily placed in this file only for debugging purpose, and will be moved to independent files in the R folder after the App development is completed.

library(shiny); library(shinydashboard); library(shinydashboardPlus); library(yaml); library(plotly); library(visNetwork); library(DT); library(shinyWidgets); library(shinyBS); library(shinyjs)
requireNamespace('DESeq2'); requireNamespace('av'); requireNamespace('BiocGenerics'); requireNamespace('distinct')
requireNamespace('dendextend'); requireNamespace('HDF5Array'); requireNamespace('magick'); requireNamespace('DT');
requireNamespace('pROC'); requireNamespace('shinyWidgets'); requireNamespace('shinyjs'); requireNamespace('htmltools');
requireNamespace('shinyBS'); requireNamespace('sortable'); requireNamespace('org.Hs.eg.db'); requireNamespace('org.Mm.eg.db'); requireNamespace('org.At.tair.db'); requireNamespace('org.Dr.eg.db'); requireNamespace('org.Dm.eg.db'); requireNamespace('AnnotationDbi'); requireNamespace('Seurat')


# A module can only have a "ui" element without the "server".

ui <- function(request) {
  dashboardPage( 
    dashboardHeader(title = NULL, titleWidth = 0),
    dashboardSidebar(collapsed = TRUE, disable = TRUE, width = 0, sidebarMenu() ),
    controlbar = dashboardControlbar(id = "right.bar", collapsed = FALSE, overlay = FALSE, width = 51),
    dashboardBody(
     useShinyjs(), 
     includeCSS("www/style.css"),
     tags$head(tags$link(rel="stylesheet", type="text/css", href="style.css")),
   
     includeScript("www/script.js"),
     # tags$head(tags$script(src = "script.js"), tags$script(HTML(js))),
     # tags$script(src="script.js"),
     # tags$head(HTML("<script src='script.js'></script>")),
     # HTML('
     #    <head>
     #    <link rel="stylesheet" type="text/css" href="style.css">
     #    <script src="script.js"></script>
     #    </head>
     #'),

   # Place title on the right of dashboard header.
   tags$script(HTML('
     $(document).ready(function() {
       $("header").find("nav").append(\'<span class="mainTitle" style="color:black;font-weight:bold;font-size:18px">spatialHeatmap Shiny App</span>\');
     })
   ')),
    tags$script(HTML('$(document).tooltip({show: {effect:"none", delay:0}})')),
    column(12, style='border-color:#3c8dbc;border-width:1px;border-style:solid;margin-bottom:5px;margin-top:-10px;text-align:center', textOutput('sumar')),
    fluidRow( 
      do.call(tabsetPanel, append(list(type = "pills", id = 'shm.sup', selected="landing", upload_ui('upl'), shm_ui('shmAll', data_ui('dat'), search_ui('sear')), deg_ui('deg'), scell_ui('scell')),
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
