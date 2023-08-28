# The Shiny modules (e.g. search_ui) are temporarily placed in this file only for debugging purpose, and will be moved to independent files in the R folder after the App development is completed.

library(shiny); library(shinydashboard); library(shinydashboardPlus); library(yaml); library(plotly); library(visNetwork); library(DT); library(shinyWidgets); library(shinyBS); library(shinyjs)

# A module can only have a "ui" element without the "server".

ui <- function(request) {
  dashboardPage( 
    dashboardHeader(title = NULL, titleWidth = 0),
    dashboardSidebar(collapsed = TRUE, disable = TRUE, width = 0, sidebarMenu() ),
    # controlbar = dashboardControlbar(id = "right.bar", collapsed = FALSE, overlay = FALSE, width = 51),
    dashboardBody(
     useShinyjs(), 
     includeCSS("www/style.css"),
     #tags$head(tags$link(rel="stylesheet", type="text/css", href="style.css")),
   
     # includeScript("www/script.js"),
     tags$head(tags$script(src = "script.js")),
     # tags$script(src="script.js"),
     # tags$head(HTML("<script src='script.js'></script>")),
     # HTML('
     #    <head>
     #    <link rel="stylesheet" type="text/css" href="style.css">
     #    <script src="script.js"></script>
     #    </head>
     #'),

    tags$script(HTML('$(document).tooltip({show: {effect:"none", delay:0}})')),
    div(id='logo', 'spatialHeatmap',
      tags$button(class="fa-solid fa-backward-step shiny-bound-input action-button", id='backBut', style="color:#337ab7;z-index:5;margin-left:5px;font-size:40px") %>% spsComps::bsTip("Previous tab", placement = "bottom"),
      tags$button(class="fa-solid fa-forward-step shiny-bound-input action-button", id="forBut", style="color:#337ab7;z-index:5;font-size:40px") %>% spsComps::bsTip("Next tab", placement = "bottom")
      ),
    absolutePanel(id = "absInf", style="z-index:5;opacity:1;background-color:white;border-color:#3c8dbc;border-width:2px;border-style:solid;box-shadow: 5px 5px 5px 0px rgba(212,201,212,1);overflow: auto;", class = "absinf", draggable=TRUE, fixed=TRUE, width = 500, height = '15%', top = '2%', left = '80%',
     p(htmlOutput('datInf')), actionButton("btnInf", "More metadata", style=paste0('padding:1px;margin-top:-10px;margin-left:10px;', run.col)), (htmlOutput('covisInf'))
    ),
    fluidRow( 
      do.call(tabsetPanel, append(list(type = "pills", id = 'tabTop', selected="ldg", 
      tabPanel("Welcome", value='ldg',
        navbarPage('', id='navLdg', selected='wel', 
        tabPanel(title="Welcome", value='wel',
        absolutePanel(id = "absLdg", style="z-index:5; opacity:1", class = "absPan", draggable=FALSE, fixed=TRUE, width = 180, height = 75, top = '20%', left = '5%',
        div(style='display: inline-flex', actionButton("btnLdg", "Getting Started!", style=run.col),
        actionButton("man", "User manual", style=run.col))
        ), htmlOutput('ldg')
        ), tabPanel(title="Gallery", value='gallery', htmlOutput('gallery')),
        tabPanel(title=span('Help', style=hp), value='help', htmlOutput('manual'))
        )
      ), upload_ui('upl'), shm_ui('shmAll', data_ui('dat'), search_ui('sear')), deg_ui('deg'), 
      tabPanel(span(id='datMTab', "Data Mining"), value='ana', 
      bsTooltip('datMTab', title = "This tab is activated only after spatial heatmaps are plotted.", placement = "bottom", trigger = "hover"), analysis_ui('ana')
      ), scell_ui('scell')
      ),
      list(tabPanel("About", value='about',
        # includeHTML("instruction/about.html")
        htmlOutput('about')
      ))
    )) # do.call append
    )
 ) # dashboardBody
)

}
