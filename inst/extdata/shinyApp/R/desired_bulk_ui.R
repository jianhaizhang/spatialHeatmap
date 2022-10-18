# The Shiny modules (e.g. search_ui) are temporarily placed in this file only for debugging purpose, and will be moved to independent files in the R folder after the App development is completed.

library(shiny); library(shinydashboard); library(shinydashboardPlus); library(DT); library(shinyWidgets); library(shinyBS); library(shinyjs)

desired_bulk_ui <- function(request) {
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
       $("header").find("nav").append(\'<span class="mainTitle">Assigning Desired Bulk for Selected Cells</span>\');
     })
   ')),

  tabsetPanel(type = "pills", id='asgBulk', selected='upl',
    tabPanel(title="Upload data", value='upl', br(),
    fileInput("sglCell", "Upload co-clustering results (.rds)", accept=c(".rds"), multiple=FALSE)
    ),
    tabPanel(title="Assign bulk", value='asg',
      column(6, id='dimCellBut', 
      fluidRow(splitLayout(cellWidths=c('1%', '50%', '1%', '20%', '1%'), '',
      selectInput('dimCell', label='Embedding plot', choices=c("UMAP", "PCA", "TSNE")), '', 
      actionButton("tailorHelp", "Help", style='margin-bottom:-60px', icon = icon('question-circle'))
      ))
      ),

      column(6, id='selBlkButCan',
      fluidRow(splitLayout(cellWidths=c('1%', '40%', '1%', '24%', '1%', '10%', '1%', '15%'), '', 
      uiOutput('selBlk'), '',
      actionButton("selBlkBut", label='Final confirmation', style='margin-bottom:-60px'), '',
      actionButton("selBlkCancel", label='Reset', style='margin-bottom:-60px'), '',
      downloadButton("dld", "Download", style = "margin-top: 24px;")
      ))
      ), uiOutput('dim.ui')
    )
  ) 
 
  ) # dashboardBody
)

}
