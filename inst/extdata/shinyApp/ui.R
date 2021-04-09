library(shiny); library(shinydashboard); library(shinydashboardPlus); library(yaml); library(plotly); library(visNetwork); library(DT); library(shinyWidgets); library(shinyBS); library(shinyjs)
source('R/function.R')
lis.cfg <- yaml.load_file('config/config.yaml')
tit <- sub('^(title|width):', '', lis.cfg$title)


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


jsFlip <- "var card = document.getElementById('card');
document.getElementById('flip').addEventListener('click', function() {
card.classList.toggle('flipped'); }, false)"


jsFlip <- "(function() {
  var card = document.querySelector('#card');
  var flipBtn = document.querySelector('#flip');
  flipBtn.addEventListener('click', function() {
    card.classList.toggle('flipped');
  });
})();"


jsFlip <- "
function flip() {
    $('.card').toggleClass('flipped');
}
"


cssFlip <- '
.container {
    width: 100%;
    height: 100%;
    position: relative;
    border: 1px solid #ccc;
    -webkit-perspective: 800px;
    -moz-perspective: 800px;
    -o-perspective: 800px;
    perspective: 800px;
}
.card {
    width: 100%; height: 100%;
    /*position: absolute;
   */
    -webkit-transition: -webkit-transform 1s;
    -moz-transition: -moz-transform 1s;
    -o-transition: -o-transform 1s;
    transition: transform 1s;
  /*  -webkit-transform-style: preserve-3d;
    -moz-transform-style: preserve-3d;
    -o-transform-style: preserve-3d;
    transform-style: preserve-3d;
    -webkit-transform-origin: 50% 50%;
  */
}
.card div {
    /*display: block; height: 100%; width: 100%;

    line-height: 260px; text-align: center; color: white; font-weight: bold; font-size: 140px; position: absolute;*/

    -webkit-backface-visibility: hidden;
    -moz-backface-visibility: hidden;
    -o-backface-visibility: hidden;
    backface-visibility: hidden;
}
.card .front {
  background: red;
}
.card .back {
    background: blue;
    -webkit-transform: rotateY( 180deg );
    -moz-transform: rotateY( 180deg );
    -o-transform: rotateY( 180deg );
    transform: rotateY( 180deg );
}
.card.flipped {
    -webkit-transform: rotateY( 180deg );
    -moz-transform: rotateY( 180deg );
    -o-transform: rotateY( 180deg );
    transform: rotateY( 180deg );
}
'


shinyUI(dashboardPage(

  
  # includeCSS("style.css"),
  dashboardHeader(title = NULL, titleWidth = 0),
  dashboardSidebar(collapsed = TRUE, disable = TRUE, width = 0, sidebarMenu() ),
  dashboardBody(
   # tags$head(HTML('<title>spatialHeatmap</title>')),
   useShinyjs(), 
   includeCSS("www/style.css"),
   tags$head(tags$link(rel="stylesheet", type="text/css", href="style.css")),
   
   tags$head(tags$script(src = "javascript.js"), tags$script(HTML(js)), tags$script(HTML(jsFlip)), tags$style(HTML(cssFlip))),
   includeScript(path = "www/javascript.js"),
   tags$script(src="javascript.js"),
   tags$head(HTML("<script type='text/javascript' src='javascript.js'></script>")),
   HTML('<head>
         <link rel="stylesheet" type="text/css" href="style.css">
         <script type="text/javascript" src="www/javascript.js"></script>
         </head>
   '),
   # To place title on the right of dashboard header.V
   tags$script(HTML('
     $(document).ready(function() {
       $("header").find("nav").append(\'<span class="myClass">spatialHeatmap</span>\');
     })
   ')),

    fluidRow(
    # tabItems(

      # tabItem(tabName="hm_net", 
     
     tabsetPanel(type = "pills", id=NULL, selected="shm",
      tabPanel("Upload", value='upload',
      box(width = 12, title = 'Upload data & aSVGs', closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      
      # menuItem("Input", icon=icon("list"),
      # menuSubItem("View", tabName="hm_net"), #br(),
      h4(strong("Step1: choose custom or default data sets")),
      fluidRow(splitLayout(cellWidths=c('1%', '30%'), '', 
      selectInput("fileIn", NULL, c('none', 'customData', 'customComputedData'), 'none')
      )), br(),
      h4(strong("Step 2: upload custom data")),
      fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '18%', '1%', '25%', '1%', '25%'), '',
      fileInput("geneInpath", "2A: upload formatted data matrix", accept=c(".txt", ".csv"), multiple=FALSE), '',
      radioButtons(inputId='dimName', label='2B: is column or row gene?', choices=c("None", "Row", "Column"), selected='None', inline=TRUE), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "columns in the data matrix corresponds with "rows" in the targets file respectively.'),
      fileInput("target", "2C (optional): upload targets file for columns", accept=c(".txt", ".csv"), multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'Ensure "rows" in the data matrix corresponds with "rows" in the row metadata file respectively.'),
      fileInput("met", "2D (optional): upload metadata file for rows", accept=c(".txt", ".csv"), multiple=FALSE))

      )),
      h4(strong("Step 3: upload custom aSVG(s)")),
      fluidRow(splitLayout(cellWidths=c('1%', '27%', '1%', '28%'), '',
      tags$div(class='tp', span(class='tpt', 'The data is matched with a single aSVG file.'),
      fileInput("svgInpath1", "3A: upload one aSVG file", accept=".svg", multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'The data is matched with multiple aSVG files (e.g. developmental stages).'),
      fileInput("svgInpath2", "3B (optional): upload multiple aSVG files", accept=".svg", multiple=TRUE))
      )),
      bsTooltip(id='svgInpath', title = "The data is matched with a single aSVG file.", placement = "bottom", trigger = "hover"),
      # h4(strong("Custom computed data")),
      # fileInput("adj.modInpath", "Upload the adjacency matrix and module definition file", accept=".txt", multiple=TRUE)
      div(style = "font-size: 10px; padding: 0px 0px; margin:0%",
      fluidRow(splitLayout(cellWidths=c('1%', '27%', '1%', '28%', '1%', '28%'), '',
      downloadButton("dld.sgl", "Example1: data matched with a single aSVG"), '',
      downloadButton("dld.mul", "Example2: data matched with multiple aSVGs"), '',
      downloadButton("dld.st", "Example3:  spatiotemporal data-aSVG")
      ))), br(), 

      #box(title="Download", status="primary", solidHeader=TRUE, collapsible=TRUE,
      #htmlOutput("dld"), downloadButton("dld.cfg", "Download config file"), 
      #width=12),

      h4(strong("Additional files")), 
      fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '35%'), '',
      tags$div(class='tp', span(class='tpt', 'Upload a config file in ".yaml" format.'),
      fileInput("config", "Upload a config file (optional)", accept=".yaml", multiple=FALSE)), '',
      tags$div(class='tp', span(class='tpt', 'The batched data sets will be listed under "Step 1".'),
      fileInput("tar", "Upload batched data, aSVGs in two separate tar files (optional)", accept=c(".tar"), multiple=TRUE))
      )),
      div(style = "font-size: 10px; padding: 0px 0px; margin:0%",
      fluidRow(splitLayout(cellWidths=c('1%', '24%', '1%', '35%'), '',
      downloadButton("dld.cfg", "Example config file"), '', downloadButton("dld.bat", "Example data/aSVGs in batch")
      )))
      )
      # menuItem("Instruction", tabName="ins", icon=icon("info")),
      # menuItem("Acknowledgement", tabName="ack", icon=icon("info")),
      # menuItem('Download', tabName='dld', icon=icon('download'))
      
      ),# tabPanel
      tabPanel("Primary Visualization", value='shm',

  box(width = 12, title = "Data (replicates aggregated)", closable = FALSE, solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE,
      navbarPage('Parameters:',
      tabPanel("Basic",
      fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%', '1%', '10%', '1%', '10%'), '',
      radioButtons(inputId='scale', label='Scale by', choices=c('No', 'Row', 'Column'), selected='No', inline=TRUE), '', 
      #log_exp_ui('data')
      radioButtons(inputId='log', label='Log/exp-transform', choices=c("No", "log2", "exp2"), selected='No', inline=TRUE)
      # dropdownButton(inputId='drdnFil', label='Filter', circle=FALSE, icon=NULL, status='primary')
      ))), # tabPanel
      tabPanel("Filter",
      fluidRow(splitLayout(cellWidths=c('1%', '14%', '1%', '30%', '1%', '24%', '1%', '24%'), '',
      numericInput(inputId="A", label="Value (A) to exceed", value=0), '',
      numericInput(inputId="P", label="Proportion (P) of samples with values >= A", value=0), '',
      numericInput(inputId="CV1", label="Min coefficient of variation (CV1)", value=-10^4), '', 
      numericInput(inputId="CV2", label="Max coefficient of variation (CV2)", value=10^4)
      )), actionButton(inputId='fil.but', label="Submit"), verbatimTextOutput("fil.par") 
      ),
      tabPanel("Re-order columns",
      column(12, uiOutput('col.order')) 
      )
      ), # navbarPage 
      fluidRow(column(1, offset=0, style='padding-left:10px; padding-right:10px; padding-top:0px; padding-bottom:5px',
      bsTooltip(id='drdnFil', title="Filter the data by a proportion of a threthold value across all samples (pOverA) and coefficient of variance (CV).", placement = "bottom", trigger = "hover"),
      bsTooltip(id='P', title="Filter the data by a proportion of a threthold value across all samples (pOverA) and coefficient of variance (CV).", placement = "bottom", trigger = "hover")
      ),
      #column(1, offset=0, style='padding-left:10px; padding-right:10px; padding-top:0px; padding-bottom:5px',
      #dropdownButton(inputId='drdn.scale', label='Transform', circle=FALSE, icon=NULL, status='primary', width=400
      #))
      # column(10, offset=0, style='padding-left:40px; padding-right:10px; padding-top:0px; padding-bottom:5px', uiOutput('col.order'))
      ),
      fluidRow(
      column(6, textInput(inputId='search', label='Search:', value='', placeholder='Muliple IDs must only be separated by space or comma.', width='100%')),
      column(1, actionButton('search.but', 'Submit'), align = "center", style = "margin-bottom: 2px;", style = "margin-top: 20px;")
      ),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput("dt"), ""))
      )
      ),
      tabPanel('Spatial Enrich', value='deg', 
      
  box(width = 12, title = "Data with Replicates", closable = FALSE,
  solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary",
  enable_dropdown = FALSE,
      
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput("dt.rep"), ""))
      ),


box(title='Spatial Enrich', status="primary", solidHeader=TRUE, collapsible=TRUE, width=12,
  # tabBox(title="", width=12, id="ssg", selected='w.meth', side='right',

  navbarPage(NULL, 
    # tabPanel(title="Across methods table", dataTableOutput("a.table"), value='a.table'), 
    # tabPanel(title="Upset Plot Down", plotOutput("upset2"), value='a.meth'), 
    tabPanel(title="Parameters",
      fluidRow(splitLayout(cellWidths=c('1%', '30%', '1%', '30%', '1%', '12%', '1%', '12%'), '',
      uiOutput('ssg.sam'), '', uiOutput('ssg.con'), '', 
      selectInput("sam.con", "Compare by", c('feature', 'factor', 'feature__factor'), 'feature'), '',
      uiOutput('ssg.tar')
      )),
      fluidRow(splitLayout(cellWidths=c('1%', '20%', '1%', '20%', '1%', '20%', '1%', '10%', '1%', '10%'), '',
      selectInput('ssg.meth', label='Select methods', choices=c('edgeR', 'limma', 'DESeq2', 'ROKU', 'distinct'), selected=c('edgeR', 'limma'), multiple=TRUE), '',
      selectInput("edg.lim.nor", "Normalize-edgeR/limma", c("CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "none"), 'CNF-TMM'), '',
      selectInput("rok.dis.nor", "Normalize-ROKU/distinct", c("CNF-TMM", "CNF-TMMwsp", "CNF-RLE", "CNF-upperquartile", "ESF", "VST", "rlog", "none"), 'CNF-TMM'), '',
      numericInput('ssg.fc', 'Log2-fold change', 1, min=0, max=Inf, step=1), '',
      numericInput('ssg.fdr', 'FDR', 0.05, min=0, max=1, step=0.01)
      )),
      actionButton("ssg.update", "Update", icon=icon("refresh"), style="color: #fff; background-color:purple;border-color: #2e6da4") 
    ), # tabPanel
    tabPanel(title="Results in plots",
      #column(12, align = "center", style = "font-weight:bold; font-size:20px;", textOutput('deg.vs')),
      strong('Over-Expressed'), br(), column(12, column(8, plotOutput("upset1")), column(4, plotOutput("ovl1"))),
      strong('Under-Expressed'), br(), column(12, column(8, plotOutput("upset2")), column(4, plotOutput("ovl2"))), 
      value='w.meth'
    ),
    tabPanel(title="Pairwise comparisons", value = 'vs', 
      column(12, align = "center", style = "font-weight:bold; font-size:20px;", textOutput('deg.vs')),
    ),
    tabPanel(title="Link with spatial heatmap", 
      downloadButton("dld.ssg.tab", "Download"), dataTableOutput("dt.deg"), value='dt.deg'
    ) 
  ) # navbarPage
)

      ), # tabPanel(value='deg',

  # box(width = 12, height = NULL, title = p("Spatial Heatmap", HTML(paste0(rep('&nbsp;', 2), collapse=''))), closable = FALSE, sidebar_icon = "info", solidHeader = TRUE, collapsible = TRUE, enable_sidebar = FALSE, status = "primary", enable_dropdown = FALSE, dropdown_icon = "wrench", dropdown_menu =NULL,

#tags$div(id = 'options', HTML('<button id="flip">Flip Card</button>')),
#tags$div(class="container", 
#div(id = "card", 
# tags$figure(class = 'front', 'sdgfg'), tags$figure(class = 'back', 'dfdfd'))

#HTML('<figure class="front">1</figure>
#<figure class="back">2</figure>'))
#), br(), br(), br(),


tabsetPanel(type = "pills", id = 'shm.sup', selected = "shm2",

   tabPanel("Spatial Heatmap", value = 'shm2', icon = icon('image'),
  
   p(tags$head(tags$style(HTML('#lgdB{background-color:#d2d6de;}'))), actionButton("lgdB", "Toggle Legend Plot", icon=icon('toggle-on'), class="btn-xs", title = "Click to toggle legend plot"), hover = TRUE),
   uiOutput('shm.ui'), uiOutput('lgd.ui')
# HTML('<button onclick="flip()">flip the card</button>'),
# tags$section(class="container",
# div(class="card",
#  div(class="front",  uiOutput('shm.ui'), uiOutput('lgd.ui')),
#  div(class="back", 'ddfdf')
#)
#)
  ),


  tabPanel("Matrix Heatmap", id = 'mhm', icon = icon('border-all'),
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '10%', '1%', '20%', '1%', '7%', '1%', '20%', '1%', '10%'), '', 
      radioButtons(inputId='measure', label="Measure:", choices=c('correlation', 'distance'), selected='correlation', inline=TRUE), '', 
      radioButtons(inputId="cor.abs", label="Cor.absolute:", choices=c('No', 'Yes'), selected='No', inline=TRUE), '', 
      radioButtons(inputId="thr", label="Select by:", choices=c('proportion'='p', 'number'='n', 'value'='v'), selected='p', inline=TRUE), '',
      numericInput(inputId='mhm.v', label='Cutoff: ', value=0.2, min=-Inf, max=Inf, step=NA, width=NULL), '',
      radioButtons(inputId="mat.scale", label="Scale by:", choices=c("No", "Column", "Row"), selected='No', inline=TRUE), '',
      # radioButtons(inputId="mhm.but", label="Show plot:", choices=c("Yes", "No"), selected='No', inline=TRUE)
      actionButton("mhm.but", "Update", icon=icon("refresh"), style="color: #fff; background-color:purple;border-color: #2e6da4")
      )),
      plotlyOutput("HMly")
      ),
  tabPanel("Interactive Network", id = 'net', icon = icon('project-diagram'), 
      fluidRow(
      # If column widths are not integers, columns are vertically aligned.
      column(2, offset=0, style='padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      textInput("color.net", "Color scheme:", 'yellow,orange,red', placeholder='Eg: yellow,orange,red', width=150), actionButton("col.but.net", "Go", icon=icon("refresh"))
      ),
      column(1, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      selectInput(inputId="net.type", label="Network type:", choices=c('signed', 'unsigned', 'signed hybrid', 'distance'), selected='signed', width=100)
      ),
      column(1, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      numericInput("min.size", "Minmum module size:", value=15, min=15, max=5000, width=150)
      ),
      column(2, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      tags$span(style="color:brown;font-weight:NULL", selectInput("gen.sel", "Select a target gene:", c("None"), selected='None', width=150))
      ),
      column(2, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      selectInput("ds","Module splitting sensitivity level:", 3:2, selected=3, width=170)
      ),
      column(1, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      selectInput("adj.in", "Adjacency threshold (the smaller, the more edges):", sort(seq(0, 1, 0.002), decreasing=TRUE), selected=1, width=200)
      ),
      column(2, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      tags$span(style="color:brown", numericInput("max.edg", "Maximun edges (too many edges may crash the app):", value=50, min=1, max=500, width=200)),
      htmlOutput("edge")
      ),
      column(1, offset=0, style='padding-left:10px; padding-right:0px; padding-top:0px; padding-bottom:5px',
      actionButton("cpt.nw", "Update", icon=icon("refresh"), style="color: #fff; background-color:purple;border-color: #2e6da4")
      )
      ),
      fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput("bar.net"), visNetworkOutput("vis"), ""))
      ) # tabPanel("Interactive Network"

    ) # tabsetPanel(id = 'shm.sup',


#   )

# )
)
      # 'height=NULL': height is automatic.
      #box(title="Matrix Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, width=12, height=NULL ), br(),
      #box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, width=12 ),

      #tabItem(tabName="ins", 
      #box(title="Summary", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("sum"), width=12),
      #box(title="Input", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("input"), width=12),
      #box(title="Data Matrix", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("matrix"), width=12), 
      #box(title="Spatial Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("shm.ins"), width=12),
      #box(title="Matrix Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("mhm.ins"), width=12),
      #box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("net.ins"), width=12)
      #),
      #tabItem(tabName="ack", 
      #htmlOutput("ack"),
      #),
      #tabItem(tabName="dld",
      #) 

      # )

    )
 )

))
