library(shiny); library(shinydashboard); library(plotly); library(visNetwork); library(DT)


shinyUI(dashboardPage(

  dashboardHeader(title="spatialHeatmap (2020-4-10)", titleWidth=300),

  dashboardSidebar(
  
    sidebarMenu(

      menuItem("Input", icon=icon("dashboard"),
      menuSubItem("View", tabName="hm_net"), br(),
      selectInput("fileIn", "Step1: data sets", c("none", "custom_data", "custom_computed_data", "brain_Prudencio", "mouse_Merkin", "chicken_Cardoso.Moreira", "shoot_Mustroph", "organ_Mustroph", "root_Mustroph", "shoot_root_Mustroph", "root_roottip_Mustroph", "map_Census"), "organ_Mustroph"),
      fileInput("svgInpath", "Step 2A: upload one aSVG file", accept=".svg", multiple=FALSE),
      fileInput("svgInpath1", "Step 2B: upload multiple aSVG files", accept=".svg", multiple=TRUE),
      fileInput("geneInpath", "Step 3: upload formatted data matrix", accept=c(".txt", ".csv"), multiple=FALSE),
      radioButtons(inputId='dimName', label='Step 4: is column or row gene?', choices=c("None", "Row", "Column"), selected="None", inline=TRUE),
      selectInput('sep', 'Step 5: separator', c("None", "Tab", "Space", "Comma", "Semicolon"), "None"),
      h4(strong("Custom computed data")),
      fileInput("adj.modInpath", "Upload the adjacency matrix and module definition file", accept=".txt", multiple=TRUE)

      ),
      menuItem("Instruction", tabName="ins", icon=icon("dashboard")),
      menuItem("Acknowledgement", tabName="ack", icon=icon("dashboard"))

     )

  ),

  dashboardBody(
 
    tags$head(tags$style(HTML(".shiny-output-error-validation { color: red; } "))),
    tags$head(tags$style(".shiny-notification {position: fixed; opacity: 1; 
    top: 35%; left: 40%; height: 90px; width: 400px}"
    )),
    tabItems(
      tabItem(tabName="hm_net", 

      box(title="Data Matrix", status="primary", solidHeader=TRUE, collapsible=TRUE, height=NULL, width=12,
      fluidRow(splitLayout(cellWidths=c('1%', '11%', '1%', '23%', '1%', '18%', '1%', '18%', '7%', '1%', '14%'), '',
      numericInput(inputId="A", label="Value (A) to exceed:", value=0), '',
      numericInput(inputId="P", label="Proportion (P) of samples with values >= A:", value=0), '',
      numericInput(inputId="CV1", label="Min coefficient of variation (CV1):", value=-Inf), '', 
      numericInput(inputId="CV2", label="Max coefficient of variation (CV2):", value=Inf), 
      actionButton(inputId='fil.but', label="Filter", icon=icon("refresh")), '',
      radioButtons(inputId='log', label='Log/exp scaling:', choices=c("No", "log2", "exp2"), selected="No", inline=TRUE)
      )), 
      fluidRow(column(2, uiOutput('linear.scl')), column(10, uiOutput('col.order'))),
      verbatimTextOutput("fil.par"),
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput("dt"), ""))
      ),
      fluidRow(column(10), column(2, radioButtons(inputId='hide.lgd', label="Hide legend:", choices=c('Yes'='Y', 'No'='N'), selected="N", inline=TRUE))),
      uiOutput('shm.ui'), uiOutput('lgd.ui'), 
      # 'height=NULL': height is automatic.
      box(title="Matrix Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, width=12, height=NULL, 
      fluidRow(splitLayout(cellWidths=c('1%', '15%', '1%', '10%', '1%', '20%', '1%', '7%', '1%', '20%', '1%', '10%'), '', 
      radioButtons(inputId='measure', label="Measure:", choices=c('correlation'='cor', 'distance'='dis'), selected="cor", inline=TRUE), '', 
      radioButtons(inputId="cor.abs", label="Cor.absolute:", choices=c('No', 'Yes'), selected='No', inline=TRUE), '', 
      radioButtons(inputId="thr", label="Select by:", choices=c('proportion'='p', 'number'='n', 'value'='v'), selected='p', inline=TRUE), '',
      numericInput(inputId='mhm.v', label='Cutoff: ', value=0.2, min=-Inf, max=Inf, step=NA, width=NULL), '',
      radioButtons(inputId="mat.scale", label="Scale:", choices=c("No", "By column", "By row"), selected="No", inline=TRUE), '',
      radioButtons(inputId="mhm.but", label="Show plot:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)
      )),

      plotlyOutput("HMly")), br(),
      box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, width=12,  
      tabBox(title="", width=12, id='inter_net', selected='inter_net', side='right', 
      tabPanel(title='Parameter', value='par', 
      selectInput("net.type", "Network type:", c('signed', 'unsigned', 'signed hybrid', 'distance'), "signed", width=210),
      numericInput("min.size", "Minmum module size:", 15, min=15, max=5000, width=210),
      selectInput("gen.sel","Select a target gene:", c("None"), selected="None", width=210),
      selectInput("ds","Module splitting sensitivity level:", 3:2, selected="3", width=210),

      selectInput("TOM.in", "Adjcency threshold:", c("None", sort(seq(0, 1, 0.002), decreasing=TRUE)), "None", width=210), 
      htmlOutput("edge"),
      textInput("color.net", "Color scheme:", "purple,yellow,blue", placeholder='Eg: "purple,yellow,blue"', width=210), actionButton("col.but.net", "Go", icon=icon("refresh")),
      radioButtons(inputId="cpt.nw", label="Show plot:", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)

      ), 
      tabPanel(title="Network", fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput("bar.net"), visNetworkOutput("vis"), "")), value='inter_net')))

      ),

      tabItem(tabName="ins", 
      box(title="Summary", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("sum"), width=12),
      box(title="Input", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("input"), downloadButton("dld.sgl", "Download single aSVG/data"), downloadButton("dld.mul", "Download multiple aSVGs/data"), htmlOutput("input1"), width=12),
      box(title="Data Matrix", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("matrix"), width=12), 
      box(title="Spatial Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("shm.ins"), width=12),
      box(title="Matrix Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("mhm.ins"), width=12),
      box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, htmlOutput("net.ins"), width=12)
      ),
      tabItem(tabName="ack", htmlOutput("ack")) 

      )

    )

))
