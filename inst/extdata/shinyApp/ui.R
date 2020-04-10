library(shiny); library(shinydashboard); library(plotly); library(visNetwork); library(DT)

shinyUI(dashboardPage(

  dashboardHeader(title="spatialHeatmap (2020-03-05)", titleWidth=300),

  dashboardSidebar(
  
    sidebarMenu(

      menuItem("Input", icon=icon("dashboard"),
      menuSubItem("View", tabName="hm_net"), br(),
      selectInput("fileIn", "Select a work mode", c("None", "Compute online", "Compute locally", "brain_Prudencio", "mouse_Merkin", "chicken_Cardoso.Moreira", "shoot_Mustroph", "organ_Mustroph", "root_Mustroph", "shoot_root_Mustroph", "root_roottip_Mustroph", "map_Census"), "organ_Mustroph"),
      fileInput("svgInpath", "Step 1: upload formatted SVG image", accept=".svg", multiple=FALSE),
      fileInput("geneInpath", "Step 2: upload formatted data matrix", accept=c(".txt", ".csv"), multiple=FALSE),
      radioButtons(inputId='dimName', label='Step 3: is column or row gene?', choices=c("None", "Row", "Column"), selected="None", inline=TRUE),
      selectInput('sep', 'Step 4: separator', c("None", "Tab", "Comma", "Semicolon"), "None"),
      radioButtons(inputId='log', label='Data transform', choices=c("No", "log2", "exp2"), selected="No", inline=TRUE),
      div(style="display:inline-block;width:75%;text-align:left;",textInput("color", "Step 5: Color scheme of Spatial Heatmap", "yellow,purple,blue", placeholder="Eg: yellow,purple,blue", width=200)),
      div(style="display:inline-block;width:25%;text-align:left;", actionButton("col.but", "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
      radioButtons(inputId='cs.v', label='Colour scale based on:', choices=c("Selected genes"="sel.gen", "Whole matrix"="w.mat"), selected="sel.gen", inline=TRUE),
      h4(strong("Compute locally")),
      fileInput("adj.modInpath", "Upload the adjacency matrix and module definition file", accept=".txt", multiple=TRUE),
      h4(strong("Compute online")),
      textInput(inputId="P", label="Filter genes: the proportion (P) of samples whose values exceed A:", value=0, width=NULL, placeholder='a numeric: 0-1'),
      textInput(inputId="A", label="Filter genes: the value (A) to be exceeded:", value='-Inf', width=NULL, placeholder='a numeric'),
      textInput(inputId="CV1", label="Filter genes: lower bound of coefficient of variation (CV1):", value='-Inf', width=NULL, placeholder='a numeric'),
      textInput(inputId="CV2", label="Filter genes: upper bound of coefficient of variation (CV2):", value='Inf', width=NULL, placeholder='a numeric'), actionButton(inputId='fil.but', label="Filter", icon=icon("refresh")),
      numericInput("min.size", "Minmum module size:", 15, min=15, max=1000),
      radioButtons(inputId="net.type", label="Network type", choices=c("Signed"="S", "Unsigned"="U"), selected="S", inline=TRUE)

      ),

      menuItem("Heatmap & Network", icon=icon("dashboard"), 
      h4(strong("Spatial Heatmap")),
      selectInput("height", "Overall canvas height:", seq(100, 15000, 20), "400", width=170), 
      selectInput("width", "Overall canvas width:", seq(100, 15000, 20), "820", width=170),
      selectInput("col.n", "No. of columns for sub-plots", seq(1, 15, 1), "2", width=210), 
      radioButtons(inputId="gen.con", label="Display by:", choices=c("Gene"="gene", "Condition"="con"), selected="gene", inline=TRUE),
      h4(strong("Matrix Heatmap & Network")), 
      selectInput("gen.sel","Select a gene to display matrix heatmap & network.", c("None"),
      selected="None", width=190),
      selectInput("mat.scale", "Scale matrix heatmap", c("No", "By column/sample", 
      "By row/gene"), "By row/gene", width=190),
      selectInput("ds","Select a module splitting sensitivity level", 3:2, selected="3", width=190),

      selectInput("TOM.in", "Input an adjcency threshold to display the adjacency network.", c("None", sort(seq(0, 1, 0.002), decreasing=TRUE)), "None", width=190), 
      htmlOutput("edge"),
      div(style="display:inline-block;width:75%;text-align:left;",textInput("color.net", "Color scheme of Interactive Network", "yellow,purple,blue", placeholder="Eg: yellow,black,blue", width=200)),
      div(style="display:inline-block;width:25%;text-align:left;", actionButton("col.but.net", "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
      radioButtons(inputId="cpt.nw", label="Display or not?", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE),
      #menuSubItem("View network", tabName="network")
      menuSubItem("View", tabName="hm_net"), br()
      ),
      menuItem("Instruction", tabName="ins", icon=icon("dashboard")),
      menuItem("Acknowledgement", tabName="ack", icon=icon("dashboard"))

     )

  ),

  dashboardBody(
 
    tabItems(
 
      tabItem(tabName="hm_net", 

      box(title="Expression Matrix", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput("dt"), "")), height=430, width=12),
      box(title="Spatial Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, 
      fluidRow(splitLayout(cellWidths=c('1%', "99%"), '', checkboxGroupInput(inputId="tis", label="Select tissues to be transparent:", choices='', selected='', inline=TRUE))),
      fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar"), plotOutput("shm"), "")), width=9),
      box(title="Original Image", status="primary", solidHeader=TRUE, collapsible=TRUE, splitLayout(cellWidths=c("1%", "98%", "1%"), "", plotOutput("lgd"), ""), width=3), br(),
      box(title="Matrix Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, plotlyOutput("HMly"), width=12, height=460), br(),
      box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput("bar.net"), visNetworkOutput("vis"), "")), width=12)
      ),

      tabItem(tabName='ins', htmlOutput("ins")),
      tabItem(tabName="ack", htmlOutput("ack")) 

      )

    )

))
