library(shiny); library(shinydashboard); library(plotly); library(visNetwork); library(DT)

ins <- "This Shiny App based Spatial Heatmap can be used for interactive visualisation as long as a data matrix and an associated SVG image are provided. In the following the instructions are given with a gene expression matrix and an associated root tissue image in SVG format. This app has three main functionalities. First, it generates spatial heatmaps where user defined tissue regions are coloured by the expression profile of a gene of interest. The gene expression information is uploaded as a table matrix and the associated tissue image is uploaded as an SVG format. Second, the app computes gene network modules based on gene expression profiles across tissue samples. It uses a matrix heatmap to visualize the expression of a chosen gene in the context of the corresponding gene network module the chosen gene belongs to. Third, the app displays the same network module in the matrix heatmap in the form of an interactive network graph. The network module identification is computationally demanding for large gene expression matrix (e.g.: > 10,000 genes), so to make this app more widely applicable the \"Compute locally\" mode is developed for processing large data matrix. If the data matrix is small (e.g.: < 10,000 genes), the \"Compute online\" mode can be used. If the app times out at some time, then users want to refresh the page in their web browser."

ins.input1 <- "At first, users need to select a mode under \"Select a work mode\" in the left \"Input\" menu. The \"Default\" is most convenient for users to test the app, since this option relies on pre-uploaded files and users do not need to upload any files at all. The \"Compute locally\" should be selected if users have a large gene expression file (e.g.: > 10,000 genes) while the \"Compute online\" can be selected if users have a small gene expression file (e.g.: < 10,000 genes). In \"Step 1: upload an svg file\" and \"Step 2: upload a gene expression file\", users are asked to upload the svg file and associated gene expression file respectively. Details about how to properly format and associate custom SVG images with expression tables are provided "

ins.input2 <- "The \"Step 3: is column or row gene?\" option specifies if column or row is gene in the gene expression table, and \"Step 4: separator\" specifies the separator among the expression values. \"Step 5: Color scheme\" allows users to input colour components to construct colour scale for gene expression levels. Colours must only be sepatated by comma, e.g. the default is \"green,blue,purple,yellow,red\"."

ins.input3 <- "In the gene matrix, the dimension names are gene IDs and sample/conditions. The sample/condition names MUST be fomatted this way: a sample name is followed by double underscore then the condition, such as \"epidermis__standard_1h\", where epidermis is the tissue and standard_1h is the condition. One column or row of meta data (e.g. gene annotation) can also be included in parallel with sample/condition. In the names of sample/condition and meta data, only letters, digits, single underscore, dots are allowed. The example SVG image and associated gene expression matrix can be downloaded in the instruction page of this app, which can be uploaded directly for testing in \"Compute online\" mode: "

ins.input4 <- "The \"Compute locally\" option is designed for large gene expression data (e.g.: > 10,000 genes), since gene network modules are identified using the R package WGCNA and the computation of topological overlap matrix (TOM) is time comsuming for large expression matrix. To maintain good performance, this 
process is expected to be performed on user's local computer. The tutorial of 
how to compute locally is provided in the documentation of the function \"filter.data\" and \"adj.mod\" in the R package \"spatialHeatmap\"."

ins.input5 <- "The \"Compute online\" option is designed for small gene expression data (e.g.: < 10,000 genes). The first two items filter genes according to a proportion that a gene's expression values exceed a threthold A across all samples. Only the genes exceeding the proportion will be maintained. The third and fourth items filter genes according to the coefficient of variation (CV). Only the genes with CV between the two specified values are maintained. The genes passing all these criteria are retained for downstream analysis. To save time, the app is designed to internally compute TOM only once when the matrix heatmap is displayed for the first time, but if the gene expression martix or its filter parameters are changed, the TOM will be re-computed."

ins.input6 <- "The \"Minmum module size\" sets the minimum module size in gene module indentification. In \"Network type\", \"Signed\" means both positive and negative adjacency between genes are used in network module identification while \"Unsigned\" takes the absolute values of negative adjacency."

ins.mat <- "The gene expression data is represented as an interactive table under \"Expression Matrix\", where row names are gene IDs and column names are samples/conditions and meta data (the dimension names of the original table are adjusted internally in \"Step 3: is column or row gene?\"). Users can sort the expression values for a sample or search for a particular gene by its ID. Users can select multiple genes by clicking IDs in the table to display corresponding spatial heatmaps. In each such spatial heatmap, the gene expression levels are represented by colours for each sample under each condition. In the \"Spatial Heatmap\" section, users can customise the dimension and layout of the spatial heatmaps."

ins.mat_hm <- "In the \"Matrix Heatmap & Network\" section, all gene IDs chosen in \"Expression Matrix\" are listed under \"Select a gene to display matrix heatmap & network.\". After a gene is selected from this list, the gene module containing the selected gene would be displayed in the form of interactive matrix heatmap, where the rows and columns are sorted by hierarchical clustering dendrograms and the chosen gene is tagged by a red rectangle. To explore the results, the matrix heatmap has several interactive features. For instance, users can zoom in and out by drawing a rectangle and by double clicking the image, respectively. Users can scale the expression values by gene or sample. The scaled values are only used for matrix heatmap display, not for downstream module identifications."

ins.net <- "The gene network modules are identified at two alternative sensitivities levels (3, 2). From 3 to 2, the sensitivity decreases and results in less modules with larger sizes. The \"Select a module splitting sensitivity level\" option allows users to choose which level to use for displaying the iteractive matrix heatmap and network."

ins.net2 <- "The module in \"Matrix Heatmap\" is displayed as an interactive network. Nodes and edges mean genes and adjacency between genes respectively. There is an interactive colour bar to denote gene connectivity (sum of a gene's adjacency with all its direct neighbours). The colour ingredients must only be separated by comma, e.g. the default are \"blue,green,red\", which means gene connectivity increases from blue to red. The edge length is inversely proportional to gene adjacency. If too many edges (e.g.: > 300) are displayed in this network, the app can possibly get stuck. So the \"Input an adjacency threshold to display the adjacency network.\" option sets a threthold to filter out some weak edges. Only edges above the threshold are displayed in the network. The app outputs the total number of remaining edges resulting from each input adjacency threshold. If it is not too large (e.g.: < 300), users can check \"Yes\" under \"Display or not?\", then the network can be responsive smoothly. To maintain acceptable performance, users are advised to choose a stringent threshold (e.g. 0.9) initially, then decrease the value gradually. The interactive feature allows users to zoom in and out, or drag a gene around. All the gene IDs in the network module are listed in \"Select by id\" in decreasing order according to gene connectivity. The selected gene ID is appended \"_selected\", which can be easily identified from the list. By clicking an ID in this list, users can identify the corresponding gene in the network."

shinyUI(dashboardPage(

  dashboardHeader(title="spatialHeatmap (updated: 2019-03-01)", 
  titleWidth=380),

  dashboardSidebar(
  
    sidebarMenu(

      menuItem("Input", icon=icon("dashboard"),
      menuSubItem("View", tabName="hm_net"), br(),
      selectInput("fileIn", "Select a work mode", c("None", "Default_organ", "Default_shoot_root", "Default_root_roottip", "Default_shoot", "Default_root", "Default_brain", "Default_map", "Compute locally", "Compute online"), "Default_shoot"),
      fileInput("svgInpath", "Step 1: upload an svg file", accept=".svg", multiple=FALSE),
      fileInput("geneInpath", "Step 2: upload a gene expression file", accept=c(".txt", ".csv"), multiple=FALSE),
      radioButtons('dimName', 'Step 3: is column or row gene?', c("None", "Row", "Column"), inline=TRUE),
      selectInput('sep', 'Step 4: separator', c("None", "Tab", "Comma", "Semicolon"), "None"),
      div(style="display:inline-block;width:75%;text-align:left;",textInput("color", "Step 5: Color scheme of Spatial Heatmap", "green,blue,purple,yellow,red", placeholder="Eg: green,yellow,red", width=200)),
      div(style="display:inline-block;width:25%;text-align:left;", actionButton("col.but", "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
      radioButtons('cs.v', 'Colour scale based on:', c("Selected genes"="sel.gen", "Whole matrix"="w.mat"), inline=TRUE),
      h4(strong("Compute locally")),
      fileInput("adj.modInpath", "Upload the adjacency matrix and module definition file", accept=".txt", multiple=TRUE),
      h4(strong("Compute online")),
      numericInput("A", "The value A to be exceeded (filter genes):", 0, min=0, max=10000), numericInput("p", "The proportion that need to exceed A (filter genes):", 0, min=0, max=1),
      numericInput("cv1", "Lower bound of coefficient of variation (CV) (filter genes):", 0, min=0, max=1),
      numericInput("cv2", "Upper bound of coefficient of variation (CV) (filter genes):", 10000, min=0, max=10000),
      numericInput("min.size", "Minmum module size:", 5, min=3, max=10000),
      radioButtons("net.type", "Network type", c("Signed"="S", "Unsigned"="U"), "S", inline=TRUE)

      ),

      menuItem("Heatmap & network", icon=icon("dashboard"), 
      h4(strong("Spatial heatmap")),
      selectInput("height", "Overall canvas height:", seq(100, 15000, 20), "400", width=170), 
      selectInput("width", "Overall canvas width:", seq(100, 15000, 20), "820", width=170),
      selectInput("col.n", "No. of columns for sub-plots", seq(1, 15, 1), "2", width=210), 
      radioButtons("gen.con", "Display by:", c("Gene"="gene", "Condition"="con"), "gene", inline=TRUE),
      h4(strong("Matrix heatmap & network")), 
      selectInput("gen.sel","Select a gene to display matrix heatmap & network.", c("None"),
      selected="None", width=190),
      selectInput("mat.scale", "Scale matrix heatmap", c("No", "By column/sample", 
      "By row/gene"), "No", width=190),
      selectInput("ds","Select a module splitting sensitivity level", 3:2, selected="3", width=190),

      selectInput("TOM.in", "Input an adjcency threshold to display the adjacency network.", c("None", sort(seq(0, 1, 0.002), decreasing=TRUE)), "None", width=190), 
      htmlOutput("edge"),
      div(style="display:inline-block;width:75%;text-align:left;",textInput("color.net", "Color scheme of Interactive Network", "green,blue,red", placeholder="Eg: green,yellow,red", width=200)),
      div(style="display:inline-block;width:25%;text-align:left;", actionButton("col.but.net", "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
      radioButtons("cpt.nw", "Display or not?", c("Yes"="Y", "No"="N"), "N", inline=TRUE),
      #menuSubItem("View network", tabName="network")

      menuSubItem("View", tabName="hm_net"), br()
      ),

      menuItem("Instruction", tabName="instruction", icon=icon("dashboard"))

     )

  ),

  dashboardBody(
 
    tabItems(

      tabItem(tabName="instruction", 
      box(title="General Instruction", status="primary", solidHeader=TRUE, collapsible=TRUE, ins, width=12),
      box(title="Input Instruction", status="primary", solidHeader=TRUE, collapsible=TRUE, p(ins.input1, HTML("<a href=
      http://biocluster.ucr.edu/~jzhan067/shiny_HM_tutorial/tutorial/shiny_heatmap_tutorial.html>here</a>"), "."), p(ins.input2), ins.input3, HTML("&nbsp"), downloadButton("dld.svg", "Download svg file"), downloadButton("dld.data", "Download gene matrix"), p(ins.input4), p(ins.input5), p(ins.input6), width=12), 
      box(title="Spatial Heatmap Instruction", status="primary", solidHeader=TRUE, 
      collapsible=TRUE, p(ins.mat), width=12),
      box(title="Matrix Heatmap & Network Instruction", status="primary", solidHeader=TRUE, collapsible=TRUE, p(ins.mat_hm), p(ins.net), p(ins.net2), width=12)
      ),
      
      tabItem(tabName="hm_net", 

      box(title="Expression Matrix", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput("dt"), "")), width=12),

      box(title="Spatial Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar"), plotOutput("tissue"), "")), width=9),

      box(title="Original image", status="primary", solidHeader=TRUE, collapsible=TRUE, splitLayout(cellWidths=c("1%", "98%", "1%"), "", plotOutput("ori.svg"), ""), width=3), br(),

      box(title="Matrix Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, plotlyOutput("HMly"), width=12, height=NULL), br(),

      box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput("bar.net"), visNetworkOutput("vis"), "")), width=12)
      )

      )

    )

))
