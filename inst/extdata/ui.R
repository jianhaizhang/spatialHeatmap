library(shiny); library(shinydashboard); library(plotly); library(visNetwork); library(DT)

ins0 <- "For a quick test, the top 7 examples under the \"Select a work mode\" menu can be selected. Alternatively, the pre-uploaded svg image and expression matrix can be downloaded from below and uploaded to the \"Compute online\" work mode."

ins <- "This Shiny-app based spatialHeatmap is an integrated application of the R/Bioconductor package \"spatialHeatmap\". It is designed for interactively visualising gene expression data in a matrix format such RNA-seq, microarray, qPCR, etc., but can be used for visualisation on other data sources as long as a data matrix and an associated svg image are provided. In the following the instructions are given with a gene expression matrix and an associated root tissue image in svg format. This app has three main functionalities. First, it generates spatial heatmaps where user defined tissue regions are coloured by the expression profile of a gene of interest underdifferent conditions. The gene expression information is uploaded as a table matrix and the associated tissue image is uploaded as an svg format. Second, the app computes gene network modules based on gene expression profiles across tissue samples. It uses a matrix heatmap to visualise the expression of a chosen gene in the context of the corresponding gene network module. Third, the app displays the same network module in the matrix heatmap in the form of an interactive network graph. The network module identification is computationally intensive for large gene expression matrix (e.g.: > 10,000 genes), so to make this app more widely applicable the \"Compute locally\" mode is developed for processing large data matrix. If the data matrix is small (e.g.: < 10,000 genes), the \"Compute online\" mode can be used. If the app times out at some time, then users want to refresh the page in their web browser."

ins.input1 <- "At first, users need to select a mode under \"Select a work mode\" in the left \"Input\" menu. The top 7 examples are most convenient for users to test the app, since this option employs on pre-uploaded files and users do not need to upload any files at all. The \"Compute locally\" should be selected if users have a large gene expression file (e.g.: > 10,000 genes) while the \"Compute online\" can be selected if users have a small gene expression file (e.g.: < 10,000 genes). In \"Step 1: upload an svg file\" and \"Step 2: upload a gene expression file\", users are asked to upload the svg image and associated gene expression file respectively. Details about how to properly format and associate custom svg images with expression tables are provided in the package vignette and "

ins.input2 <- "The \"Step 3: is column or row gene?\" option specifies if column or row is gene in the expression table, and \"Step 4: separator\" specifies the separator in the table. \"Step 5: Color scheme\" allows users to input colour components to construct a colour scale for gene expression levels. Colours MUST only be sepatated by comma, e.g. the default is \"green,blue,purple,yellow,red\"."

ins.input3.1 <- "In the expression matrix, the row and column names are usually gene IDs and sample/conditions, respectively. " 
ins.input3.2 <- "The sample/condition names MUST be fomatted this way: a sample name is followed by double underscore then the condition, such as \"stele__140mM_48h\", where stele is the sample and 140mM_48h is the condition. "
ins.input3.3 <- "One column of meta data (e.g. gene annotation) can also be included in parallel with sample/condition at the end. "
ins.input3.4 <- "In the column names of sample/condition and meta data, only letters, digits, single underscore, dots are allowed. "
ins.input3.5 <- "Not all samples in the matrix necessarily need to be present in the svg image, and vice versa. Only samples present in the svg image are recognised and coloured. The example svg image and associated gene expression matrix can be downloaded in this page and uploaded directly for testing in \"Compute online\" mode: "

ins.input4 <- "The \"Compute locally\" option is designed for large gene expression data (e.g.: > 10,000 genes), since gene network modules are identified using the R package WGCNA and the computation is intensive for large expression matrix. To maintain good performance on large matrix, this process is expected to be performed on user's local computer. The instruction on how to compute locally is provided in the documentation of the function \"filter.data\" and \"adj.mod\" in the R package \"spatialHeatmap\"."

ins.input5 <- "The \"Compute online\" option is designed for small gene expression matrix (e.g.: < 10,000 genes). The first two items filter genes according to a proportion of samples where a gene's expression values exceed a threthold A. Only the genes exceeding the proportion will be maintained. The third and fourth items filter genes according to the coefficient of variation (CV). Only the genes with CV between the two specified values are maintained. The genes passing all these criteria are retained for downstream analysis. To save time, the app is designed to compute network modules only once when the matrix heatmap is displayed for the first time, and the module results are saved, so displaying matrix heatmap of another gene will not trigger the computation again. But if the gene expression martix or filter parameters are changed, the modules are computated again."

ins.input6 <- "The \"Minmum module size\" sets the minimum module size in gene module indentification. In \"Network type\", \"Signed\" means both positive and negative adjacency between genes are used in network module identification while \"Unsigned\" takes the absolute values of negative adjacency."

ins.mat <- "The gene expression data is represented as an interactive table under \"Expression Matrix\", where row names are gene IDs and column names are samples/conditions and meta data (the dimension names of the original table are adjusted internally in \"Step 3: is column or row gene?\"). Users can sort the expression values for a sample or search for a particular gene by its ID. Users can select multiple genes by clicking IDs in the table to display corresponding spatial heatmaps. In each such spatial heatmap, the gene expression profiles are represented by colours for each sample under each condition. In the \"Spatial Heatmap\" section, users can customise the dimension and layout of the spatial heatmaps."

ins.mat_hm <- "In the \"Matrix Heatmap & Network\" section, all gene IDs chosen in \"Expression Matrix\" are listed under \"Select a gene to display matrix heatmap & network.\". Once selected from this list, the gene is displayed in the context of its module in the form of interactive matrix heatmap, where the rows and columns are sorted by hierarchical clustering dendrograms and the chosen gene is tagged by a red rectangle. To explore the results, users can zoom in and out by drawing a rectangle and by double clicking the image, respectively. Users can scale the expression values by gene or sample. The scaled values are only used for matrix heatmap display, not for module identification."

ins.net <- "The network modules are identified at two alternative sensitivities levels (3, 2). From 3 to 2, the sensitivity decreases and results in less modules with larger sizes. The \"Select a module splitting sensitivity level\" option allows users to choose which level to use for displaying the iteractive matrix heatmap and network."

ins.net2 <- "The same module in \"Matrix Heatmap\" is displayed as an interactive network. Nodes and edges are genes and adjacency between genes respectively. There is an interactive colour bar to denote gene connectivity (sum of a gene's adjacency with all its direct neighbours). The colour ingredients MUST only be separated by comma, e.g. the default are \"green,blue,red\", which means gene connectivity increases from blue to red. The edge length is inversely proportional to gene adjacency. If too many edges (e.g.: > 300) are displayed in this network, the app can possibly get stuck. So the \"Input an adjacency threshold to display the adjacency network.\" option sets a threthold to filter out some weak edges. Only edges above the threshold are displayed. The app outputs the total number of remaining edges resulting from each input adjacency threshold. If it is not too large (e.g.: < 300), users can check \"Yes\" under \"Display or not?\", then the network will be displayed and can be responsive smoothly. To maintain acceptable performance, users are advised to choose a stringent threshold (e.g. 0.9) initially, then decrease the value gradually. The interactive feature allows users to zoom in and out, or drag a gene around. All the gene IDs in the network module are listed in \"Select by id\" in decreasing order according to gene connectivity. The selected gene ID is appended \"_selected\" for easy identification. By clicking an ID in this list, users can identify the corresponding gene in the network."
ack1 <- "<b>Fund source:</b> NSF (https://www.nsf.gov/) award PGRP-1810468. <p/>"
ack2 <- "<b>Authors:</b> Jianhai Zhang, PhD student at Girke Lab, University of California, Riverside; PI: Prof. Dr. Thomas Girke at University of California, Riverside. <p/>"
ack3 <- "<b>Data and image source:</b> the source of expression data and SVG templates of \"organ_Mustroph\", \"shoot_root_Mustroph\", \"root_roottip_Mustroph\", \"shoot_Mustroph\" was from Mustroph et al. (2009) and Dr. Julia Bailey-Serres Lab at University of California, Riverside, respectively. The exression data of \"root_Geng\" and SVG template was from Geng et al. (2013) and Dr. Julia Bailey-Serres lab at University of California, Riverside, respectively. The expression data and SVG templates of \"brain_Chen\" was from Chen-Plotkin et al. (2008), epilepsyresearch (2017) and anatomybodysystem.com (2017), respectively. The data matrix and SVG template of \"map_Census\" was from Bureau (2018) and Trip8.Co (2018), respectively. <p/>"

ack4 <- "<p/> Mustroph, Angelika, M Eugenia Zanetti, Charles J H Jang, Hans E Holtan, Peter P Repetti, David W Galbraith, Thomas Girke, and Julia Bailey-Serres. 2009. \"Profiling Translatomes of Discrete Cell Populations Resolves Altered Cellular Priorities During Hypoxia in Arabidopsis.\" Proc Natl Acad Sci U S A 106 (44): 18843–8. <br/> Geng, Yu, Rui Wu, Choon Wei Wee, Fei Xie, Xueliang Wei, Penny Mei Yeen Chan, Cliff Tham, Lina Duan, and José R Dinneny. 2013. \"A Spatio-Temporal Understanding of Growth Regulation During the Salt Stress Response in Arabidopsis.\" Plant Cell 25 (6): 2132–54. <br/> Chen-Plotkin, Alice S, Felix Geser, Joshua B Plotkin, Chris M Clark, Linda K Kwong, Wuxing Yuan, Murray Grossman, Vivianna M Van Deerlin, John Q Trojanowski, and Virginia M-Y Lee. 2008. \"Variations in the Progranulin Gene Affect Global Gene Expression in Frontotemporal Lobar Degeneration.\" Hum. Mol. Genet. 17 (10): 1349–62. <br/> epilepsyresearch. 2017. \"The Hippocampus: What Is It?\" https://www.epilepsyresearch.org.uk/the-hippocampus-what-is-it/. <br/> anatomybodysystem.com. 2017. \"HEAD ANATOMY.\" http://anatomybodysystem.com/lateral-view-of-the-brain-labeled/lateral-view-of-the-brain-labeled-brain-diagram-and-label-anatomy-body-list/. <br/> Bureau, U.S. Census. 2018. \"Annual Population Estimates, Estimated Components of Resident Population Change, and Rates of the Components of Resident Population Change for the United States, States, and Puerto Rico: April 1, 2010 to July 1, 2018.\" https://www.census.gov/data/datasets/time-series/demo/popest/2010s-state-total.html. <br/> Trip8.Co. 2018. \"Us Map States.\" http://trip8.co/interactive-map-of-us-states-us-map-states-interactive-us-map/interactive-map-of-us-states-us-map-states-interactive-us-map-new-united-states-map-interactive-detail-color-usa-with-name-inside/."


shinyUI(dashboardPage(

  dashboardHeader(title="spatialHeatmap (updated: 2019-04-25)", 
  titleWidth=380),

  dashboardSidebar(
  
    sidebarMenu(

      menuItem("Input", icon=icon("dashboard"),
      menuSubItem("View", tabName="hm_net"), br(),
      selectInput("fileIn", "Select a work mode", c("None", "organ_Mustroph", "shoot_root_Mustroph", "root_roottip_Mustroph", "shoot_Mustroph", "root_Geng", "brain_Chen", "map_Census", "Compute locally", "Compute online"), "shoot_Mustroph"),
      fileInput("svgInpath", "Step 1: upload an svg file", accept=".svg", multiple=FALSE),
      fileInput("geneInpath", "Step 2: upload a gene expression file", accept=c(".txt", ".csv"), multiple=FALSE),
      radioButtons(inputId='dimName', label='Step 3: is column or row gene?', choices=c("None", "Row", "Column"), selected="None", inline=TRUE),
      selectInput('sep', 'Step 4: separator', c("None", "Tab", "Comma", "Semicolon"), "None"),
      div(style="display:inline-block;width:75%;text-align:left;",textInput("color", "Step 5: Color scheme of Spatial Heatmap", "yellow,purple,blue", placeholder="Eg: yellow,purple,blue", width=200)),
      div(style="display:inline-block;width:25%;text-align:left;", actionButton("col.but", "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
      radioButtons(inputId='cs.v', label='Colour scale based on:', choices=c("Selected genes"="sel.gen", "Whole matrix"="w.mat"), selected="sel.gen", inline=TRUE),
      h4(strong("Compute locally")),
      fileInput("adj.modInpath", "Upload the adjacency matrix and module definition file", accept=".txt", multiple=TRUE),
      h4(strong("Compute online")),
      numericInput("p", "The proportion of samples whose values exceed A (filter genes):", 0, min=0, max=1), numericInput("A", "The value A to be exceeded (filter genes):", 0, min=0, max=10000), 
      numericInput("cv1", "Lower bound of coefficient of variation (CV) (filter genes):", 0, min=0, max=1),
      numericInput("cv2", "Upper bound of coefficient of variation (CV) (filter genes):", 10000, min=0, max=10000),
      numericInput("min.size", "Minmum module size:", 5, min=3, max=10000),
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
      div(style="display:inline-block;width:75%;text-align:left;",textInput("color.net", "Color scheme of Interactive Network", "green,blue,red", placeholder="Eg: yellow,black,blue", width=200)),
      div(style="display:inline-block;width:25%;text-align:left;", actionButton("col.but.net", "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
      radioButtons(inputId="cpt.nw", label="Display or not?", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE),
      #menuSubItem("View network", tabName="network")

      menuSubItem("View", tabName="hm_net"), br()
      ),

      menuItem("Instruction", tabName="instruction", icon=icon("dashboard")),
      menuItem("Acknowledgement", tabName="ack", icon=icon("dashboard"))

     )

  ),

  dashboardBody(
 
    tabItems(

      tabItem(tabName="instruction", 
      box(title="General Instruction", status="primary", solidHeader=TRUE, collapsible=TRUE, p(strong(ins0), align="justify"), p(ins, align="justify"), width=12),
      box(title="Input Instruction", status="primary", solidHeader=TRUE, collapsible=TRUE, p(ins.input1, HTML("<a href=
      https://jianhaizhang.github.io/SVG_tutorial_file/SVG_tutorial.html>here</a>"), ".", align="justify"), p(ins.input2, align="justify"), p(ins.input3.1, strong(ins.input3.2), ins.input3.3, strong(ins.input3.4), ins.input3.5, HTML("&nbsp"), downloadButton("dld.svg", "Download svg file"), downloadButton("dld.data", "Download gene matrix"), align="justify"), p(ins.input4, align="justify"), p(ins.input5, align="justify"), p(ins.input6, align="justify"), width=12), 
      box(title="Spatial Heatmap Instruction", status="primary", solidHeader=TRUE, 
      collapsible=TRUE, p(ins.mat, align="justify"), width=12),
      box(title="Matrix Heatmap & Network Instruction", status="primary", solidHeader=TRUE, collapsible=TRUE, p(ins.mat_hm, align="justify"), p(ins.net, align="justify"), p(ins.net2, align="justify"), width=12)
      ),
      
      tabItem(tabName="hm_net", 

      box(title="Expression Matrix", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput("dt"), "")), height=430, width=12),

      box(title="Spatial Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "7%", "91%", "1%"), "", plotOutput("bar"), plotOutput("tissue"), "")), width=9),

      box(title="Original image", status="primary", solidHeader=TRUE, collapsible=TRUE, splitLayout(cellWidths=c("1%", "98%", "1%"), "", plotOutput("ori.svg"), ""), width=3), br(),

      box(title="Matrix Heatmap", status="primary", solidHeader=TRUE, collapsible=TRUE, plotlyOutput("HMly"), width=12, height=460), br(),

      box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput("bar.net"), visNetworkOutput("vis"), "")), width=12)
      ),

      tabItem(tabName="ack", HTML(ack1), HTML(ack2), HTML(ack3), HTML(ack4))

      )

    )

))
