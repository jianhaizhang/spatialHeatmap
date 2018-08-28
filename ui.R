library(shiny); library(shinydashboard); library(plotly); library(visNetwork); library(DT)

ins <- "This shiny app has three main functionalities. First, it generates spatial tissue 
       heatmaps where user defined tissue regions are colored by the expression profile 
       of a gene of interest. The gene expression information is uploaded as a table matrix
       and the tissue image is uploaded as an SVG format. Second, the app allows to 
       visualize the expression of a chosen gene in the context of many other genes in the 
       form of a matrix heatmap where the rows and columns have been sorted by pre-computed 
       hierarchical clustering dendrograms. To explore the results, the matrix heatmap has 
       several interactive features. For instance, users can zoom in and out by drawing a 
       rectangle or by double clicking the image, respectively. To identify a query gene 
       in the matrix heatmap, it is tagged by a red rectangle. Third, the app can compute 
       gene network modules across tissue samples, and display them as interactive network 
       graphs." 

ins.input1 <- "For testing purposes, the app includes pre-uploaded SVG image files and gene
              expression tables. They can be viewed by selecting \"Default\" in \"Use 
              default or your own files?\" in the left \"Input\" menu and clicking a gene 
              ID of interest in the matrix table. Subsequently, the spatial tissue heatmap 
              will show the expression of the chosen gene and it is also highlighted in the
              matrix heatmap. Alternatively, users can upload custom SVG image files and 
              expression tables. Details about how to properly format and associate custom 
              SVG images with expression tables are provided here: "

ins.input2 <- "The test SVG image, the expression matrix, and the file of subsetting gene 
              IDs can be downloaded here:"

ins.input3 <- "In the gene matrix, if there are multiple conditions, the sample names MUST 
              be fomatted this way: a sample name is followed by double underscore and 
              then the condition, such as \"S1__con1\". In the sample/condition description,
              only letters, digits, single underscore, dots are allowed. When upload gene 
              expression table, users need to specify the dimension names (Step 3: is 
              column or row gene?) and the seperator. Users can also subset the gene table 
              by uploading a file of gene IDs, which can be separated by comma, space, tab 
              or semicolon."

ins.tis <- "The app displays the gene expression matrix as an interactive table where gene 
           IDs are in the first column. Users can sort the gene values for a sample or 
           search for a particular gene by its ID. Once users click on a specific gene ID, 
           the corresponding spatial tissue heatmap and matrix heatmap will be displayed."

ins.tis1 <- "For the tissue heatmap, users can customize the width/height, image 
            organization, color scheme. In the color scheme, colors must be sepatated by 
            comma, such as \"green,yellow,red\". The \"Scale\" option allows to scale the 
            data by gene or sample. In the matrix heatmap, the chosen gene is shown in the
            context of its related genes. Their dendrogram clustering is exactly exrtacted 
            from the original clustering of the complete gene table rather than 
            re-clustering among themself. The chosen gene is labeled with a red rectangle."

ins.tis2 <- "If the complete gene table contains a whole genome such as 30 thousand genes,
            it may take several mins to display heatmaps for the first time when a gene is 
            chosen in the table, since the app needs to cluster all the genes, which is time
            consuming. From the second time onward, it only takes 6-8 seconds to display
            all the heatmaps."

ins.tis3 <- "The upper limit for matrix heatmap controls the gene group size included. If 
            users are not interested in the matrix heatmap, it can be disabled by selecting 
            0 for the upper limit, and this can seed up the generation of the tissue heatmap." 

ins.net1 <- "The network feature is computationally demanding, especially for large data 
            matrices. If the app times out then users want to refresh the page in their web
            browser. The network is computed among the same group of genes in the matrix 
            heatmap. If no modules are identified, a message will pop up, then users can 
            either increase the upper limit for gene group or change the gene ID."

ins.net2 <- "Gene network modules are computed with the R package WGCNA. In the 
            visualization result each module is assigned a unique color. The similarity 
	    among gene expression profiles is calculated and represented by similarity 
	    values internally. A similarity threshold can be used to restrict the network 
            edges to genes with similar expression profiles. The length of the edges is 
	    proportional to the expression similarity among genes, while the size of the 
	    nodes is proportional to the average expressions of the corresponding genes."

ins.net4 <- "The number of edges to display in the network is shown. It tells how possible 
            the network is going to get stuck. If it exceeds 300, then the network may not 
            be well responsive and more possibly get stuck, in this case, users may not 
            choose \"Yes\"(Compute or not?). To maintain acceptable performance, users are 
            advised to choose a stringent similarity threshold (e.g. 0.9) initially, then 
            decrease the value gradually."

ins.net3 <- "After choosing a similarity threshold, the app will compute the network. The 
            network result can be accessed by clicking \"View network\". \"Select by id\" 
	    allows users to highlight a particular gene/node along with its associated 
	    genes, while \"Select by group\" highlights a particular module."

shinyUI(dashboardPage(

  dashboardHeader(title="Integrated Tissue Heatmap (updated: 2017-12-30)", 
  titleWidth=480),

  dashboardSidebar(
  
    sidebarMenu(

      menuItem("Instruction", tabName="instruction", icon=icon("dashboard")),
      menuItem("Input", icon=icon("dashboard"),
      selectInput("fileIn", "Select a mode", c("None", "Default", 
      "Compute locally", "Compute online"), "Compute locally"),
      fileInput("svgInpath", "Step 1: upload an svg file", accept=".svg", multiple=F),
      fileInput("geneInpath", "Step 2: upload a gene expression file", accept=c(".txt", 
      ".csv"), multiple=F),
      selectInput('dimName', 'Step 3: is column or row gene?', c("None", "Row", 
      "Column")),
      selectInput('sep', 'Step 4: separator', c("None", "Tab", "Comma", "Semicolon"), 
      "None"),
      h4(strong("Compute locally")),
      fileInput("adj.modInpath", "Upload the adjacency matrix and module definition file.",
      accept=".txt", multiple=T),
      h4(strong("Compute online")),
      numericInput("A", "The value A to be exceeded (filter genes):", 0, min=0, max=10000),
      numericInput("p", "The proportion that need to exceed A (filter genes):", 0, min=0, max=1),
      numericInput("cv1", "Lower bound of coefficient of variation (CV) (filter genes):", 
      0, min=0, max=1),
      numericInput("cv2", "Upper bound of coefficient of variation (CV) (filter genes):", 
      10000, min=0, max=10000),
      numericInput("min.size", "Minmum module size:", 15, min=15, max=10000),
      radioButtons("net.type", "Network type", c("Signed"="S", "Unsigned"="U"), "S", inline=T)

      ),

      menuItem("Heatmap & network", icon=icon("dashboard"), 
      h4(strong("Tissue heatmap")),
      selectInput("height", "Overall height of tissue heatmap:", seq(100, 1500, 50), "800",
      width=150), 
      selectInput("width", "Overall width of tissue heatmap:", seq(100, 1500, 50), "800",
      width=150),
      selectInput("col.n", "No. of columns for sub-plots", seq(1, 15, 1), "3", width=150), 
      div(style="display:inline-block;width:65%;text-align:left;",textInput("color", 
      "Color scheme:", "green,blue,purple,yellow,red", placeholder="Eg: green,yellow,red",
      width=150)),
      div(style="display:inline-block;width:35%;text-align:left;", actionButton("col.but", 
      "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
      selectInput("mat.scale", "Scale:", c("No", "By column/sample", "By row/gene"), 
      "No", width=167),
      textOutput("but"),
      radioButtons("gen.con", "Display by:", c("Gene"="gene", "Condition"="con"), "gene", 
      inline=T),
      h4(strong("Matrix heatmap & network")), 
      selectInput("gen.sel","Select a gene to display matrix heatmap & network.", c("None"),
      selected="None"),
      selectInput("ds","Select a module splitting sensitivity level", 0:3, selected="0", 
      width=190),

      selectInput("TOM.in", "Input a similarity threshold to display the similarity 
      network.", c("none", sort(seq(0, 1, 0.002), decreasing=T)), "none"), 
      htmlOutput("edge"),
      radioButtons("cpt.nw", "Compute or not?", c("Yes"="Y", "No"="N"), "N", inline=T),
      #menuSubItem("View network", tabName="network")

      menuSubItem("Display", tabName="hm_net"), br()
      )

     )

  ),

  dashboardBody(
 
    tabItems(

      tabItem(tabName="instruction", 
      box(title="General instruction", status="primary", solidHeader=T, collapsible=T, ins,
      width=12),
      box(title="Input instruction", status="primary", solidHeader=T, collapsible=T, 
      p(ins.input1, HTML("<a href=
      http://biocluster.ucr.edu/~jzhan067/shiny_HM_tutorial/shiny_heatmap_tutorial.html>
      Shiny Heatmap Tutorial</a>"), "."), p(ins.input3), ins.input2, HTML("&nbsp"), 
      downloadButton("dld.svg", "Download svg file"), downloadButton("dld.data", "Download 
      gene matrix"), downloadButton("dld.sub", "Download subsetting IDs"), width=12), 
      box(title="Heatmap instruction", status="primary", solidHeader=T, collapsible=T, 
      p(ins.tis), p(ins.tis1), p(ins.tis2), p(ins.tis3), width=12),
      box(title="Network instruction", status="primary", solidHeader=T, collapsible=T, 
      p(ins.net1), p(ins.net2), p(ins.net4), p(ins.net3), width=12)
      ),
      
      tabItem(tabName="hm_net", 

      box(title="Expr matrix", status="primary", solidHeader=T, collapsible=T, 
      fluidRow(splitLayout(cellWidths=c("1%", "98%", "1%"), "", dataTableOutput("dt"), 
      "")), width=12),

      box(title="Tissue heatmap", status="primary", solidHeader=T, collapsible=T, 
      fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput("bar"), 
      plotOutput("tissue"), "")), width=9),
      
      box(title="Original image", status="primary", solidHeader=T, collapsible=T,
      splitLayout(cellWidths=c("1%", "98%", "1%"), "", plotOutput("ori.svg"), ""), 
      width=3), br(),

      box(title="Matrix heatmap", status="primary", solidHeader=T, collapsible=T,
      plotlyOutput("HMly"), width=12), br(),

      box(title="Interactive network", status="primary", solidHeader=T, collapsible=T, 
      visNetworkOutput("vis"), width=12)
      )

      )

    )

))
