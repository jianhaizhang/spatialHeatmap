#' Network
#'
#' The gene network modules are identified at two alternative sensitivities levels (3, 2). From 3 to 2, the sensitivity decreases and results in less modules with larger sizes. The "Select a module splitting sensitivity level" option allows users to choose which level to use for displaying the iteractive network. \cr Nodes and edges mean genes and adjacency between genes respectively. There is an interactive colour bar to denote gene connectivity (sum of a gene's adjacency with all its direct neighbours). The colour ingredients must only be separated by comma, e.g. the default are "blue,green,red", which means gene connectivity increases from blue to red. The edge length is inversely proportional to gene adjacency. If too many edges (e.g.: > 300) are displayed in this network, the app can possibly get stuck. So the "Input an adjacency threshold to display the adjacency network." option sets a threthold to filter out some weak edges. Only edges above the threshold are displayed in the network. The app outputs the total number of remaining edges resulting from each input adjacency threshold. If it is not too large (e.g.: < 300), users can check "Yes" under "Display or not?", then the network can be responsive smoothly. To maintain acceptable performance, users are advised to choose a stringent threshold (e.g. 0.9) initially, then decrease the value gradually. The interactive feature allows users to zoom in and out, or drag a gene around. All the gene IDs in the network module are listed in "Select by id" in decreasing order according to gene connectivity. The selected gene ID is appended "_selected", which can be easily identified from the list. By clicking an ID in this list, users can identify the corresponding gene in the network. If the input expression matrix has an annotation column, then the annotation can be seen by hovering the cursor over a node. \cr The same module can also be displayed in the form of an interactive matrix heatmap with the function "matrix.heatmap". \cr The network can seen by clicking the "View" button at the bottom of the side menu.

#' @param geneID A gene ID from the expression matrix. 
#' @param expr.table The gene expression table, where rows are gene IDs and columns are samples. If annotation is included, it must be in the last column.

#' @param adj.mod The "adjacency matrix" and "modules" definitions resulting from the function "adj.mod".

#' @return An interactive network interface lauched on the web browser. 

#' @examples

#' data.path <- system.file("extdata/example", "root_expr_ann_row_gen.txt", package = "spatialHeatmap")
#' exp <- filter.data(data=data.path, sep="\t", isRowGen=TRUE, c(0, 0), c(0.1, 10000), "processed_data")
#' adj_mod <- adj.mod(data=exp, type="signed", minSize=15)
#' \donttest{ network("PSAC", exp, adj_mod) }

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2018). shiny: Web Application Framework for R. R package version 1.1.0. https://CRAN.R-project.org/package=shiny \cr Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard \cr Almende B.V., Benoit Thieurmel and Titouan Robert (2018). visNetwork: Network Visualization using 'vis.js' Library. R package version 2.0.4. https://CRAN.R-project.org/package=visNetwork

#' @export
#' @importFrom shinydashboard dashboardSidebar
#' @importFrom shiny shinyApp shinyUI selectInput htmlOutput div textInput icon actionButton radioButtons fluidRow splitLayout plotOutput shinyServer reactive reactiveValues observeEvent showModal modalDialog withProgress incProgress renderPlot renderUI HTML 
#' @importFrom shinydashboard dashboardPage dashboardHeader sidebarMenu menuItem menuSubItem dashboardBody tabItems tabItem box
#' @importFrom visNetwork visNetworkOutput visNetwork visOptions renderVisNetwork


network <- function(geneID, expr.table, adj.mod) {

ui <- shinyUI(dashboardPage(

  dashboardHeader(),

  dashboardSidebar(
  
    sidebarMenu(

      menuItem("Network", icon=icon("dashboard"), 
      selectInput("ds","Select a module splitting sensitivity level", 3:2, selected="3", width=190),
      selectInput("TOM.in", "Input an adjcency threshold to display the adjacency network.", c("None", sort(seq(0, 1, 0.002), decreasing=TRUE)), "None", width=190), 
      htmlOutput("edge"),
      div(style="display:inline-block;width:75%;text-align:left;",textInput("color.net", "Color scheme of Interactive Network", "green,blue,red", placeholder="Eg: green,yellow,red", width=200)),
      div(style="display:inline-block;width:25%;text-align:left;", actionButton("col.but.net", "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
      radioButtons("cpt.nw", "Display or not?", c("Yes"="Y", "No"="N"), "N", inline=TRUE),
      menuSubItem("View", tabName="net")

      )

     )

  ),

  dashboardBody(
 
    tabItems(

      tabItem(tabName="net", 
      box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput("bar.net"), visNetworkOutput("vis"), "")), width=12)
      )

      )

    )

))

server <- shinyServer(function(input, output, session) {

  color_scale <- y <- NULL 
  col.sch.net <- reactive({ if(input$color.net=="") { return(NULL) }
  unlist(strsplit(input$color.net, ",")) }); color.net <- reactiveValues(col.net="none")

  len.cs.net <- 500
  observeEvent(input$col.but.net, {

    if (is.null(col.sch.net())) return (NULL)

      color.net$col.net <- colorRampPalette(col.sch.net())(len.cs.net)

  })

  visNet <- reactive({

    if (input$TOM.in=="None") return(NULL)

    adj <- adj.mod[["adj"]]; mods <- adj.mod[["mod"]]; col.len <- ncol(expr.table)
  if (!is.numeric(expr.table[, col.len])) { gene <- expr.table[, -col.len]; ann <- expr.table[, col.len, drop=FALSE] } else { gene <- expr.table; ann <- NULL }
  lab <- mods[, input$ds][rownames(gene)==geneID]

    if (lab=="0") { showModal(modalDialog(title="Module", "The selected gene is not assigned to any module. Please select a different gene.")); return() }
    idx.m <- mods[, input$ds]==lab; adj.m <- adj[idx.m, idx.m]
    withProgress(message="Computing network:", value=0, {
   
      incProgress(0.8, detail="making network data frame")
      idx = adj.m > as.numeric(input$TOM.in)
      link <- data.frame(from=rownames(adj.m)[row(adj.m)[idx]], 
      to=colnames(adj.m)[col(adj.m)[idx]], length=adj.m[idx])
      # Should not exclude duplicate rows by "length".
      link1 <- subset(link, length!=1 & !duplicated(link[, "length"]))
      node <- data.frame(id=colnames(adj.m),  
      value=colSums(adj.m), title=ann[colnames(adj.m), ], borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      
      idx.sel <- grep(paste0("^", geneID, "$"), node$id)
      rownames(node)[idx.sel] <- node$id[idx.sel] <- paste0(geneID, "_selected")

      # Match colours with gene connectivity by approximation.
      node <- node[order(-node$value), ]; col <- color.net$col.net; col.nod <- NULL
      node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
      for (i in node$value) {

        ab <- abs(i-v.net); col.nod <- c(col.nod, col[which(ab==min(ab))[1]])

      }; node$color <- col.nod; net.lis <- list(node=node, link=link1, v.net=v.net)

    }); net.lis

  })

  output$bar.net <- renderPlot({  

    if (input$TOM.in=="None"|input$cpt.nw=="N") return(NULL)
    if (length(color.net$col.net=="none")==0) return(NULL)
    v.net <- visNet()[["v.net"]]
    if(input$col.but.net==0) color.net$col.net <- colorRampPalette(c("green", "blue", "red"))(length(v.net))
     
      withProgress(message="Color scale: ", value = 0, {

        incProgress(0.25, detail="Preparing data. Please wait.")
        cs.df.net <- data.frame(color_scale=v.net, y=1); #save(cs.df, file="cs.df")
        incProgress(0.75, detail="Plotting. Please wait.")
        cs.g.net <- ggplot()+geom_bar(data=cs.df.net, aes(x=color_scale, y=y), fill=color.net$col.net, stat="identity", width=0.2)+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=margin(3, 0.1, 3, 0.1, "cm"), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"))+coord_flip()+scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand = c(0,0)); return(cs.g.net)

      })

  })

  output$edge <- renderUI({ 

    if (input$TOM.in=="None") return(NULL)
    HTML(paste0("&nbsp&nbsp&nbsp&nbsp Total edges to display (If > 300, the <br/> 
    &nbsp&nbsp&nbsp App can possibly get stuck.): ", dim((visNet()[["link"]]))[1]))

  })

  vis.net <- reactive({ 

    if (input$TOM.in=="None"|input$cpt.nw=="N") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="prepare for plotting.")
      visNetwork(visNet()[["node"]], visNet()[["link"]], height="300px", width="100%", background="", main=paste0("Network Module Containing ", geneID), submain="", footer= "") %>% 
      visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE)

    })
    
  })

  output$vis <- renderVisNetwork({

    if (input$TOM.in=="None"|input$cpt.nw=="N") return(NULL)
    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="plotting.")
      vis.net()

    })

  })

})

shinyApp(ui, server) 

}







