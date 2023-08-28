#' Network graphs
#'
#' @description
#' This function presents a network module returned by \code{\link{adj_mod}} in form of a static or interactive network graph. In the network graph, nodes are biomolecules (genes, proteins, etc) and edges are adjacencies (coexpression similarities) between biomolecules. The thicker edge denotes higher adjacency between nodes while larger node indicates higher connectivity (sum of a node's adjacencies with all its direct neighbours).
#'
#' In the interactive mode, there is an interactive color scale to denote node connectivity. The color ingredients can only be separated by comma, semicolon, single space, dot, hypen, or underscore. \emph{E.g.} `yellow,orange,red`, which means node connectivity increases from yellow to red. If too many edges (\emph{e.g.}: > 500) are displayed, the app may get crashed, depending on the computer RAM. So the "Adjacency threshold" option sets a threthold to filter out weak edges. Meanwhile, the "Maximun edges" limits the total edges to show. In case a very low adjacency threshold is chosen and introduces too many edges that exceed the "Maximun edges", the App will internally increase the adjacency threshold until the edge total is within the "Maximun edges". The adjacency threshold of 1 leads to no edges. In this case the App wil internally decrease this threshold until the number of edges reaches the "Maximun edges". If adjacency threshold of 0.998 is selected and no edge remains, this App will also internally update the edges to 1 or 2. To maintain acceptable performance, users are advised to choose a stringent threshold (\emph{e.g.} 0.9) initially, then decrease the value gradually. The interactive feature allows users to zoom in and out, or drag a node around. All the node IDs in the network module are listed in "Select by id" in decreasing order according to node connectivity. The ID of a chosen target biomolecule is appended "_target" as a label. By clicking an ID in this list, users can identify the corresponding node in the network. If the input data matrix has annotations for biomolecules, then the annotation can be seen by hovering the cursor over a node. 

#' @inheritParams matrix_hm
#' @param ID A biomolecule of interest.
#' @param adj.mod A two-component list returned by \code{\link{adj_mod}}, consisting of the adjacency matrix and module assignments.
#' @param ds One of `0`, `1`, `2`, or `3` (default), the module identification sensitivity level. The smaller `ds` indicates larger but less modules while the larger `ds` denotes smaller but more modules. See function \code{\link{adj_mod}} for details.
#' @param adj.min A cutoff of adjacency between nodes. Edges with adjacency below this cutoff will not be shown. Default is 0. Applicable to static network.
#' @param con.min A cutoff of node connectivity. Nodes with connectivity below which will not be shown. Default is 0. Applicable to static network. 
#' @param desc A column name in the \code{rowData} slot of \code{SummarizedExperiment}. If provided, when mousing over the nodes in the interactive network, the corresponding description will show up.
#' @param node.col A vector of color ingredients for node color scale in the static image. The default is `c("turquoise", "violet")`, where node connectivity increases from "turquoise" to "violet".
#' @param edge.col A vector of color ingredients for edge color scale in the static image. The default is `c("yellow", "blue")`, where edge adjacency increases from "yellow" to "blue".
#' @param vertex.label.cex The size of node label in the static and interactive networks. The default is 1.
#' @param vertex.cex The size of nodes in the static image. The default is 3.
#' @param edge.cex The size of edges in the static image. The default is 10.
#' @param layout The layout of the network in static images, either `circle` (default) or `fr`. The `fr` stands for force-directed layout algorithm developed by Fruchterman and Reingold.
#' @param color.key.lab.size,color.key.text.size The size of color key label and text respectively.
#' @param main The title in the static image. Default is NULL.
#' @param static Logical. `TRUE` and `FALSE` return a static and interactive network graphs respectively.  
#' @param ... Other arguments passed to the generic function \code{\link[graphics]{plot.default}}, \emph{e.g.}: \code{asp=1}. 
#' @return A static or interactive network graph.

#' @inherit submatrix examples

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1
#' Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. http://igraph.org
#' R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/
#' Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2018). shiny: Web Application Framework for R. R package version 1.1.0. https://CRAN.R-project.org/package=shiny
#' Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard
#' Almende B.V., Benoit Thieurmel and Titouan Robert (2018). visNetwork: Network Visualization using 'vis.js' Library. R package version 2.0.4. https://CRAN.R-project.org/package=visNetwork
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. “Gene Expression Across Mammalian Organ Development.” Nature 571 (7766): 505–9

#' @export
#' @importFrom SummarizedExperiment assays
#' @importFrom igraph V E graph_from_data_frame delete_edges delete_vertices as_data_frame layout_in_circle layout_with_fr
#' @importFrom shiny shinyApp shinyUI selectInput htmlOutput div textInput icon actionButton radioButtons fluidRow splitLayout plotOutput shinyServer reactive reactiveValues observeEvent withProgress incProgress renderPlot renderUI HTML observe updateSelectInput updateRadioButtons numericInput validate need span tags
#' @importFrom shinydashboard dashboardSidebar dashboardPage dashboardHeader sidebarMenu menuItem menuSubItem dashboardBody tabItems tabItem box

network <- function(ID, data, assay.na=NULL, adj.mod, ds="3", adj.min=0, con.min=0, desc=NULL, node.col=c("turquoise", "violet"), edge.col=c("yellow", "blue"), vertex.label.cex=1, vertex.cex=3, edge.cex=10, layout="circle", color.key.lab.size=1.5, color.key.text.size=1, main=NULL, static=TRUE, ...) {
 # save(ID, data, assay.na, adj.mod, ds, adj.min, con.min, node.col, edge.col, vertex.label.cex, vertex.cex, edge.cex, layout, color.key.lab.size, color.key.text.size, main, static, file='network.arg')
  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(data, 'matrix')|is(data, 'DFrame')|is(data, 'dgCMatrix')) {
    dat.lis <- check_data(data=data, assay.na=assay.na); gene <- dat.lis$dat; ann <- dat.lis$row.meta
    if (ncol(ann)>0) { 
      if ('desc' %in% colnames(ann)) ann <- ann[, 'desc', drop=FALSE] 
    } else ann <- NULL
  } else if (is(data, 'SummarizedExperiment') | is(data, 'SingleCellExperiment')) { 
    if (is.null(assay.na)) {
      if (length(assays(data)) > 1) stop("Please specify which assay to use by assigning the assay name to 'assay.na'!") else if (length(assays(data)) == 1) assay.na <- 1
    }; gene <- assays(data)[[assay.na]]
    if (!is.null(rowData(data)) & !is.null(desc)) { ann <- rowData(data)[, desc, drop=FALSE]; rownames(gene) <- rownames(ann) <- make.names(rownames(gene)) } else ann <- NULL
  } else { stop('Accepted data classes are "data.frame", "matrix", "DFrame", "dgCMatrix", "SummarizedExperiment", or "SummarizedExperiment", except that "spatial_hm" also accepts a "vector".') } 
  from <- to <- width <- size <- NULL 
  adj <- adj.mod[["adj"]]; mods <- adj.mod[["mod"]]
  if (length(ID)!=1) return('Only one ID is required!')
  if (!ID %in% rownames(gene)) return('ID is not in data!')
  ds <- as.character(ds); lab <- mods[, ds][rownames(gene)==ID]
  msg <- 'The selected ID is not assigned to any module. Please select a different one!'
  if (!check_obj(lab)) {message(msg); return(msg) }
  if (lab=="0") { message(msg); return(msg) }
  
  if (static==TRUE) { 

  nod.lin <- nod_lin(ds=ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adj.min)
  node <- nod.lin[['node']]; link1 <- nod.lin[['link']]
  if (nrow(link1)==0) stop('All edges are filtered out, please reduce \'adj.min\'!')
  net <- graph_from_data_frame(d=link1, vertices=node, directed=FALSE)
  # Delete edges and nodes.
  edg.del <- delete_edges(net, E(net)[width <= adj.min])
  net <- delete_vertices(edg.del, igraph::V(edg.del)[size <= con.min])
  if (length(V(net))==0) stop('All nodes are filtered out, please reduce \'con.min\'!')
  # Remaining nodes and edges.
  node <- as_data_frame(net, what="vertices"); link1 <- as_data_frame(net, what="edges")
  
  # Match colours with gene connectivity by approximation.
  col.len <- 350; col.nod <- colorRampPalette(node.col)(col.len)
  # Assign node colours.
  col.node <- NULL; node.v <- node$size; v.nod <- seq(min(node.v), max(node.v), len=col.len)
      for (i in node$size) {

        ab <- abs(i-v.nod); col.node <- c(col.node, col.nod[which(ab==min(ab))[1]])

      }; igraph::V(net)$color <- col.node
  # Match colours with gene adjacency by approximation.
  col.lin <- colorRampPalette(edge.col)(col.len)
  # Assign edge colours.
  col.link <- NULL; link.v <- link1$width; v.link <- seq(min(link.v), max(link.v), len=col.len)
      for (i in link1$width) {
        ab <- abs(i-v.link); col.link <- c(col.link, col.lin[which(ab==min(ab))[1]])
      }; igraph::E(net)$color <- col.link
  if (layout=="circle") lay <- layout_in_circle(net); if (layout=="fr") lay <- layout_with_fr(net)
  if (vertex.label.cex==0) igraph::V(net)$name <- NA 
  # Network.
  graphics::layout(mat=matrix(c(1, 2, 1, 3), nrow=2, ncol=2), height=c(6, 1))
  par(mar=c(2, 2.5, 2, 2.5), new=FALSE); plot(net, edge.width=igraph::E(net)$width*edge.cex, vertex.size=igraph::V(net)$size*vertex.cex, vertex.label.cex=vertex.label.cex, layout=lay, ...)
  # Node colour bar.
  mat1 <- matrix(v.nod, ncol=1, nrow=length(v.nod)) 
  par(mar=c(2.2, 3, 1, 2.5), new=FALSE); image(x=seq_len(length(v.nod)), y=1, mat1, col=col.nod, xlab="", ylab="", axes=FALSE)
  title(xlab="Node colour scale", line=1, cex.lab=color.key.lab.size)
  mtext(text=round(seq(min(node.v), max(node.v), len=5), 1), side=3, line=0.3, at=seq(1, col.len, len=5), las=0, cex=color.key.text.size)
  mat2 <- matrix(v.link, ncol=1, nrow=length(v.link)) 
  par(mar=c(2.2, 2.5, 1, 3), new=FALSE); image(x=seq_len(length(v.link)), y=1, mat2, col=col.lin, xlab="", ylab="", axes=FALSE)
  title(xlab="Edge colour scale", line=1, cex.lab=color.key.lab.size)
  mtext(text=round(seq(min(link.v), max(link.v), len=5), 1), side=3, line=0.3, at=seq(1, col.len, len=5), las=0, cex=color.key.text.size)

  node.df <- node
  colnames(node.df) <- c('id', 'connectivity')
  rownames(node.df) <- sub('_target$', '', rownames(node.df))
  node.df <- node.df[, 'connectivity', drop=FALSE]; invisible(node.df)
  } else if (static==FALSE) {
    pkg <- check_pkg('visNetwork'); if (is(pkg, 'character')) { warning(pkg); return(pkg) }
    ui <- shinyUI(dashboardPage(
      dashboardHeader(title='Interactive Network'),
      dashboardSidebar(
      sidebarMenu(
        menuItem("Network", icon=icon("list"), 
        div(style="display:inline-block;width:75%;text-align:left;",textInput("color.net", "Color scheme:", "yellow,orange,red", placeholder="Eg: yellow,orange,red", width=200)),
        div(style="display:inline-block;width:25%;text-align:left;", actionButton("col.but.net", "Go", icon=icon("refresh"), style="padding:7px; font-size:90%; margin-left: 0px")),
        selectInput("ds","Module splitting sensitivity level:", 3:2, selected="3", width=190),
        selectInput("adj.in", "Adjacency threshold (the smaller, the more edges):", sort(seq(0, 1, 0.002), decreasing=TRUE), 1, width=190),
        tags$span(style="color:yellow", numericInput(inputId="max.edg", label="Maximum edges (too many edges may crash the app):", value=10, min=1, max=500, width=200)),
        htmlOutput("edge"),
        menuSubItem("View graph", tabName="net")
        ),
        menuItem("Instruction", tabName="ins", icon=icon("info"))

     )

     ),

     dashboardBody(
       
       tags$head(tags$style('.shiny-output-error-validation { color: red }')),
       tabItems(

        tabItem(tabName="net", 
        box(title="Interactive Network", status="primary", solidHeader=TRUE, collapsible=TRUE, fluidRow(splitLayout(cellWidths=c("1%", "6%", "91%", "2%"), "", plotOutput("bar.net"), visNetwork::visNetworkOutput("vis"), "")), width=12)
        ),
        tabItem(tabName="ins", 
        box(title=NULL, status="primary", solidHeader=TRUE, collapsible=TRUE, HTML('By default, only top edges are shown, since too many edges might crash the session. To display more edges, the adjcency threshold should be decreased gradually. <br/> <br/> The "Maximun edges" limits the total of shown edges. In case a very low adjacency threshold is choosen and introduces too many edges that exceed the Maximun edges, the app will internally increase the adjacency threshold until the edge total is within the Maximun edges, which is a protection against too many edges. The adjacency threshold of 1 produces no edges, in this case the app wil internally decrease this threshold until the number of edges reaches the Maximun edges. If adjacency threshold of 0.998 is selected and no edge is left, this app will also internally update the edges to 1 or 2.'), width=12),
        )
        )

     )

  ))

  server <- shinyServer(function(input, output, session) {

    color_scale <- y <- NULL 
    col.sch.net <- reactive({ if(input$color.net=="") { return(NULL) }
    unlist(strsplit(input$color.net, ",")) }); color.net <- reactiveValues(col.net="none")

    len <- 350
    observeEvent(input$col.but.net, {

      if (is.null(col.sch.net())) return (NULL)
      color.net$col.net <- colorRampPalette(col.sch.net())(len)

    })
    
    visNet <- reactive({

      withProgress(message="Computing network:", value=0, {
   
        incProgress(0.8, detail="making network data frame")
        adjs <- 1; lin <- 0; adj.lin.vec <- NULL
        validate(need(try(as.integer(input$max.edg)==input$max.edg), 'The number of edges should be an integer!'))
        # Compute the min adj.
        while (lin<input$max.edg) {
          
          adjs <- adjs-0.002; if (adjs<=10^-15) adjs <- 0
          nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adjs)
          lin <- nrow(nod.lin[['link']])
          vec0 <- adjs; names(vec0) <- lin
          adj.lin.vec <- c(adj.lin.vec, vec0)
          if (adjs==0) break

        }
        # The first version of links computed from the min adj or the input adj, which is subjected to the following check.
        nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=ifelse(input$adj.in==1, adjs, input$adj.in))

        link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'
        # If the links are 0 due to the input adj, change the "adjs" to the value bringing 1 or 2 links.
        lins <- NULL; if (nrow(link1)==0) {

          adjs <- adj.lin.vec[names(adj.lin.vec)>=1][1]
          nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adjs)
          link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)
 
        } else if (nrow(link1)>input$max.edg) {
       
          # If the links are larger than max links due to the input adj, change the "adjs" to the value producing max links.
          adjs <- adjs+0.002
          nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adjs)
          link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)

        } else if (nrow(link1)<=input$max.edg & input$adj.in!=1) {
        
          # If 0<link total<max links, use the input adj.
          adjs <- input$adj.in
          nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=ID, adj.min=adjs)
          link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'; lins <- nrow(link1)
 
        }; node <- nod.lin[['node']]; colnames(node) <- c('id', 'value')
        if (nrow(link1)!=0) { 
        
          link1$title <- link1$value # 'length' is not well indicative of adjacency value, so replaced by 'value'.
          link1$color <- 'lightblue'
        
        }; node <- cbind(node, label=node$id, font.size=vertex.label.cex*20, borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
        if (!is.null(ann)) node <- cbind(node, title=sub('_target$', '', ann[node$id, ]), stringsAsFactors=FALSE)
        net.lis <- list(node=node, link=link1, adjs=adjs, lins=lins)

      }); net.lis

    })

    observe({
      input$ds; lins <- visNet()[["lins"]]
      if (input$adj.in==1|is.null(lins)|is.numeric(lins)) updateSelectInput(session, "adj.in", "Adjacency threshold (the  smaller, the more edges):", sort(seq(0, 1, 0.002), decreasing=TRUE), as.numeric(visNet()[["adjs"]]))
    })

    output$bar.net <- renderPlot({  

      if (length(color.net$col.net=="none")==0) return(NULL)
      if(input$col.but.net==0) color.net$col.net <- colorRampPalette(c('yellow', 'orange', 'red'))(len)
     
        withProgress(message="Color scale: ", value = 0, {

          incProgress(0.25, detail="Preparing data. Please wait.")
          incProgress(0.75, detail="Plotting. Please wait.")
          node.v <- visNet()[["node"]]$value; v.net <- seq(min(node.v), max(node.v), len=len)
          cs.net <- col_bar(geneV=v.net, cols=color.net$col.net, width=1); return(cs.net)

        })

    })

    output$edge <- renderUI({ 
      span(style="color:yellow;font-weight:bold;", HTML(paste0("Remaining edges for display: ", nrow((visNet()[["link"]])))))
    })

    vis.net <- reactive({ 
      withProgress(message="Network:", value=0.5, {
        incProgress(0.3, detail="prepare for plotting.")
        # Match colours with gene connectivity by approximation.
        node <- visNet()[["node"]]; node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len)
        col.nod <- NULL; for (i in node$value) {
          ab <- abs(i-v.net); col.nod <- c(col.nod, color.net$col.net[which(ab==min(ab))[1]])
        }; node$color <- col.nod
        visNetwork::visNetwork(node, visNet()[["link"]], height="300px", width="100%", background="", main=paste0("Network Module Containing ", ID), submain="", footer= "") %>% visNetwork::visIgraphLayout(physics=FALSE, smooth=TRUE) %>%
        visNetwork::visOptions(highlightNearest=list(enabled=TRUE, hover=TRUE), nodesIdSelection=TRUE)
      })
    })

    output$vis <- visNetwork::renderVisNetwork({

      withProgress(message="Network:", value=0.5, {

        incProgress(0.3, detail="plotting.")
        vis.net()

      })

    })

  }); shinyApp(ui, server) 

  }

}






