#' Network
#'
#' This function exhibits the input gene in the context of corresponding gene network module, where nodes are genes and edges are adjacencies between genes. The network can be dispayed in static or interactive mode. \cr The gene modules are identified at two alternative sensitivity levels (ds=3 or 2). See function \code{\link{adj_mod}} for details. The thicker edge denotes higher adjacency (coexpression similarity) between genes while larger node indicates higher gene connectivity (sum of a gene's adjacency with all its direct neighbours). \cr In the interactive mode, there is an interactive colour bar to denote gene connectivity. The colour ingredients must only be separated by comma, e.g. "yellow,purple,blue", which means gene connectivity increases from yellow to blue. If too many edges (e.g.: > 300) are displayed, the network could get stuck. So the "Input an adjacency threshold to display the adjacency network." option sets a threthold to filter out weak edges and all remaining edges are displayed. If not too many (e.g.: < 300), users can check "Yes" under "Display or not?", then the network will be displayed and would be responsive smoothly. To maintain acceptable performance, users are advised to choose a stringent threshold (e.g. 0.9) initially, then decrease the value gradually. The interactive feature allows users to zoom in and out, or drag a gene around. All the gene IDs in the network module are listed in "Select by id" in decreasing order according to gene connectivity. The input gene ID is appended "_selected" as a label. By clicking an ID in this list, users can identify the corresponding gene in the network. If the input has gene annotations, then the annotation can be seen by hovering the cursor over a node. \cr The same module can also be displayed in the form of a matrix heatmap with the function \code{\link{matrix_hm}}. 

#' @param geneID A gene ID from the expression matrix. 
#' @param data A "SummarizedExperiment" object containing the processed data matrix and metadata returned by the function \code{\link{filter_data}}.
#' @param ann A character. The column name corresponding to row (gene) annotation in the "rowData" of "se" parameter.
#' @param adj.mod A list of adjacency matrix and module definition returned by the function \code{\link{adj_mod}}.
#' @param ds The module identification sensitivity, either 2 or 3. See function \code{\link{adj_mod}} for details. Used for static network.
#' @param adj.min Minimum adjacency between genes, edges with adjacency below which will be removed. Used for static network.
#' @param con.min Minimun connectivity of a gene, genes with connectivity below which will be removed. Used for static network.
#' @param node.col A vector of colour ingredients for constructing node colour scale in the static image. The default is c("mediumorchid1", "chocolate4").
#' @param edge.col A vector of colour ingredients for constructing edge colour scale in the static image. The default is c("yellow", "blue").
#' @param vertex.label.cex The size of node label. The default is 1.
#' @param vertex.cex The size of node in the static image. The default is 3.
#' @param edge.cex The size of edge in the static image. The default is 10.
#' @param layout The layout of the network in static image, either "circle" or "fr". The "fr" stands for force-directed layout algorithm by Fruchterman and Reingold. The default is "circle".
#' @param main The title in the static image.
#' @param static Either "TRUE" or "FALSE", which displays static or interactive network respectively.
#' @param ... Other arguments passed to the generic function \code{\link[graphics]{plot}}, e.g: "asp=1". 
#' @return A static or interactive network display.

#' @examples

#' # The example data (E-GEOD-67196) is an RNA-seq data measured in cerebellum and frontal cortex of human brain across normal and amyotrophic lateral sclerosis (ALS) subjects (Prudencio et al. 2015). 
#' library(ExpressionAtlas); library(SummarizedExperiment)
#' rse.hum <- getAtlasData('E-GEOD-67196')[[1]][[1]]; assay(rse.hum)[1:3, 1:3]
#'
#' # A targets file describing replicates of samples and conditions is required, which is made based on the "colData" slot in the downloaded "RangedSummarizedExperiment" and available in spatialHeatmap. See the "se" parameter for details. 
#' brain.pa <- system.file('extdata/shinyApp/example/target_brain.txt', package='spatialHeatmap')
#' target.hum <- read.table(brain.pa, header=TRUE, row.names=1, sep='\t')
#' # The "organism_part" and "disease" column describes tissue and condition replicates respectively.  
#' target.hum[c(1:3, 41:42), 4:5]
#' # Place the targets file into "colData" slot as a DataFrame class. 
#' colData(rse.hum) <- DataFrame(target.hum)
#' 
#' # For users with little R expertise, if the gene expression matrix comes as a data frame, it should be placed into "SummarizedExperiment" before proceeding to next step. An example is shown below by borrowing a data matrix from the brain data.
#' # Borrow a data matrix.
#' df <- assay(rse.hum); df[1:2, 1:3]
#' # Place the data matrix and targets file (target.hum) into "SummarizedExperiment".
#' rse.hum <- SummarizedExperiment(assay=df, colData=target.hum, rowData=NULL)
#' 
#' # The count matrix is normalised with estimateSizeFactors (type=‘ratio’).
#' se.nor.hum <- norm_data(data=rse.hum, method.norm='CNF', data.trans='log2')
#'
#' # Average replicates of concatenated sample__condition.
#' se.aggr.hum <- aggr_rep(data=se.nor.hum, sam.factor='organism_part', con.factor='disease', aggr='mean')
#' assay(se.aggr.hum)[49939:49942, ] # The concatenated tissue__conditions are the column names of the output data matrix.
#' 
#' # Genes with low expression level and low variantion are always filtered. 
#' se.fil.hum <- filter_data(data=se.aggr.hum, sam.factor='organism_part', con.factor='disease', pOA=c(0.01, 5), CV=c(0.3, 100), dir=NULL)

#' # Detect modules. 
#' adj.mod <- adj_mod(data=se.fil.hum, type="signed", minSize=15, dir=NULL)
#' # The first column is ds=2 while the second is ds=3. The numbers in each column are module labels with "0" meaning genes not assigned to any modules.
#' adj.mod[['mod']][1:3, ]

#' # Plot network on gene ENSG00000008196 with ds='3'. Set "static=TRUE" to launch the interactive mode. 
#' network(geneID="ENSG00000008196", data=se.fil.hum, ann=NULL, adj.mod=adj.mod, ds="3", adj.min=0.999, con.min=0, vertex.label.cex=1, vertex.cex=0.1, static=TRUE)


#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container. R package version 1.10.1 \cr Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. http://igraph.org \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ \cr Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2018). shiny: Web Application Framework for R. R package version 1.1.0. https://CRAN.R-project.org/package=shiny \cr Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard \cr Almende B.V., Benoit Thieurmel and Titouan Robert (2018). visNetwork: Network Visualization using 'vis.js' Library. R package version 2.0.4. https://CRAN.R-project.org/package=visNetwork
#' Prudencio, Mercedes, Veronique V Belzil, Ranjan Batra, Christian A Ross, Tania F Gendron, Luc J Pregent, Melissa E Murray, et al. 2015. "Distinct Brain Transcriptome Profiles in C9orf72-Associated and Sporadic ALS." Nat. Neurosci. 18 (8): 1175–82
#' Keays, Maria. 2019. ExpressionAtlas: Download Datasets from EMBL-EBI Expression Atlas
#' Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. "Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2." Genome Biology 15 (12): 550. doi:10.1186/s13059-014-0550-8


#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom igraph V E graph_from_data_frame delete_edges delete_vertices as_data_frame layout_in_circle layout_with_fr
#' @importFrom shiny shinyApp shinyUI selectInput htmlOutput div textInput icon actionButton radioButtons fluidRow splitLayout plotOutput shinyServer reactive reactiveValues observeEvent showModal modalDialog withProgress incProgress renderPlot renderUI HTML observe updateSelectInput updateRadioButtons 
#' @importFrom shinydashboard dashboardSidebar dashboardPage dashboardHeader sidebarMenu menuItem menuSubItem dashboardBody tabItems tabItem box
#' @importFrom visNetwork visNetworkOutput visNetwork visOptions renderVisNetwork visIgraphLayout

network <- function(geneID, data, ann, adj.mod, ds="3", adj.min=0, con.min=0, node.col=c("mediumorchid1", "chocolate4"), edge.col=c("yellow", "blue"), vertex.label.cex=1, vertex.cex=3, edge.cex=10, layout="circle", main=NULL, static=TRUE, ...) {

  options(stringsAsFactors=FALSE)
  if (is(data, 'data.frame')|is(se, 'matrix')) {

    data <- as.data.frame(data); rna <- rownames(data); cna <- colnames(data) 
    na <- vapply(seq_len(ncol(data)), function(i) { tryCatch({ as.numeric(data[, i]) }, warning=function(w) { return(rep(NA, nrow(data)))
    }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(data)) )
    na <- as.data.frame(na); rownames(na) <- rna
    idx <- colSums(apply(na, 2, is.na))!=0
    gene <- na[!idx]; colnames(gene) <- cna[!idx]
    if (any(idx)) ann <- data[which(idx)[1]] else ann <- NULL

  } else if (is(data, 'SummarizedExperiment')) { 

    gene <- assay(data); if (!is.null(rowData(data)) & !is.null(ann)) { ann <- rowData(data)[, ann, drop=FALSE]; rownames(gene) <- rownames(ann) <- make.names(rownames(gene)) } else ann <- NULL

  }; from <- to <- width <- size <- NULL 
  adj <- adj.mod[["adj"]]; mods <- adj.mod[["mod"]]
  ds <- as.character(ds); lab <- mods[, ds][rownames(gene)==geneID]
  if (lab=="0") { return('The selected gene is not assigned to any module. Please select a different one') }
  
  if (static==TRUE) { 

  nod.lin <- nod_lin(ds=ds, lab=lab, mods=mods, adj=adj, geneID=geneID, adj.min=adj.min)
  node <- nod.lin[['node']]; link1 <- nod.lin[['link']]
  net <- graph_from_data_frame(d=link1, vertices=node, directed=FALSE)
  # Delete edges and nodes.
  edg.del <- delete_edges(net, E(net)[width <= adj.min])
  net <- delete_vertices(edg.del, igraph::V(edg.del)[size <= con.min])
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
  
  # Network.
  graphics::layout(mat=matrix(c(1, 2, 1, 3), nrow=2, ncol=2), height=c(6, 1))
  par(mar=c(2, 2.5, 2, 2.5), new=FALSE); plot(net, edge.width=igraph::E(net)$width*edge.cex, vertex.size=igraph::V(net)$size*vertex.cex, vertex.label.cex=vertex.label.cex, layout=lay, ...)
  # Node colour bar.
  mat1 <- matrix(v.nod, ncol=1, nrow=length(v.nod)) 
  par(mar=c(2.2, 3, 1, 2.5), new=FALSE); image(x=seq_len(length(v.nod)), y=1, mat1, col=col.nod, xlab="", ylab="", axes=FALSE)
  title(xlab="Node colour scale", line=1, cex.lab=1)
  mtext(text=round(seq(min(node.v), max(node.v), len=5), 1), side=1, line=0.3, at=seq(1, col.len, len=5), las=0, cex=0.8)
  mat2 <- matrix(v.link, ncol=1, nrow=length(v.link)) 
  par(mar=c(2.2, 2.5, 1, 3), new=FALSE); image(x=seq_len(length(v.link)), y=1, mat2, col=col.lin, xlab="", ylab="", axes=FALSE)
  title(xlab="Edge colour scale", line=1, cex.lab=1)
  mtext(text=round(seq(min(link.v), max(link.v), len=5), 1), side=1, line=0.3, at=seq(1, col.len, len=5), las=0, cex=0.8)

  } else if (static==FALSE) {


    ui <- shinyUI(dashboardPage(

      dashboardHeader(),

      dashboardSidebar(
  
      sidebarMenu(

        menuItem("Network", icon=icon("dashboard"), 
        selectInput("ds","Select a module splitting sensitivity level", 3:2, selected="3", width=190),
        selectInput("adj.in", "Input an adjcency threshold to display the adjacency network.", c("None", sort(seq(0, 1, 0.002), decreasing=TRUE)), "None", width=190), 
        htmlOutput("edge"),
        div(style="display:inline-block;width:75%;text-align:left;",textInput("color.net", "Color scheme of Interactive Network", "yellow,purple,blue", placeholder="Eg: yellow,purple,blue", width=200)),
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

    observe({

      input$ds
      updateSelectInput(session, "adj.in", "Input an adjcency threshold to display the adjacency network.", c("None", sort(seq(0, 1, 0.002), decreasing=TRUE)), "None")
      updateRadioButtons(session, "cpt.nw", "Display or not?", c("Yes"="Y", "No"="N"), "N", inline=TRUE)

    })

    observe({

      input$adj.in; updateRadioButtons(session, "cpt.nw", "Display or not?", c("Yes"="Y", "No"="N"), "N", inline=TRUE)

    })

    color_scale <- y <- NULL 
    col.sch.net <- reactive({ if(input$color.net=="") { return(NULL) }
    unlist(strsplit(input$color.net, ",")) }); color.net <- reactiveValues(col.net="none")

    len <- 350
    observeEvent(input$col.but.net, {

      if (is.null(col.sch.net())) return (NULL)
      color.net$col.net <- colorRampPalette(col.sch.net())(len)

    })
    
    visNet <- reactive({

      if (input$adj.in=="None") return(NULL)
      withProgress(message="Computing network:", value=0, {
   
        incProgress(0.8, detail="making network data frame")
        nod.lin <- nod_lin(ds=ds, lab=lab, mods=mods, adj=adj, geneID=geneID, adj.min=input$adj.in)
        node <- nod.lin[['node']]; colnames(node) <- c('id', 'value')
        link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'
        if (nrow(link1)!=0) { 
        
          link1$title <- link1$value # 'length' is not well indicative of adjacency value, so replaced by 'value'.
          link1$color <- 'lightblue'
        
        }; node <- cbind(node, label=node$id, font.size=vertex.label.cex*20, borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
        if (!is.null(ann)) node <- cbind(node, title=ann[node$id, ], stringsAsFactors=FALSE)
        net.lis <- list(node=node, link=link1)

      }); net.lis

    })

    output$bar.net <- renderPlot({  

      if (input$adj.in=="None"|input$cpt.nw=="N") return(NULL)
      if (length(color.net$col.net=="none")==0) return(NULL)
      if(input$col.but.net==0) color.net$col.net <- colorRampPalette(c("yellow", "purple", "blue"))(len)
     
        withProgress(message="Color scale: ", value = 0, {

          incProgress(0.25, detail="Preparing data. Please wait.")
          incProgress(0.75, detail="Plotting. Please wait.")
          node.v <- visNet()[["node"]]$value; v.net <- seq(min(node.v), max(node.v), len=len)
          cs.net <- col_bar(geneV=v.net, cols=color.net$col.net, width=1, mar=c(3, 0.1, 3, 0.1)); return(cs.net)

        })

    })

    output$edge <- renderUI({ 

      if (input$adj.in=="None") return(NULL)
      HTML(paste0("&nbsp&nbsp&nbsp&nbsp Total edges to display (If > 300, the <br/> 
      &nbsp&nbsp&nbsp App can possibly get stuck.): ", nrow((visNet()[["link"]]))))

    })

    vis.net <- reactive({ 

      if (input$adj.in=="None"|input$cpt.nw=="N") return(NULL)
      withProgress(message="Network:", value=0.5, {

        incProgress(0.3, detail="prepare for plotting.")
        # Match colours with gene connectivity by approximation.
        node <- visNet()[["node"]]; node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len)
        col.nod <- NULL; for (i in node$value) {

          ab <- abs(i-v.net); col.nod <- c(col.nod, color.net$col.net[which(ab==min(ab))[1]])

        }; node$color <- col.nod

        visNetwork(node, visNet()[["link"]], height="300px", width="100%", background="", main=paste0("Network Module Containing ", geneID), submain="", footer= "") %>% visIgraphLayout(physics=FALSE, smooth=TRUE) %>%
        visOptions(highlightNearest=list(enabled=TRUE, hover=TRUE), nodesIdSelection=TRUE)

      })
    
    })

    output$vis <- renderVisNetwork({

      if (input$adj.in=="None"|input$cpt.nw=="N") return(NULL)
      withProgress(message="Network:", value=0.5, {

        incProgress(0.3, detail="plotting.")
        vis.net()

      })

    })

  }); shinyApp(ui, server) 

  }

}






