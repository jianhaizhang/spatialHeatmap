library(shiny); library(shinydashboard); library(Cairo); library(visNetwork)

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

