source('~/tissue_specific_gene/function/fun.R')
options(shiny.maxRequestSize=7*1024^3, stringsAsFactors=FALSE) 

# Import internal functions.
filter_data <- get('filter_data', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
svg_df <- get('svg_df', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
nod_lin <- get('nod_lin', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
grob_list <- get('grob_list', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
col_bar <- get('col_bar', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
lay_shm <- get('lay_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
adj_mod <- get('adj_mod', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
matrix_hm <- get('matrix_hm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
network <- get('network', envir=asNamespace('spatialHeatmap'), inherits=FALSE)






library(SummarizedExperiment); library(shiny); library(shinydashboard); library(grImport); library(rsvg); library(ggplot2); library(DT); library(gridExtra); library(ggdendro); library(WGCNA); library(grid); library(XML); library(plotly); library(data.table); library(genefilter); library(flashClust); library(visNetwork); library(reshape2); library(igraph)

# Import input matrix.
fread.df <- function(input, isRowGene, header, sep, fill, rep.aggr='mean') {
        
  df0 <- fread(input=input, header=header, sep=sep, fill=fill)
  cna <- colnames(df0)[-ncol(df0)]
  df1 <- as.data.frame(df0); rownames(df1) <- df1[, 1]
  # Subsetting identical column names in a matrix will not trigger appending numbers.
  df1 <- as.matrix(df1[, -1]); colnames(df1) <- cna
  if(isRowGene==FALSE) df1 <- t(df1)
  cna <- colnames(df1)
  idx <- grep("__", cna); idx1 <- setdiff(seq_len(length(cna)), idx)
  gene2 <- df1[, idx, drop=FALSE]; gene3 <- df1[, idx1, drop=FALSE]; if (ncol(gene3)>0) colnames(gene3) <- 'ann'
  if(sum(is.na(as.numeric(gene2)))>=1) return('Make sure all values in data matrix are numeric.')
  gen.rep <- gene2; rna <- rownames(gen.rep); gen.rep <-apply(gen.rep, 2, as.numeric); rownames(gen.rep) <- rna
  # Aggregate replicates.
  if (any(duplicated(cna)) & !is.null(rep.aggr)) {

    # To keep colnames, "X" should be a character, not a factor.
    if (rep.aggr=='mean') gene2 <- sapply(X=unique(cna), function(x) rowMeans(gene2[, cna==x, drop=FALSE]))
    if (aggr=='median') {

      gene2 <- sapply(X=unique(cna), function(x) Biobase::rowMedians(gene2[, cna==x, drop=FALSE]))
      rownames(gene2) <- rna

    }

  }; gene2 <-apply(gene2, 2, as.numeric); rownames(gene2) <- rna 
  return(list(gene2=as.data.frame(gene2), gene3=as.data.frame(gene3), gen.rep=as.data.frame(gen.rep)))

}

# enableWGCNAThreads()
inter.svg <- readLines("example/root_cross_final.svg")
inter.data <- read.table("example/expr_arab.txt", header=TRUE, row.names=1, sep="\t")

shinyServer(function(input, output, session) {

  output$dld.svg <- downloadHandler(

    filename=function(){ "root_cross_final.svg"}, 
    content=function(file){ writeLines(inter.svg, file) }

  )

  output$dld.data <- downloadHandler(

    filename=function(){ "expr_arab.txt" },
    content=function(file){ write.table(inter.data, file, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t") }

  )
  # Instruction.
  output$sum <-renderUI({ includeHTML("file/summary.html") })
  output$input <-renderUI({ includeHTML("file/input.html") })
  output$input1 <-renderUI({ includeHTML("file/input1.html") })
  output$shm.ins <-renderUI({ includeHTML("file/shm.html") })
  output$matrix_net <-renderUI({ includeHTML("file/matrix_net.html") })
  # Acknowledgement.
  output$ack <-renderUI({ includeHTML("file/acknowledgement.html") })

  # Filter parameters.
  fil <- reactiveValues(P=0, A=-Inf, CV1=-Inf, CV2=Inf)

  observe({

    input$fileIn; input$geneInpath
    updateRadioButtons(session, inputId="dimName", label="Step 3: is column or row gene?", 
    inline=TRUE, choices=c("None", "Row", "Column"), selected="None")
    updateSelectInput(session, 'sep', 'Step 4: separator', c("None", "Tab", "Comma", "Semicolon"), "None")
    updateRadioButtons(session, inputId='log', label='Data transform', choices=c("No", "log2", "exp2"), selected="No", inline=TRUE)
    updateRadioButtons(session, inputId='cs.v', label='Colour scale based on:', choices=c("Selected genes"="sel.gen", "Whole matrix"="w.mat"), selected="sel.gen", inline=TRUE)
    updateSelectInput(session, "height", "Overall canvas height:", seq(100, 15000, 20), "400")
    updateSelectInput(session, "width", "Overall canvas width:", seq(100, 15000, 20), "820")
    updateSelectInput(session, "col.n", "No. of columns for sub-plots", seq(1, 15, 1), "2")
    updateTextInput(session, inputId="P", label="Filter genes: the proportion (P) of samples whose values exceed A:", value=0, placeholder='a numeric: 0-1')                                                                                    
    updateTextInput(session, inputId="A", label="Filter genes: the value (A) to be exceeded:", value='-Inf', placeholder='a numeric')                                                                                                           
    updateTextInput(session, inputId="CV1", label="Filter genes: lower bound of coefficient of variation (CV1):", value='-Inf', placeholder='a numeric')                                                                                        
    updateTextInput(session, inputId="CV2", label="Filter genes: lower bound of coefficient of variation (CV2):", value='Inf', placeholder='a numeric') 

  })

  observe({

    input$fileIn; input$geneInpath; input$log
    updateTextInput(session, inputId="P", label="Filter genes: the proportion (P) of samples whose values exceed A:", value=0, placeholder='a numeric: 0-1')                                                                                    
    updateTextInput(session, inputId="A", label="Filter genes: the value (A) to be exceeded:", value='-Inf', placeholder='a numeric')                                                                                                           
    updateTextInput(session, inputId="CV1", label="Filter genes: lower bound of coefficient of variation (CV1):", value='-Inf', placeholder='a numeric')                                                                                        
    updateTextInput(session, inputId="CV2", label="Filter genes: lower bound of coefficient of variation (CV2):", value='Inf', placeholder='a numeric') 
    fil$P <- 0; fil$A=-Inf; fil$CV1 <- -Inf; fil$CV2 <- Inf

  })

  # As long as a button is used, observeEvent should be used.
  observeEvent(input$fil.but, { 
                 
    fil$P <- as.numeric(input$P); fil$A <- as.numeric(input$A); fil$CV1 <- as.numeric(input$CV1); fil$CV2 <- as.numeric(input$CV2) 
    if (is.na(fil$P) | !(fil$P>=0) | !(fil$P<=1)) { showModal(modalDialog(title="Filtering", "P should be a numeric between 0 and 1.")); return() }
    if (is.na(fil$A)) { showModal(modalDialog(title="Filtering", "A should be a numeric.")); return() }
    if (is.na(fil$CV1)) { showModal(modalDialog(title="Filtering", "CV1 should be a numeric.")); return() }
    if (is.na(fil$CV2)) { showModal(modalDialog(title="Filtering", "CV2 should be a numeric.")); return() }

  })

  geneIn0 <- reactive({

    if (input$fileIn=="None") return(NULL)  
    withProgress(message="Loading data: ", value = 0, {
    if (grepl("_Mustroph$|_Prudencio$|_Merkin$|_Cardoso.Moreira$|_Census$", input$fileIn)) {

        incProgress(0.5, detail="Loading matrix. Please wait.")
        if (input$fileIn=="brain_Prudencio") { df.te <- fread.df(input="example/expr_human.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE) }
        if (input$fileIn=="mouse_Merkin") df.te <- fread.df(input="example/expr_mouse.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="chicken_Cardoso.Moreira") df.te <- fread.df(input="example/expr_chicken.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="shoot_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="organ_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="root_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="shoot_root_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="root_roottip_Mustroph") df.te <- fread.df(input="example/expr_arab.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="map_Census") df.te <- fread.df(input="example/us_population2018.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        gen.rep <- df.te[['gen.rep']]
        gene2 <- df.te[['gene2']]; if (input$log=='log2') { 
          
          g.min <- min(gene2) 
          if (g.min<0) gene2 <- gene2-g.min+1; if (g.min==0) gene2 <- gene2+1; gene2 <- log2(gene2)  

        }; if (input$log=='exp2') gene2 <- 2^gene2
        gene3 <- df.te[['gene3']][, , drop=FALSE]

    }

    if ((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & 
    (is.null(input$geneInpath)|input$dimName=="None"|input$sep=="None")) return(NULL)
    if ((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & 
    !is.null(input$geneInpath) & input$dimName!="None" & input$sep!="None") {

      incProgress(0.25, detail="Importing matrix. Please wait.")
      geneInpath <- input$geneInpath; if (input$sep=="Tab") sep <- "\t" else if (
      input$sep=="Comma") sep <- "," else if (input$sep=="Semicolon") sep <- ";"
      df.upl <- fread.df(input=geneInpath$datapath, isRowGene=(input$dimName=='Row'), header=TRUE, sep=sep, fill=TRUE, rep.aggr='mean'); gen.rep <- df.upl[['gen.rep']]
      gene2 <- df.upl[['gene2']]; if (input$log=='log2') {
             
        g.min <- min(gene2)
        if (g.min<0) gene2 <- gene2-g.min+1; if (g.min==0) gene2 <- gene2+1; gene2 <- log2(gene2)
      
      }; if (input$log=='exp2') gene2 <- 2^gene2 
      gene3 <- df.upl[['gene3']]

      }; return(list(gene2=gene2, gene3=gene3, gen.rep=gen.rep))

    })

  })

  geneIn <- reactive({
    
    if (is.null(geneIn0())) return(NULL)    
    gene2 <- geneIn0()[['gene2']]; gene3 <- geneIn0()[['gene3']]; input$fil.but
    # Input variables in "isolate" will not triger re-excution, but if the whole reactive object is trigered by "input$fil.but" then code inside "isolate" will re-excute.
    isolate({
      
      se <- SummarizedExperiment(assays=list(expr=as.matrix(gene2)), rowData=gene3)
      if (ncol(gene3)>0) ann.col <- colnames(gene3)[1] else ann.col <- NULL
      se <- filter_data(se=se, ann=ann.col, sam.factor=NULL, con.factor=NULL, pOA=c(fil$P, fil$A), CV=c(fil$CV1, fil$CV2), dir=NULL)
      if (nrow(se)==0) { showModal(modalDialog(title="Filtering", "All genes are filtered out. Please refresh this app to restart.")); return() }
      gene2 <- as.data.frame(assay(se), stringsAsfactors=FALSE); colnames(gene2) <- make.names(colnames(gene2))
      gene3 <- as.data.frame(rowData(se))[, , drop=FALSE]
      
    }); return(list(gene2=gene2, gene3=gene3))

  })

  output$dt <- renderDataTable({

    if (is.null(geneIn())) return()
    if (((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & is.null(geneIn()))|input$fileIn=="None") return(NULL)

    withProgress(message="Data table: ", value = 0, {

      incProgress(0.5, detail="Displaying. Please wait.")
      if (input$fileIn!="None") {

      gene <- geneIn(); gene.dt <- cbind.data.frame(gene[["gene2"]][, , drop=FALSE], gene[["gene3"]][, , drop=FALSE], stringsAsFactors=FALSE) 

   }

    datatable(gene.dt, selection=list(mode="multiple", target="row", selected=c(1)),
    filter="top", extensions='Scroller', options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
    formatRound(colnames(geneIn()[["gene2"]]), 2)

    })

  })

  gID <- reactiveValues(geneID="none", new=NULL, all=NULL)
  observe({ input$geneInpath; input$fileIn; gID$geneID <- "none" })
  observe({ if (is.null(geneIn())) gID$geneID <- "none" })
  observeEvent(input$fileIn, { gID$all <- gID$new <- NULL })

  geneV <- reactive({

    if (is.null(geneIn())) return(NULL)
    if (input$cs.v=="sel.gen" & is.null(input$dt_rows_selected)) return(NULL)
    if (input$fileIn!="None") { if (input$cs.v=="sel.gen" & !is.null(input$dt_rows_selected)) gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
    if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]] } 
    seq(min(gene), max(gene), len=1000) # len must be same with that from the function "spatial_hm()". Otherwise the mapping of a gene value to the colour bar is not accurate. 

  })


  col.sch <- reactive({ 

    if(input$color=="") return(NULL); unlist(strsplit(input$color, ","))

  }); color <- reactiveValues(col="none")

  observeEvent(input$col.but, {

    if (is.null(col.sch())) return (NULL)
    if (input$fileIn!="None") { color$col <- colorRampPalette(col.sch())(length(geneV())) }

  })

  observeEvent(input$dt_rows_selected, {
    
    if (is.null(input$dt_rows_selected)) return()
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    gID$new <- setdiff(gID$geneID, gID$all); gID$all <- c(gID$all, gID$new)

    })

  # To make the "gID$new" and "gID$all" updated with the new "input$fileIn", since the selected row is fixed (3rd row), the "gID$new" is not updated when "input$fileIn" is changed, and the downstream is not updated either. The shoot/root examples use the same data matrix, so the "gID$all" is the same (pre-selected 3rd row) when change from the default "shoot" to others like "organ". As a result, the "gene$new" is null and downstream is not updated. Also the "gene$new" is the same when change from shoot to organ, and downstream is not updated, thus "gene$new" and "gene$all" are both set null above upon new "input$fileIn".
  observeEvent(input$fileIn, {

    if (!grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn)) return()
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    gID$new <- setdiff(gID$geneID, gID$all); gID$all <- c(gID$all, gID$new)

    })

  output$bar <- renderPlot({

    if ((grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn) & !is.null(geneIn()))|((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & !is.null(input$svgInpath) & !is.null(geneIn()))) {

      if (length(color$col=="none")==0|input$color==""|is.null(geneV())) return(NULL)

      if(input$col.but==0) color$col <- colorRampPalette(c("yellow", "purple", "blue"))(length(geneV()))

      withProgress(message="Color scale: ", value = 0, {

        incProgress(0.75, detail="Plotting. Please wait.")
        cs.g <- col_bar(geneV=geneV(), cols=color$col, width=1, mar=c(3, 0.1, 3, 0.1)); return(cs.g)

      })

    }

  })


  observe({

    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    r.na <- rownames(geneIn()[["gene2"]]); gen.sel <- r.na[input$dt_rows_selected]
    updateSelectInput(session, "gen.sel", choices=c("None", gen.sel), selected="None")

  })

  svg.path <- reactive({

    if (input$fileIn=="Compute locally"|input$fileIn=="Compute online") { svg.path <- input$svgInpath$datapath; svg.na <- input$svgInpath$name } else if (input$fileIn=="brain_Prudencio") { svg.path <- "example/homo_sapiens.brain.svg"; svg.na <- "homo_sapiens.brain.svg" } else if (input$fileIn=="mouse_Merkin") { svg.path <- "example/mus_musculus.male.svg"; svg.na <- "mus_musculus.male.svg" } else if (input$fileIn=="chicken_Cardoso.Moreira") { svg.path <- "example/gallus_gallus.svg"; svg.na <- "gallus_gallus.svg" } else if (input$fileIn=="shoot_Mustroph") { svg.path <- "example/shoot_final.svg"; svg.na <- "shoot_final.svg" } else if (input$fileIn=="organ_Mustroph") { svg.path <- "example/organ_final.svg"; svg.na <- "organ_final.svg" } else if (input$fileIn=="root_Mustroph") { svg.path <- "example/root_cross_final.svg"; svg.na <- "root_cross_final.svg" } else if (input$fileIn=="shoot_root_Mustroph") { svg.path <- "example/shoot_root_final.svg"; svg.na <- "shoot_root_final.svg" } else if (input$fileIn=="root_roottip_Mustroph") { svg.path <- "example/root_roottip_final.svg"; svg.na <- "root_roottip_final.svg" } else if (input$fileIn=="map_Census") { svg.path <- "example/us_map_final.svg"; svg.na <- "us_map_final.svg" } else return(NULL)
    return(list(svg.path=svg.path, svg.na=svg.na))

  })

  svg.df <- reactive({ 

    if (((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & 
    !is.null(input$svgInpath))|(grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn) & is.null(input$svgInpath))) {

      withProgress(message="Tissue heatmap: ", value=0, {
    
        incProgress(0.5, detail="Extracting coordinates. Please wait.") 
        df_tis <- svg_df(svg.path=svg.path()[['svg.path']])
        validate(need(!is.character(df_tis), df_tis))
        return(df_tis)

      })

    }

  })

  sam <- reactive({ 

    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\1", c.na) else return(NULL) 

  })

  observe({

    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; svg.df()
    updateCheckboxGroupInput(session, inputId="tis", label='Select tissues to be transparent:', choices=intersect(unique(sam()), unique(svg.df()[['tis.path']])), selected='', inline=TRUE)

  })

  con <- reactive({ 

    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(.*$)", "\\3", c.na) else return(NULL) 

  })

  grob <- reactiveValues(all=NULL); observeEvent(input$fileIn, { grob$all <- NULL })
  gs <- reactive({ 

    if (is.null(svg.df())|is.null(gID$new)) return(NULL); 
    withProgress(message="Tissue heatmap: ", value=0, {

      incProgress(0.25, detail="preparing data.")
      if (input$cs.v=="sel.gen") gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
      if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]]
      g.df <- svg.df()[["df"]]; tis.path <- svg.df()[["tis.path"]]
      grob.lis <- grob_list(gene=gene, geneV=geneV(), coord=g.df, ID=gID$new, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=21) # Only gID$new is used.
      return(grob.lis)

    })

  })

  observe({
    
    input$log; input$tis; input$col.but; input$cs.v # input$tis as an argument in "grob_list" will not cause evaluation of all code, thus it is listed here.
    grob$all <- NULL; gs <- reactive({ 

      if (is.null(svg.df())) return(NULL); if (input$cs.v=="sel.gen" & is.null(input$dt_rows_selected)) return(NULL)
      withProgress(message="Spatial heatmap: ", value=0, {
        incProgress(0.25, detail="preparing data.")

        g.df <- svg.df()[["df"]]; tis.path <- svg.df()[["tis.path"]]
        if (input$cs.v=="sel.gen") gene <- geneIn()[["gene2"]][input$dt_rows_selected, ]
        if (input$cs.v=="w.mat") gene <- geneIn()[["gene2"]]
        g.df <- svg.df()[["df"]]; tis.path <- svg.df()[["tis.path"]]
        grob.lis <- grob_list(gene=gene, geneV=geneV(), coord=g.df, ID=gID$all, cols=color$col, tis.path=tis.path, tis.trans=input$tis, sub.title.size=20) # All gene IDs are used.
        return(grob.lis)

      })

    }); grob$all <- gs()[['grob.lis']]

  })

  observeEvent(gID$new, { grob.all <- c(grob$all, gs()[['grob.lis']]); grob$all <- grob.all[unique(names(grob.all))] })
  # In "observe" and "observeEvent", if one code return (NULL), then all the following code stops. If one code changes, all the code renews.
  observe({
    
    if (is.null(geneIn())|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)

    output$shm <- renderPlot(width=as.numeric(input$width)/2*as.numeric(input$col.n), height=as.numeric(input$height)*length(input$dt_rows_selected), {

    if (is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)
    if (length(color$col=="none")==0|input$color=="") return(NULL)

    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    grob.na <- names(grob$all); con <- unique(con())
    idx <- NULL; for (i in gID$geneID) { idx <- c(idx, grob.na[grob.na %in% paste0(i, '_', con)]) } 
    grob.lis.p <- grob$all[sort(idx)] #grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    lay_shm(lay.shm=input$gen.con, con=con, ncol=input$col.n, ID.sel=gID$geneID, grob.list=grob.lis.p, width=input$width, height=input$height, shiny=TRUE) 

    })

  })

  output$lgd <- renderPlot(width=260, {

    if (is.null(svg.path())|is.null(gs())) return(ggplot())

      svg.path <- svg.path()[['svg.path']]
      # Width and height in original SVG.
      na <- c('width', 'height'); lis <- xmlToList(svg.path)
      for (i in seq_len(length(lis))) {

        if (sum(na %in% names(lis[[i]]))==2) w.h <- as.character(xmlToList(svg.path)[[i]][na])

      }; w.h <- as.numeric(gsub("^(\\d+\\.\\d+|\\d+).*", "\\1", w.h)); r <- w.h[1]/w.h[2]
      g.lgd <- gs()[['g.lgd']]; g.lgd <- g.lgd+coord_fixed(ratio =r); return(g.lgd)

  })

  adj.mod <- reactive({ 

    if (input$fileIn=="Compute locally") {

      name <- input$adj.modInpath$name; path <- input$adj.modInpath$datapath
      path1 <- path[name=="adj.txt"]; path2 <- path[name=="mod.txt"]

      withProgress(message="Loading: ", value = 0, {
        incProgress(0.5, detail="adjacency matrix and module definition.")
        adj <- fread(path1, sep="\t", header=TRUE, fill=TRUE); c.na <- colnames(adj)[-ncol(adj)]
        r.na <- as.data.frame(adj[, 1])[, 1];  adj <- as.data.frame(adj)[, -1] 
        rownames(adj) <- r.na; colnames(adj) <- c.na

        mcol <- fread(path2, sep="\t", header=TRUE, fill=TRUE); c.na <- colnames(mcol)[-ncol(mcol)]
        r.na <- as.data.frame(mcol[, 1])[, 1]; mcol <- as.data.frame(mcol)[, -1] 
        rownames(mcol) <- r.na; colnames(mcol) <- c.na

      }); return(list(adj=adj, mcol=mcol))

    }

  })
  
  adj.tree <- reactive({ 
    #gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    if (input$fileIn=="Compute online"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn)) {

      gene <- geneIn()[["gene2"]]; if (is.null(gene)) return()
      if (input$net.type=="S") { sft <- 12; type <- "signed" } else if (input$net.type=="U") { sft <- 6; type <- "unsigned" }

      withProgress(message="Computing: ", value = 0, {
        incProgress(0.3, detail="adjacency matrix.")
        incProgress(0.5, detail="topological overlap matrix.")
        incProgress(0.1, detail="dynamic tree cutting.")
        se <- SummarizedExperiment(assays=list(expr=as.matrix(gene)))
        adjMod <- adj_mod(se=se, type=type, minSize=input$min.size, dir=NULL)
        adj <- adjMod[['adj']]; mod4 <- adjMod[['mod']]

      }); return(list(adj=adj, mod4=mod4))

    }

  })

  observe({
    
    input$gen.sel; updateSelectInput(session, 'ds', "Select a module splitting sensitivity level", 3:2, selected="3")

  })

  mcol <- reactive({

    if (!is.null(adj.tree()) & (input$fileIn=="Compute online"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn))) { 
      
    withProgress(message="Computing dendrogram:", value=0, {
      incProgress(0.7, detail="hierarchical clustering.")
      mod4 <- adj.tree()[['mod4']]
      
    }); return(mod4) 
    
  }

  })

  observe({

    geneIn(); input$adj.modInpath; input$A; input$p; input$cv1
    input$cv2; input$min.size; input$net.type
    updateSelectInput(session, "mat.scale", "Scale matrix heatmap", c("No", "By column/sample", "By row/gene"), "By row/gene")

  })

  output$HMly <- renderPlotly({
    
    gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's matrix heatmap to another example.
    if (input$fileIn=="Compute locally") { 

      if (is.null(adj.mod())|input$gen.sel=="None") return(NULL); adj <- adj.mod()[[1]]; mods <- adj.mod()[[2]] 

    } else if (input$fileIn=="Compute online"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn)) { 

      if (is.null(adj.tree())|input$gen.sel=="None") return(NULL); adj <- adj.tree()[['adj']]; mods <- mcol() 

    }
    lab <- mods[, input$ds][rownames(gene)==input$gen.sel]
    if (!is.null(lab)) if (lab=="0") { showModal(modalDialog(title="Module", "The selected gene is not assigned to any module. Please select a different one.")); return() } # All the arguments in 'if' statement are evaluated, regardless of the order. Therefore a 'NULL' object should not be compared with others (e.g. >, <) in 'if' statement. But it can be evaluated exclusively in a separate 'if' statement, e.g. 'if (!is.null(lab))'.

    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="Plotting...")
      se <- SummarizedExperiment(assays=list(expr=as.matrix(gene))); mods <- list(mod=mods); if (input$mat.scale=="By column/sample") scale.hm <- 'column' else if (input$mat.scale=="By row/gene") scale.hm <- 'row'
      matrix_hm(geneID=input$gen.sel, se=se, adj.mod=mods, ds=input$ds, scale=scale.hm, main=paste0('Network module containing ', input$gen.sel), title.size=10, static=FALSE)

    })

  })

  observe({
  
    geneIn(); gID$geneID; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    updateSelectInput(session, "TOM.in", label="Input a similarity threshold to display the similarity network.", choices=c("None", sort(seq(0, 1, 0.002), decreasing=TRUE)), selected="None")

  })

  observe({
  
    geneIn(); gID$geneID; input$TOM.in; input$gen.sel; input$ds; input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; input$min.size; input$net.type
    updateRadioButtons(session, inputId="cpt.nw", label="Display or not?", choices=c("Yes"="Y", "No"="N"), selected="N", inline=TRUE)

  })

  col.sch.net <- reactive({ if(input$color.net=="") { return(NULL) }
  unlist(strsplit(input$color.net, ",")) }); color.net <- reactiveValues(col.net="none")

  len.cs.net <- 350
  observeEvent(input$col.but.net, {

    if (is.null(col.sch.net())) return (NULL)

      color.net$col.net <- colorRampPalette(col.sch.net())(len.cs.net)

  })

  visNet <- reactive({

    if (input$TOM.in=="None") return(NULL)
    if (input$fileIn=="Compute locally") { adj <- adj.mod()[[1]]; mods <- adj.mod()[[2]] } else if (input$fileIn=="Compute online"|grepl("_Mustroph$|_Merkin$|_Cardoso.Moreira$|_Prudencio$|_Census$", input$fileIn)) { adj <- adj.tree()[[1]]; mods <- mcol() }
    gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    lab <- mods[, input$ds][rownames(gene)==input$gen.sel]
    if (lab=="0") { showModal(modalDialog(title="Module", "The selected gene is not assigned to any module. Please select a different gene.")); return() }
    idx.m <- mods[, input$ds]==lab; adj.m <- adj[idx.m, idx.m]; gen.na <- colnames(adj.m) 
    idx.sel <- grep(paste0("^", input$gen.sel, "$"), gen.na); gen.na[idx.sel] <- paste0(input$gen.sel, "_selected")
    colnames(adj.m) <- rownames(adj.m) <- gen.na
    withProgress(message="Computing network:", value=0, { 
      incProgress(0.8, detail="making network data frame")
      nod.lin <- nod_lin(ds=input$ds, lab=lab, mods=mods, adj=adj, geneID=input$gen.sel, adj.min=input$TOM.in)
      node <- nod.lin[['node']]; colnames(node) <- c('id', 'value')
      link1 <- nod.lin[['link']]; colnames(link1)[3] <- 'value'
      if (nrow(link1)!=0) { 
        
        link1$title <- link1$value # 'length' is not well indicative of adjacency value, so replaced by 'value'.
        link1$color <- 'lightblue'
        
      }; ann <- geneIn()[[2]]
      if (!is.null(ann)) node <- cbind(node, title=ann[node$id, ], borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      if (is.null(ann)) node <- cbind(node, borderWidth=2, color.border="black", color.highlight.background="orange", color.highlight.border="darkred", color=NA, stringsAsFactors=FALSE)
      net.lis <- list(node=node, link=link1)

    }); net.lis

  })

  output$bar.net <- renderPlot({  

    if (input$TOM.in=="None"|input$cpt.nw=="N") return(NULL)
    if (length(color.net$col.net=="none")==0) return(NULL)
    gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.
    if(input$col.but.net==0) color.net$col.net <- colorRampPalette(c("yellow", "purple", "blue"))(len.cs.net) # color.net$col.net is changed alse outside renderPlot, since it is a reactive value.
   
      withProgress(message="Color scale: ", value = 0, {
      incProgress(0.25, detail="Preparing data. Please wait.")
      incProgress(0.75, detail="Plotting. Please wait.")
      node <- visNet()[["node"]]; node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
        cs.net <- col_bar(geneV=v.net, cols=color.net$col.net, width=1, mar=c(3, 0.1, 3, 0.1)); return(cs.net) # '((max(v.net)-min(v.net))/len.cs.net)*0.7' avoids bar overlap.

      })

  })

  output$edge <- renderUI({ 

    if (input$TOM.in=="None") return(NULL)
    if (input$fileIn=="None"|(input$fileIn=="Your own" & is.null(geneIn()))|
    input$gen.sel=="None") return(NULL)
    HTML(paste0("&nbsp&nbsp&nbsp&nbsp Total edges to display (If > 300, the <br/> 
    &nbsp&nbsp&nbsp App can possibly get stuck.): ", dim((visNet()[["link"]]))[1]))

  })

  vis.net <- reactive({ 
    
    if (input$TOM.in=="None"|input$cpt.nw=="N") return(NULL)
    gene <- geneIn()[["gene2"]]; if (!(input$gen.sel %in% rownames(gene))) return() # Avoid unnecessary computing of 'adj', since 'input$gen.sel' is a cerain true gene id of an irrelevant expression matrix, not 'None', when switching from one defaul example's network to another example.

    withProgress(message="Network:", value=0.5, {
    incProgress(0.3, detail="prepare for plotting.")
    # Match colours with gene connectivity by approximation.
    node <- visNet()[["node"]]; node.v <- node$value; v.net <- seq(min(node.v), max(node.v), len=len.cs.net)
    col.nod <- NULL; for (i in node$value) {

      ab <- abs(i-v.net); col.nod <- c(col.nod, color.net$col.net[which(ab==min(ab))[1]])

    }; node$color <- col.nod
    visNetwork(node, visNet()[["link"]], height="300px", width="100%", background="", main=paste0("Network Module Containing ", input$gen.sel), submain="", footer= "") %>% visIgraphLayout(physics=FALSE, smooth=TRUE) %>% visOptions(highlightNearest=list(enabled=TRUE, hover=TRUE), nodesIdSelection=TRUE)

    })
    
  })

  output$vis <- renderVisNetwork({

    if (input$fileIn=="None"|(input$fileIn=="Your own" & is.null(geneIn()))|input$TOM.in=="None"|input$gen.sel=="None") return(NULL)
    if (input$cpt.nw=="N") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="plotting.")
      vis.net()

    })

  })

})



