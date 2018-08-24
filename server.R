# Gene value is mapped to color bar in an approximate way.

# The order of paths is not necessarily to be the same with tissues in expr matrix. And 
# multiple polygons can stand for the same tissue. 

# Show tissue heatmap of different conditions.

# Organize multiple condition plots in a single image not through renderUI.

# Append annotation columns to matrix table.

# Display plots for multiple genes. 


library(shiny); library(shinydashboard); library(grImport); library(rsvg); library(ggplot2)
library(DT); library(gridExtra); library(ggdendro); library(WGCNA); library(Cairo)
library(grid); library(XML); library(plotly); library(data.table)

options(shiny.maxRequestSize=50*1024^2) 
# enableWGCNAThreads()
inter.svg <- readLines("interdata/test_final.svg")
inter.data <- read.table("interdata/gene_expr_test.txt", header=T, row.names=1, sep="\t")
inter.sub <- readLines("interdata/subset.txt")

shinyServer(function(input, output, session) {

  output$dld.svg <- downloadHandler(

    filename=function(){ "test_final.svg"}, 
    content=function(file){ writeLines(inter.svg, file) }

  )

  output$dld.data <- downloadHandler(

    filename=function(){ "gene_expr_test.txt" },
    content=function(file){ write.table(inter.data, file, col.names=T, row.names=T, 
    quote=F, sep="\t") }

  )

 output$dld.sub <- downloadHandler(

    filename=function(){ "subset.txt"}, 
    content=function(file){ writeLines(inter.sub, file) }

  )

  observe({

    input$geneInpath
    updateSelectInput(session, "dimName", label="Step 3: is column or row gene?", 
    c("None", "Row", "Column"), "None")
    updateSelectInput(session, 'sep', 'Step 4: separator', c("None", "Tab", "Comma", 
    "Semicolon"), "None")

  })

  geneIn <- reactive({
    
    if (input$fileIn=="None") return(NULL) 
    if (input$fileIn=="Default") { 
      
      withProgress(message="Loading data: ", value = 0, {
 
        incProgress(0.5, detail="Loading matrix. Please wait.")

        if (is.null(input$subset)) {

          if (input$mat.scale=="No") { load("precompute/gene") } 
          else if (input$mat.scale=="By column/gene") { load("precompute/sc.gen")
          gene <- sc.gen } else if (input$mat.scale=="By row/sample") { 
          load("precompute/sc.sam"); gene <- sc.sam }

        } else { 
        
          load("precompute/gene")
          sub <- unlist(strsplit(readLines(input$subset$datapath), ",|;|\\s+|\t"))
          idx <- colnames(gene) %in% sub; gene <- gene[ , idx] 

        }; return(gene)

      }); return(gene)

    }

    if (input$fileIn=="Your own" & (is.null(input$geneInpath)|
    input$dimName=="None")) return(NULL)
    if (input$fileIn=="Your own" & !is.null(input$geneInpath) & input$dimName!="None" &
    input$sep!="None") {

      withProgress(message="Reading data: ", value = 0, {

        incProgress(0.25, detail="Reading matrix. Please wait.")
        geneInpath <- input$geneInpath; if (input$sep=="Tab") sep <- "\t" else if (
        input$sep=="Comma") sep <- "," else if (input$sep=="Semicolon") sep <- ";"
        gene.f <- fread(geneInpath$datapath, header=T, sep=sep, fill=T)
        if (input$dimName=="Row") { 
 
          c.na <- colnames(gene.f)[1:(ncol(gene.f)-1)]
          gene1 <- data.frame(gene.f[, -1], stringsAsFactors=F); idx <- grep("__", c.na)
          idx1 <- setdiff(1:length(c.na), idx)  
          gene2 <- gene1[, idx, drop=F]; gene3 <- gene1[, idx1, drop=F]
          rownames(gene2) <- rownames(gene3) <- as.data.frame(gene.f[, 1], 
          stringsAsFactors=F)[, 1]; colnames(gene2) <- c.na[idx] 
          colnames(gene3) <- c.na[idx1]

        } else if (input$dimName=="Column") {

          col.na <- strsplit(readLines(geneInpath$datapath, n=1), sep)[[1]]
          gene.f1 <- as.data.frame(gene.f); gene <- gene.f1[, -1]
          r.na <- data.frame(gene.f1[, 1])[, 1]
          colnames(gene) <- col.na; rownames(gene) <- r.na
          idx <- grep("__", r.na); gene1 <- t(gene[idx, , drop=F])
          gene2 <- apply(gene1, 2, as.numeric)
          idx1 <- setdiff(1:length(r.na), idx); gene3 <- t(gene[idx1, , drop=F])
          rownames(gene2) <- col.na
          
        }

        if (!is.null(input$subset)) { 
        
          sub <- unlist(strsplit(readLines(input$subset$datapath), ",|;|\\s+|\t"))
          idx <- rownames(gene2) %in% sub; gene2 <- gene2[idx, , drop=F]
          gene3 <- gene3[idx, , drop=F] 

        }; return(list(gene2=gene2, gene3=gene3)) 

      })

    }
  
  })

  geneSub.sc.c <- reactive({

    if (input$fileIn=="None"|input$fileIn=="Default"|is.null(geneIn())) return(NULL) 
    gene <- geneIn()
    if (input$mat.scale=="By column/sample") {

      withProgress(message="Scaling data: ", value = 0, {

        incProgress(0.5, detail="column/sample scaling.")
        gene2 <- scale(gene[["gene2"]])
        return(list(gene2=gene2, gene3=gene[["gene3"]]))
  
      })

    }

  })

  geneSub.sc.r <- reactive({

    if (input$fileIn=="None"|input$fileIn=="Default"|is.null(geneIn())) return(NULL) 
    gene <- geneIn()
    if (input$mat.scale=="By row/gene") {

      withProgress(message="Scaling data: ", value = 0, {

        incProgress(0.25, detail="row/gene scaling.")
        gene2 <- t(scale(t(gene[["gene2"]])))
        return(list(gene2=gene2, gene3=gene[["gene3"]])) 

      })

    }

  })

  geneIn.sc.c <- reactive({

    if (input$fileIn=="None"|input$fileIn=="Default"|is.null(geneIn())) return(NULL) 
    gene <- geneIn()
    if (input$mat.scale=="By column/sample") {

      withProgress(message="Scaling data: ", value = 0, {

        incProgress(0.5, detail="column/sample scaling.")
        gene2 <- scale(gene[["gene2"]])
        return(list(gene2=gene2, gene3=gene[["gene3"]]))
  
     })

    }

  })

  geneIn.sc.r <- reactive({

    if (input$fileIn=="None"|input$fileIn=="Default"|is.null(geneIn())) return(NULL) 
    gene <- geneIn()
    if (input$mat.scale=="By row/gene") {

      withProgress(message="Scaling data: ", value = 0, {

        incProgress(0.25, detail="row/gene scaling.")
        gene2 <- t(scale(t(gene[["gene2"]])))
        return(list(gene2=gene2, gene3=gene[["gene3"]])) 

      })

    }

  })


  output$dt <- renderDataTable({
   
    if ((input$fileIn=="Your own" & is.null(geneIn()))|input$fileIn=="None") return(NULL)
    withProgress(message="Data table: ", value = 0, {
 
      incProgress(0.5, detail="Ploting. Please wait.")
 
    if (input$fileIn=="Default") { 

      if (is.null(input$subset)) gene <- geneIn() else {

        if (input$mat.scale=="By row/gene") gene <- geneSub.sc.r() else if 
        (input$mat.scale=="By column/sample") gene <- geneSub.sc.c() else if 
        (input$mat.scale=="No") gene <- geneIn()

      }

    } else if (input$fileIn=="Your own") {

      if (is.null(input$subset)) {
 
        if (input$mat.scale=="By row/gene") gene <- geneIn.sc.r() else if 
        (input$mat.scale=="By column/sample") gene <- geneIn.sc.c() else if 
        (input$mat.scale=="No") gene <- geneIn()

        gene.dt <- cbind.data.frame(gene[["gene2"]], gene[["gene3"]], 
        stringsAsFactors=F) 

      } else {

        if (input$mat.scale=="By row/gene") gene <- geneSub.sc.r() else if 
        (input$mat.scale=="By column/sample") gene <- geneSub.sc.c() else if 
        (input$mat.scale=="No") gene <- geneIn()

        gene.dt <- cbind.data.frame(gene[["gene2"]], gene[["gene3"]], 
        stringsAsFactors=F) 

      } 

   }
    
    datatable(gene.dt, selection=list(mode="multiple", target="row"),
    filter="top", extensions='Scroller', options=list(autoWidth=T, scrollCollapse=T, 
    deferRender=T, scrollX=T, scrollY=200, scroller=T), class='cell-border strip hover') %>% 
    formatStyle(0, backgroundColor="orange", cursor='pointer') %>% 
    formatRound(colnames(geneIn()[["gene2"]]), 3)

    })
   
  })

  gID <- reactiveValues(geneID="none", new=NULL, all=NULL)
  observe({ input$geneInpath; input$subset; input$fileIn; input$mat.scale; 
  gID$geneID <- "none" })

  geneV <- reactive({
 
    if (is.null(geneIn())) return(NULL)
    if (input$fileIn=="Your own") {

      if (is.null(input$subset)) {
 
        if (input$mat.scale=="By row/gene") gene <- geneIn.sc.r()[["gene2"]] else if 
        (input$mat.scale=="By column/sample") gene <- geneIn.sc.c()[["gene2"]] else if 
        (input$mat.scale=="No") gene <- geneIn()[["gene2"]]

      } else {

        if (input$mat.scale=="By row/gene") gene <- geneSub.sc.r()[["gene2"]] else if 
        (input$mat.scale=="By column/sample") gene <- geneSub.sc.c()[["gene2"]] else if 
        (input$mat.scale=="No") gene <- geneIn()[["gene2"]]

      } 

    } else if (input$fileIn=="Default" & !is.null(input$subset)) {
 
      if (input$mat.scale=="By row/sample") gene <- geneSub.sc.r() else if 
      (input$mat.scale=="By column/gene") gene <- geneSub.sc.c() else if 
      (input$mat.scale=="No") gene <- geneIn() 

    } else if (input$fileIn=="Default" & is.null(input$subset)) { gene <- geneIn() } 

    seq(min(gene), max(gene), len=1000)

  })

  col.sch <- reactive({ if(input$color=="") { return(NULL) }
  unlist(strsplit(input$color, ","))})
  color <- reactiveValues(col="none")

  observeEvent(input$col.but, {

    if (is.null(col.sch())) return (NULL)
    if (input$fileIn=="Your own"|(input$fileIn=="Default" & !is.null(input$subset))) {

      color$col <- colorRampPalette(col.sch())(length(geneV())) 

    }

  })


  observeEvent(input$dt_rows_selected, {
    
    if (is.null(input$dt_rows_selected)) return()
    r.na <- rownames(geneIn()[["gene3"]]); gID$geneID <- r.na[input$dt_rows_selected]
    gID$new <- setdiff(gID$geneID, gID$all)
    
    gID$all <- c(gID$all, gID$new)

    })

  output$bar <- renderPlot({  

    if (input$fileIn=="Default" & is.null(input$subset)) {

      if (input$mat.scale=="No") { load("precompute/cs.g"); return(cs.g) } 
      else if (input$mat.scale=="By column/gene") { load("precompute/cs.scale.gen")
      return(cs.scale.gen) } else if (input$mat.scale=="By row/sample") { 
      load("precompute/cs.scale.sam"); return(cs.scale.sam) }

    } 

    if ((input$fileIn=="Your own" & !is.null(geneIn()) & !is.null(input$svgInpath))|
    (input$fileIn=="Default" & !is.null(input$subset))) {
    if (length(color$col=="none")==0) return(NULL)

    if(input$col.but==0) color$col <- colorRampPalette(c("green", "blue", "purple", 
    "yellow", "red"))(length(geneV()))
     
      withProgress(message="Color scale: ", value = 0, {

        incProgress(0.25, detail="Fetching data. Please wait.")
        cs.df <- data.frame(color_scale=geneV(), y=1)

        incProgress(0.75, detail="Plotting. Please wait.")
        cs.g <- ggplot()+geom_bar(data=cs.df, aes(x=color_scale, y=y), fill=color$col, 
        stat="identity", width=0.2)+theme(axis.title.x=element_blank(), axis.text.x=
        element_blank(), axis.ticks.x=element_blank(), plot.margin=margin(3, 0.1, 3, 0.1, 
        "cm"), panel.grid=element_blank(), panel.background=element_rect(fill="white", 
        colour="grey80"))+coord_flip()+scale_y_continuous(expand=c(0,0))+scale_x_continuous(
        expand = c(0,0)); return(cs.g)

      })

    }

  })

  observe({

    r.na <- rownames(geneIn()[["gene2"]]); gen.sel <- r.na[input$dt_rows_selected]
    updateSelectInput(session, "gen_sel", choices=c("None", gen.sel), selected="None")

  })

  svg.df <- reactive({ 

    if (input$fileIn=="Default") { load("precompute/df.path"); return(df.path) }
    if (input$fileIn=="Your own" & !is.null(input$svgInpath)) {

      withProgress(message="Tissue heatmap: ", value=0, {
    
        incProgress(0.5, detail="Extracting coordinates. Please wait.") 
        svg.path <- input$svgInpath$datapath
        ps.path <- paste0("tmp/", sub(".svg$", ".ps", input$svgInpath$name))
        xml.path <- paste0(ps.path, ".xml"); rsvg_ps(svg.path, ps.path)
        PostScriptTrace(ps.path, xml.path); grml <- xmlParse(xml.path); top <- xmlRoot(grml)

        # Match paths in svg with samples in matrix
        xml <- xmlParse(svg.path); xmltop <- xmlRoot(xml); size <- xmlSize(xmltop)
  
        lis.ma <- xmlSApply(xmltop[[size]], xmlAttrs)
        if (class(lis.ma)=="matrix") { id.xml <- lis.ma["id", ] } 
        else if (class(lis.ma)=="list") {

          id.xml <- NULL
          for (i in 1:length(lis.ma)) { id.xml <- c(id.xml, lis.ma[[i]][["id"]]) }

        }

        xml.na <- NULL
        for (i in 1:xmlSize(xmltop[[size]])) { xml.na <- c(xml.na, 
        xmlName(xmltop[[size]][[i]])) } 

        for (j in 1:length(xml.na)) {

          if (xml.na[j]=="g") {
      
            len.dif <- length(id.xml)-length(xml.na); g.size <- xmlSize(xmltop[[size]][[j]])
            if ((j+1+len.dif) <= length(id.xml)) {

              if (j==1) { 

                id.xml <- c(paste0(id.xml[j+len.dif], "_", 1:g.size), 
                id.xml[(j+1+len.dif):length(id.xml)]) 

              } else if (j>1) {

                id.xml <- c(id.xml[1:(j-1+len.dif)], paste0(id.xml[j+len.dif], "_", 
                1:g.size), id.xml[(j+1+len.dif):length(id.xml)])

              }

            } else if ((j+1+len.dif) >= length(id.xml)) { 

              id.xml <- c(id.xml[1:(j-1+len.dif)], paste0(id.xml[j+len.dif], "_", 1:g.size)) 
 
            }

          }

        }

        tis.path <- gsub("_\\d+$", "", id.xml)

        k <-0; df <- NULL; 
        for (i in 1:((xmlSize(top)-1))) {
  
          if (xmlAttrs(top[[i]])[1]=="fill") {

            k <- k+1
            chil <- xmlChildren(top[[i]])
            xy <- chil[grep("move|line", names(chil))]
            coor <- matrix(NA, nrow=length(xy), ncol=2, dimnames=list(NULL, c("x", "y")))

            for (j in 1:length(xy)) {

              coor[j, "x"] <- as.numeric(xmlAttrs(xy[[j]])["x"])
              coor[j, "y"] <- as.numeric(xmlAttrs(xy[[j]])["y"])

            }

            df0 <- cbind(tissue=id.xml[k], data.frame(coor), stringsAsFactors=T)

          } 

          df <- rbind(df, df0)

        }; return(list(df=df, tis.path=tis.path))

      })

    }

  })


  con <- reactive({ 

    cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
    if (length(grep("__", c.na))>=1) gsub("(.*)(__)(\\w+$)", "\\3", c.na) else 
    return(NULL) 
       
  })

grob <- reactiveValues(all=NULL)

  gs <- reactive({ 

    if (is.null(svg.df())|is.null(gID$new)) return(NULL)
    withProgress(message="Tissue heatmap: ", value=0, {

      incProgress(0.25, detail="preparing data.")

      if (input$fileIn=="Default" & is.null(input$subset)) {

        if (input$mat.scale=="By row/sample") { load("precompute/color.scale.r")
        color$col <- color.scale.r } else if (input$mat.scale=="By column/gene") { 
        load("precompute/color.scale.c"); color$col <- color.scale.c } else if 
        (input$mat.scale=="No") { load("precompute/color"); color$col <- color}
        gene <- geneIn()

      }

      if (input$fileIn=="Your own") {
 
        if (is.null(input$subset)) {
 
          if (input$mat.scale=="By row/gene") gene <- geneIn.sc.r()[["gene2"]] else if 
          (input$mat.scale=="By column/sample") gene <- geneIn.sc.c()[["gene2"]] else if 
          (input$mat.scale=="No") gene <- geneIn()[["gene2"]]

        } else {

          if (input$mat.scale=="By row/gene") gene <- geneSub.sc.r()[["gene2"]] else if 
          (input$mat.scale=="By column/sample") gene <- geneSub.sc.c()[["gene2"]] else if 
          (input$mat.scale=="No") gene <- geneIn()[["gene2"]]

        } 

      } else if (input$fileIn=="Default" & !is.null(input$subset)) {
 
        if (input$mat.scale=="By row/sample") gene <- geneSub.sc.r() else if 
        (input$mat.scale=="By column/gene") gene <- geneSub.sc.c() else if 
        (input$mat.scale=="No") gene <- geneIn() 

     }

      g.df <- svg.df()[["df"]]; tis.path <- svg.df()[["tis.path"]]
      # Assign colors to paths in svg.

      g.list <- function(j) { 

        withProgress(message="Tissue heatmap: ", value=0, {

        incProgress(0.25, detail=paste0("plotting ", j))
        g.col <- NULL; con.idx <- grep(paste0("^", j, "$"), con)
        tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]

      if (k=="none") g.col <- rep("white", length(unique(g.df[, "tissue"]))) else {

        for (i in tis.path) {

          tis.idx <- which(tis.col1 %in% i)
          if (length(tis.idx)==1) { g.col <- c(g.col, scol1[tis.idx]) 
          } else if (length(tis.idx)==0) { g.col <- c(g.col, "white") }

        }

      }

     g <- ggplot()+geom_polygon(data=g.df, aes(x=x, y=y, fill=tissue), color="black")+
     scale_fill_manual(values=g.col, guide=F)+theme(axis.text=element_blank(), 
     axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=
     element_rect(fill="white", colour="grey80"), plot.margin=
     margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), 
     plot.title=element_text(hjust=0.5, size=20))+labs(x="", y="")+
     scale_y_continuous(expand=c(0.01,0.01))+scale_x_continuous(expand=c(0.01,0.01))+
     ggtitle(paste0(k, "_", j)); return(g)
   
      })

    }

      grob.na <- grob.lis <- NULL
      if (!is.null(gID$new)) {

      for (k in gID$new) {

      if (k=="none") scol <- NULL else {

        scol <- NULL
        for (i in gene[k, ]) { 

          ab <- abs(i-geneV()); col.ind <- which(ab==min(ab))[1]
          scol <- c(scol, color$col[col.ind])

        }

      }

      cname <- colnames(geneIn()[["gene2"]]); idx <- grep("__", cname); c.na <- cname[idx]
      tis.col <- gsub("(.*)(__)(\\w+$)", "\\1", c.na); g.lis <- NULL

     if (!is.null(con())) {

     con <- con(); con.uni <- unique(con)
     grob.na0 <- paste0(k, "_", con.uni)
     g.lis <- lapply(con.uni, g.list); grob <- lapply(g.lis, ggplotGrob)
     names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob) 

    } else {

      withProgress(message="Tissue heatmap: ", value=0, {

        incProgress(0.25, detail="plotting ")
        g.col <- NULL
        if (k=="none") g.col <- rep("white", length(unique(g.df[, "tissue"]))) 
        else {

        for (i in tis.path) {

          tis.idx <- which(tis.col %in% i)
          if (length(tis.idx)==1) { g.col <- c(g.col, scol[tis.idx]) 
          } else if (length(tis.idx)==0) { g.col <- c(g.col, "white") }

        }

      }

     g <- ggplot()+geom_polygon(data=g.df, aes(x=x, y=y, fill=tissue), color="black")+
     scale_fill_manual(values=g.col, guide=F)+theme(axis.text=element_blank(), 
     axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=
     element_rect(fill="white", colour="grey80"), plot.margin=
     margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), 
     plot.title=element_text(hjust=0.5, size=20))+labs(x="", y="")+
     scale_y_continuous(expand=c(0.01,0.01))+scale_x_continuous(expand=c(0.01,0.01))+
     ggtitle(paste0(k, "_Tissue Heatmap")); g.grob <- ggplotGrob(g); return(list(g.grob))
   
      })

   }

  }; return(grob.lis)

 }

    })

  })

  observeEvent(gID$new, { grob$all <- c(grob$all, gs()) })

  observe({

  output$tissue <- renderPlot(width=as.numeric(input$width), 
                   height=as.numeric(input$height), {

    if (is.null(svg.df())|gID$geneID[1]=="none") return(NULL)
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    idx <- NULL; for (i in gID$geneID) idx <- c(idx, grep(i, names(grob$all)))
    grob.lis.p <- grob$all[idx]

    if (input$gen.con=="gene") {

      all.cell <- ceiling(length(unique(con()))/as.numeric(input$col.n)) * 
      as.numeric(input$col.n)
      cell.idx <- c(1:length(unique(con())), rep(NA, all.cell-length(unique(con()))))
      m <- matrix(cell.idx, ncol=as.numeric(input$col.n), byrow=T)
      lay <- NULL
      for (i in 1:length(gID$geneID)) { lay <- rbind(lay, m+(i-1)*length(unique(con()))) }
      grid.arrange(grobs=grob.lis.p, layout_matrix=lay)

    } else if (input$gen.con=="con") {
     
      grob.p.na <- names(grob.lis.p); grob.p.idx <- NULL
      for (i in unique(con())) {

        grob.p.idx <- c(grob.p.idx, grep(paste0(i, "$"), grob.p.na))

      } 

      grob.lis.p.con <- grob.lis.p[grob.p.idx]
      all.cell <- ceiling(length(gID$geneID)/as.numeric(input$col.n)) * 
      as.numeric(input$col.n)
      cell.idx <- c(1:length(gID$geneID), rep(NA, all.cell-length(gID$geneID)))
      m <- matrix(cell.idx, ncol=as.numeric(input$col.n), byrow=T)
      lay <- NULL
      for (i in 1:length(unique(con()))) { lay <- rbind(lay, m+(i-1)*length(gID$geneID)) }
      grid.arrange(grobs=grob.lis.p.con, layout_matrix=lay)
      
    }

    do.call(file.remove, list(list.files(".", "capture.*.ps")))
    do.call(file.remove, list(list.files("./tmp", ".ps$", full.name=T)))
    do.call(file.remove, list(list.files("./tmp", ".ps.xml$", full.name=T)))

  })

  })

  output$ori.svg <- renderImage({

    if ((is.null(input$svgInpath) & input$fileIn!="Default")|input$fileIn==
    "None") return(list(src="precompute/blank.png", contentType="image/png"))

    w <- as.numeric(input$width); h <- as.numeric(input$height); con.n <- length(con())
    W <- w/as.numeric(input$col.n); H <- h/(ceiling(con.n/as.numeric(input$col.n)))
    if (input$fileIn=="Default") {

      list(src="precompute/default.png", contentType="image/png", width=W, height=H*3.7, 
      alt=NULL)

    } else if (input$fileIn=="Your own"|!is.null(input$svgInpath)) {

       svg.path <- input$svgInpath$datapath; rsvg_png(svg.path, "tmp/user.png")
       list(src="tmp/user.png", contentType="image/png", width=W, height=H*3.5, 
       alt=NULL)

    }

  }, deleteFile=T)


  gene.dend <- reactive({

    if (is.null(geneIn())) return(NULL)

    withProgress(message="Matrix heatmap: ", value = 0, {

      incProgress(0.2, detail="Computing dendrogram. Please wait.")

      if (input$fileIn=="Default" & is.null(input$subset)) {
   
        gene <- geneIn(); load("precompute/dd.col"); load("precompute/dd.row")
        return(list(gene=gene, dd.col=dd.col, dd.row=dd.row))
 
      } else if (input$fileIn=="Your own"|(input$fileIn=="Default" & 
      !is.null(input$subset))) {

        geneIn <- geneIn()$gene2; gene.mat <- as.matrix(geneIn)
        dd.row <- as.dendrogram(hclust(dist(gene.mat)))
    
        incProgress(0.4, detail="Computing dendrogram. Please wait.")
        dd.col <- as.dendrogram(hclust(dist(t(gene.mat))))

        if (input$fileIn=="Your own") {
 
          if (input$mat.scale=="By row/sample") gene <- geneIn.sc.r() else if 
          (input$mat.scale=="By column/gene") gene <- geneIn.sc.c() else if 
          (input$mat.scale=="No") gene <- geneIn() 

        } else if (input$fileIn=="Default" & !is.null(input$subset)) {

          if (input$mat.scale=="By row/sample") gene <- geneSub.sc.r() else if 
          (input$mat.scale=="By column/gene") gene <- geneSub.sc.c() else if 
          (input$mat.scale=="No") gene <- geneIn() 

        }; return(list(gene=gene, dd.col=dd.col, dd.row=dd.row))

     } 

    })

  })


  output$HMly <- renderPlotly({

    if (input$gen_sel=="None"|as.numeric(input$clus.lim)==0) return(NULL)
    gene <- gene.dend()$gene[[1]]; dd.col <- gene.dend()$dd.col; dd.row <- gene.dend()$dd.row
    lab.len <- length(labels(dd.row))

    if (lab.len > as.numeric(input$clus.lim)) {

      while(lab.len > as.numeric(input$clus.lim)) {

        if (length(grep(input$gen_sel, labels(dd.row[[1]])))) {

          dd.row <- dd.row[[1]]; lab.len <- length(labels(dd.row)) 

        } else if (length(grep(input$gen_sel, labels(dd.row[[2]])))) { 

          dd.row <- dd.row[[2]]; lab.len <- length(labels(dd.row)) 

        }

      }

    }

    dx <- dendro_data(dd.row); dy <- dendro_data(dd.col)
  
    g.dengra <- function(df) {
    
      ggplot()+geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend))+labs(x = "", 
      y = "") + theme_minimal()+ theme(axis.text = element_blank(), axis.ticks=
      element_blank(), panel.grid = element_blank())
  
    }

    px <- g.dengra(dx$segments); py <- g.dengra(dy$segments)+coord_flip()
  
    df.row <- data.frame(x=1:length(labels(dd.row)), y=0, lab=labels(dd.row))
    gen.idx <- grep(input$gen_sel, labels(dd.row))
    df.rec <- data.frame(x1=gen.idx-0.5, x2=gen.idx+0.5, y1=-1, y2=8)
    df.col <- data.frame(x=1:length(labels(dd.col)), y=0, lab=labels(dd.col)) 
    px1 <- px+geom_text(data=df.row, aes(x=x, y=y, label=lab), position=position_dodge(0.9),
    vjust=1, hjust=0, size=2, angle=90)+geom_rect(data=df.rec, aes(xmin=x1, xmax=x2, 
    ymin=y1, ymax=y2), fill=NA, color="red", size=1, alpha=1)
    py1 <- py+geom_text(data=df.col, aes(x=x, y=y, label=lab), position=position_dodge(0.9),    vjust=1, hjust=0, size=2, angle=0) 

    col.ord <- order.dendrogram(dd.col); row.ord <- order.dendrogram(dd.row)
 
    gene.clus <- rbind(Y=0, cbind(X=0, gene[row.ord, col.ord]))

    ply <- plot_ly(z=as.matrix(t(gene.clus)), type="heatmap") %>% layout(yaxis=
    list(domain=c(0, 1), showticklabels=F, showgrid=F, ticks="", zeroline=F), 
    xaxis=list(domain=c(0, 1), showticklabels=F, showgrid=F, ticks="", zeroline=F))

    subplot(px1, plot_ly(), ply, py1, nrows=2, shareX=T, shareY=T, margin=0, heights=
    c(0.2, 0.8), widths=c(0.8, 0.2))

  })


  observe({
  
    geneIn(); gID$geneID; geneIn.sc.r(); geneIn.sc.c()
    updateSelectInput(session, "TOM.in", label="Input a similarity threshold to display the 
    similarity network.", choices=c("none", sort(seq(0, 1, 0.002), decreasing=T)), 
    selected="none")

  })

  observe({
  
    geneIn(); gID$geneID; input$TOM.in; geneIn.sc.r(); geneIn.sc.c()
    updateRadioButtons(session, "cpt.nw", label="Compute or not?", c("Yes"="Y", "No"="N"),
    "N", inline=T, selected="N")

  })

  gID.idx <- reactive({ if (input$gen_sel=="None") return(0) else return(1) })
  tom <- reactiveValues(idx=0)
  observeEvent(input$TOM.in, {

    if (tom$idx==1) return(NULL)
    if (input$TOM.in!="none") tom$idx <- 1

  })

  observeEvent(input$gen_sel, { tom$idx <- 0 })

  TOM <- reactive({

    if (is.null(geneIn())) return(NULL)
    if (gID.idx()==0) return(NULL); if (tom$idx==0) return(NULL)
    geneID <- input$gen_sel
    withProgress(message="Computing network:", value = 0, {
 
      incProgress(0.8, detail="computing Topological Overlap Matrix")

      if (input$fileIn=="Default" & is.null(input$subset)) { 

        load("precompute/dd.row"); load("precompute/gene")

      }

      if (input$fileIn=="Your own" & !is.null(geneIn())|(input$fileIn=="Default" & 
      !is.null(input$subset))) { dd.row <- gene.dend()$dd.row; gene <- geneIn()[[1]] }

      lab.len <- length(labels(dd.row)); clus.lim <- as.numeric(input$clus.lim)

      if (lab.len > clus.lim) {

        while(lab.len > clus.lim) {

	  if (length(grep(geneID, labels(dd.row[[1]])))) {

	    dd.row <- dd.row[[1]]; lab.len <- length(labels(dd.row)) 

	  } else if (length(grep(geneID, labels(dd.row[[2]])))) { 

	    dd.row <- dd.row[[2]]; lab.len <- length(labels(dd.row)) 

	  }

	}

      }

      g.mat <- t(gene)[, labels(dd.row)]
      sft <- pickSoftThreshold(g.mat, dataIsExpr=T, 
      powerVector=c(seq(1, 10, by=1), seq(12, 20, by=2)))
      adj <- adjacency(g.mat, power=6); TOM <- TOMsimilarity(adj)
      colnames(TOM) <- rownames(TOM) <- colnames(g.mat)
    
    }); list(TOM=TOM, g.mat=g.mat)

  })


  dy.col <- reactive({

    if (is.null(TOM())) return
    TOM <- TOM()$TOM
    withProgress(message="Computing network:", value=0, {

      incProgress(0.8, detail="extracting colors.")
      dissTOM <- 1-TOM
      geneTree <- hclust(as.dist(dissTOM), method = "average") 
      mod <- cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=2, pamRespectsDendro=F,
      minClusterSize=30)
      dy.col <- labels2colors(mod)
     
      return(dy.col)

      })

    })

  observe({

    if (!is.null(TOM()) & input$TOM.in!="none") { 
   
      if (tom$idx==0) return 
      dy.col <- dy.col()
      if (unique(dy.col)=="grey") { showModal(modalDialog(title="Network",
      "No modules are identified. Please increase the \"upper limit\" for displaying genes
      in matrix heatmap or change to another gene.")) }
    
    }

  })

   merge <- reactive({
 
     if (is.null(TOM())|unique(dy.col())=="grey") return
     dynamicColors <- dy.col(); g.mat <- TOM()$g.mat
     withProgress(message="Computing network:", value=0, {

      incProgress(0.8, detail="merging similar modules.")
      MEList <- moduleEigengenes(g.mat, colors = dynamicColors)
      MEDiss <- 1-cor(MEList$eigengenes)
      mer <- mergeCloseModules(g.mat, dynamicColors, cutHeight=0.25, verbose=3)

    }); mer

  })

  

  visNet <- reactive({

    if (input$TOM.in=="none"|unique(dy.col())=="grey") return(NULL)
    if (input$TOM.in=="none") return(NULL)
    TOM <- TOM()$TOM; g.mat <- TOM()$g.mat

    withProgress(message="Computing network:", value=0, {
   
      incProgress(0.8, detail="making network data frame")
      tom <- TOM; idx = tom > input$TOM.in
      link <- data.frame(from=rownames(tom)[row(tom)[idx]], 
      to=colnames(tom)[col(tom)[idx]], length=tom[idx])
      # Should not exclude duplicate rows by "length".
      link1 <- subset(link, length!=1 & !duplicated(link[, "length"]))

      grp <- integer(length(merge()[["colors"]]))
      for (i in unique(merge()[["colors"]])) { grp[merge()[["colors"]]==i] <- i } 
      node <- data.frame(id=colnames(tom), group=paste0("Module_", grp), 
      value=colMeans(g.mat), color=grp, stringsAsFactors=F)	
      node[, "group"][grep("grey", node[,"group"])] <- "unassigned genes" 
      net.lis <- list(node=node, link=link1)

    }); net.lis

  })


  output$edge <- renderUI({ 

    if (input$TOM.in=="none") return(NULL)
    if (input$fileIn=="None"|(input$fileIn=="Your own" & is.null(geneIn()))|
    input$gen_sel=="None") return(NULL)
    if (unique(dy.col())=="grey") return(NULL)

    HTML(paste0("&nbsp&nbsp&nbsp&nbsp Edges (If > 300, the App can get <br/> &nbsp&nbsp
    &nbsp stuck easily.): ", dim((visNet()[["link"]]))[1]))

  })


  vis.net <- reactive({ 

    if (input$TOM.in=="none"|input$cpt.nw=="N") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="prepare for plotting.")
      visNetwork(visNet()[["node"]], visNet()[["link"]]) %>% 
      visOptions(highlightNearest=T, nodesIdSelection=T, selectedBy="group")

    })
    
  })

  output$vis <- renderVisNetwork({
  
    if (input$fileIn=="None"|(input$fileIn=="Your own" & is.null(geneIn()))|input$TOM.in
    =="none"|input$gen_sel=="None") return(NULL)
    if (input$cpt.nw=="N"|unique(dy.col())=="grey") return(NULL)

    withProgress(message="Network:", value=0.5, {

      incProgress(0.3, detail="plotting.")
      vis.net()

    })

  })

})



