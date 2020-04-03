source('~/tissue_specific_gene/function/fun.R')
options(shiny.maxRequestSize=7*1024^3, stringsAsFactors=FALSE) 

# Import internal functions.
# svg_df <- get('svg_df', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
# nod_lin <- get('nod_lin', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
# grob_list <- get('grob_list', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
# col_bar <- get('col_bar', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
# lay_shm <- get('lay_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

filter_data <- function(data, pOA=c(0, 0), CV=c(-Inf, Inf), ann=NULL, samples, conditions,  dir=NULL) {

    if (!is.null(dir)) { path <- paste0(dir, "/local_mode_result/"); if (!dir.exists(path)) dir.create(path) }
    df <- assay(data); col.met <- as.data.frame(colData(data), stringsAsFactors=FALSE)
    if (!is.null(samples) & !is.null(conditions)) { colnames(df) <- paste(col.met[, samples], col.met[, conditions], sep='__') }
    ffun <- filterfun(pOverA(pOA[1], pOA[2]), cv(CV[1], CV[2]))
    filtered <- genefilter(df, ffun); df <- df[filtered, ]
    row.met <- as.data.frame(rowData(data), stringsAsFactors=FALSE)[filtered, , drop=FALSE]

    df1 <- NULL; if (!is.null(dir) & !is.null(ann) & ncol(row.met)>0) { 

      df1 <- cbind.data.frame(df, row.met[ann, ], stringsAsFactors=FALSE)
      colnames(df1)[ncol(df1)] <- ann

    }

    if (!is.null(dir)) {
      
      if (!is.null(df1)) write.table(df1, paste0(path, "processed_data.txt"), sep="\t", row.names=TRUE, col.names=TRUE) else write.table(df, paste0(path, "processed_data.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
      
    }; rownames(col.met) <- NULL
    expr <- SummarizedExperiment(assays=list(expr=df), rowData=row.met, colData=col.met); return(expr)

}


svg_df <- function(svg.path) {

  # Make sure the style is correct. If the stroke width is not the same across polygons such as '0.0002px', '0.216px', some stroke outlines cannot be recognised by 'PostScriptTrace'. Then some polygons are missing. Since the ggplot is based on 'stroke' not 'fill'.
  tmp <- tempdir()
  xmlfile <- xmlParse(svg.path); xmltop <- xmlRoot(xmlfile); ply <- xmltop[[xmlSize(xmltop)]]
  style <- 'stroke:#000000;stroke-width:5.216;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1' # 'fill' is not necessary. In Inkscape, "group" or move an object adds transforms (relative positions), and this can lead to related polygons uncolored in the spatial heatmaps. Solution: ungroup and regroup to get rid of transforms and get absolute positions.
  # Change 'style' of all polygons.
  for (i in seq_len(xmlSize(ply))) {                      
         
    addAttributes(ply[[i]], style=style) 
    if (xmlSize(ply[[i]])>=1) for (j in seq_len(xmlSize(ply[[i]]))) { addAttributes(ply[[i]][[j]], style=style) }
        
  }; svg.inter <- paste0(tmp, '/internal.svg'); saveXML(doc=xmlfile, file=svg.inter)
  # SVG file conversion. 
  rsvg_ps(svg.inter, file=sub("svg$", "ps", svg.inter))
  p1 <- sub("svg$", "ps", svg.inter); p2 <- paste0(sub("svg$", "ps", svg.inter), ".xml")  
  if (length(grep("~", svg.inter))) {

    wd1 <- getwd(); setwd("~"); hm <- getwd(); setwd(wd1)
    p1 <- sub("~", hm, p1); p2 <- sub("~", hm, p2)

  }; PostScriptTrace(p1, p2) 
  grml <- xmlParse(p2); top <- xmlRoot(grml) # Use internal svg to get coordinates.
  xml <- xmlParse(svg.path); xmltop <- xmlRoot(xml); size <- xmlSize(xmltop) # Use original not internal svg to get path ids. Otherwise, errors can come up.
  do.call(file.remove, list(svg.inter, p1, p2))

  lis.ma <- xmlApply(xmltop[[size]], xmlAttrs)
  if (is(lis.ma, "matrix")) { id.xml <- lis.ma["id", ] } else if (is(lis.ma, "list")) {

    id.xml <- NULL; for (i in seq_len(length(lis.ma))) { id.xml <- c(id.xml, lis.ma[[i]][["id"]]) }

  }

  xml.na <- NULL; for (i in seq_len(xmlSize(xmltop[[size]]))) { xml.na <- c(xml.na, xmlName(xmltop[[size]][[i]])) } 

  for (j in seq_len(length(xml.na))) {

    if (xml.na[j]=="g") {

      len.dif <- length(id.xml)-length(xml.na); g.size <- xmlSize(xmltop[[size]][[j]])
      if ((j+1+len.dif) <= length(id.xml)) {

        if (j==1) { 

          id.xml <- c(paste0(id.xml[j+len.dif], "_", seq_len(g.size)), id.xml[(j+1+len.dif):length(id.xml)]) 

        } else if (j>1) {

          id.xml <- c(id.xml[seq_len(j-1+len.dif)], paste0(id.xml[j+len.dif], "_", 
          seq_len(g.size)), id.xml[(j+1+len.dif):length(id.xml)])

       } } else if ((j+1+len.dif) >= length(id.xml)) { 

        id.xml <- c(id.xml[seq_len(j-1+len.dif)], paste0(id.xml[j+len.dif], "_", seq_len(g.size))) 

       }

    }

  }; tis.path <- gsub("_\\d+$", "", id.xml)

  # Detect groups that use relative coordinates ("transform", "matrix" in Inkscape.), which leads to some plygons missed in ".ps.xml" file.
  fil.stk <- sapply(seq_len(xmlSize(top)-1), function (i) xmlAttrs(top[[i]])['type']); tab <- table(fil.stk)
  w <- which(fil.stk=='fill')%%2==0
  if (any(w)) { 

    # All path ids in original SVG.
    id.svg <- NULL; for (i in seq_len(xmlSize(xmltop[[size]]))) {

     node <- xmltop[[size]][[i]] 
     if (xmlName(node)=='g') for (j in seq_len(xmlSize(node))) { id.svg <- c(id.svg, xmlAttrs(node[[j]])[['id']]) } else id.svg <- c(id.svg, xmlAttrs(node)[['id']])  

    }
    
    # Index of wrong path.
    w1 <- which(w)[1]
    # Wrong path and related group. 
    tis.wrg <- paste0(id.svg[c(w1-1, w1)], tis.wrg <- paste0(" (", tis.path[c(w1-1, w1)], ")"), collapse='; ')
    mg1 <- paste0("Error detected in ",  "'", tis.wrg, "'", " in SVG image", collapse='')
    mg2 <- "Possible solutions: 1. If they belong to a group, ungroup and regroup it; 2. If they are tiny polygons, remove them."
    return(paste0(mg1, '. ', mg2)) 

  }

  k <-0; df <- NULL; for (i in seq_len(xmlSize(top)-1)) {

    if (xmlAttrs(top[[i]])['type']=='stroke') {

      k <- k+1; chil <- xmlChildren(top[[i]]); xy <- chil[grep("move|line", names(chil))]
      coor <- matrix(NA, nrow=length(xy), ncol=2, dimnames=list(NULL, c("x", "y")))

      for (j in seq_len(length(xy))) {

        coor[j, "x"] <- as.numeric(xmlAttrs(xy[[j]])["x"])
        coor[j, "y"] <- as.numeric(xmlAttrs(xy[[j]])["y"])

      }

      df0 <- cbind(tissue=id.xml[k], data.frame(coor, stringsAsFactors=FALSE), stringsAsFactors=TRUE)
      df <- rbind(df, df0, stringsAsFactors=TRUE)

     }

    }; g.df <- df; return(list(df=df, tis.path=tis.path))

}


lay_shm <- function(lay.shm, con, ncol, ID.sel, grob.list, width, height, shiny) {

    width <- as.numeric(width); height <- as.numeric(height); ncol <- as.numeric(ncol); con <- unique(con)
    if (lay.shm=="gene") {

        all.cell <- ceiling(length(con)/ncol)*ncol
        cell.idx <- c(seq_len(length(con)), rep(NA, all.cell-length(con)))
        m <- matrix(cell.idx, ncol=as.numeric(ncol), byrow=TRUE)
        lay <- NULL; for (i in seq_len(length(ID.sel))) { lay <- rbind(lay, m+(i-1)*length(con)) }
      if (shiny==TRUE & length(grob.list)>=1) return(grid.arrange(grobs=grob.list, layout_matrix=lay, newpage=TRUE))
       
      g.tr <- lapply(grob.list[seq_len(length(grob.list))], grobTree)
      n.col <- ncol(lay); n.row <- nrow(lay)
      g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc"))

    } else if (lay.shm=="con") {

      grob.all.na <- names(grob.list); na.rev <- NULL
      # Reverse the "gene_condition" names, and re-oder them.
      for (i in grob.all.na) { na.rev <- c(na.rev, paste0(rev(strsplit(i, NULL)[[1]]), collapse='')) }
      grob.all.con <- grob.list[order(na.rev)]
      all.cell <- ceiling(length(ID.sel)/ncol)*ncol
      cell.idx <- c(seq_len(length(ID.sel)), rep(NA, all.cell-length(ID.sel)))
      m <- matrix(cell.idx, ncol=ncol, byrow=TRUE)
      lay <- NULL; for (i in seq_len(length(con))) { lay <- rbind(lay, m+(i-1)*length(ID.sel)) }
 
      if (shiny==TRUE & length(grob.all.con)>=1) return(grid.arrange(grobs=grob.all.con, layout_matrix=lay, newpage=TRUE))
      g.tr <- lapply(grob.all.con, grobTree); g.tr <- g.tr[names(grob.all.con)]
      n.col <- ncol(lay); n.row <- nrow(lay)
      g.arr <- arrangeGrob(grobs=g.tr, layout_matrix=lay, widths=unit(rep(width/n.col, n.col), "npc"), heights=unit(rep(height/n.row, n.row), "npc")) 

    }; return(g.arr)

}


grob_list <- function(gene, geneV, coord, ID, cols, tis.path, tis.trans=NULL, sub.title.size, sam.legend='identical', legend.title=NULL, legend.ncol=NULL, legend.nrow=NULL, legend.position='right', legend.direction='vertical', legend.key.size=0.5, legend.label.size=8, legend.title.size=8, line.size=0.2, line.color='grey70', ...) {
  
  x <- y <- tissue <- NULL
  # Map colours to samples according to expression level.
  g.list <- function(j) {

    g.col <- NULL; con.idx <- grep(paste0("^", j, "$"), con)
    tis.col1 <- tis.col[con.idx]; scol1 <- scol[con.idx]

    for (i in tis.path) {

      tis.idx <- which(tis.col1 %in% i); if (length(tis.idx)==1) { g.col <- c(g.col, scol1[tis.idx])
      } else if (length(tis.idx)==0) { g.col <- c(g.col, NA) }

    }
    names(g.col) <- tis.df <- unique(coord[, 'tissue']) # The colors might be internally re-ordered alphabetically during mapping, so give them names to fix the match with tissues. E.g. c('yellow', 'blue') can be re-ordered to c('blue', 'yellow'), which makes tissue mapping wrong. Correct: colours are not re-ordered. The 'tissue' in 'data=coord' are internally re-ordered according to a factor. Therfore, 'tissue' should be a factor with the right order. Otherwise, disordered mapping can happen.
    # Make selected tissues transparent by setting their colours as NA.
    if (!is.null(tis.trans)) for (i in tis.df) { if (sub('_\\d+$', '', i) %in% tis.trans) g.col[i] <- NA }
    # Show selected or all samples in legend.
    if (length(sam.legend)==1) if (sam.legend=='identical') sam.legend <- unique(tis.path[!is.na(g.col)]) else if (sam.legend=='all') sam.legend <- unique(tis.path)
    leg.idx <- !duplicated(tis.path) & (tis.path %in% sam.legend)
    g <- ggplot()+geom_polygon(data=coord, aes(x=x, y=y, fill=tissue), color=line.color, size=line.size, linetype='solid')+scale_fill_manual(values=g.col, breaks=tis.df[leg.idx], labels=tis.path[leg.idx], guide=guide_legend(title=legend.title, ncol=legend.ncol, nrow=legend.nrow))+theme(axis.text=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_rect(fill="white", colour="grey80"), plot.margin=margin(0.1, 0.1, 0.1, 0.3, "cm"), axis.title.x=element_text(size=16,face="bold"), plot.title=element_text(hjust=0.5, size=sub.title.size), legend.position=legend.position, legend.background = element_rect(fill=alpha(NA, 0)), legend.direction=legend.direction, legend.key.size=unit(legend.key.size, "cm"), legend.text=element_text(size=legend.label.size), legend.title=element_text(size=legend.title.size))+labs(x="", y="")+scale_y_continuous(expand=c(0.01, 0.01))+scale_x_continuous(expand=c(0.01, 0.01))+ggtitle(paste0(k, "_", j)); return(g)


  }
  cname <- colnames(gene); con <- gsub("(.*)(__)(.*)", "\\3", cname); con.uni <- unique(con)
  grob.na <- grob.lis <- NULL; for (k in ID) {

    scol <- NULL; for (i in gene[k, ]) { 
      ab <- abs(i-geneV); col.ind <- which(ab==min(ab))[1]; scol <- c(scol, cols[col.ind])
    }

    idx <- grep("__", cname); c.na <- cname[idx]
    tis.col <- gsub("(.*)(__)(.*)", "\\1", c.na); g.lis <- NULL
    grob.na0 <- paste0(k, "_", con.uni); g.lis <- lapply(con.uni, g.list)
    # Repress popups by saving it to a png file, then delete it.
    tmp <- tempdir(); pa <- paste0(tmp, '/delete.png')
    png(pa); grob <- lapply(g.lis, ggplotGrob); dev.off(); if (file.exists(pa)) do.call(file.remove, list(pa))
    names(grob) <- grob.na0; grob.lis <- c(grob.lis, grob) 

  }; return(grob.lis)

}

nod_lin <- function(ds, lab, mods, adj, geneID, adj.min) {

  from <- to <- NULL
  idx.m <- mods[, ds]==lab; adj.m <- adj[idx.m, idx.m]; gen.na <- colnames(adj.m) 
  idx.sel <- grep(paste0("^", geneID, "$"), gen.na); gen.na[idx.sel] <- paste0(geneID, "_selected")
  colnames(adj.m) <- rownames(adj.m) <- gen.na; idx = adj.m > as.numeric(adj.min)
  link <- data.frame(from=rownames(adj.m)[row(adj.m)[idx]], to=colnames(adj.m)[col(adj.m)[idx]], width=adj.m[idx], stringsAsFactors=FALSE)
  # Should not exclude duplicate rows by "length".
  node.pas <- NULL; for (i in seq_len(nrow(link))) { node.pas <- c(node.pas, paste0(sort(c(link[i, 'from'], link[i, 'to'])), collapse='')) }
  w <- which(duplicated(node.pas)); link <- link[-w, ]
  link1 <- subset(link, from!=to, stringsAsFactors=FALSE); link1 <- link1[order(-link1$width), ]
  node <- data.frame(label=colnames(adj.m), size=colSums(adj.m), stringsAsFactors=FALSE)
  node <- node[order(-node$size), ]; return(list(node=node, link=link1))

}

col_bar <- function(geneV, cols, width, mar=c(3, 0.1, 3, 0.1)) {        

  color_scale <- y <- NULL
  cs.df <- data.frame(color_scale=geneV, y=1)
  cs.g <- ggplot()+geom_bar(data=cs.df, aes(x=color_scale, y=y), fill=cols, stat="identity", width=((max(geneV)-min(geneV))/length(geneV))*width)+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=margin(t=mar[1], r=mar[2], b=mar[3], l=mar[4], "cm"), panel.grid=element_blank(), panel.background=element_blank())+coord_flip(); # save(cs.g, file='cs.g')
  if (max(geneV)<=10000) cs.g <- cs.g+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))
  if (max(geneV)>10000) cs.g <- cs.g+scale_x_continuous(expand=c(0,0), labels=function(x) format(x, scientific=TRUE))+scale_y_continuous(expand=c(0,0))
  return(cs.g)

}


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

    gene2 <- aggregate(x=t(gene2), by=list(sam.var=cna), FUN=rep.aggr)
    sam.var <- gene2[, 1]; gene2 <- t(gene2[, -1]); colnames(gene2) <- sam.var

  }; rna <- rownames(gene2); gene2 <-apply(gene2, 2, as.numeric); rownames(gene2) <- rna 
  return(list(gene2=as.data.frame(gene2), gene3=as.data.frame(gene3), gen.rep=as.data.frame(gen.rep)))

}

# enableWGCNAThreads()
inter.svg <- readLines("example/root_cross_final.svg")
inter.data <- read.table("example/root_expr_ann_row_gen.txt", header=TRUE, row.names=1, sep="\t")

shinyServer(function(input, output, session) {

  output$dld.svg <- downloadHandler(

    filename=function(){ "root_cross_final.svg"}, 
    content=function(file){ writeLines(inter.svg, file) }

  )

  output$dld.data <- downloadHandler(

    filename=function(){ "root_expr_ann_row_gen.txt" },
    content=function(file){ write.table(inter.data, file, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t") }

  )
  # Filter parameters.
  fil <- reactiveValues(P=0, A=-Inf, CV1=-Inf, CV2=Inf)

  observe({

    input$fileIn; input$geneInpath
    updateRadioButtons(session, inputId="dimName", label="Step 3: is column or row gene?", 
    inline=TRUE, choices=c("None", "Row", "Column"), selected="None")
    updateSelectInput(session, 'sep', 'Step 4: separator', c("None", "Tab", "Comma", "Semicolon"), "None")
    updateRadioButtons(session, inputId='log', label='Data transform', choices=c("No", "log2", "exp.2"), selected="No", inline=TRUE)
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
    if (grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn)) {

        incProgress(0.5, detail="Loading matrix. Please wait.")
        if (input$fileIn=="organ_Mustroph") df.te <- fread.df(input="example/ucr_efp_expr_ann_row_gen.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="shoot_root_Mustroph") df.te <- fread.df(input="example/ucr_efp_expr_ann_row_gen.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="root_roottip_Mustroph") df.te <- fread.df(input="example/ucr_efp_expr_ann_row_gen.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="shoot_Mustroph") df.te <- fread.df(input="example/ucr_efp_expr_ann_row_gen.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="root_Geng") df.te <- fread.df(input="example/root_expr_ann_row_gen.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="brain_Chen") df.te <- fread.df(input="example/brain_expr_ann_row_gen.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        if (input$fileIn=="map_Census") df.te <- fread.df(input="example/us_population2018.txt", isRowGene=TRUE, header=TRUE, sep="\t", fill=TRUE)
        gen.rep <- df.te[['gen.rep']]
        gene2 <- df.te[['gene2']]; if (input$log=='log2') { 
          
          g.min <- min(gene2) 
          if (g.min<0) gene2 <- gene2-g.min+1; if (g.min==0) gene2 <- gene2+1; gene2 <- log2(gene2)  

        }; if (input$log=='exp.2') gene2 <- 2^gene2
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
      
      }; if (input$log=='exp.2') gene2 <- 2^gene2 
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
      se <- filter_data(data=se, ann=ann.col, samples=NULL, conditions=NULL, pOA=c(fil$P, fil$A), CV=c(fil$CV1, fil$CV2), dir=NULL)
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

    if (!grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn)) return()
    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    gID$new <- setdiff(gID$geneID, gID$all); gID$all <- c(gID$all, gID$new)

    })

  output$bar <- renderPlot({

    if ((grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn) & !is.null(geneIn()))|((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & !is.null(input$svgInpath) & !is.null(geneIn()))) {

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


  svg.df <- reactive({ 

    if (((input$fileIn=="Compute locally"|input$fileIn=="Compute online") & 
    !is.null(input$svgInpath))|(grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn) & is.null(input$svgInpath))) {

      withProgress(message="Tissue heatmap: ", value=0, {
    
        incProgress(0.5, detail="Extracting coordinates. Please wait.") 
	if (input$fileIn=="Compute locally"|input$fileIn=="Compute online") { svg.path <- input$svgInpath$datapath; svg.na <- input$svgInpath$name } else if     
(input$fileIn=="organ_Mustroph") { svg.path <- "example/organ_final.svg"; svg.na <- "organ_final.svg" } else if (input$fileIn=="shoot_root_Mustroph") { svg.path <- "example/shoot_root_final.svg"; svg.na <- "shoot_root_final.svg" } else if (input$fileIn=="root_roottip_Mustroph") { svg.path <- "example/root_roottip_final.svg"; svg.na <- "root_roottip_final.svg" } else if (input$fileIn=="shoot_Mustroph") { svg.path <- "example/shoot_final.svg"; svg.na <- "shoot_final.svg" } else if (input$fileIn=="root_Geng") { svg.path <- "example/root_cross_final.svg"; svg.na <- "root_cross_final.svg" } else if (input$fileIn=="brain_Chen") { svg.path <- "example/brain_final.svg"; svg.na <- "brain_final.svg" } else if (input$fileIn=="map_Census") { svg.path <- "example/us_map_final.svg"; svg.na <- "us_map_final.svg" }   
       df_tis <- svg_df(svg.path=svg.path)
       validate(need(!is.character(df_tis), df_tis))
       return(df_tis)

      })

    }

  })

  observe({

    input$fileIn; geneIn(); input$adj.modInpath; input$A; input$p; input$cv1; input$cv2; svg.df()
    updateCheckboxGroupInput(session, inputId="tis", label='Select tissues to be transparent:', choices=unique(svg.df()[['tis.path']]), selected='', inline=TRUE)

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

    }); grob$all <- gs()

  })

  observeEvent(gID$new, { grob.all <- c(grob$all, gs()); grob$all <- grob.all[unique(names(grob.all))] })
  # In "observe" and "observeEvent", if one code return (NULL), then all the following code stops. If one code changes, all the code renews.
  observe({
    
    if (is.null(geneIn())|is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)

    output$tissue <- renderPlot(width=as.numeric(input$width)/2*as.numeric(input$col.n), height=as.numeric(input$height)*length(input$dt_rows_selected), {

    if (is.null(input$dt_rows_selected)|is.null(svg.df())|gID$geneID[1]=="none"|is.null(grob$all)) return(NULL)
    if (length(color$col=="none")==0|input$color=="") return(NULL)

    r.na <- rownames(geneIn()[["gene2"]]); gID$geneID <- r.na[input$dt_rows_selected]
    grob.na <- names(grob$all); con <- unique(con())
    idx <- NULL; for (i in gID$geneID) { idx <- c(idx, grob.na[grob.na %in% paste0(i, '_', con)]) } 
    grob.lis.p <- grob$all[sort(idx)] #grob.lis.p <- grob.lis.p[unique(names(grob.lis.p))]
    lay_shm(lay.shm=input$gen.con, con=con, ncol=input$col.n, ID.sel=gID$geneID, grob.list=grob.lis.p, width=input$width, height=input$height, shiny=TRUE) 

    })

  })

  output$ori.svg <- renderImage({

    if ((is.null(input$svgInpath) & !grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn))|input$fileIn==
    "None") return(list(src="precompute/blank.png", contentType="image/png"))
    w <- as.numeric(input$width); h <- as.numeric(input$height); con.n <- length(con())
    W <- w/as.numeric(input$col.n); H <- h/(ceiling(con.n/as.numeric(input$col.n)))

    if (((input$fileIn=="Compute locally"|input$fileIn=="Compute online")|
    !is.null(input$svgInpath))|grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn)) {

      if (input$fileIn=="Compute locally"|input$fileIn=="Compute online") { svg.path <- input$svgInpath$datapath } else if (input$fileIn=="organ_Mustroph") { svg.path <- "example/organ_final.svg" } else if (input$fileIn=="shoot_root_Mustroph") { svg.path <- "example/shoot_root_final.svg" } else if (input$fileIn=="root_roottip_Mustroph") { svg.path <- "example/root_roottip_final.svg" } else if (input$fileIn=="shoot_Mustroph") { svg.path <- "example/shoot_final.svg" } else if (input$fileIn=="root_Geng") { svg.path <- "example/root_cross_final.svg" } else if (input$fileIn=="brain_Chen") { svg.path <- "example/brain_final.svg" } else if (input$fileIn=="map_Census") { svg.path <- "example/us_map_final.svg" }
      # svg.ln <- readLines(svg.path, 200); w.h <- svg.ln[grep(" width| height", svg.ln)]
      xmlfile <- xmlParse(svg.path); saveXML(doc=xmlfile, file='tmp/original.svg')
      na <- c('width', 'height'); lis <- xmlToList(svg.path)
      for (i in seq_len(length(lis))) {

        if (sum(na %in% names(lis[[i]]))==2) w.h <- as.character(xmlToList(svg.path)[[i]][na])

      }; w.h <- as.numeric(gsub("^(\\d+\\.\\d+|\\d+).*", "\\1", w.h)); r <- w.h[1]/w.h[2]
      list(src="tmp/original.svg", contentType="image/svg+xml", width=250, height=250/r, alt=NULL)

    }

  }, deleteFile=FALSE)

  edg <- reactive({

    if (input$fileIn!="Compute online") return(NULL)
    if (is.null(geneIn0())) return(NULL)
    gen.rep <- geneIn0()[['gen.rep']]
    log2.fc=1; fdr=0.5
    se <- SummarizedExperiment(assays=list(expr=as.matrix(gen.rep)), colData=data.frame(fct=colnames(gen.rep)))
    edgeR(se=se, method.norm='TMM', sample.factor='fct', method.adjust='BH', log2.fc=log2.fc, fdr=fdr)

  })

  dsq <- reactive({

    if (input$fileIn!="Compute online") return(NULL)
    if (is.null(geneIn0())) return(NULL)
    gen.rep <- geneIn0()[['gen.rep']]
    log2.fc=1; fdr=0.5
    se <- SummarizedExperiment(assays=list(expr=as.matrix(gen.rep)), colData=data.frame(fct=colnames(gen.rep)))
    deseq2(se=se, sample.factor='fct', method.adjust='BH', log2.fc=log2.fc, fdr=fdr)

  })

  output$ssg.sep <- renderPlot({ 

    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
    width=0.85
    lis.all <- list(edg=edg(), dsq=dsq()); sam.all <- names(lis.all[[1]])
    sam <- sam.all[1:3]
    ssg_sep(lis.all=lis.all, sam=sam, width=width)

  })
 
  
  output$w.table <- renderDataTable({
    
    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
    if (is.null(geneIn0())) return(NULL)
    sam.con <- unique(colnames(geneIn()))[1:3]
    lis.all <- list(edg=edg(), dsq=dsq()); sam.all <- names(lis.all[[1]])
    sam <- sam.all[1:3]; sam.tar <- sam.all[2]; sam.vs <- sam.all[c(1, 3)]
    meth <- 'edg'
    df.up <- lis.all[[meth]][[sam[2]]][[1]]
    df.dn <- lis.all[[meth]][[sam[2]]][[2]]
    df.up.dn <- rbind(df.up, df.dn)
    if (nrow(df.up.dn)==0) return("No SSGs detected.")
    print(df.up.dn)
    pat <- c(paste0('^', sam.tar, '_VS_', sam.vs, '_'), paste0('^', sam.vs, '_VS_', sam.tar, '_'))
    df.up.dn <- df.up.dn[, c(1, grep(paste0(pat, collapse='|'), colnames(df.up.dn)))]
    datatable(df.up.dn, selection=list(mode="multiple", target="row", selected=NULL), filter="top", extensions='Scroller', options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer')
 
  })
 
  output$ssg.sum <- renderPlot({ 

    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
    if (input$fileIn!="Compute online") return(NULL)
    if (is.null(geneIn0())) return(NULL)
    per=0.5; width=0.85
    lis.all <- list(edg(), dsq()); sam.all <- names(lis.all[[1]])
    sam <- sam.all[1:3]
    lis.aggr <- sig_frq(lis.all=lis.all, per=per, sam=sam)
    ssg_sum(lis.aggr=lis.aggr, width=width)

  })


  output$a.table <- renderDataTable({
    
    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
    if (is.null(geneIn0())) return(NULL)
    per=0.5
    sam.con <- unique(colnames(geneIn()))[1:3]
    lis.all <- list(edg=edg(), dsq=dsq()); sam.all <- names(lis.all[[1]])
    lis.aggr.r <- sig_frq(lis.all=lis.all, per=per)
    df.sum <- as.data.frame(lis.aggr.r[[1]][[1]][[2]])
    datatable(df.sum, selection=list(mode="multiple", target="row", selected=NULL), filter="top", extensions='Scroller', options=list(pageLength=5, lengthMenu=c(5, 15, 20), autoWidth=TRUE, scrollCollapse=TRUE, deferRender=TRUE, scrollX=TRUE, scrollY=200, scroller=TRUE), class='cell-border strip hover') %>% formatStyle(0, backgroundColor="orange", cursor='pointer')

  })

    expr.nor <- reactive({

    if (is.null(geneIn0())|is.null(edg())|is.null(dsq())) return(NULL)
      if (is.null(geneIn0())) return(NULL)
      gen.rep <- geneIn0()[['gen.rep']]
      se <- SummarizedExperiment(assays=list(expr=as.matrix(gen.rep)), colData=data.frame(fct=colnames(gen.rep)))
      norm_aggr(se=se, method.norm='TMM', data.trans='log2', sample.factor='fct', rep.aggr='mean')

    })


    #id.r <- rownames(se)[1]

    #plot_gen(se=se.nor, id=id.r)


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
    if (input$fileIn=="Compute online"|grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn)) {

      gene <- geneIn()[["gene2"]]; if (is.null(gene)) return()
      if (input$net.type=="S") { sft <- 12; type <- "signed" } else if (input$net.type=="U") { sft <- 6; type <- "unsigned" }

      withProgress(message="Computing: ", value = 0, {
        incProgress(0.3, detail="adjacency matrix.")
        incProgress(0.5, detail="topological overlap matrix.")
        incProgress(0.1, detail="dynamic tree cutting.")
        se <- SummarizedExperiment(assays=list(expr=as.matrix(gene)))
        adjMod <- adj_mod(data=se, type=type, minSize=input$min.size, dir=NULL)
        adj <- adjMod[['adj']]; mod4 <- adjMod[['mod']]

      }); return(list(adj=adj, mod4=mod4))

    }

  })

  observe({
    
    input$gen.sel; updateSelectInput(session, 'ds', "Select a module splitting sensitivity level", 3:2, selected="3")

  })

  mcol <- reactive({

    if (!is.null(adj.tree()) & (input$fileIn=="Compute online"|grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn))) { 
      
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

    } else if (input$fileIn=="Compute online"|grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn)) { 

      if (is.null(adj.tree())|input$gen.sel=="None") return(NULL); adj <- adj.tree()[['adj']]; mods <- mcol() 

    }
    lab <- mods[, input$ds][rownames(gene)==input$gen.sel]
    if (!is.null(lab)) if (lab=="0") { showModal(modalDialog(title="Module", "The selected gene is not assigned to any module. Please select a different one.")); return() } # All the arguments in 'if' statement are evaluated, regardless of the order. Therefore a 'NULL' object should not be compared with others (e.g. >, <) in 'if' statement. But it can be evaluated exclusively in a separate 'if' statement, e.g. 'if (!is.null(lab))'.

    withProgress(message="Matrix heatmap:", value=0, {
      incProgress(0.7, detail="Plotting...")
      se <- SummarizedExperiment(assays=list(expr=as.matrix(gene))); mods <- list(mod=mods); if (input$mat.scale=="By column/sample") scale.hm <- 'column' else if (input$mat.scale=="By row/gene") scale.hm <- 'row'
      matrix_heatmap(geneID=input$gen.sel, data=se, adj.mod=mods, ds=input$ds, scale=scale.hm, main=paste0('Network module containing ', input$gen.sel), title.size=10, static=FALSE)

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
    if (input$fileIn=="Compute locally") { adj <- adj.mod()[[1]]; mods <- adj.mod()[[2]] } else if (input$fileIn=="Compute online"|grepl("_Mustroph$|_Geng$|_Chen$|_Census$", input$fileIn)) { adj <- adj.tree()[[1]]; mods <- mcol() }
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



