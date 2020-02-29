#' Extract coordinates and tissue/path names from the SVG file
#'

#' @param svg.path The path of an SVG file.

#' @return A 2-length list, the first component is a data frame of the coordinates and the second is a vector of all tissue/path names.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' https://www.gimp.org/tutorials/ \cr https://inkscape.org/en/doc/tutorials/advanced/tutorial-advanced.en.html \cr http://www.microugly.com/inkscape-quickguide/
#' Jeroen Ooms (2018). rsvg: Render SVG Images into PDF, PNG, PostScript, or Bitmap Arrays. R package version 1.3. https://CRAN.R-project.org/package=rsvg \cr Paul Murrell (2009). Importing Vector Graphics: The grImport Package for R. Journal of Statistical Software, 30(4), 1-37. URL http://www.jstatsoft.org/v30/i04/ \cr Duncan Temple Lang and the CRAN Team (2018). XML: Tools for Parsing and Generating XML Within R and S-Plus. R package version 3.98-1.16. https://CRAN.R-project.org/package=XML \cr R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. RL https://www.R-project.org/ \cr    

#' @importFrom rsvg rsvg_ps 
#' @importFrom grImport PostScriptTrace 
#' @importFrom XML addAttributes xmlParse xmlRoot xmlSize xmlSApply xmlAttrs xmlName xmlChildren saveXML xmlApply

svg_df <- function(svg.path) {

  # Make sure the style is correct. If the stroke width is not the same across polygons such as '0.0002px', '0.216px', some stroke outlines cannot be recognised by 'PostScriptTrace'. Then some polygons are missing. Since the ggplot is based on 'stroke' not 'fill'.
  tmp <- system.file("extdata/shinyApp/tmp", package="spatialHeatmap")
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


