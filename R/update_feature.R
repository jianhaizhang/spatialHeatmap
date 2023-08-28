#' Update aSVG Spatial Features
#' 
#' Successful spatial heatmap plotting requires the aSVG features of interest have matching samples (cells, tissues, \emph{etc}) in the data. If this requirement is not fulfiled, either the sample identifiers in the data or the spatial feature identifiers in the aSVG should be changed. This function is designed to replace existing feature identifiers, stroke (outline) widths, and/or feature colors in aSVG files with user-provided entries.

#' @param df.new The custom feature identifiers, stroke (outline) widths, and/or feature colors, should be included in the data frame returned by \link{return_feature} as independent columns, and the corresponding column names should be "featureNew", "strokeNew", and "colorNew" respectively in order to be recognized. \cr To color the corresponding features, the identifiers in "featureNew" should be the same with matching sample identifiers. The numeric values in "strokeNew" would be the outline widths of corresponding features. The colors in "colorNew" would be the default colors for highlighting target features in the legend plot.  

#' @param dir The directory path where the aSVG files to update. It should be the same with \code{dir} in \code{\link{return_feature}}.
#' @return Nothing is returned. The aSVG files of interest in \code{dir} are updated with provided attributes, and are ready to use in function \code{\link{spatial_hm}}.

#' @examples

#' # The following shows how to download a chicken aSVG containing spatial features of 'brain'
#' # and 'heart' from the EBI aSVG repository directly 
#' # (https://github.com/ebi-gene-expression-group/anatomogram/tree/master/src/svg). An empty
#' # directory is recommended so as to avoid overwriting existing SVG files with the same names.
#' # Here "~/test" is used. 
#'
#' # Remote aSVG repos.
#' data(aSVG.remote.repo)
#' tmp.dir <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
#' tmp.dir.ebi <- paste0(tmp.dir, '/ebi.zip')
#' tmp.dir.shm <- paste0(tmp.dir, '/shm.zip')
#' \donttest{
#' # Download the remote aSVG repos as zip files. According to Bioconductor's 
#' # requirements, downloadings are not allowed inside functions, so the repos are 
#' # downloaded before calling "return_feature".  
#' download.file(aSVG.remote.repo$ebi, tmp.dir.ebi)
#' download.file(aSVG.remote.repo$shm, tmp.dir.shm)
#' remote <- list(tmp.dir.ebi, tmp.dir.shm)
#' # Make an empty directory "~/test" if not exist.
#' if (!dir.exists('~/test')) dir.create('~/test')
#' # Query the remote aSVG repos.
#' feature.df <- return_feature(feature=c('heart', 'brain'), species=c('gallus'), dir='~/test',
#' match.only=TRUE, remote=remote)
#' feature.df
#'
#' # New features, stroke widths, colors.
#' ft.new <- c('BRAIN', 'HEART')
#' stroke.new <- c(0.05, 0.1)
#' col.new <- c('green', 'red')
#' # Include new features, stroke widths, colors to the feature data frame.
#' feature.df.new <- cbind(featureNew=ft.new, strokeNew=stroke.new, colorNew=col.new, feature.df)
#' feature.df.new
#'
#' # Update features.
#' update_feature(df.new=feature.df.new, dir='~/test')
#' }
#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Hadley Wickham, Jim Hester and Jeroen Ooms (2019). xml2: Parse XML. R package version 1.2.2. https://CRAN.R-project.org/package=xml2
#' Cardoso-Moreira, Margarida, Jean Halbert, Delphine Valloton, Britta Velten, Chunyan Chen, Yi Shao, Angélica Liechti, et al. 2019. "Gene Expression Across Mammalian Organ Development." Nature 571 (7766): 505-9
#' Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber, Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2020). gplots: Various R Programming Tools for Plotting Data. R package version 3.0.3. https://CRAN.R-project.org/package=gplots

#' @export update_feature
#' @importFrom gplots col2hex
#' @importFrom xml2 read_xml xml_children xml_length xml_attr xml_set_attr xml_name xml_add_child xml_set_text write_xml

update_feature <- function(df.new, dir) {

  options(stringsAsFactors=FALSE); SVG <- parent <- feature <- NULL
  dir <- normalizePath(dir, winslash = "/", mustWork=FALSE)
  if (!is(df.new, 'data.frame')) stop('"df.new" must be a "data.frame"!')
  idx.upd <- grep('featureNew|colorNew|strokeNew', colnames(df.new), ignore.case=TRUE)
  if (length(idx.upd)==0) return('No new columns are detected!\n')
  df.upd <- df.new[idx.upd]; df.old <- df.new[-idx.upd]
  colnames(df.upd) <- cna.new <- tolower(colnames(df.upd))
  # Convert factors to vectors.
  df.upd <- as.data.frame(apply(df.upd, 2, as.vector))
  if ('strokenew' %in% cna.new) df.upd$strokenew <- as.numeric(df.upd$strokenew)
  df.new <- cbind(df.upd, df.old)
  # df.new$index <- as.numeric(df.new$index)
  # df.new$index1 <- as.numeric(df.new$index1)
  svgs.na <- unique(df.new$SVG)
  # Update each SVG file in order.
  for (i in svgs.na) {
    df0 <- subset(df.new, SVG==i)
    if ('featurenew' %in% cna.new) {
      dup <- duplicated(df0[, 'featurenew'])
      if (any(dup)) stop(paste0("Duplicated feature \'", paste0(df0[, 'featurenew'][dup], collapse=', '), "\' detected in ", i, "!"))
    } 
    path.in <- paste0(dir, '/', i); doc <- read_xml(path.in)
    lis <- svg_attr(doc, NULL, FALSE); df.attr <- lis$df.attr
    out <- lis$out; ply <- lis$ply
    # Combine new columns in df.new to df.attr.
    int.ft <- intersect(df.attr$feature, df0$feature)
    if (length(int.ft)==0) stop('No matching spatial features detected for ', i, '!')
    ft.no <- setdiff(df0$feature, int.ft)
    if (length(ft.no)>0) cat(ft.no, 'No matching spatial features detected for ', ft.no, '!\n')
    df.attr <- subset(df.attr, feature %in% int.ft)
    df0 <- subset(df0, feature %in% int.ft)
    df.attr <- df.attr[order(df.attr$feature), ]
    df0 <- df0[order(df0$feature), ]
    df.sub <- cbind(df0[cna.new], df.attr)
    # len <- xml_length(doc) out <- xml_children(doc)[[len-1]]; ply <- xml_children(doc)[[len]]
    # Update features for two parent layers.
    for (k in list(out, ply)) { 

      # Check if the parent needs to be updated. 
      p.id <- make.names(xml_attr(k, 'id')); if (is.na(p.id)) next
      df1 <- subset(df.sub, parent==p.id); if (nrow(df1)==0) next
      cat(paste0('Updating attributes', ' in ', path.in, '\n'))
      # Update features by order.
      for (j in seq_len(nrow(df1))) {
 
        nod0 <- xml_children(k)[[df1$index.sub[j]]]
        if ('featurenew' %in% cna.new) {
          chil0 <- xml_children(nod0); nas <- xml_name(chil0)
        # If there is no "title" node, new feature is added to "id", and no need to add a "title" node. Otherwise, added to the text in title.
          if (!('title' %in% nas)) {

            if (make.names(xml_attr(nod0, 'id'))!=df1$feature[j]) {
              cat(df1$feature[j], 'is not updated! \n'); next
            } else { xml_set_attr(nod0, 'id', df1[j, 'featurenew']) }
          
          } else xml_set_text(chil0[[which(nas=='title')]], df1[j, 'featurenew'])
        # Add a 'title' node if the feature to update does not have it, and update 'chil0' and 'nas'.
        # if (!('title' %in% nas)) { xml_add_child(nod0, 'title', .where=0); xml_set_attr(xml_children(nod0)[[1]], 'id', df1[j, 1]); chil0 <- xml_children(nod0); nas <- c('title', nas)
        # }; xml_set_text(chil0[[which(nas=='title')]], df1[j, 1])

        }; sty0 <- xml_attr(nod0, 'style')
        if ('strokenew' %in% cna.new) {
          stk <- paste0('stroke-width:', df1$strokenew[j])
          sty0 <- col_str(sty0, '^stroke-width:', stk)
          #xml_set_attr(nod0, 'style', sty.new)
        } # Update stroke width.
        if ('colornew' %in% cna.new) {
          col <- paste0('fill:', col2hex(df1$colornew[j]))
          sty0 <- col_str(sty0, '^fill:', col)
          # xml_set_attr(nod0, 'style', sty.new)
        } # Update fill colors.
        xml_set_attr(nod0, 'style', sty0)
      }

    }; write_xml(doc, file=path.in)

  }

}

#' # Update fill and stroke width in a node.
#' @keywords Internal
#' @noRd
col_str <- function(style, at.pat, attr.val) {
  sty <- strsplit(style, ';')[[1]]
  attr.idx <- grepl(at.pat, sty)
  if (any(attr.idx)) sty[attr.idx][1] <- attr.val else sty <- c(attr.val, sty)
  sty <- paste0(sty, collapse=';'); return(sty)
}



