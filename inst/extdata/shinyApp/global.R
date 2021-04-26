# Import internal functions.
deter_core <- get('deter_core', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

rela_size <- get('rela_size', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

norm_data <- get('norm_data', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

cord_parent <- get('cord_parent', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

use <- get('use', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

cord <- get('cord', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

xy0 <- get('xy0', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

xy <- get('xy', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

tit_id <- get('tit_id', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

out_ply <- get('out_ply', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

sort_gen_con <- get('sort_gen_con', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

test_ffm <- get('test_ffm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

read_hdf5 <- get('read_hdf5', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

matrix_hm <- get('matrix_hm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Function to extract nearest genes.
sub_na <- get('sub_na', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

adj_mod <- get('adj_mod', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

filter_data <- get('filter_data', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

col_bar <- get('col_bar', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

lay_shm <- get('lay_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

nod_lin <- get('nod_lin', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Break combined path to a group (g=TRUE) or siblings (g=FALSE).
path_br <- get('path_br', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# The outline or tissue nodes are checked for combines paths. If combined paths are detected, those outside a group are broken to a group while those inside a group are broken as siblings.  
path_br_all <- get('path_br_all', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# 'a' nodes are not removed.
svg_attr <- get('svg_attr', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

svg_df <- get('svg_df', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Separate SHMs of grob and ggplot. Different SHMs (under different SVGs) of same 'gene_condition' are indexed with suffixed of '_1', '_2', ...
gg_shm <- get('gg_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
grob_shm <- get('grob_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Subset data matrix by correlation or distance measure.
submatrix <- get('submatrix', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Adjust legend key size and rows in ggplot.
gg_lgd <- get('gg_lgd', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Add value keys SHMs.
gg_2lgd <- get('gg_2lgd', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Prepare interactive SHMs in html.
html_ly <- get('html_ly', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Make videos.
video <- get('video', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
# Shown popup window. 
modal <- function(title = NULL, msg) {
  modalDialog(title = title, span(msg),
    footer = tagList(modalButton("Dismiss")), size = c("m")
  )
}
# Import input matrix, able to deal with/separate numeric matrix, character matrix, and mixture of both.
fread_df <- function(input, isRowGene=TRUE, header=TRUE, sep='auto', fill=TRUE, rep.aggr='mean', check.names=FALSE) {
  
  if (!is(input, 'data.frame') & !is(input, 'matrix')) { 
  df0 <- tryCatch({
    fread(input=input, header=header, sep=sep, fill=fill, check.names=check.names)
  }, error = function(error_condition) {
    # Deals with only one column with row names.
    fread(input = input, header = FALSE, sep = sep, fill = FALSE, check.names = check.names)
  }) 
    cna <- make.names(colnames(df0))
    if (cna[1]=='V1') cna <- cna[-1] else cna <- cna[-ncol(df0)] 
    df1 <- as.data.frame(df0); rownames(df1) <- df1[, 1]
    df1 <- df1[, -1, drop = FALSE]; colnames(df1) <- cna
    if(isRowGene==FALSE) df1 <- t(df1)
    cna <- colnames(df1); rna <- rownames(df1) 
  } else { df1 <- input; rna <- rownames(df1); cna <- colnames(df1) }
  # Covert factors to character. Only data.frame works not matrix.
  fct.idx <- vapply(df1, is.factor, logical(1))
  df1[fct.idx] <- lapply(df1[fct.idx], as.character) 
  # Subsetting identical column names in a matrix will not trigger appending numbers.
  df1 <- as.matrix(df1)
  # Isolate data and row metadata.
  na <- vapply(seq_len(ncol(df1)), function(i) { tryCatch({ as.numeric(df1[, i]) }, warning=function(w) { return(rep(NA, nrow(df1))) }, error=function(e) { stop("Please make sure input data are numeric!") }) }, FUN.VALUE=numeric(nrow(df1)) )
  na <- as.data.frame(na); rownames(na) <- rna
  idx <- colSums(apply(na, 2, is.na))!=0
  df.num <- na[!idx]; colnames(df.num) <- cna <- cna[!idx]
  df.met <- as.data.frame(df1)[idx]
  # Only one column is selected for row metadata.
  if (ncol(df.met) > 0) {
    met.idx <- grep('metadata', colnames(df.met), ignore.case = TRUE)[1]
    if (!is.na(met.idx)) df.met <- df.met[, met.idx, drop = FALSE] else { colnames(df.met)[1] <- 'metadata'; df.met <- df.met[, 1, drop = FALSE] }
  }
  # Only row metadata.
  if (ncol(df.num) == 0) {
    return(list(df.aggr = NULL, df.met=as.data.frame(df.met), df.rep = NULL, con.na = FALSE))
  }
  form <- grepl("__", cna); if (sum(form)==0) { colnames(df.num) <- paste0(cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE
  if(sum(is.na(as.numeric(as.matrix(df.num))))>=1) return('Make sure all values in data matrix are numeric.')
  
  df.rep <- df.num; rna <- rownames(df.rep); df.rep <- apply(df.rep, 2, as.numeric); rownames(df.rep) <- rna
  # Aggregate replicates.
  if (any(duplicated(cna)) & !is.null(rep.aggr)) {

    # To keep colnames, "X" should be a character, not a factor.
    if (rep.aggr=='mean') df.num <- sapply(X=unique(cna), function(x) rowMeans(df.num[, cna==x, drop=FALSE]))
    if (rep.aggr=='median') {
      df.num <- sapply(X=unique(cna), function(x) Biobase::rowMedians(df.num[, cna==x, drop=FALSE]))
      rownames(df.num) <- rna
    }

  }; df.aggr <-apply(df.num, 2, as.numeric); rownames(df.aggr) <- rna 
  return(list(df.aggr=as.data.frame(df.aggr), df.met=as.data.frame(df.met), df.rep=as.data.frame(df.rep), con.na=con.na))

}

# Separate colour ingredients.
col_sep <- function(color) {

  color <- gsub(' |\\.|-|;|,|/', '_', color)
  color <- strsplit(color, '_')[[1]]
  color <- color[color!='']; return(color)

}

# Check/process "sample__condition" in the se extracted from tar.
se_from_db <- function(se) {
  dat <- assay(se); cold <- colData(se)
  form <- grepl("__", colnames(dat))
  if (sum(form)==0) {
    if (all(c('sample', 'condition') %in% colnames(cold))) { 
      cna <- colnames(dat) <- paste0(cold$sample, '__', cold$condition)
      if (any(duplicated(cna))) return('Duplicated "sample__condition" replicates are detected in the selected dataset!')
    } else if ('sample' %in% colnames(cold)) {
      if (any(duplicated(cold$sample))) return('The "sample" should not be duplicated in the absence of "condition"!')
      colnames(dat) <- cold$sample
    }
  }; return(dat)
}

# Extract target svgs in tar into tmp folder, and return the paths. 
extr_svg <- function(file, name) {
  dir <- paste0(tempdir(check=TRUE), '/svg_shm')
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  untar(file, exdir=dir, tar='tar')
  pa <- paste0(dir, '/', name) 
  if (file.exists(pa)) return(pa) else return()
}

# Extract svg path/na from uploaded or internal tar files.
svg_pa_na <- function(svg.path, pa.svg.upl) {
  svg.na <- NULL; for (i in seq_along(svg.path)) {
    # Extract svg names. 
    str <- strsplit(svg.path[[i]], '/')[[1]]
    na0 <- str[length(str)]
    if (!grepl('\\.svg$', na0)) return('No aSVG file is detected! Solution: 1) select another aSVG and rematch it to data; 2) add an aSVG file for the selected data in the backend aSVG tar file or uploaded aSVG tar file.')
    svg.na <- c(svg.na, na0)
    # Complete uploaded svg paths.
    if (!grepl('example/', svg.path[[i]])) {
      # The data/svg precedence: uploaded tar > internal tar > default examples. The duplicated data/svgs are removed according to this precedence when processing data/svg upstream.
      pa0 <- NULL; if (!is.null(pa.svg.upl)) pa0 <- extr_svg(file=pa.svg.upl, name=na0)
      if (is.null(pa.svg.upl)|is.null(pa0)) {
        tar.all <- list.files('example', pattern='\\.tar$', full.names=TRUE)
        tar.svg <- tar.all[!grepl('data_shm.tar$', tar.all)][1]
        pa0 <- extr_svg(file=tar.svg, name=na0)
      }; if (is.null(pa0)) return(paste0("This aSVG file is not detected: ", na0, "!")) else svg.path[[i]] <- pa0
    }
  }; return(list(svg.path=svg.path, svg.na=svg.na))
}

# Check suffixes if multiple svgs.
svg_suffix <- function(svg.path, svg.na) {
  if (length(svg.na)>1) {
    if (!all(grepl('_shm\\d+\\.svg$', svg.na, perl=TRUE))) return("Suffixes of aSVGs should be indexed as '_shm1.svg', '_shm2.svg', '_shm3.svg', ...")
    ord <- order(gsub('.*_(shm.*)$', '\\1', svg.na))
    svg.path <- svg.path[ord]; svg.na <- svg.na[ord]  
  }; return(list(svg.path=svg.path, svg.na=svg.na))
}

## Rematch features.
# Create a panel for each data feature, where aSVG features can be dropped.
# ns() is the namespace in shiny modules.
ft_dat <- function(x, ns) {
 span(class = "panel panel-default",
   div(class = "panel-heading", x), 
   div(class = "panel-body", id = ns(x))
  )
}

# Allow features are draggable across panels.
# ns() is the namespace in shiny modules.
ft_js <- function(x, ns) {
  sortable_js(css_id = ns(x),
    options = sortable_options(
      multiDrag = NULL, sort = FALSE, animation = 1000, direction = NULL, 
      group = list(name = "sortGroup1", put = TRUE),
      onSort = sortable_js_capture_input(ns(x))
    )
  )
}

# Convert aSVG features to draggable items.
# ns() is the namespace in shiny modules.
ft2tag <- function(ft){
  lapply(ft, function(i) { tag("span", list(class = class(i), tags$span(class = "glyphicon glyphicon-move"), i)) }
  )
}

# Clean trash files in animation and video.
ggly_rm <- function() {
  if (dir.exists('www/ggly/')) {
    cat("Removing animation files in 'www/ggly/' ... \n")
    unlink('www/ggly/lib', recursive=TRUE)
    file.remove(list.files('www/ggly/', '*.html$', full.names=TRUE))
  } else dir.create('www/ggly', recursive=TRUE)
}
vdo_rm <- function() {
  if (dir.exists('www/video/')) {
    cat("Removing video file in 'www/video/' ... \n")
    file.remove(list.files('www/video/', '*.mp4$', full.names=TRUE))
  } else dir.create('www/video/', recursive=TRUE)
}


## Spatial enrichment
# Translate overlap up/down genes from different methods in a list to a data frame.
venn_inter <- function(lis.all) {
  
  gen.all <- unique(unlist(lis.all))
  # Create an empty data frame, where rows are all genes and columns are methods.
  zero <- rep(0, length(gen.all))
  df.all <- as.data.frame(matrix(rep(0, length(gen.all)*length(lis.all)), ncol=length(lis.all)))
  cna <- names(lis.all)
  suf <- unique(gsub('(.*)\\.(.*)', '\\2', cna))
  meth <- unique(gsub('(.*)\\.(.*)', '\\1', cna))
  names(lis.all) <- colnames(df.all) <- meth
  rownames(df.all) <- gen.all
  # Retrieve all overlaps.
  lis <- venn(lis.all, show.plot=FALSE) 
  inter <- attributes(lis)$intersections
  # Translate all overlaps into a data frame.
  for (i in seq_along(inter)) {
    lis0 <- inter[[i]]; na0 <- strsplit(names(inter[i]), ':')[[1]]
    w1 <- which(gen.all %in% lis0); w2 <- which(meth %in% na0)  
    df.all[w1, w2] <- 1
  }

  df.all <- cbind(type=suf, total=rowSums(df.all), df.all)
  # Matrix accepts dupliated rows. Some genes might be up in one method while down in other methods, so when combine up and down table in a single table, there could be duplicated row names. In data frame, duplicated row names are appended 1.
  df.all <- as.matrix(df.all[order(df.all$total, decreasing=TRUE), ])
  return(df.all)
}

# Given a DEG list of different methods, plot the overlap matrix.
deg_olp <- function(deg.lis) {
  mat <- vapply(names(deg.lis), function(x) vapply(names(deg.lis), function(y) length(intersect(deg.lis[[x]], deg.lis[[y]])), numeric(1)), numeric(length(deg.lis)))
  mel <- reshape2::melt(mat)
  g <- ggplot(data=mel, aes(x=Var1, y=Var2, fill=value))+geom_tile(colour="white")+scale_fill_gradient(low="lightcyan3", high="darkorange")+theme_minimal()+theme(axis.text=element_text(angle=45, vjust=1, size=10, hjust=1))+coord_fixed()+geom_text(aes(Var2, Var1, label=value), color="black", size=4)+theme(axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.ticks=element_blank()); return(g)
}





