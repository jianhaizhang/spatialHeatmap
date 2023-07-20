if (interactive()) {
  requireNamespace('DESeq2'); requireNamespace('av'); requireNamespace('BiocGenerics'); requireNamespace('distinct') 
  requireNamespace('dendextend'); requireNamespace('HDF5Array'); requireNamespace('magick'); requireNamespace('DT'); 
  requireNamespace('pROC'); requireNamespace('shinyWidgets'); requireNamespace('shinyjs'); requireNamespace('htmltools');
  requireNamespace('shinyBS'); requireNamespace('sortable'); requireNamespace('org.Hs.eg.db'); requireNamespace('org.Mm.eg.db')
  requireNamespace('org.At.tair.db'); requireNamespace('org.Dr.eg.db'); requireNamespace('org.Dm.eg.db');
  requireNamespace('AnnotationDbi'); requireNamespace('sparkline'); requireNamespace('spsComps'); requireNamespace('spsUtil')
}

# Accessing local html files by iframes.
# addResourcePath("tmpuser", getwd())
# Application-level cache, max size is 1G.
shiny::shinyOptions(cache = cachem::cache_disk(dir="./app_cache/cache/", max_size = 1024 * 1024^2))

na.sgl.def <- '^covis_'
na.cus <- c('customBulkData', 'customCovisData')
na.cus.dis <- c('customBulkData', 'customCovisData')
na.sgl <- c('^covis_|^customCovisData$')

# Temporary directory.
tmp.dir <- normalizePath(tempdir(check=TRUE), winslash="/", mustWork=FALSE)
tmp.file <- normalizePath(tempfile(), winslash='/', mustWork=FALSE)

# Tooltip of colData table in covis.
msg.meta.ann <- 'Metadata of single-cell (see "bulkCell") data. "label", "label1": cell group labels obtained from annotation labels, marker genes, etc.'
msg.meta.coclus <- "Metadata of bulk and single-cell data after co-clustering. <br/> 1. 'cluster': clusters containing only cells or cells and tissues (co-clusters); <br/> 2. 'bulkCell': 'bulk' and 'cell' indicate tissues and cells respectively; <br/> 3. 'assignedBulk': tissue lables assinged to cells as group labels, 'none' indicates no assignment; <br/> 4. 'similarity': Spearman correlation efficient used for tissue-cell assignment through a nearest-neighbor approach."

# Run button colors.
hp <- 'background-color:#e7e7e7'
hp.txt <- 'color:black'
run.msg <- 'Always click the <span style="color:white;font-weight:bold;background-color:#369ef7;font-size:20px">"Run"</span> button to see the latest results!'
run.col <- "color:white;background-color:#369ef7;border-color:#ddd"
run.top <- "margin-top:24px;color:white;background-color:#369ef7;border-color:#ddd"
# Confirm buttons.
conf.col <- 'margin-top:2px;margin-bottom:-10px;margin-left:20px;padding-top:2px;padding-bottom:2px;color:white;background-color:#369ef7;border-color:#ddd'
# Rectangle.
rec <- 'border-color:#3c8dbc;border-width:1px;border-style:solid'

# Bitmap extensions accepted in uploaded images.
raster.ext <- c('.jpg', '.JPG', '.png', '.PNG')

# Confirm button labels.
lab.sgl <- 'Search by single gene ID (e.g. Cav2)'
lab.mul <- 'Search by single or multiple gene IDs (e.g. Cav2,Apoh)'

# Search portion of URLs on the landing page.

# Extract parameter values from url.
url_val <- function(na, lis.url, def=NULL) {
  # if (!exists('lis.url')) return('null')
  if (!na %in% names(lis.url$par)) val <- 'null' else if (length(lis.url$par)==0) val <- 'null' else val <- lis.url$par[[na]]
    # In "ifelse", the length of returned value is same with the first argument.
    # val <- ifelse(length(lis.url$par)==0, 'null', lis.url$par[[na]])
   val <- gsub('\\"', '', val)
   if (!is.null(def)) ifelse('null' %in% val[1], def, val) else val
}

# Import internal functions.
read_svg <- get('read_svg', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
read_svg_m <- memoise::memoise(read_svg, cache = getShinyOption("cache"))

reduce_dim <- get('reduce_dim', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
reduce_dim_m <- memoise::memoise(reduce_dim, cache = getShinyOption("cache"))

norm_cell <- get('norm_cell', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
norm_cell_m <- memoise::memoise(norm_cell, cache = getShinyOption("cache"))


ovl_dat_db <- get('ovl_dat_db', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
check_se <- get('check_se', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
check_sce <- get('check_sce', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
cut_dendro <- get('cut_dendro', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
com_roc <- get('com_roc', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
qc_cell <- get('qc_cell', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
check_exp <- get('check_exp', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
check_obj <- get('check_obj', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#img_pa_na <- get('img_pa_na', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
covis_trans <- get('covis_trans', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
svg_separ <- get('svg_separ', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
sc_qc_plot <- get('sc_qc_plot', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
detect_cluster <- get('detect_cluster', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
detect_cluster_m <- memoise::memoise(detect_cluster, cache = getShinyOption("cache"))

dim_color_idp <- get('dim_color_idp', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
dim_color_coclus <- get('dim_color_coclus', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
dim_color2cell <- get('dim_color2cell', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
dim_color <- get('dim_color', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
nn_graph <- get('nn_graph', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
svg_raster <- get('svg_raster', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

raster_path <- get('raster_path', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# dat_fun <- get('dat_fun', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

df_is_as <- get('df_is_as', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

thrsd <- get('thrsd', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
scale_all <- get('scale_all', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

edgeR <- get('edgeR', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

limma <- get('limma', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

deseq2 <- get('deseq2', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

distt <- get('distt', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

sf_var <- get('sf_var', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
up_dn <- get('up_dn', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

deg_ovl_mat  <- get('deg_ovl_mat', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

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
# Application-level cache.
gg_shm_m <- memoise::memoise(gg_shm, cache = getShinyOption("cache"))

grob_shm <- get('grob_shm', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
#grob_shm_m <- memoise::memoise(grob_shm, cache = getShinyOption("cache"))

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
show_mod <- get('show_mod', envir=asNamespace('spatialHeatmap'), inherits=FALSE)
modal <- get('modal', envir=asNamespace('spatialHeatmap'), inherits=FALSE)

# Extract a 1-column data frame of URLs. If no column of URL is present, the default google-search URLs are composed.
link_dat <- function(df.met, link.only=TRUE) {
  # save(df.met, link.only, file='link.dat.arg')
  cna <- colnames(df.met); rna <- rownames(df.met)
  link.idx <- grep('link|links', cna, ignore.case=TRUE)[1] 
  if (is.na(link.idx)) {   
    # Iterative operation on data frame: vectorization is faster than for/lapply loop.  
    link <- paste0('<a href=\"https://www.google.com/search?q=', rna, '" target="_blank">link</a>')  
    # link <- lapply(rownames(df.met), function(x) a("link", href=paste0('https://www.google.com/search?q=', x), target="_blank"))
    # link <- unlist(lapply(link, as.character))  
  } else { link <- df.met[, link.idx]
    link <- paste0('<a href="', gsub('\\\\"', '\'', link), '" target="_blank">link</a>')  
  }
  df.lk <- data.frame(link=link, row.names=rownames(df.met))
  # df.met: 0 or 1 column, return(df.lk)
  if (link.only==TRUE | ncol(df.met)<=1) return(df.lk) else {
    cna.other <- setdiff(colnames(df.met), colnames(df.lk))
    df.met <- cbind(df.met[, cna.other[1], drop=FALSE], df.lk, df.met[, cna.other[setdiff(seq_along(cna.other), 1)], drop=FALSE]); return(df.met)
  }
} 

# Import input matrix, able to deal with separate numeric matrix, character matrix, and mixture of both.
fread_df <- function(input, isRowGene=TRUE, header=TRUE, sep='auto', fill=TRUE, rep.aggr='mean', check.names=FALSE) { 
  if (is(input, 'dgCMatrix')|is(input, 'matrix')|is(input, 'data.frame')|is(input, 'DFrame')|is(input, 'DelayedMatrix')) input <- as.data.frame(as.matrix(input))
  if (!is(input, 'data.frame')) {
  df0 <- tryCatch({
    fread(input=input, header=header, sep=sep, fill=fill, check.names=check.names)
  }, error = function(error_condition) {
    # Deals with only one column with row names.
    fread(input = input, header = FALSE, sep = sep, fill = FALSE, check.names = check.names)
  }) 
    cna <- make.names(colnames(df0))
    if (cna[1]=='V1') cna <- cna[-1] else cna <- cna[-ncol(df0)] 
    df1 <- as.data.frame(df0); rownames(df1) <- make.names(df1[, 1])
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
  if (nrow(df1)==1) na <- matrix(na, byrow=TRUE, ncol=ncol(df1))
  na <- as.data.frame(na); rownames(na) <- rna; colnames(na) <- cna
  vap <- df_is_as(na, is.na); idx <- colSums(vap)!=0
  df.num <- na[!idx]; colnames(df.num) <- cna <- cna[!idx]
  df.met.all <- as.data.frame(df1)[idx]
  cat('Preparing URLs .. \n')
  df.link <- link_dat(df.met.all) # Works if ncol(df.met.all) is 0.
  if (ncol(df.met.all) > 0) {
  cat('Preparing metadata .. \n')
    met.idx <- grep('^metadata$', colnames(df.met.all), ignore.case = TRUE)[1]
    if (!is.na(met.idx)) { 
      df.met <- df.met.all[, met.idx, drop = FALSE] 
      colnames(df.met) <- 'metadata'
      df.met <- cbind(df.met, df.link)
    } else df.met <- df.link
  } else df.met <- df.link

  # Only row metadata.
  if (ncol(df.num) == 0) {
    return(list(df.aggr = NULL, df.met=as.data.frame(df.met), df.rep = NULL, con.na = FALSE))
  }
  form <- grepl("__", cna); if (sum(form)==0) { cna <- colnames(df.num) <- paste0(cna, '__', 'con'); con.na <- FALSE } else con.na <- TRUE
  sf <- gsub('(.*)__(.*)', '\\1', cna)
  vari <- gsub('(.*)__(.*)', '\\2', cna)
  df.cdat <- DataFrame(spFeature=sf, variable=vari)
  if(sum(is.na(as.numeric(as.matrix(df.num))))>=1) return('Make sure all values in data matrix are numeric.')
  
  df.rep <- df.num; df.rep <- df_is_as(df.rep, as.numeric)
  #if (TRUE %in% rdat) {
  #  rdat <- cbind(DataFrame(df.met[, !colnames(df.met) %in% colnames(rdat)]), DataFrame(rdat))
  # } else rdat <- NULL
  se.rep <- SummarizedExperiment(assays=list(rep=as.matrix(df.rep)), colData=df.cdat, rowData=df.met)
  if (!is.null(rep.aggr)) { 
    se.aggr <- aggr_rep(data=se.rep, assay.na=NULL, sam.factor=NULL, con.factor=NULL, aggr=rep.aggr)
    assayNames(se.aggr) <- 'aggr'
  } else se.aggr <- NULL
  return(list(se.rep=se.rep, se.aggr=se.aggr, con.na=con.na))
}

# Separate colour ingredients.
col_sep <- function(color) {
  color <- gsub(' |\\.|-|;|,|/', '_', color)
  color <- strsplit(color, '_')[[1]]
  color <- color[color!='']; return(color)
}

# Extract target svgs in tar into tmp folder, and return the paths. 
extr_svg <- function(file, name) {
  dir <- paste0(tempdir(check=TRUE), '/svg_shm')
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  untar(file, exdir=dir, tar='tar')
  pa <- paste0(dir, '/', name) 
  if (file.exists(pa)) return(pa) else return()
}

# Extract svg path/na from uploaded or internal tar files if not found in 'data' folder.
svg_pa_na <- function(svg.path, pa.svg.upl, raster.ext) {
  svg.na <- NULL; for (i in seq_along(svg.path)) {
    # Extract svg names. 
    str <- strsplit(svg.path[[i]], '/')[[1]]
    na0 <- str[length(str)]
    if (!grepl(paste0('\\', c('.svg', raster.ext), '$', collapse='|'), na0)) return('No aSVG/template file is detected! Solution: 1) select another aSVG and rematch it to data; 2) add an aSVG/template file for the selected data in the backend aSVG tar file or uploaded aSVG tar file.')
    svg.na <- c(svg.na, na0)
    # Complete uploaded svg paths.
    if (!grepl('data/', svg.path[[i]])) {
      # The data/svg precedence: uploaded tar > internal tar > default examples. The duplicated data/svgs are removed according to this precedence when processing data/svg upstream.
      pa0 <- NULL; if (!is.null(pa.svg.upl)) pa0 <- extr_svg(file=pa.svg.upl, name=na0)
      if (is.null(pa.svg.upl)|is.null(pa0)) {
        tar.all <- list.files('data', pattern='\\.tar$', full.names=TRUE)
        tar.svg <- tar.all[!grepl('data_shm.tar$', tar.all)][1]
        pa0 <- extr_svg(file=tar.svg, name=na0)
      }; if (is.null(pa0)) return(paste0("This aSVG/template file is not detected: ", na0, "!")) else svg.path[[i]] <- pa0
    }
  }; return(list(svg.path=svg.path, svg.na=svg.na))
}

# Check suffixes if multiple svgs.
svg_suffix <- function(svg.path, svg.na, raster.ext) {
  if (length(svg.na)>1) {
    ext <- paste0('_shm\\d+\\', c('.svg', raster.ext), '$', collapse='|')
    if (!all(grepl(ext, svg.na, perl=TRUE))) return("Suffixes of aSVGs and templates should be indexed as '_shm1.svg', '_shm1.png', '_shm2.svg', '_shm2.png', '_shm3.svg', '_shm3.png', ...")
    ord <- order(gsub('.*_(shm.*)$', '\\1', svg.na))
    svg.path <- svg.path[ord]; svg.na <- svg.na[ord]  
  }; return(list(svg.path=svg.path, svg.na=svg.na))
}

# Convert aSVG features to draggable items.
# ns() is the namespace in shiny modules.
ft2tag <- function(ft){
  lapply(ft, function(i) { tag("span", list(class = class(i), tags$span(class = "glyphicon glyphicon-move"), i)) }
  )
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

ft_lis_dat <- function(x, ns) {
 span(class = "panel panel-default",
   div(class = "panel-heading", names(x)), 
   div(class = "panel-body", id = ns(names(x)), ft2tag(x[[1]]))
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

# Complete matching interface.
match_interface <- function(to.ft, to.div.id='ftSVG', to.div.tit='Features in aSVG', from.ft, from.div.tit='Features in data',   ns) {
  to.ft <- sort(to.ft); from.ft <- sort(from.ft) 
  frow <- fluidRow(
    span(class = "panel panel-default", style = 'margin-left:0px',
      div(class = "panel-heading", strong(to.div.tit)), 
      div(class = "panel-body", id = ns(to.div.id), ft2tag(to.ft)) 
      ),
    div(class = "panel panel-default", 
      div(class = "panel-heading", strong(from.div.tit)),  
      lapply(from.ft, ft_dat, ns = ns)  
    ), lapply(c(to.div.id, from.ft), ft_js, ns = ns) # Items are interchangeable across ftSVG and sam.all.
  ); return(frow) 
} 

# Clean temporary files in data mining.
data_mining_rm <- function() {
  if (dir.exists('www/tmp/')) {
    message("Removing files in 'www/tmp/' ... \n")
    files <- list.files('www/tmp/', '*.txt$', full.names=TRUE)
    file.remove(grep('README', list.files('www/tmp/', '*.txt$', full.names=TRUE), invert=TRUE, value=TRUE))   
  } 
}

# Clean trash files in animation and video.
ggly_rm <- function() {
  if (dir.exists('www/html_shm')) {
    cat("Removing animation files in 'www/html_shm' ... \n")
    unlink('www/html_shm/lib', recursive=TRUE)
    file.remove(list.files('www/html_shm/', '*.html$', full.names=TRUE))
  } else dir.create('www/html_shm', recursive=TRUE)
  if (dir.exists('R/www')) {
    cat("Removing animation files in 'R/www/html_shm' ... \n")
    unlink('R/www', recursive=TRUE)
  }
}
vdo_rm <- function() {
  if (dir.exists('www/video/')) {
    cat("Removing video file in 'www/video/' ... \n")
    file.remove(list.files('www/video/', '*.mp4$', full.names=TRUE))
  } else dir.create('www/video/', recursive=TRUE)
 if (dir.exists('R/www')) {
    cat("Removing video file in 'R/www/video' ... \n")
    unlink('R/www', recursive=TRUE)
 }
}






