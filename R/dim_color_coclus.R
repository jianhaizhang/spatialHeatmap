#' Map colours in SHMs to embedding plots 
#'
#' True matching between cells and bulk tissues is not included in \code{sce}.
#'
#' @param sce A \code{SingleCellExperiment} containing the (un-)aggregated cells that have source tissue assignments. The \code{lis.match} will be built internally.
#' @param row.sel The selected row indexes in the table of \code{colData} slot in Shiny app.
#' @param gg.dim The ggplot of embedding plot.
#' @param gg.shm.all The list of SHM ggplot.
#' @param grob.shm.all The list of SHM grob.
#' @param col.shm.all The list of SHM colours.
#' @param color.by A column name in the \code{colData} slot such as \code{label}.  
#' @param lgd.all.dim Logical. The default is \code{FALSE}, and only cells with source tissue assignment have legends.
#' @param con.na Logical, TRUE or FALSE. Default is TRUE, meaning conditions are available.
#' @param lis.match The remaching list of spatial features in data and aSVGs.
#' @param sub.title.size The title size of embedding plots. The default is 11.
#' @param dim.lgd.pos The legend position. The default is \code{bottom}.
#' @param dim.lgd.nrow The number of legend rows. The default is \code{1}.
#' @param dim.lgd.text.size The size of legend text. The default is \code{10}.

#' @return A nested list of embedding and SHM plots.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Morgan M, Obenchain V, Hester J, Pagès H (2021). SummarizedExperiment: SummarizedExperiment container. R package version 1.24.0, https://bioconductor.org/packages/SummarizedExperiment.
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016
#' R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

#' @importFrom SummarizedExperiment colData
#' @importFrom methods as
#' @importFrom ggplot2 layer_data ggplot geom_point theme_classic theme element_text element_blank labs scale_shape_manual scale_color_manual margin guide_legend

dim_color_coclus <- function(sce=NULL, row.sel=NULL, targ=NULL, profile=FALSE, gg.dim, gg.shm.all, grob.shm.all, gg.lgd.all, col.shm.all, col.lgd.all, grob.lgd.all, con.na=TRUE, lis.match=NULL, sub.title.size=11, dim.lgd.pos='bottom', dim.lgd.nrow=1, dim.lgd.key.size=4, dim.lgd.text.size=13, dim.capt.size=13) {
  save(sce, row.sel, targ, profile, gg.dim, gg.shm.all, grob.shm.all, gg.lgd.all, col.shm.all, col.lgd.all, grob.lgd.all, con.na, lis.match, sub.title.size, dim.lgd.pos, dim.lgd.nrow, dim.lgd.text.size, dim.lgd.key.size, dim.capt.size, file='dim.color.coclus.arg')
  response <- feature <- idx <- x <- y <- NULL
  if (!is.null(sce) & !is.null(lis.match)) stop("Only one of 'sce' and 'lis.match' is required!")
  cdat <- colData(sce)
  if (!is.null(sce)) {
    # The matching list between aggregated cells and aSVG spatial features. The former are cells with a source tissue assignment in co-clustering.   
    blk.uni <- unique(cdat$SVGBulk)
    blk.uni <- setdiff(blk.uni, 'none')
    lis.match <- as(blk.uni, 'list'); names(lis.match) <- blk.uni
    if (targ[1]=='matched') targ <- blk.uni
    if (any(!targ %in% blk.uni)) stop("Make sure all entries in 'targ' are in 'SVGBulk'!")
  }
  lis.match <- lis.match[!unlist(lapply(lis.match, is.null))]
  
  # Ggplots of all reduced dim.
  if (profile==TRUE) {
    n <- length(gg.shm.all); gg.dim.all <- rep(list(gg.dim), n)
    names(gg.dim.all) <- paste0('dim_', names(gg.shm.all))
  } else if (profile==FALSE) {
    n <- length(gg.lgd.all); gg.dim.all <- rep(list(gg.dim), n)
    names(gg.dim.all) <- paste0('dim_', names(gg.lgd.all))
  }
   for (i in seq_along(gg.dim.all)) {  

  # Map colours in SHMs to dim plots.
  col_dim_shm <- function(gg.dim, gcol.all, lis.match) {
    na <- sub('^dim_', '', names(gg.dim))
    g.col <- gcol.all[[paste0('col_', na)]]
    dat.ft.na <- names(lis.match)
    # 'gray80' is a reserved color.
    # In the case of multiple development stages, colours in all SHM/SVGs are mapped to dim.col together. dim.col is similar to the data table, where all SVG features are included. 
    dim.col <- rep('gray80', length(lis.match))
    names(dim.col) <- dat.ft.na
    for (j in dat.ft.na) {
      # Matched svg fts have the same color, so only the 1st is taken. 
      ft.svg <- lis.match[[j]][1]; matched <- g.col[ft.svg]
      if (length(matched)==0) next else if (is.na(matched)) {
        matched <- g.col[sub('__\\d+', '', names(g.col))==ft.svg][1]
        if (!is.na(matched)) if (matched!='NA') dim.col[j] <- matched
      } else if (length(matched)>0) { if (matched!='NA') dim.col[j] <- matched }
    }; return(dim.col)
  }
  gg.dim <- gg.dim.all[i]
  if (profile==TRUE) dim.col <- col_dim_shm(gg.dim=gg.dim, gcol.all=col.shm.all, lis.match=lis.match)
  if (profile==FALSE) dim.col <- col_dim_shm(gg.dim=gg.dim, gcol.all=col.lgd.all, lis.match=lis.match)
    # Target cells without assignments ('none') are assigned 'gray50'.
    true.tar <- unique(subset(cdat, SVGBulk!='none')$SVGBulk)
    dim.col[setdiff(targ, true.tar)] <- 'gray50'
    
    gg.dim0 <- gg.dim[[1]]
    dat <- gg.dim0$data; lay.dat <- layer_data(gg.dim0) 
    df.all <- cbind(lay.dat, colour_by=dat$colour_by)
    df.all$fill <- df.all$colour <- dim.col[df.all$colour_by]
    df.all <- df.all[, !colnames(df.all) %in% 'colour']

    # The row order is the same between cdat and df.all.
    df.all$feature <- cdat$SVGBulk
    df.all <- cbind(df.all, cdat[, c('cell', 'assignedBulk', 'dataBulk')])
    df.all$idx <- seq_len(nrow(df.all))
    na.sel <- paste0(unique(df.all$feature[row.sel]), '.selected')

    # Target cells are place on top of non-target cells.
    idx.top <- which(df.all$feature != 'none')
    df.all <- rbind(df.all[setdiff(seq_len(nrow(df.all)), idx.top), ], df.all[idx.top, , drop=FALSE])
    # Selected cells are placed on top.
    if (!is.null(row.sel)) {
      # The cells are divided into "selected" and non-selected. The "other" cells in the coclustering results are treated as ordinary cells.
      df.all$feature <- sub('\\.other$', '', df.all$feature)
      un.sel.idx <- setdiff(df.all$idx, row.sel)
      # "row.sel" referrs to original row index, i.e. df.all$idx.
      df.all <- rbind(subset(df.all, idx %in% un.sel.idx), subset(df.all, idx %in% row.sel))
      idx.sel <- df.all$idx %in% row.sel
      df.all$feature[idx.sel] <- paste0(df.all$feature[idx.sel], '.selected')
    }

    dim.col.final <- dim.col
    if (is.null(row.sel)) {
      # Non-target cells are assigned colour gray80.
      dim.col.final[!names(dim.col.final) %in% targ] <- 'gray80'
    } else { # Colours of selected cells.
      col.sel <- rep('gray50', length(na.sel))
      names(col.sel) <- na.sel
      dim.col.final[names(dim.col.final)] <- 'gray80'
      dim.col.final <- c(dim.col.final, col.sel)
    }

    # Non-target cells are assigned shape 1.
    sp <- rep(1, length(dim.col.final))
    names(sp) <- names(dim.col.final)

    # Legal shapes: c(0:25, 32:127)
    sp.sel <- c(15:18, 7:14)
    sp.all <- c(0, 2:25, 32:127)
    sp.all <- c(sp.sel, setdiff(sp.all, sp.sel))
    if (is.null(row.sel)) {
      sp[targ] <- sp.all[seq_along(targ)]
    } else {
      sp[na.sel] <- sp.all[seq_along(na.sel)]
    }
    # Cells without true bulk in the matching data frame.
    if ('none' %in% df.all$dataBulk) {
      dim.col.final <- c('gray80', dim.col.final)
      names(dim.col.final)[1] <- 'none'
      sp <- c(1, sp); names(sp)[1] <- 'none'
    }
  if (con.na==TRUE) tit <- gsub('^dim_(.*)_\\d+$', '\\1', names(gg.dim)) else tit <- gsub('^dim_(.*)_con_\\d+$', '\\1', names(gg.dim))
  if (profile==FALSE) tit <- NULL
  if (is.null(row.sel)) {
    lgd.tit <- ''; lgd.show <- targ
    br <- targ 
  } else {
    lgd.tit <- 'Selected'; lgd.show <- sub('\\.selected$', '', na.sel) 
    br <- na.sel
  }
  # Combine color and shape: use identical name and labels values for both shape and colour scale.
  gg <- ggplot(df.all, aes(x=x, y=y, text=df.all$feature)) + geom_point(size=2, alpha=1, aes(shape=feature, colour=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank(), legend.position=dim.lgd.pos, legend.text=element_text(size=dim.lgd.text.size), legend.margin = margin(t=-0.02, l=0.05, r=0.1, unit='npc'), plot.caption = element_text(hjust = 0, size=dim.capt.size), legend.title=element_text(size=dim.lgd.text.size+1)) + labs(title=tit, x=gg.dim0$labels$x, y=gg.dim0$labels$y, caption = NULL) + scale_colour_manual(name=lgd.tit, values=dim.col.final, breaks=br, labels=lgd.show, guide=guide_legend(title=lgd.tit, nrow=dim.lgd.nrow)) + scale_shape_manual(name=lgd.tit, values=sp, breaks=br, labels=lgd.show, guide=guide_legend(title=lgd.tit, nrow=dim.lgd.nrow, override.aes = list(size=dim.lgd.key.size)))
  gg.dim.all[[i]] <- gg
  } 
  
  # Convert all reduced dim of ggplots to grobs.  
  grob.dim.all <- grob_shm(gg.dim.all, lgd.pos=NULL)
  # Empty list of all reduced dim and SHMs. 
  dim.shm.gg.lis <- dim.shm.grob.lis <- rep(list(NULL), 2*n)

  # Assign all reduced dims to the empty list.
  dim.shm.gg.lis[seq(1, 2*n, 2)] <- gg.dim.all
  dim.shm.grob.lis[seq(1, 2*n, 2)] <- grob.dim.all
  names(dim.shm.gg.lis)[seq(1, 2*n, 2)] <- names(dim.shm.grob.lis)[seq(1, 2*n, 2)] <- names(grob.dim.all)
  # Assign all SHMs to the empty list.
  if (profile==TRUE) {
    dim.shm.gg.lis[seq(2, 2*n, 2)] <- gg.shm.all
    dim.shm.grob.lis[seq(2, 2*n, 2)] <- grob.shm.all
    names(dim.shm.gg.lis)[seq(2, 2*n, 2)] <- names(dim.shm.grob.lis)[seq(2, 2*n, 2)] <- names(grob.shm.all)
  } else {
    dim.shm.gg.lis[seq(2, 2*n, 2)] <- gg.lgd.all
    dim.shm.grob.lis[seq(2, 2*n, 2)] <- grob.lgd.all
    names(dim.shm.gg.lis)[seq(2, 2*n, 2)] <- names(dim.shm.grob.lis)[seq(2, 2*n, 2)] <- names(grob.lgd.all)

  }
  return(list(dim.shm.gg.lis=dim.shm.gg.lis, dim.shm.grob.lis=dim.shm.grob.lis))
}
