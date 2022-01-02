dim_color_coclus <- function(sce=NULL, row.sel=NULL, tar.cell=NULL, profile=FALSE, gg.dim, gg.shm.all, grob.shm.all, gg.lgd.all, col.shm.all, col.lgd.all, grob.lgd.all, con.na=TRUE, lis.match=NULL, sub.title.size=11) {
 # save(sce, row.sel, tar.cell, profile, gg.dim, gg.shm.all, grob.shm.all, gg.lgd.all, col.shm.all, col.lgd.all, grob.lgd.all, con.na, lis.match, sub.title.size, file='dim.color.arg')
  if (!is.null(sce) & !is.null(lis.match)) stop("Only one of 'sce' and 'lis.match' is required!")
  cdat <- colData(sce)
  if (!is.null(sce)) {
    # The matching list between aggregated cells and aSVG spatial features. The former are cells with a source tissue assignment in co-clustering.   
    blk.uni <- unique(cdat$SVGBulk)
    blk.uni <- setdiff(blk.uni, 'none')
    lis.match <- as(blk.uni, 'list'); names(lis.match) <- blk.uni
    if (tar.cell[1]=='all') tar.cell <- blk.uni
    if (any(!tar.cell %in% blk.uni)) stop("Make sure all items in 'tar.cell' are in 'SVGBulk'!")
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
    # Target cells without true assignments are assigned 'gray50'.
    true.tar <- unique(subset(cdat, response==TRUE)$SVGBulk)
    dim.col[setdiff(tar.cell, true.tar)] <- 'gray50'
    
    gg.dim0 <- gg.dim[[1]]
    dat <- gg.dim0$data; lay.dat <- layer_data(gg.dim0) 
    df.all <- cbind(lay.dat, colour_by=dat$colour_by)
    df.all$fill <- df.all$colour <- dim.col[df.all$colour_by]
    df.all <- df.all[, !colnames(df.all) %in% 'colour']

    # The row order is the same between cdat and df.all.
    df.all$tissue <- df.all$feature <- cdat$SVGBulk
    df.all <- cbind(df.all, cdat[, c('cell', 'assignedBulk', 'response', 'trueBulk')])
    df.all$idx <- seq_len(nrow(df.all))
    na.sel <- paste0(unique(df.all$feature[row.sel]), '.selected')

# Separate true/false/other bulk assignments for each target cell type in the data frame for embedding plot.
tar_cell <- function(tar.cell, df.all, row.sel=NULL) {  
  lis <- NULL; for (i in tar.cell) {
    # True bulk items.
    tar.all <- subset(df.all, feature == i)
    true.bulk <- strsplit(tar.all$trueBulk[1], ',|;| ')[[1]]
    true.bulk <- true.bulk[true.bulk!='']
    # True/false/other index.
    idx.true <- which(df.all$assignedBulk %in% true.bulk & df.all$response == TRUE)
    idx.false <- which(df.all$assignedBulk %in% true.bulk & df.all$response == FALSE)
    idx.other.tar <- which(df.all$trueBulk == tar.all$trueBulk[1] & df.all$response %in% c(FALSE, 'none'))
    idx.all <- seq_len(nrow(df.all))

    # all.other <- df.all[setdiff(idx.all, c(idx.true, idx.false, idx.other.tar)), , drop=FALSE]
    # all.other$border.col <- 'NA'
    tar.other <- df.all[idx.other.tar, , drop=FALSE]
    tar.other$feature <- tar.other$tissue <- paste0(tar.other$feature, '.other')
    # tar.other$border.col <- 'NA'
    tar.other.col <- 'gray50'; na.other <- paste0(i, '.other')
    names(tar.other.col) <- na.other
    
    tar.false <- df.all[idx.false, , drop=FALSE] 
    # tar.false$border.col <- 'NA'
    
    tar.true <- df.all[idx.true, , drop=FALSE]
    # tar.true$border.col <- 'black'
    rate.true1 <- round(nrow(tar.true)/(nrow(tar.true) + nrow(tar.other)), 2)
    rate.true2 <- round(nrow(tar.true)/(nrow(tar.true) + nrow(tar.false)), 2)
    if (is.null(row.sel)) capt <- paste0('Black: cells with false or no bulk assignments, gray: cells corresponding to other bulk', '\n', 'True/(true + false + cells without assignments): ', rate.true1, ', ', 'True/(true + false): ', rate.true2) else capt <- NULL 
    lis0 <- list(list(tar.true=tar.true, tar.false=tar.false, tar.other=tar.other, idx.true=idx.true, idx.false=idx.false, idx.other.tar=idx.other.tar, capt=capt))
    names(lis0) <- i; lis <- c(lis0, lis) 
  }; return(lis)
}

    tar.lis <- tar_cell(tar.cell, df.all, row.sel=row.sel)
    # All true/other target cells are combined and treated as top layer in the embedding plot.
    idx.top <- df.top <- NULL; for (k in tar.cell) {
      lis0 <- tar.lis[[k]]
      idx.top <- c(idx.top, lis0$idx.other.tar, lis0$idx.true)
      df.top <- rbind(df.top, lis0$tar.other, lis0$tar.true)
    }
    # Legend key text of each target cell.
    capt <- NULL; for (k in tar.cell) {
      capt <- c(capt, tar.lis[[k]]$capt)
    }

    # Target cells are place on top of non-target cells.
    df.all <- rbind(df.all[setdiff(seq_len(nrow(df.all)), idx.top), ], df.top)
    # Selected cells are placed on top.
    if (!is.null(row.sel)) {
      # The cells are divided into "selected" and non-selected. The "other" cells in the coclustering results are treated as ordinary cells.
      df.all$feature <- sub('\\.other$', '', df.all$feature)
      un.sel.idx <- setdiff(df.all$idx, row.sel)
      # "row.sel" referrs to original row index, i.e. df.all$idx.
      df.all <- rbind(subset(df.all, idx %in% un.sel.idx), subset(df.all, idx %in% row.sel))
      idx.sel <- df.all$idx %in% row.sel
      df.all$feature[idx.sel] <- paste0(df.all$feature[idx.sel], '.selected')
    }; df.all$tissue <- df.all$feature

    dim.col.final <- dim.col
    if (is.null(row.sel)) {
      # Non-target cells are assigned colour gray80.
      dim.col.final[!names(dim.col.final) %in% tar.cell] <- 'gray80'
      # Other target cells are assigned colour gray50.
      tar.other.col <- rep('gray50', length(tar.cell))
      na.other <- paste0(tar.cell, '.other')
      names(tar.other.col) <- na.other
      dim.col.final <- c(tar.other.col, dim.col.final) 
    } else { # Colours of selected cells.
      col.sel <- rep('gray50', length(na.sel))
      names(col.sel) <- na.sel
      dim.col.final[names(dim.col.final)] <- 'gray80'
      dim.col.final <- c(dim.col.final, col.sel)
    }

    # Non-target cells are assigned shape 1.
    sp <- rep(1, length(dim.col.final))
    names(sp) <- names(dim.col.final)

    if (is.null(row.sel)) {
      # Legal shapes: c(0:25, 32:127)
      sp[tar.cell] <- c(0, 2:25, 32:127)[seq_along(tar.cell)]
      sp[na.other] <- sp[tar.cell]
    } else {
      sp[na.sel] <- c(0, 2:25, 32:127)[seq_along(na.sel)]
    }

    # gg.dim.all[[i]] <- ggplot(dat.all, aes(x=x, y=y, shape=feature, text=dat.all$feature)) + geom_point(size=2, alpha=1, aes(colour=feature)) + scale_color_manual(values=dim.col) + scale_shape_manual(values=sp) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank()) + labs(title=gsub('^dim_(.*)_\\d+$', '\\1', names(gg.dim)), x=gg.dim0$labels$x, y=gg.dim0$labels$y, colour=color.by, shape=color.by)
  if (con.na==TRUE) tit <- gsub('^dim_(.*)_\\d+$', '\\1', names(gg.dim)) else tit <- gsub('^dim_(.*)_con_\\d+$', '\\1', names(gg.dim))
  if (profile==FALSE) tit <- NULL
  # Re-plot dimensionlaity plot.
  # gg <- ggplot(dat.all, aes(x=x, y=y, text=dat.all$feature)) + geom_point(size=2, alpha=1, aes(colour=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank()) + labs(title=tit, x=gg.dim0$labels$x, y=gg.dim0$labels$y, colour=color.by, shape=color.by)
  # if (lgd.all.dim==FALSE) gg <- gg + scale_color_manual(values=dim.col.all, breaks=names(dim.col)) else gg <- gg + scale_color_manual(values=dim.col.all, breaks=names(dim.col.all))


  # gg <- ggplot(dat.all, aes(x=x, y=y, shape=feature, text=dat.all$feature)) + geom_point(size=2, alpha=1, aes(colour=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank()) + labs(title=tit, x=gg.dim0$labels$x, y=gg.dim0$labels$y, colour=color.by, shape=color.by) + scale_color_manual(values=dim.col.all.uni, breaks=names(dim.col)) + scale_shape_manual(values=sp.all.uni, breaks=names(dim.col)) 
  if (is.null(row.sel)) {
    lgd.tit <- 'True assignment'; lgd.show <- tar.cell
    br <- tar.cell 
  } else {
    lgd.tit <- 'Selected'; lgd.show <- sub('\\.selected$', '', na.sel) 
    br <- na.sel
  }
  gg <- ggplot(df.all, aes(x=x, y=y, text=df.all$feature)) + geom_point(size=2, alpha=1, aes(shape=feature, colour=feature)) + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=sub.title.size), axis.text = element_blank(), axis.ticks = element_blank(), plot.caption = element_text(hjust = 0)) + labs(title=tit, x=gg.dim0$labels$x, y=gg.dim0$labels$y, caption =capt) + scale_color_manual(name=lgd.tit, values=dim.col.final, breaks=br, labels=lgd.show) + scale_shape_manual(name=lgd.tit, values=sp, breaks=br, labels=lgd.show)
  gg.dim.all[[i]] <- gg

  }
  grob_shm <- spatialHeatmap:::grob_shm
  # Convert all reduced dim of ggplots to grobs.
  grob.dim.all <- grob_shm(gg.dim.all, lgd.pos='right')
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
