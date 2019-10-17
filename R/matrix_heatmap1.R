# geneID='PSAC'; data=exp; adj.mod=adj_mod; ds="3"; scale="cl"; angleCol=45; angleRow=45; margin=c(8, 8)
matrix.heatmap1 <- function(geneID, data, adj.mod, ds="3", scale="row", angleCol=45, angleRow=45, margin=c(8, 8)) {

  library(gplots)
  mods <- adj.mod[["mod"]]

  # The last column of data is expected to be annotation.
  col.len <- ncol(data)
  if (!is.numeric(data[, col.len])) gene <- data[, -col.len] else gene <- data
  lab <- mods[, as.character(ds)][rownames(gene)==geneID]
  if (lab=="0") stop("The input gene is not assigned to any module. Please input a different one.")

  mod <- as.matrix(gene[mods[, ds]==lab, ])
  hm <- heatmap.2(mod, main=paste0('Netowrk module containing ', geneID), trace="none")
  # Logical matrix with the same dimensions as module matrix.
  selection <- matrix(rep(FALSE, nrow(mod)*ncol(mod)), nrow=nrow(mod))
  # Select the row of target gene.  
  idx <- which(rev(colnames(hm$carpet) %in% geneID))
  selection[idx, ] <- TRUE
  
  nx <- ncol(mod); ny <- nrow(mod)
  makeRects <- function(cells, nx, ny){

    coords = expand.grid(rev(seq_len(ny)), seq_len(nx))[cells,]
    xl=coords[, 2]-0.49; yb=coords[, 1]-0.49
    xr=coords[, 2]+0.49; yt=coords[, 1]+0.49
    rect(xl, yb, xr, yt, border="black", lwd=1)

  }
  # Option 1 to label a row.
  heatmap.2(mod, scale=scale, main=paste0('Gene module containing ', geneID), key=TRUE, trace="none", density.info="none", dendrogram="both", Rowv=TRUE, Colv=TRUE, srtCol=angleCol, srtRow=angleRow, margins=margin, rowsep=c(idx-1, idx), sepcolor="black", sepwidth=c(0.02,0.02))
  # Option 2 to label a row.
  heatmap.2(mod, scale=scale, main=paste0('Gene module containing ', geneID), key=TRUE, trace="none", density.info="none", dendrogram="both", Rowv=TRUE, Colv=TRUE, srtCol=angleCol, srtRow=angleRow, margins=margin, add.expr={makeRects(selection, nx, ny)})

}

# library(spatialHeatmap)
# data.path <- system.file("extdata/example", "root_expr_ann_row_gen.txt", package = "spatialHeatmap")
# exp <- filter.data(data=data.path, sep="\t", isRowGen=TRUE, pOA=c(0, 0), CV=c(0.1, 10000), dir="./")
# adj_mod <- adj.mod(data=exp, type="signed", minSize=15, dir="./")
# source("matrix_heatmap1.R")
# matrix.heatmap1(geneID='PSAC', data=exp, adj.mod=adj_mod, ds="3", scale="cl", angleCol=45, angleRow=45, margin=c(8, 8))

