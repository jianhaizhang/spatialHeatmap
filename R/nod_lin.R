#' Calculate Gene Connectivity and Edge Adjacency
#'
#' @param ds Module detecting sensativity, "2" or "3".
#' @param lab An interger, which is the module the target gene belongs to.
#' @param mods A 2-column data frame, where each coloumn is a module assignmet.
#' @param adj The complete adjacency matrix.
#' @param geneID The target gene.
#' @param adj.min The minimum adjacency to keep for the edges.
#' @return A 2-length list, containing the gene connectiviy and edge adjaceny respectively.
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

nod_lin <- function(ds, lab, mods, adj, geneID, adj.min) {

  from <- to <- NULL
  idx.m <- mods[, ds]==lab; adj.m <- adj[idx.m, idx.m]; gen.na <- colnames(adj.m) 
  idx.sel <- grep(paste0("^", geneID, "$"), gen.na); gen.na[idx.sel] <- paste0(geneID, "_target")
  colnames(adj.m) <- rownames(adj.m) <- gen.na; idx = adj.m > as.numeric(adj.min)
  link <- data.frame(from=rownames(adj.m)[row(adj.m)[idx]], to=colnames(adj.m)[col(adj.m)[idx]], width=adj.m[idx], stringsAsFactors=FALSE)
  # Should not exclude duplicate rows by "length".
  node.pas <- NULL; for (i in seq_len(nrow(link))) { node.pas <- c(node.pas, paste0(sort(c(link[i, 'from'], link[i, 'to'])), collapse='')) }
  w <- which(duplicated(node.pas)); link <- link[-w, ]
  link1 <- subset(link, from!=to, stringsAsFactors=FALSE); link1 <- link1[order(-link1$width), ]
  node <- data.frame(label=colnames(adj.m), size=colSums(adj.m), stringsAsFactors=FALSE)
  node <- node[order(-node$size), ]; return(list(node=node, link=link1))

}

