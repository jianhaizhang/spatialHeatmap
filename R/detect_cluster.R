#' Detect subgraphs to define clusters
#'
#' @param graph A \code{igraph} object returned by \code{nn_graph}.
#' @param clustering The clustering method. One of \code{wt} (cluster_walktrap, default), \code{fg} (cluster_fast_greedy), \code{le} (cluster_leading_eigen), \code{sl} (cluster_spinglass), \code{ed} (cluster_edge_betweenness). 
#' @param wt.arg,fg.arg,sl.arg,le.arg,eb.arg A list of arguments passed to \code{wt}, \code{fg}, \code{le}, \code{sl}, \code{ed} respectively.

#' @return A class of \code{communities}.
 
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. https://igraph.org

#' @importFrom igraph cluster_walktrap cluster_fast_greedy cluster_leading_eigen cluster_spinglass cluster_edge_betweenness 

detect_cluster <- function(graph, clustering='wt', wt.arg=list(steps = 4), fg.arg=list(), sl.arg=list(spins = 25), le.arg=list(), eb.arg=list()) {
  cat('Scell: clustering ... \n')
  if (clustering=='wt') clus <- do.call('cluster_walktrap', c(list(graph=graph), wt.arg)) else if (clustering=='fg') clus <- do.call('cluster_fast_greedy', c(list(graph=graph),fg.arg)) else if (clustering=='le') clus <- do.call('cluster_leading_eigen', c(list(graph=graph), le.arg)) else if (clustering=='sl') clus <- do.call('cluster_spinglass', c(list(graph=graph), sl.arg)) else if (clustering=='eb') clus <- do.call('cluster_edge_betweenness', c(list(graph=graph), eb.arg))
  return(clus)
}
