#' Detect subgraphs to define clusters
#'
#' @param graph A \code{igraph} object returned by \code{nn_graph}.
#' @param clustering The clustering method. One of \code{wt} (cluster_walktrap, default), \code{fg} (cluster_fast_greedy), \code{le} (cluster_leading_eigen), \code{sl} (cluster_spinglass), \code{eb} (cluster_edge_betweenness). 
#' @param wt.arg,fg.arg,sl.arg,le.arg,eb.arg A named list of arguments passed to \code{wt}, \code{fg}, \code{le}, \code{sl}, \code{eb} respectively.

#' @return A class of \code{communities}.
 
#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references
#' Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. https://igraph.org

#' @importFrom igraph cluster_walktrap cluster_fast_greedy cluster_leading_eigen cluster_spinglass cluster_edge_betweenness 

detect_cluster <- function(graph, clustering='fg', wt.arg=list(steps = 4), fg.arg=list(), sl.arg=list(spins = 25), le.arg=list(), eb.arg=list()) {
  cat('Scell: clustering ... \n') 
  clus <- tryCatch( 
  if (clustering=='wt') { 
    do.call('cluster_walktrap', c(list(graph=graph), wt.arg)) 
  } else if (clustering=='fg') { 
    do.call('cluster_fast_greedy', c(list(graph=graph), fg.arg)) 
  } else if (clustering=='le') { 
    do.call('cluster_leading_eigen', c(list(graph=graph), le.arg)) 
  } else if (clustering=='sl') { 
    do.call('cluster_spinglass', c(list(graph=graph), sl.arg)) 
  } else if (clustering=='eb') { 
    do.call('cluster_edge_betweenness', c(list(graph=graph), eb.arg))
  }, warning = function(w){ 'w' }, error = function(e){ 'e' } 
  )
  if (!is(clus, 'communities')) {
    message('Failed to detect clusters!'); return() 
  } else return(clus)
}
