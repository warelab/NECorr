#' NetTopology
#' @title Network topology statistics
#' @description Calculates centrality and connectivity metrics for each gene in the network.
#' @param network Data frame with 'source' and 'target' columns
#' @return Data frame of topology metrics for each gene
#' @export
NetTopology <- function(network) {
  requireNamespace("igraph")

  g <- graph.data.frame(network, directed = TRUE)

  data.frame(
    EdgeCount = degree(g, v = V(g), mode = "total", normalized = TRUE),
    BetweennessCentrality = betweenness(g, directed = FALSE, normalized = TRUE),
    eigenvector_centrality = evcent(g, scale = TRUE)$vector,
    ClusteringCoefficient = transitivity(g, type = "local", isolates = "zero"),
    pagerank = page.rank(g)$vector,
    row.names = V(g)$name
  )
}
