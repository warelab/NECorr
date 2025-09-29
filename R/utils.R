# =========================
#' hub_edge_significance
#'
#' @param network.int network file
#' @param gene.rank.vec gene ranking hash
#' @description define the hub edge ranking
#' @return res
hub_edge_significance <- function(network.int, gene.rank.vec) {
  # Ensure gene.rank.vec is named
  stopifnot(!is.null(names(gene.rank.vec)))
  # Ensure all genes in network are in gene.rank.vec
  if (!all(c(network.int[, 1], network.int[, 2]) %in% names(gene.rank.vec))) {
    stop("All genes in network.int must be present in names of gene.rank.vec")
  }
  # Calculate sum of ranks for each edge
  gene.rank.vec <- as.numeric(gene.rank.vec) |> setNames(names(gene.rank.vec))
  edge_scores <- gene.rank.vec[network.int[, 1]] + gene.rank.vec[network.int[, 2]]
  break.points <- c(-Inf, sort(unique(edge_scores)), Inf)
  p2 <- 1 - cut(edge_scores, breaks = break.points, labels = FALSE) / length(break.points)

  data.frame(sourceIDs = network.int[, 1],
             targetIDs = network.int[, 2],
             ranks.sum = edge_scores,
             p2 = as.numeric(p2))
}

#' activator_significant
#'
#' @param hub.int.significant significance of the hub genes
#' @param network.int network edges
#' @param Desc description file genes and gene names
#' @description define the activator significance
#' @return res
activator_significant <- function(hub.int.significant, network.int, Desc) {
  requireNamespace("dplyr")
  sig.hub <- unique(c(hub.int.significant$sourceIDs, hub.int.significant$targetIDs))
  net.extension <- network.int[network.int[, 1] %in% sig.hub | network.int[, 2] %in% sig.hub, ]

  g <- igraph::graph.data.frame(net.extension, directed = TRUE)

  datagrid <- data.frame(
    gene_name = igraph::V(g)$name,
    pagerank = igraph::page_rank(g)$vector,
    EdgeCount = igraph::degree(g, mode = "total", normalized = TRUE),
    degree_out = igraph::degree(g, mode = "out")
  )

  activators <- ScalN(setNames(datagrid$pagerank, datagrid$gene_name))
  colnames(activators)[2] <- "pagerank_scaled"
  # order

  dplyr::left_join(
    tibble::rownames_to_column(activators, "rowname"),
    Desc,
    by = c("rowname"="gene_id")
  ) %>%
    dplyr::select(-rowname)
}


#' linked_act_hub_net
#'
#' @param hub.int.significant significant network of the hub genes
#' @param gene.rank.act.significant the ranking of the activator genes
#' @param network.int full gene network
#' @return act.net
linked_act_hub_net <- function(hub.int.significant, gene.rank.act.significant, network.int) {
  sig.hub <- unique(c(hub.int.significant[, 1], hub.int.significant[, 2]))

  net.extension <- network.int[network.int[, 1] %in% sig.hub | network.int[, 2] %in% sig.hub, ]

  act.genes <- as.character(gene.rank.act.significant[, 1])

  act.net <- net.extension[net.extension[, 1] %in% act.genes | net.extension[, 2] %in% act.genes, ]
  meanSig <- mean(as.numeric(hub.int.significant[, 3]))

  data.frame(source = act.net[, 1], target = act.net[, 2],
             score = meanSig, node.type = "act")
}

#' linked_eff_hub_net
#'
#' @param hub.int.significant significant network of the hub genes
#' @param eff.int.significant significant network of the effector genes
#'
#' @return eff.net
linked_eff_hub_net <- function(hub.int.significant, eff.int.significant) {
  sig.hub <- unique(c(hub.int.significant[, 1], hub.int.significant[, 2]))

  eff.net <- eff.int.significant[eff.int.significant[, 1] %in% sig.hub |
                                   eff.int.significant[, 2] %in% sig.hub, ]

  data.frame(source = eff.net[, 1],
             target = eff.net[, 2],
             score = eff.net[, 3],
             node.type = "eff")
}


#' linked_eff_hub_net
#' @param hub.int.significant significant network of the hub genes
#' @param eff.int.significant significant network of the effector genes
#' @return netindex
linked_eff_hub_net <- function(hub.int.significant, eff.int.significant) {
  sig.hub <- unique(c(hub.int.significant[, 1], hub.int.significant[, 2]))

  eff.net <- eff.int.significant[eff.int.significant[, 1] %in% sig.hub |
                                   eff.int.significant[, 2] %in% sig.hub, ]

  data.frame(source = eff.net[, 1],
             target = eff.net[, 2],
             score = eff.net[, 3],
             node.type = "eff")
}

# =========================
# NECorr Utility Functions
# =========================

#' ScalN
#' @description Scales a numeric vector to a range from 0 to 1 using min-max scaling.
#' This is useful for normalising data so that all values fall within the same range.
#' @param x Numeric vector to be scaled.
#' @return A numeric vector of the same length as `x`, with values rescaled to fall between 0 and 1.
#' If all values in `x` are the same, returns a vector of zeros.
#' @export
ScalN <- function(x) {
  # Handle empty or NULL input
  if (is.null(x) || length(x) == 0) return(NULL)
  # Force numeric conversion
  x <- setNames(as.numeric(x), names(x))
  # Scale between 0 and 1
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    scaled <- rep(0, length(x))  # avoid division by zero
  } else {
    scaled <- (x - rng[1]) / diff(rng)
  }
  # Preserve names if present
  if (!is.null(names(x))) {
    df <- data.frame(Gene = names(x), value = scaled, stringsAsFactors = FALSE)
  } else {
    stop("ScalN requires named vector so it can merge by Gene")
  }
  return(df)
}


#' scaling_param
#' @description Scale named numeric vector from 0 to 100
#' @param x Named numeric vector
#' @return Data frame with scaled values
#' @export
scaling_param <- function(x) {
  rng <- range(x[!is.infinite(x)], na.rm = TRUE)
  x[is.infinite(x) & sign(x) < 0] <- rng[1]
  x[is.infinite(x)] <- rng[2]
  data.frame(val = scales::rescale(x, to = c(0, 100)))
}

#' fishersMethod
#' @description Combine p-values using Fisher's method.
#' @param x Numeric vector of p-values
#' @return Combined p-value
#' @export
# Fisher's method to combine p-values
fishersMethod <- function(x) {
  x <- x[!is.na(x) & x > 0]  # remove NAs and zeros to avoid log(0)
  if (length(x) == 0) return(NA_real_)
  pchisq(-2 * sum(log(x)), df = 2 * length(x), lower.tail = FALSE)
}

#' indexing.network
#' @description Replace gene names in network with their row indices in expression matrix
#' @param tab Expression matrix (rownames = genes)
#' @param network Two-column data frame of edges (source, target)
#' @return Matrix of same size as network with integer indices
#' @export
indexing.network <- function(tab, network) {
  cbind(match(network[, 1], rownames(tab)),
        match(network[, 2], rownames(tab)))
}

#' mergeMetricsByGene
#' @description Merge multiple gene metric data frames by 'Gene' column.
#' Handles NULL inputs and missing genes.
#' @param ... Data frames with 'Gene' column and metric columns
#' @return Merged data frame with all genes and metrics
#' @export
mergeMetricsByGene <- function(...) {
  dfs <- list(...)
  # Remove NULLs
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) {
    return(data.frame(Gene = character(0)))
  }
  # Check all have Gene column
  for (df in dfs) {
    if (!"Gene" %in% colnames(df)) {
      stop("All metric data frames must have a 'Gene' column")
    }
  }
  merged <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfs)
  return(merged)
}

#' safe_as_data_table
#' @description Safely convert a list to a data.table, padding shorter elements with NA if needed.
#' @param x List to convert
#' @param pad_with_NA Logical, if TRUE pads shorter elements with NA to match longest length
#' @return data.table
#' @export
#' @examples
#' #' # Example with equal lengths
#' lst_equal <- list(a = 1:3, b = letters[1:3])
#' safe_as_data_table(lst_equal)
#' #' # Example with differing lengths, padding with NA
#' lst_diff <- list(a = 1:3, b = letters[1:2])

safe_as_data_table <- function(x, pad_with_NA = TRUE) {
  # Check that input is a list
  if (!is.list(x)) {
    stop("Input must be a list.")
  }
  # Get lengths of each item
  lens <- lengths(x)
  # If all lengths are equal, convert directly
  if (length(unique(lens)) == 1) {
    return(as.data.table(x))
  }
  # If lengths differ
  message("Column lengths differ: ", paste(lens, collapse = ", "))
  if (!pad_with_NA) {
    stop("Lengths differ. Set pad_with_NA = TRUE to pad shorter columns with NA.")
  }
  # Pad shorter elements with NA to match the max length
  max_len <- max(lens)
  x <- lapply(x, function(col) {
    length(col) <- max_len
    col
  })
  return(as.data.table(x))
}
